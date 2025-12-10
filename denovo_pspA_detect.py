#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
De-novo pspA detection with BLASTX + 200-bp upstream validation.
Outputs under {--work-dir}/01_denovo/ :
 - filtered_pspA_hits.fasta
 - summary_fixed_200bp_upstream.csv
 - denovo_summary.csv
Optional: de-novo NT tree (off by default).

Parameters (defaults from user decisions):
 BLASTX: -evalue 1e-6 -seg yes -max_target_seqs 10 -outfmt 6
 Filters: pident >= 75, AA length >= 200
 Upstream: 200 nt window; blastn; keep if pident >= 75 and aln_len >= 100;
 HighConfidence if pident >= 90 and aln_len >= 150

 Copy cap: max-hits-per-sample 3, auto-raise for separate loci:
 - hits on different contigs don't count toward the cap
 - same-contig hits allowed if non-overlapping and >= min-inter-hit-distance (default 2000 nt)
 - hard safety cap max-hits-hard-cap 6
"""
import os
import re
import tempfile
import argparse
import multiprocessing as mp
from typing import Dict, List, Tuple, Optional, Any
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pspA_utils import (
    run_cmd, blast_db_exists, blast_nucl_db_exists,
    safe_id, norm_interval, overlap_len
)

# ---------- helpers ----------
def enumerate_assemblies(assemblies_dir: str) -> List[str]:
    return sorted(os.path.splitext(f)[0] for f in os.listdir(assemblies_dir) if f.endswith(".fasta"))

def assemble_blastx_cmd(
    query_fa: str,
    protein_db: Optional[str],
    protein_subject_fasta: Optional[str],
    evalue: str,
    max_target_seqs: int,
    threads: int,
    seg: str = "yes"
) -> List[str]:
    outfmt = "6 qseqid sseqid pident length evalue bitscore qstart qend sstart send qframe"
    if protein_db and blast_db_exists(protein_db):
        return [
            "blastx", "-query", query_fa, "-db", protein_db,
            "-evalue", evalue, "-seg", seg, "-max_target_seqs", str(max_target_seqs),
            "-outfmt", outfmt, "-num_threads", str(threads)
        ]
    elif protein_subject_fasta and os.path.exists(protein_subject_fasta):
        return [
            "blastx", "-query", query_fa, "-subject", protein_subject_fasta,
            "-evalue", evalue, "-seg", seg, "-max_target_seqs", str(max_target_seqs),
            "-outfmt", outfmt, "-num_threads", str(threads)
        ]
    else:
        raise FileNotFoundError("No valid protein DB or subject FASTA provided for BLASTX.")

def assemble_blastn_cmd(
    query_fa: str,
    nucl_db: Optional[str],
    nucl_subject_fasta: Optional[str],
    evalue: str = "1e-10",
    task: str = "blastn",
    threads: int = 1
) -> List[str]:
    outfmt = "6 qseqid sseqid pident length evalue bitscore qstart qend sstart send"
    if nucl_db and blast_nucl_db_exists(nucl_db):
        return [
            "blastn", "-query", query_fa, "-db", nucl_db,
            "-task", task, "-evalue", evalue, "-outfmt", outfmt, "-num_threads", str(threads)
        ]
    elif nucl_subject_fasta and os.path.exists(nucl_subject_fasta):
        return [
            "blastn", "-query", query_fa, "-subject", nucl_subject_fasta,
            "-task", task, "-evalue", evalue, "-outfmt", outfmt, "-num_threads", str(threads)
        ]
    else:
        raise FileNotFoundError("No valid nucleotide DB or subject FASTA provided for BLASTN.")

def hit_row_to_dict(row: pd.Series) -> Dict[str, Any]:
    s = dict(row)
    (qs, qe) = norm_interval(int(s["qstart"]), int(s["qend"]))
    strand = "+" if int(s["qframe"]) > 0 else "-"
    return {
        "contig": str(s["qseqid"]),
        "sseqid": str(s["sseqid"]),
        "pident": float(s["pident"]),
        "length_aa": int(s["length"]),
        "bitscore": float(s["bitscore"]),
        "evalue": str(s["evalue"]),
        "qstart": qs,
        "qend": qe,
        "strand": strand,
        "qframe": int(s["qframe"])
    }

def select_nonoverlapping_hits(
    hit_rows: List[Dict[str, Any]],
    max_hits_per_sample: int = 3,
    min_inter_hit_distance: int = 2000,
    max_hits_hard_cap: int = 6
) -> List[Dict[str, Any]]:
    """
    Keep best-scoring hits; avoid overlaps on same contig/strand.
    Auto-raise cap for clearly separate loci (diff contigs or >= min_inter_hit_distance on same contig).
    """
    kept: List[Dict[str, Any]] = []
    per_contig: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    for h in sorted(hit_rows, key=lambda x: (x["bitscore"], x["length_aa"], x["pident"]), reverse=True):
        key = (h["contig"], h["strand"])
        s1, s2 = h["qstart"], h["qend"]
        # overlap check
        conflict = False
        min_sep = None
        for (k1, k2) in per_contig.get(key, []):
            ov = overlap_len(s1, s2, k1, k2)
            if ov > 0:
                conflict = True
                break
            sep = min(abs(s1 - k2), abs(k1 - s2))
            min_sep = sep if (min_sep is None or sep < min_sep) else min_sep
        if conflict:
            continue
        # cap logic
        total_kept = len(kept)
        diff_contig = key not in per_contig
        allow = False
        if diff_contig:
            allow = True
        else:
            if (min_sep is not None and min_sep >= min_inter_hit_distance) or (total_kept < max_hits_per_sample):
                allow = True
        if allow and total_kept < max_hits_hard_cap:
            kept.append(h)
            per_contig.setdefault(key, []).append((s1, s2))
    return kept

def extract_upstream_200(contig_seq: SeqRecord, hit: Dict[str, Any]) -> SeqRecord:
    """
    Return 200-nt upstream sequence oriented 5'->3' relative to the gene.
    """
    L = len(contig_seq)
    if hit["strand"] == "+":
        start = max(1, hit["qstart"] - 200)
        end = hit["qstart"] - 1
        if end < start:
            seq = Seq("")
        else:
            seq = contig_seq.seq[start - 1:end]
        return SeqRecord(seq, id=f"{contig_seq.id}:{start}-{end}(up200)+", description="")
    else:
        start = hit["qend"] + 1
        end = min(L, hit["qend"] + 200)
        if end < start:
            seq = Seq("")
        else:
            seq = contig_seq.seq[start - 1:end].reverse_complement()
        return SeqRecord(seq, id=f"{contig_seq.id}:{start}-{end}(up200)-", description="")

def extract_coding_nt(contig_seq: SeqRecord, hit: Dict[str, Any]) -> SeqRecord:
    """
    Extract coding nucleotide for the BLASTX HSP, oriented 5'->3' relative to the gene.
    """
    start, end = hit["qstart"], hit["qend"]
    seq = contig_seq.seq[start - 1:end]
    if hit["strand"] == "-":
        seq = seq.reverse_complement()
    rid = f"{contig_seq.id}:{start}-{end}({hit['strand']})"
    return SeqRecord(seq, id=rid, description="")

def validate_upstream_blastn(
    up_200: SeqRecord,
    upstream_db: Optional[str],
    upstream_subject_fasta: Optional[str],
    threads: int = 1,
    min_pident: float = 75.0,
    min_length: int = 100
) -> Dict[str, Any]:
    """
    Run blastn for the 200-bp upstream and return best HSP metrics.
    """
    if len(up_200.seq) == 0:
        return {"kept": False, "reason": "NoUpstream", "pident": 0.0, "length": 0, "bitscore": 0.0, "evalue": ""}
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".fa") as tmpq:
        SeqIO.write(up_200, tmpq.name, "fasta")
    tmp_out = tmpq.name + ".blastn.txt"
    cmd = assemble_blastn_cmd(tmpq.name, upstream_db, upstream_subject_fasta, evalue="1e-10", task="blastn", threads=threads)
    run_cmd(cmd, stdout_path=tmp_out, check=True)
    cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qstart", "qend", "sstart", "send"]
    try:
        df = pd.read_csv(tmp_out, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        return {"kept": False, "reason": "NoHSP", "pident": 0.0, "length": 0, "bitscore": 0.0, "evalue": ""}
    if df.empty:
        return {"kept": False, "reason": "NoHSP", "pident": 0.0, "length": 0, "bitscore": 0.0, "evalue": ""}
    df.columns = cols
    hit = df.sort_values(by="bitscore", ascending=False).iloc[0].to_dict()
    pident = float(hit["pident"])
    alen = int(hit["length"])
    bits = float(hit["bitscore"])
    ev = str(hit["evalue"])
    kept = (pident >= float(min_pident) and alen >= int(min_length))
    hi_conf = (pident >= max(90.0, float(min_pident)) and alen >= max(150, int(min_length)))
    return {"kept": kept, "HighConfidence": hi_conf, "pident": pident, "length": alen, "bitscore": bits, "evalue": ev}

def _process_one_sample(
    sample: str,
    assemblies_dir: str,
    protein_db: Optional[str],
    protein_subject_fasta: Optional[str],
    upstream_db: Optional[str],
    upstream_subject_fasta: Optional[str],
    out_dir_denovo: str,
    evalue: str, seg: str, max_target_seqs: int,
    min_pident: float, min_aa_len: int,
    max_hits_per_sample: int, min_inter_hit_distance: int, max_hits_hard_cap: int,
    blastx_threads: int, blastn_threads: int,
    upstream_min_pident: float, upstream_min_length: int,
    only_high_confidence: bool
) -> Dict[str, Any]:
    """
    Worker: run BLASTX on one assembly, select hits, validate upstream, write per-sample temp outputs.
    Returns a dict with counts and the kept SeqRecords (coding NT).
    """
    sample_fa = os.path.join(assemblies_dir, f"{sample}.fasta")
    if not os.path.exists(sample_fa):
        return {"sample": sample, "status": f"AssemblyNotFound: {sample_fa}"}
    # 1) BLASTX
    blastx_cmd = assemble_blastx_cmd(
        sample_fa, protein_db, protein_subject_fasta,
        evalue=evalue, max_target_seqs=max_target_seqs, threads=blastx_threads, seg=seg
    )
    tmp_blastx = os.path.join(out_dir_denovo, "logs", f"{sample}.blastx.txt")
    os.makedirs(os.path.dirname(tmp_blastx), exist_ok=True)
    run_cmd(blastx_cmd, stdout_path=tmp_blastx, check=True)

    # Parse BLASTX
    try:
        df = pd.read_csv(tmp_blastx, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()
    cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qstart", "qend", "sstart", "send", "qframe"]
    if not df.empty:
        df.columns = cols
        # Filter by AA identity and AA length
        df = df[(df["pident"] >= min_pident) & (df["length"] >= min_aa_len)]
        df = df[cols]
    else:
        df = pd.DataFrame(columns=cols)

    total_hsps = len(df)
    if total_hsps == 0:
        return {"sample": sample, "total_hsps": 0, "kept": 0, "upstream_confirmed": 0, "high_conf": 0, "seqs": []}

    # 2) Select non-overlapping hits
    hit_dicts = [hit_row_to_dict(row) for _, row in df.iterrows()]
    selected = select_nonoverlapping_hits(
        hit_dicts,
        max_hits_per_sample=max_hits_per_sample,
        min_inter_hit_distance=min_inter_hit_distance,
        max_hits_hard_cap=max_hits_hard_cap
    )

    # Load contigs
    contigs = {rec.id: rec for rec in SeqIO.parse(sample_fa, "fasta")}
    kept_records: List[Tuple[SeqRecord, Dict[str, Any], Dict[str, Any]]] = []  # (coding_nt, hit, upstream_eval)
    upstream_csv_rows: List[Dict[str, Any]] = []
    upstream_confirmed = 0
    high_conf = 0

    # 3) Upstream validation and coding extraction
    for h in selected:
        contig = contigs.get(h["contig"])
        if contig is None:
            continue
        up200 = extract_upstream_200(contig, h)
        up_eval = validate_upstream_blastn(up200, upstream_db, upstream_subject_fasta,
                                           threads=blastn_threads, min_pident=upstream_min_pident, min_length=upstream_min_length)
        row = {
            "Sample": sample,
            "Contig": h["contig"],
            "Strand": h["strand"],
            "qstart": h["qstart"],
            "qend": h["qend"],
            "BLASTX_pident": h["pident"],
            "BLASTX_lengthAA": h["length_aa"],
            "BLASTX_bitscore": h["bitscore"],
            "Up_pident": up_eval.get("pident", 0.0),
            "Up_length": up_eval.get("length", 0),
            "Up_bitscore": up_eval.get("bitscore", 0.0),
            "Up_evalue": up_eval.get("evalue", ""),
            "Up_kept": up_eval.get("kept", False),
            "Up_HighConfidence": up_eval.get("HighConfidence", False)
        }
        upstream_csv_rows.append(row)

        # *** NEW BEHAVIOR: keep only hits whose upstream passes the cutoff,
        # and if only_high_confidence=True, require HighConfidence=True ***
        if up_eval.get("kept", False) and (not only_high_confidence or up_eval.get("HighConfidence", False)):
            upstream_confirmed += 1
            if up_eval.get("HighConfidence", False):
                high_conf += 1
            coding = extract_coding_nt(contig, h)
            kept_records.append((coding, h, up_eval))
        # else: drop this hit (still recorded in the upstream CSV)

    # Write per-sample upstream CSV (append; header if new)
    up_csv = os.path.join(out_dir_denovo, "summary_fixed_200bp_upstream.csv")
    hdr_needed = not os.path.exists(up_csv)
    with open(up_csv, "a") as f:
        if hdr_needed:
            f.write(",".join(list(upstream_csv_rows[0].keys()) if upstream_csv_rows else
                             ["Sample", "Contig", "Strand", "qstart", "qend", "BLASTX_pident", "BLASTX_lengthAA",
                              "BLASTX_bitscore", "Up_pident", "Up_length", "Up_bitscore", "Up_evalue",
                              "Up_kept", "Up_HighConfidence"]) + "\n")
        for r in upstream_csv_rows:
            f.write(",".join([str(r[k]) for k in r.keys()]) + "\n")

    # Order of kept_records respects selection (best first)
    return {
        "sample": sample,
        "total_hsps": int(total_hsps),
        "kept": int(len(selected)),
        "upstream_confirmed": int(upstream_confirmed),
        "high_conf": int(high_conf),
        "seqs": [coding for (coding, _, _) in kept_records]
    }

def detect_denovo_pspA(
    assemblies_dir: str,
    protein_db: Optional[str],
    protein_subject_fasta: Optional[str],
    upstream_db: Optional[str],
    upstream_subject_fasta: Optional[str],
    work_dir: str,
    threads: int = 12,
    jobs: Optional[int] = None,
    blastx_threads: Optional[int] = None,
    min_pident: float = 75.0,
    min_aa_len: int = 200,
    evalue: str = "1e-4",
    seg: str = "yes",
    max_target_seqs: int = 10,
    max_hits_per_sample: int = 3,
    min_inter_hit_distance: int = 2000,
    max_hits_hard_cap: int = 6,
    build_denovo_nt_tree: bool = False,
    denovo_nt_ref_fasta: Optional[str] = None,
    upstream_min_pident: float = 75.0,
    upstream_min_length: int = 100,
    only_high_confidence: bool = False
) -> str:
    """
    Run de-novo detection across all assemblies and write outputs.
    Returns path to filtered_pspA_hits.fasta.

    only_high_confidence:
        If True, only hits whose upstream BLASTN evaluation meets the HighConfidence criterion
        (pident >= 90 and length >= 150) are retained and reported.
    """
    # cap threads to <= 12
    threads = max(1, min(12, int(threads)))
    out_denovo = os.path.join(work_dir, "01_denovo")
    os.makedirs(out_denovo, exist_ok=True)
    os.makedirs(os.path.join(out_denovo, "logs"), exist_ok=True)
    assemblies = enumerate_assemblies(assemblies_dir)
    if len(assemblies) == 0:
        raise FileNotFoundError(f"No .fasta assemblies found in {assemblies_dir}")

    # Auto concurrency split
    if jobs is None or jobs <= 0:
        jobs = max(1, min(len(assemblies), threads // 2 if threads >= 4 else 1))
    if blastx_threads is None or blastx_threads <= 0:
        blastx_threads = max(1, threads // jobs) if jobs > 0 else threads
    blastn_threads = max(1, min(threads, 2))

    # Parallel per-sample
    worker_args = []
    for sample in assemblies:
        worker_args.append((
            sample, assemblies_dir, protein_db, protein_subject_fasta,
            upstream_db, upstream_subject_fasta, out_denovo,
            evalue, seg, max_target_seqs,
            min_pident, min_aa_len,
            max_hits_per_sample, min_inter_hit_distance, max_hits_hard_cap,
            blastx_threads, blastn_threads,
            upstream_min_pident, upstream_min_length,
            only_high_confidence
        ))
    results: List[Dict[str, Any]] = []
    with mp.Pool(processes=max(1, int(jobs))) as pool:
        for rec in pool.starmap(_process_one_sample, worker_args):
            results.append(rec)

    # Aggregate filtered sequences
    filtered_fa = os.path.join(out_denovo, "filtered_pspA_hits.fasta")
    total_written = 0
    with open(filtered_fa, "w") as fout:
        for rec in results:
            sid = rec.get("sample")
            seqs: List[SeqRecord] = rec.get("seqs", [])
            # write all kept (now HighConfidence-only if enabled)
            for i, srec in enumerate(seqs, start=1):
                out_id = f"{sid}__hit{i}"
                srec.id = out_id
                srec.description = ""
                SeqIO.write(srec, fout, "fasta")
                total_written += 1

    # Write per-assembly summary
    denovo_summary_csv = os.path.join(out_denovo, "denovo_summary.csv")
    with open(denovo_summary_csv, "w") as f:
        f.write("Assembly,TotalHSPs,SelectedNonOverlap,UpstreamConfirmed,HighConfidence,SequencesWritten\n")
        for rec in sorted(results, key=lambda r: r.get("sample", "")):
            f.write(",".join([
                str(rec.get("sample", "")),
                str(rec.get("total_hsps", 0)),
                str(rec.get("kept", 0)),
                str(rec.get("upstream_confirmed", 0)),
                str(rec.get("high_conf", 0)),
                str(len(rec.get("seqs", [])))
            ]) + "\n")

    # Optional: de-novo NT tree (off by default)
    if build_denovo_nt_tree and os.path.exists(filtered_fa):
        if not denovo_nt_ref_fasta or not os.path.exists(denovo_nt_ref_fasta):
            print("[WARN] --denovo-nt-tree enabled but reference fasta not found; skipping tree.")
        else:
            _build_denovo_nt_tree(
                out_denovo=out_denovo,
                filtered_fa=filtered_fa,
                ref_fa=denovo_nt_ref_fasta,
                threads=max(1, threads // 2),
                best_only=True  # unchanged: tree uses best per sample
            )

    return filtered_fa

def _build_denovo_nt_tree(
    out_denovo: str,
    filtered_fa: str,
    ref_fa: str,
    threads: int = 4,
    best_only: bool = True
) -> None:
    """
    Build MAFFT+FastTree NT tree from filtered de-novo hits + provided NT references.
    If best_only=True, include only the best hit per sample (assumed to be '__hit1').
    """
    import subprocess
    from ete3 import Tree, TreeStyle, TextFace, NodeStyle
    threads = max(1, min(12, int(threads)))
    d_tree = os.path.join(out_denovo, "nt_tree")
    os.makedirs(d_tree, exist_ok=True)
    combined = os.path.join(d_tree, "denovo_plus_refs.fna")
    aln = os.path.join(d_tree, "denovo_plus_refs.aln.fna")
    nwk = os.path.join(d_tree, "denovo_plus_refs.nwk")
    png = os.path.join(d_tree, "denovo_nt_tree.png")

    # combine: references + de-novo (optionally best per sample)
    with open(combined, "w") as f:
        # references (prefix for styling)
        for rec in SeqIO.parse(ref_fa, "fasta"):
            rec.id = "Ref\n" + rec.id.split()[0]
            rec.description = ""
            SeqIO.write(rec, f, "fasta")
        # de-novo sequences
        if best_only:
            for rec in SeqIO.parse(filtered_fa, "fasta"):
                rid = rec.id.split()[0]
                if rid.endswith("__hit1"):
                    rec.id = safe_id(rec.id)
                    rec.description = ""
                    SeqIO.write(rec, f, "fasta")
        else:
            for rec in SeqIO.parse(filtered_fa, "fasta"):
                rec.id = safe_id(rec.id)
                rec.description = ""
                SeqIO.write(rec, f, "fasta")

    # align + tree
    subprocess.run(["mafft", "--thread", str(threads), "--auto", combined], stdout=open(aln, "w"), check=True)
    subprocess.run(["fasttree", "-nt", "-gtr", aln], stdout=open(nwk, "w"), check=True)

    # simple coloring: refs in gray; samples in blue
    t = Tree(nwk)
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.title.add_face(TextFace("De-novo pspA (NT) tree — refs gray, samples blue", fsize=12, bold=True), column=0)
    for node in t.traverse():
        if node.is_leaf():
            col = "#888888" if node.name.startswith("Ref\n") else "#377eb8"
            ns = NodeStyle()
            ns["fgcolor"] = col
            ns["size"] = 8
            node.set_style(ns)
    t.render(png, w=1600, dpi=200, tree_style=ts)

def main():
    ap = argparse.ArgumentParser(description="De-novo pspA detection (BLASTX + 200-bp upstream validation)")
    ap.add_argument("--assemblies-dir", required=True, help="Folder with *.fasta assemblies (basenames are assembly IDs)")
    ap.add_argument("--protein-db", required=False, default=None, help="Protein BLAST DB basename for Hollingshead pspA")
    ap.add_argument("--protein-subject-fasta", required=False, default=None, help="Fallback protein FASTA for BLASTX")
    ap.add_argument("--upstream-db", required=False, default=None, help="Nucleotide BLAST DB basename (200-bp Rx1 upstream)")
    ap.add_argument("--upstream-subject-fasta", required=False, default=None, help="Fallback upstream FASTA for BLASTN")
    ap.add_argument("--work-dir", required=True, help="Base working directory")
    ap.add_argument("--threads", type=int, default=12, help="Global threads (capped at 12)")
    ap.add_argument("--jobs", type=int, default=0, help="Parallel samples (0=auto)")
    ap.add_argument("--blastx-threads", type=int, default=0, help="Threads per BLASTX (0=auto)")
    ap.add_argument("--min-pident", type=float, default=75.0)
    ap.add_argument("--min-aa-len", type=int, default=200)
    ap.add_argument("--evalue", default="1e-6")
    ap.add_argument("--seg", default="yes", choices=["yes", "no"])
    ap.add_argument("--max-target-seqs", type=int, default=10)
    ap.add_argument("--max-hits-per-sample", type=int, default=3)
    ap.add_argument("--min-inter-hit-distance", type=int, default=2000)
    ap.add_argument("--max-hits-hard-cap", type=int, default=6)
    ap.add_argument("--denovo-nt-tree", action="store_true", help="Build optional NT tree for de-novo hits (+ refs)")
    ap.add_argument("--denovo-nt-ref", default=None, help="Reference NT FASTA for de-novo NT tree")
    # NEW: upstream gating thresholds and HighConfidence filter
    ap.add_argument("--upstream-min-pident", type=float, default=75.0, help="Upstream BLASTN gating: minimum %identity")
    ap.add_argument("--upstream-min-length", type=int, default=100, help="Upstream BLASTN gating: minimum alignment length (nt)")
    ap.add_argument("--upstream-only-highconf", action="store_true",
                    help="If set, only retain hits whose upstream BLASTN meets HighConfidence (>=90% and >=150 nt)")
    args = ap.parse_args()

    detect_denovo_pspA(
        assemblies_dir=args.assemblies_dir,
        protein_db=args.protein_db,
        protein_subject_fasta=args.protein_subject_fasta,
        upstream_db=args.upstream_db,
        upstream_subject_fasta=args.upstream_subject_fasta,
        work_dir=args.work_dir,
        threads=args.threads,
        jobs=args.jobs,
        blastx_threads=args.blastx_threads,
        min_pident=args.min_pident,
        min_aa_len=args.min_aa_len,
        evalue=args.evalue,
        seg=args.seg,
        max_target_seqs=args.max_target_seqs,
        max_hits_per_sample=args.max_hits_per_sample,
        min_inter_hit_distance=args.min_inter_hit_distance,
        max_hits_hard_cap=args.max_hits_hard_cap,
        build_denovo_nt_tree=args.denovo_nt_tree,
        denovo_nt_ref_fasta=args.denovo_nt_ref,
        upstream_min_pident=args.upstream_min_pident,
        upstream_min_length=args.upstream_min_length,
        only_high_confidence=bool(args.upstream_only_highconf)
    )

if __name__ == "__main__":
        main()