#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified orchestrator for pspA de-novo detection, rescue, BoxB classification, and optional tree.
(Full description omitted for brevity; behavior as documented previously.)
"""
import os
import argparse
import tempfile
import pandas as pd
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from pspA_utils import (
    dump_db_or_subject_fasta,
    parse_cluster_from_protein_id,
    safe_id,
    normalize_sample_id,
)
from core_boxb_classifier import classify_boxb
from denovo_pspA_detect import detect_denovo_pspA
from rescue_refguided import rescue_for_strains

SUMMARY_HEADER = [
    "Sample", "CodingID", "Status",
    "Matching_pspA_reference", "Hollingshead_family",
    "pident", "bitscore", "AA_Length",
    "AA_BoxB", "evalue", "qstart_nt", "qend_nt", "qframe", "Source"
]
COMPACT_HEADER = [
    "Assembly", "Status", "pspA_copy_number",
    "Matching_pspA_reference", "Hollingshead_family",
    "pident", "bitscore", "AA_Length", "All_Matching_pspA_references"
]

def enumerate_assemblies(assemblies_dir: str) -> List[str]:
    return sorted(os.path.splitext(f)[0] for f in os.listdir(assemblies_dir) if f.endswith(".fasta"))

def fasta_record_count(path: str) -> int:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return 0
    try:
        return sum(1 for _ in SeqIO.parse(path, "fasta"))
    except Exception:
        return 0

def load_filtered_coding(nt_fasta: str) -> Tuple[Dict[str, List[SeqRecord]], Dict[str, int]]:
    if not os.path.exists(nt_fasta):
        return {}, {}
    by_sample: Dict[str, List[SeqRecord]] = {}
    for rec in SeqIO.parse(nt_fasta, "fasta"):
        sid = normalize_sample_id(rec.id)
        by_sample.setdefault(sid, []).append(rec)
    copy_number = {sid: len(lst) for sid, lst in by_sample.items()}
    return by_sample, copy_number

def is_nonempty_file(p: str) -> bool:
    return os.path.exists(p) and os.path.getsize(p) > 0

def extract_hollingshead_family(ref_id: str) -> str:
    import re
    match = re.search(r'F\d+_C\d+', ref_id)
    return match.group(0) if match else "UNK"

def build_tree(
    work_dir: str,
    boxb_db: str,
    boxb_subject_fasta: str,
    aa_records: List[Tuple[str, str]],
    cluster_map: Dict[str, str],
    threads: int = 12
) -> None:
    import subprocess
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as _SeqRecord
    from ete3 import Tree, TreeStyle, TextFace, NodeStyle, faces
    import csv

    threads = max(1, min(12, int(threads)))
    d_tree = os.path.join(work_dir, "07_tree")
    os.makedirs(d_tree, exist_ok=True)

    refs_fa = os.path.join(d_tree, "boxB_refs.faa")
    dump_db_or_subject_fasta(boxb_db, boxb_subject_fasta, refs_fa)

    ref_records = []
    ref_clusters: Dict[str, str] = {}
    for rec in SeqIO.parse(refs_fa, "fasta"):
        rid = rec.id.split()[0]
        rec.id = f"Ref_{safe_id(rid)}"
        rec.description = ""
        ref_records.append(rec)
        ref_clusters[rec.id] = extract_hollingshead_family(rid)

    combined_faa = os.path.join(d_tree, "boxB_combined.faa")
    with open(combined_faa, "w") as f:
        for rec in ref_records:
            SeqIO.write(rec, f, "fasta")
        for leaf_id, aa in aa_records:
            r = _SeqRecord(Seq(aa), id=safe_id(leaf_id), description="")
            SeqIO.write(r, f, "fasta")

    aln_faa = os.path.join(d_tree, "boxB_aligned.faa")
    tree_nwk = os.path.join(d_tree, "boxB_tree.nwk")
    tree_png = os.path.join(d_tree, "boxB_tree_colored.png")
    legend_csv = os.path.join(d_tree, "boxB_cluster_legend.csv")

    subprocess.run(["mafft", "--thread", str(threads), "--auto", combined_faa], stdout=open(aln_faa, "w"), check=True)
    subprocess.run(["fasttree", "-gamma", "-wag", aln_faa], stdout=open(tree_nwk, "w"), check=True)

    tree = Tree(tree_nwk)
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.title.add_face(TextFace("PspA BoxB tree (protein) — colored by Hollingshead family", fsize=14, bold=True), column=0)

    palette = [
        "#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
        "#ff7f00", "#a65628", "#f781bf", "#999999",
        "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"
    ]
    cluster_colors: Dict[str, str] = {}
    ci = 0

    def get_cluster(leaf: str) -> str:
        if leaf.startswith("Ref_"):
            return ref_clusters.get(leaf, "UNK")
        return cluster_map.get(leaf, "UNK")

    for node in tree.traverse():
        if node.is_leaf():
            leaf = node.name.split()[0]
            cl = get_cluster(leaf)
            if cl not in cluster_colors:
                cluster_colors[cl] = palette[ci % len(palette)]
                ci += 1
            ns = NodeStyle()
            ns["fgcolor"] = cluster_colors[cl]
            ns["size"] = 10
            node.set_style(ns)

    legend_rows = sorted(cluster_colors.items())
    try:
        ts.legend_position = 1
        for cl, col in legend_rows:
            ts.legend.add_face(faces.CircleFace(6, col), column=0)
            ts.legend.add_face(TextFace(f" {cl}", fsize=10), column=1)
            ts.legend.add_face(TextFace(" "), column=0)
    except Exception:
        pass

    with open(legend_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Cluster", "Color"])
        for cl, col in legend_rows:
            w.writerow([cl, col])

    tree.render(tree_png, w=2000, dpi=300, tree_style=ts)

def main():
    ap = argparse.ArgumentParser(description="Orchestrate de-novo + rescue pspA BoxB analysis")
    ap.add_argument("--assemblies-dir", required=True)
    ap.add_argument("--reads-dir", required=True)
    ap.add_argument("--rescue-ref-fasta", required=True)
    ap.add_argument("--boxb-db", required=True)
    ap.add_argument("--boxb-subject-fasta", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--build-tree", action="store_true")
    ap.add_argument("--max-per-sample", type=int, default=1)

    # De-novo
    ap.add_argument("--protein-db", default=None)
    ap.add_argument("--protein-subject-fasta", default=None)
    ap.add_argument("--upstream-db", default=None)
    ap.add_argument("--upstream-subject-fasta", default=None)
    ap.add_argument("--upstream-min-pident", type=float, default=75.0)
    ap.add_argument("--upstream-min-length", type=int, default=100)
    ap.add_argument("--denovo-only-highconf", action="store_true")
    ap.add_argument("--filtered-coding-fasta", default=None)
    ap.add_argument("--jobs", type=int, default=0)
    ap.add_argument("--blastx-threads", type=int, default=0)
    ap.add_argument("--min-pident", type=float, default=75.0)
    ap.add_argument("--min-aa-len", type=int, default=200)
    ap.add_argument("--evalue", default="1e-6")
    ap.add_argument("--seg", default="yes", choices=["yes", "no"])
    ap.add_argument("--max-target-seqs", type=int, default=10)
    ap.add_argument("--max-hits-per-sample", type=int, default=3)
    ap.add_argument("--min-inter-hit-distance", type=int, default=2000)
    ap.add_argument("--max-hits-hard-cap", type=int, default=6)
    ap.add_argument("--denovo-nt-tree", action="store_true")
    ap.add_argument("--denovo-nt-ref", default=None)

    # BoxB thresholds
    ap.add_argument("--boxb-min-aa-len", type=int, default=60)
    ap.add_argument("--boxb-min-pident", type=float, default=80.0)

    # NEW: rescue multi-candidate reporting
    ap.add_argument("--rescue-report-all-qualified", action="store_true",
                    help="Classify all rescued references (breadth-qualified) per strain")

    args = ap.parse_args()
    threads = max(1, min(12, int(args.threads)))

    os.makedirs(args.work_dir, exist_ok=True)
    d01 = os.path.join(args.work_dir, "01_denovo")
    d05 = os.path.join(args.work_dir, "05_boxb_per_sample")
    d06 = os.path.join(args.work_dir, "06_boxb_reports")
    os.makedirs(d05, exist_ok=True)
    os.makedirs(d06, exist_ok=True)

    assemblies = enumerate_assemblies(args.assemblies_dir)

    # De-novo run (if needed)
    if args.filtered_coding_fasta:
        filtered_coding_fa = args.filtered_coding_fasta
    else:
        filtered_coding_fa = os.path.join(d01, "filtered_pspA_hits.fasta")
    need_denovo = fasta_record_count(filtered_coding_fa) == 0

    if need_denovo:
        os.makedirs(d01, exist_ok=True)
        detect_denovo_pspA(
            upstream_min_pident=args.upstream_min_pident,
            upstream_min_length=args.upstream_min_length,
            assemblies_dir=args.assemblies_dir,
            protein_db=args.protein_db,
            protein_subject_fasta=args.protein_subject_fasta,
            upstream_db=args.upstream_db,
            upstream_subject_fasta=args.upstream_subject_fasta,
            work_dir=args.work_dir,
            threads=threads,
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
            only_high_confidence=bool(args.denovo_only_highconf)
        )

    if fasta_record_count(filtered_coding_fa) == 0:
        print("[WARN] No de-novo pspA coding sequences found after detection.")

    by_sample, copy_number = load_filtered_coding(filtered_coding_fa)
    total_records = sum(len(v) for v in by_sample.values())
    print(f"[INFO] De-novo file: {filtered_coding_fa}")
    print(f"[INFO] Assemblies in dir: {len(assemblies)}")
    print(f"[INFO] De-novo unique samples found: {len(by_sample)} (total coding records: {total_records})")

    missing = sorted([a for a in assemblies if a not in by_sample])
    present = sorted(set(assemblies) - set(missing))

    # Rescue
    if missing:
        rescue_results = rescue_for_strains(
            strains=missing,
            reads_dir=args.reads_dir,
            ref_fasta_with_200bp=args.rescue_ref_fasta,
            work_dir=args.work_dir,
            threads=threads,
            breadth_threshold=0.80,
            min_depth=10,
            report_all_qualified=bool(args.rescue_report_all_qualified)
        )
    else:
        rescue_results = {}

    # Classification
    summary_rows: List[Dict] = []
    sample_best: Dict[str, Dict] = {}
    all_refs_per_sample: Dict[str, List[str]] = {}

    # De-novo records — classify each coding record
    for sid in present:
        for rec in by_sample[sid]:
            with tempfile.NamedTemporaryFile("w", delete=False, suffix=".fa") as tmpq:
                SeqIO.write(rec, tmpq.name, "fasta")
            out_rec = classify_boxb(
                nt_fasta=tmpq.name,
                sample_id=sid,
                boxb_db=args.boxb_db,
                boxb_subject_fasta=args.boxb_subject_fasta,
                out_dir=d05,
                min_aa_len=args.boxb_min_aa_len,
                min_pident=args.boxb_min_pident,
                sort_field="bitscore",
                threads=threads
            )
            out_rec["Source"] = "DeNovo"
            summary_rows.append(out_rec)
            if out_rec.get("Status") == "Matched":
                if sid not in sample_best or float(out_rec["bitscore"]) > float(sample_best[sid]["bitscore"]):
                    sample_best[sid] = out_rec
                all_refs_per_sample.setdefault(sid, []).append(str(out_rec["Matching_pspA_reference"]))

    # Rescued — classify single or multiple qualified per strain
    for sid in missing:
        r = rescue_results.get(sid, {})
        qualified_list = r.get("qualified_list", []) if isinstance(r, dict) else []
        if qualified_list:
            # Multi-candidate path
            for q in qualified_list:
                rescued_fa = q.get("trimmed_fasta")
                if not rescued_fa or not os.path.exists(rescued_fa):
                    continue
                out_rec = classify_boxb(
                    nt_fasta=rescued_fa,
                    sample_id=sid,
                    boxb_db=args.boxb_db,
                    boxb_subject_fasta=args.boxb_subject_fasta,
                    out_dir=d05,
                    min_aa_len=args.boxb_min_aa_len,
                    min_pident=args.boxb_min_pident,
                    sort_field="bitscore",
                    threads=threads
                )
                out_rec["Source"] = q.get("status", "Rescue_Qualified")
                summary_rows.append(out_rec)
                if out_rec.get("Status") == "Matched":
                    if sid not in sample_best or float(out_rec["bitscore"]) > float(sample_best[sid]["bitscore"]):
                        sample_best[sid] = out_rec
                    all_refs_per_sample.setdefault(sid, []).append(str(out_rec["Matching_pspA_reference"]))
        else:
            # Backward-compatible single-best path
            rescued_fa = r.get("out_fasta")
            if not rescued_fa or not os.path.exists(rescued_fa):
                continue
            out_rec = classify_boxb(
                nt_fasta=rescued_fa,
                sample_id=sid,
                boxb_db=args.boxb_db,
                boxb_subject_fasta=args.boxb_subject_fasta,
                out_dir=d05,
                min_aa_len=args.boxb_min_aa_len,
                min_pident=args.boxb_min_pident,
                sort_field="bitscore",
                threads=threads
            )
            out_rec["Source"] = r.get("status", "Rescue")
            summary_rows.append(out_rec)
            if out_rec.get("Status") == "Matched":
                if sid not in sample_best or float(out_rec["bitscore"]) > float(sample_best[sid]["bitscore"]):
                    sample_best[sid] = out_rec
                all_refs_per_sample.setdefault(sid, []).append(str(out_rec["Matching_pspA_reference"]))

    # Reports
    d06 = os.path.join(args.work_dir, "06_boxb_reports")
    os.makedirs(d06, exist_ok=True)
    summary_csv = os.path.join(d06, "boxB_summary.csv")
    with open(summary_csv, "w") as f:
        f.write(",".join(SUMMARY_HEADER) + "\n")
        for r in summary_rows:
            row = [str(r.get(k, "")) for k in SUMMARY_HEADER]
            f.write(",".join(row) + "\n")

    compact_csv = os.path.join(d06, "pspA_boxB_compact.csv")
    with open(compact_csv, "w") as f:
        f.write(",".join(COMPACT_HEADER) + "\n")
        for asm in assemblies:
            copies = int(copy_number.get(asm, 0))
            g = [r for r in summary_rows if r.get("Sample") == asm and r.get("Status") == "Matched"]
            if copies > 0:
                if g:
                    best = max(g, key=lambda x: float(x["bitscore"]))
                    all_refs = ";".join(sorted(set(all_refs_per_sample.get(asm, []))))
                    status = "Matched"
                    f.write(",".join([
                        asm, status, str(copies),
                        str(best["Matching_pspA_reference"]),
                        str(best["Hollingshead_family"]),
                        str(best["pident"]), str(best["bitscore"]),
                        str(best["AA_Length"]), all_refs
                    ]) + "\n")
                else:
                    f.write(",".join([asm, "No_BoxB_match", str(copies), "", "", "", "", "", ""]) + "\n")
            else:
                # No de-novo pspA; consider rescue results (still best-per-assembly)
                rescued_rows = [x for x in summary_rows if x.get("Sample") == asm and str(x.get("Source", "")).startswith("Rescue")]
                if rescued_rows and any(x.get("Status") == "Matched" for x in rescued_rows):
                    best = max([x for x in rescued_rows if x.get("Status") == "Matched"], key=lambda y: float(y["bitscore"]))
                    all_refs = ";".join(sorted(set(all_refs_per_sample.get(asm, []))))
                    f.write(",".join([
                        asm, "Rescue_Matched", "1",
                        str(best["Matching_pspA_reference"]),
                        str(best["Hollingshead_family"]),
                        str(best["pident"]), str(best["bitscore"]),
                        str(best["AA_Length"]), all_refs
                    ]) + "\n")
                else:
                    # Reflect rescue status if available, else "No_pspA_detected"
                    r = rescue_results.get(asm, {})
                    status = r.get("status", "No_pspA_detected") if r else "No_pspA_detected"
                    f.write(",".join([asm, status, "0", "", "", "", "", "", ""]) + "\n")

    # Optional tree (best per sample)
    if args.build_tree:
        aa_records: List[Tuple[str, str]] = []
        cluster_map: Dict[str, str] = {}
        for sid, rec in sample_best.items():
            aa = rec.get("AA_BoxB", "")
            if not aa:
                continue
            aa_records.append((sid, aa))
            cluster_map[sid] = rec.get("Hollingshead_family", "UNK")
        if aa_records:
            build_tree(
                work_dir=args.work_dir,
                boxb_db=args.boxb_db,
                boxb_subject_fasta=args.boxb_subject_fasta,
                aa_records=aa_records,
                cluster_map=cluster_map,
                threads=threads
            )

    n_total = len(assemblies)
    df_compact = pd.read_csv(compact_csv)
    n_matched = (df_compact["Status"].isin(["Matched", "Rescue_Matched"])).sum()
    n_no_boxb = (df_compact["Status"].isin(["No_BoxB_match", "Rescue_No_BoxB_match"])).sum()
    n_no_pspa = (df_compact["Status"] == "No_pspA_detected").sum()
    print("Done.")
    print(f"Per-record summary: {summary_csv}")
    print(f"Compact report: {compact_csv}")
    print(f"Assemblies: {n_total} \n Matched: {n_matched} \n No_BoxB_match: {n_no_boxb} \n No_pspA_detected: {n_no_pspa}")

if __name__ == "__main__":
    main()