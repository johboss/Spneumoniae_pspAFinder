#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reference-guided pspA rescue using minimap2 + samtools.

For each strain (missing de-novo pspA), align cleaned paired-end reads against a
combined reference FASTA (pspA with ~200-bp upstream). Select rescued candidates
by breadth (>= breadth_threshold at >= min_depth used for breadth) and mean depth.

Emit trimmed (200-nt upstream removed) rescued FASTA(s):
 - single best when report_all_qualified=False (default, backward-compatible)
 - one per qualified reference when report_all_qualified=True

Outputs under {--work-dir}/04_best/ :
  - rescue_refs.mmi                 (minimap2 index of provided reference FASTA)
  - align/{strain}.sorted.bam (+ .bai)          (sorted & indexed BAM)
  - depth/{strain}__{ref}.depth.txt             (samtools depth per reference)
  - fasta/{strain}[__{ref}].rescued.trimmed.fa  (trimmed coding-only)
  - rescue_summary.csv              (rows per strain; if reporting-all, one row per qualified)

Function (for orchestrator):
  rescue_for_strains(
      strains: List[str],
      reads_dir: str,
      ref_fasta_with_200bp: str,
      work_dir: str,
      threads: int = 12,
      breadth_threshold: float = 0.80,
      min_depth: int = 10,
      report_all_qualified: bool = False
  ) -> Dict[str, Dict[str, Any]]

Return dict per strain:
{
  "status": "Rescue_Qualified" | "Rescue_LowBreadth" | "Rescue_NoReads" | "Rescue_Failed",
  "best_reference": "<ref_id>" | "",
  "breadth": <float>,              # metrics for best reference
  "mean_depth": <float>,           # metrics for best reference
  "out_fasta": "<path>" | None,    # single best (for backward compatibility)
  "bam": "<path-to-sorted.bam>" | None,
  "qualified_list": [              # present only when report_all_qualified=True and any qualify
     {
       "reference": "<ref_id>",
       "breadth": <float>,
       "mean_depth": <float>,
       "trimmed_fasta": "<path>",
       "status": "Rescue_Qualified"
     }, ...
  ]
}
"""
import os
import argparse
from typing import Dict, List, Tuple, Optional, Any
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pspA_utils import (
    run_cmd,
    compute_breadth_from_samtools_depth,
    trim_fixed_upstream_200,  # remove the first 200 nt
    safe_id,
    normalize_sample_id
)

# -------------------- helpers --------------------
def _cap_threads(t: int) -> int:
    return max(1, min(12, int(t)))

def _paths(work_dir: str) -> Dict[str, str]:
    d04 = os.path.join(work_dir, "04_best")
    d_align = os.path.join(d04, "align")
    d_depth = os.path.join(d04, "depth")
    d_fasta = os.path.join(d04, "fasta")
    d_logs = os.path.join(d04, "logs")
    os.makedirs(d04, exist_ok=True)
    os.makedirs(d_align, exist_ok=True)
    os.makedirs(d_depth, exist_ok=True)
    os.makedirs(d_fasta, exist_ok=True)
    os.makedirs(d_logs, exist_ok=True)
    return {
        "d04": d04,
        "align": d_align,
        "depth": d_depth,
        "fasta": d_fasta,
        "logs": d_logs,
        "mmi": os.path.join(d04, "rescue_refs.mmi"),
        "summary_csv": os.path.join(d04, "rescue_summary.csv"),
    }

def _reads_for_strain(reads_dir: str, strain: str) -> Tuple[str, str]:
    r1 = os.path.join(reads_dir, strain, "R1_001.fastq.gz")
    r2 = os.path.join(reads_dir, strain, "R2_001.fastq.gz")
    return r1, r2

def _file_nonempty(p: str) -> bool:
    return os.path.exists(p) and os.path.getsize(p) > 0

def _write_summary_header_if_needed(csv_path: str) -> None:
    hdr_needed = not os.path.exists(csv_path)
    if hdr_needed:
        with open(csv_path, "w") as f:
            f.write(",".join([
                "Strain", "Reference", "Breadth>=10x", "MeanDepth",
                "MinDepthThreshold", "Qualified(>=0.80)", "Status",
                "BAM", "TrimmedFASTA"
            ]) + "\n")

def _append_summary_row(csv_path: str, row: Dict[str, Any]) -> None:
    _write_summary_header_if_needed(csv_path)
    with open(csv_path, "a") as f:
        f.write(",".join([
            str(row.get("Strain", "")),
            str(row.get("Reference", "")),
            f"{row.get('Breadth', 0.0):.6f}",
            f"{row.get('MeanDepth', 0.0):.3f}",
            str(row.get("MinDepthThreshold", 10)),
            "yes" if bool(row.get("Qualified", False)) else "no",
            str(row.get("Status", "")),
            str(row.get("BAM", "")),
            str(row.get("TrimmedFASTA", "")),
        ]) + "\n")

# -------------------- core API --------------------
def rescue_for_strains(
    strains: List[str],
    reads_dir: str,
    ref_fasta_with_200bp: str,
    work_dir: str,
    threads: int = 12,
    breadth_threshold: float = 0.80,
    min_depth: int = 10,
    report_all_qualified: bool = False
) -> Dict[str, Dict[str, Any]]:
    """
    Run reference-guided rescue for provided strains.
    Returns a mapping {strain: {...}} as documented in module docstring.
    """
    threads = _cap_threads(threads)
    ps = _paths(work_dir)

    # Build minimap2 index (refresh if missing or older than the fasta)
    if (not _file_nonempty(ps["mmi"])) or (os.path.getmtime(ps["mmi"]) < os.path.getmtime(ref_fasta_with_200bp)):
        run_cmd(["minimap2", "-d", ps["mmi"], ref_fasta_with_200bp], check=True)

    # Load references into dict for trimming later
    ref_records: Dict[str, SeqRecord] = {rec.id.split()[0]: rec for rec in SeqIO.parse(ref_fasta_with_200bp, "fasta")}
    ref_ids = list(ref_records.keys())

    results: Dict[str, Dict[str, Any]] = {}
    _write_summary_header_if_needed(ps["summary_csv"])

    for strain in strains:
        sid = normalize_sample_id(strain)
        r1, r2 = _reads_for_strain(reads_dir, sid)

        # Validate reads
        if not (_file_nonempty(r1) and _file_nonempty(r2)):
            status = "Rescue_NoReads"
            results[sid] = {
                "status": status,
                "best_reference": "",
                "breadth": 0.0,
                "mean_depth": 0.0,
                "out_fasta": None,
                "bam": None
            }
            _append_summary_row(ps["summary_csv"], {
                "Strain": sid, "Reference": "", "Breadth": 0.0, "MeanDepth": 0.0,
                "MinDepthThreshold": min_depth, "Qualified": False, "Status": status,
                "BAM": "", "TrimmedFASTA": ""
            })
            continue

        # Paths
        sam_path = os.path.join(ps["align"], f"{safe_id(sid)}.sam")
        bam_raw  = os.path.join(ps["align"], f"{safe_id(sid)}.bam")
        bam_sorted = os.path.join(ps["align"], f"{safe_id(sid)}.sorted.bam")

        # minimap2 -> SAM
        run_cmd([
            "minimap2", "-t", str(threads), "-ax", "sr", ps["mmi"], r1, r2, "-o", sam_path
        ], check=True)

        # samtools view -> BAM
        run_cmd(["samtools", "view", "-b", "-o", bam_raw, sam_path], check=True)

        # samtools sort -> sorted BAM
        sort_threads = max(1, min(threads, 4))
        run_cmd(["samtools", "sort", "-@", str(sort_threads), "-o", bam_sorted, bam_raw], check=True)
        run_cmd(["samtools", "index", bam_sorted], check=True)

        # Compute breadth and mean depth per reference, track best and qualified set
        best_ref = ""
        best_breadth = -1.0
        best_mean_depth = -1.0
        per_ref_metrics: Dict[str, Tuple[float, float]] = {}

        for rid in ref_ids:
            depth_txt = os.path.join(ps["depth"], f"{safe_id(sid)}__{safe_id(rid)}.depth.txt")
            run_cmd(["samtools", "depth", "-a", "-r", rid, bam_sorted], stdout_path=depth_txt, check=True)
            breadth, mean_depth = compute_breadth_from_samtools_depth(depth_txt, min_depth=min_depth)
            per_ref_metrics[rid] = (breadth, mean_depth)
            # Select best by highest breadth, break ties by mean depth
            is_better = (breadth > best_breadth) or (abs(breadth - best_breadth) < 1e-9 and mean_depth > best_mean_depth)
            if is_better:
                best_ref = rid
                best_breadth = breadth
                best_mean_depth = mean_depth

        # Decide qualification for "best"
        qualified_best = (best_breadth >= float(breadth_threshold))
        status = "Rescue_Qualified" if qualified_best else "Rescue_LowBreadth"

        # Backward-compatible single best trimmed FASTA (only if qualified)
        trimmed_best = None
        candidate_best = os.path.join(ps["fasta"], f"{safe_id(sid)}.rescued.trimmed.fa")
        if qualified_best and best_ref and best_ref in ref_records:
            rec = ref_records[best_ref]
            trimmed = trim_fixed_upstream_200(rec)  # remove first 200 nt upstream
            trimmed.id = sid
            trimmed.description = f"rescued_from={best_ref}"
            with open(candidate_best, "w") as f:
                SeqIO.write(trimmed, f, "fasta")
            trimmed_best = candidate_best
        else:
            try:
                if os.path.exists(candidate_best):
                    os.remove(candidate_best)
            except Exception:
                pass

        # If reporting all qualified, emit one trimmed FASTA per qualified reference
        qualified_list: List[Dict[str, Any]] = []
        if report_all_qualified:
            for rid, (br, md) in per_ref_metrics.items():
                if br >= float(breadth_threshold):
                    # Per-reference trimmed fasta
                    candidate = os.path.join(ps["fasta"], f"{safe_id(sid)}__{safe_id(rid)}.rescued.trimmed.fa")
                    rec = ref_records[rid]
                    trimmed = trim_fixed_upstream_200(rec)
                    trimmed.id = sid
                    trimmed.description = f"rescued_from={rid}"
                    with open(candidate, "w") as f:
                        SeqIO.write(trimmed, f, "fasta")
                    qualified_list.append({
                        "reference": rid,
                        "breadth": float(br),
                        "mean_depth": float(md),
                        "trimmed_fasta": candidate,
                        "status": "Rescue_Qualified"
                    })
                    # Per-qualified row
                    _append_summary_row(ps["summary_csv"], {
                        "Strain": sid,
                        "Reference": rid,
                        "Breadth": br,
                        "MeanDepth": md,
                        "MinDepthThreshold": min_depth,
                        "Qualified": True,
                        "Status": "Rescue_Qualified",
                        "BAM": bam_sorted,
                        "TrimmedFASTA": candidate
                    })

            # If none qualified, still append one summary row with best metrics
            if len(qualified_list) == 0:
                _append_summary_row(ps["summary_csv"], {
                    "Strain": sid,
                    "Reference": best_ref,
                    "Breadth": best_breadth if best_breadth >= 0 else 0.0,
                    "MeanDepth": best_mean_depth if best_mean_depth >= 0 else 0.0,
                    "MinDepthThreshold": min_depth,
                    "Qualified": False,
                    "Status": status,
                    "BAM": bam_sorted,
                    "TrimmedFASTA": ""
                })
        else:
            # Single-row summary (best only)
            _append_summary_row(ps["summary_csv"], {
                "Strain": sid,
                "Reference": best_ref,
                "Breadth": best_breadth if best_breadth >= 0 else 0.0,
                "MeanDepth": best_mean_depth if best_mean_depth >= 0 else 0.0,
                "MinDepthThreshold": min_depth,
                "Qualified": qualified_best,
                "Status": status,
                "BAM": bam_sorted,
                "TrimmedFASTA": trimmed_best if trimmed_best else ""
            })

        # Assemble result record
        res = {
            "status": status,
            "best_reference": best_ref,
            "breadth": float(best_breadth if best_breadth >= 0 else 0.0),
            "mean_depth": float(best_mean_depth if best_mean_depth >= 0 else 0.0),
            "out_fasta": trimmed_best,   # backward-compatible single file (may be None)
            "bam": bam_sorted
        }
        if report_all_qualified and qualified_list:
            res["qualified_list"] = qualified_list

        results[sid] = res

        # Clean up large intermediates (keep sorted BAM; remove SAM and raw BAM)
        try:
            if os.path.exists(sam_path):
                os.remove(sam_path)
            if os.path.exists(bam_raw):
                os.remove(bam_raw)
        except Exception:
            pass

    return results

# -------------------- CLI --------------------
def main():
    ap = argparse.ArgumentParser(description="Reference-guided pspA rescue (minimap2 + samtools)")
    ap.add_argument("--strains", required=True, nargs="+",
                    help="List of strain IDs to rescue (space-separated)")
    ap.add_argument("--reads-dir", required=True,
                    help="Folder with per-strain subdirs containing R1_001.fastq.gz and R2_001.fastq.gz")
    ap.add_argument("--ref-fasta-with-200bp", required=True,
                    help="Combined references with upstream + pspA (trim 200 nt upstream for output)")
    ap.add_argument("--work-dir", required=True, help="Base working directory")
    ap.add_argument("--threads", type=int, default=12, help="Threads (capped at 12)")
    ap.add_argument("--breadth-threshold", type=float, default=0.80,
                    help="Breadth cutoff for 'Rescue_Qualified'")
    ap.add_argument("--min-depth", type=int, default=10,
                    help="Depth threshold used when computing breadth (>= this depth)")
    ap.add_argument("--report-all-qualified", action="store_true",
                    help="Emit and report all rescued references that meet breadth threshold (not just the best)")
    args = ap.parse_args()

    _ = rescue_for_strains(
        strains=args.strains,
        reads_dir=args.reads_dir,
        ref_fasta_with_200bp=args.ref_fasta_with_200bp,
        work_dir=args.work_dir,
        threads=args.threads,
        breadth_threshold=args.breadth_threshold,
        min_depth=args.min_depth,
        report_all_qualified=bool(args.report_all_qualified)
    )

if __name__ == "__main__":
    main()