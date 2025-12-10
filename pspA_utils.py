#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared helpers for the pspA workflow.

- No external I/O is performed automatically; functions are reusable/importable.
- Centralizes ID normalization to avoid circular imports between modules.
"""

import os
import re
import json
import subprocess
from typing import Dict, Any, List, Tuple, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# ---------------- shell utilities ----------------
def run_cmd(cmd: List[str], stdout_path: Optional[str] = None, check: bool = True) -> subprocess.CompletedProcess:
    """
    Run a command with optional redirection; raise if fails when check=True.
    Returns subprocess.CompletedProcess.
    """
    if stdout_path:
        with open(stdout_path, "w") as out:
            return subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, check=check)
    else:
        return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=check)


# ---------------- BLAST helpers ----------------
def blast_db_exists(db_basename: str) -> bool:
    """Return True if protein BLAST DB with given basename exists (.pin/.phr/.psq)."""
    return all(os.path.exists(db_basename + ext) for ext in (".pin", ".phr", ".psq"))


def blast_nucl_db_exists(db_basename: str) -> bool:
    """Return True if nucleotide BLAST DB with given basename exists (.nin/.nhr/.nsq)."""
    return all(os.path.exists(db_basename + ext) for ext in (".nin", ".nhr", ".nsq"))


def dump_db_or_subject_fasta(db_basename: str, subject_fasta: str, out_fasta: str) -> None:
    """
    Dump all sequences from DB if available; else sanitize and copy subject FASTA.
    """
    if blast_db_exists(db_basename) or blast_nucl_db_exists(db_basename):
        run_cmd(["blastdbcmd", "-db", db_basename, "-entry", "all", "-outfmt", "%f"], stdout_path=out_fasta, check=True)
    else:
        with open(out_fasta, "w") as fout:
            for rec in SeqIO.parse(subject_fasta, "fasta"):
                rec.id = rec.id.split()[0]
                rec.description = ""
                SeqIO.write(rec, fout, "fasta")


# ---------------- IDs and clusters ----------------
def safe_id(x: str) -> str:
    """Make an identifier shell/FASTA-safe enough for our use."""
    return x.replace(" ", "_")


def normalize_sample_id(raw_id: str, split_marker: str = "__hit") -> str:
    """
    Convert a FASTA header into a sample ID that should match assembly basenames.
    Rules:
    - If 'split_marker' is present, take the part before it; else take the first token up to whitespace.
    - Remove common file-like suffixes.
    - Replace spaces with underscores.
    """
    hdr = raw_id.strip()
    if split_marker in hdr:
        sid = hdr.split(split_marker)[0]
    else:
        sid = hdr.split()[0]

    # FIX: remove common FASTA-like extensions using a single-line regex
    sid = re.sub(r'\.(?:fa|fasta|fna|fas)$', '', sid, flags=re.IGNORECASE)

    sid = sid.replace(' ', '_')
    return sid




def parse_cluster_from_protein_id(protein_id: str) -> str:
    """
    Keep the current rule exactly as requested: characters [2:7] if length >= 7, else 'UNK'.
    """
    if isinstance(protein_id, str) and len(protein_id) >= 7:
        return protein_id[2:7]
    return "UNK"


# ---------------- JSON/CSV one-row writers ----------------
def write_one_row_json(record: Dict[str, Any], out_json_path: str) -> None:
    with open(out_json_path, "w") as f:
        json.dump(record, f, indent=2, ensure_ascii=False)


def write_one_row_csv(record: Dict[str, Any], out_csv_path: str, header_order: List[str]) -> None:
    """Create or append with header if missing; write values in header_order."""
    header_needed = not os.path.exists(out_csv_path)
    with open(out_csv_path, "a") as f:
        if header_needed:
            f.write(",".join(header_order) + "\n")
        row = [str(record.get(k, "")) for k in header_order]
        f.write(",".join(row) + "\n")


# ---------------- FASTA sequence trimming ----------------
def trim_fixed_upstream_200(nt_record: SeqRecord) -> SeqRecord:
    """
    Return a copy of the record with the first 200 nt removed (remove upstream).
    """
    trimmed = SeqRecord(nt_record.seq[200:], id=nt_record.id, description="")
    return trimmed


# ---------------- depth / breadth parsing ----------------
def compute_breadth_from_samtools_depth(depth_txt_path: str, min_depth: int = 10) -> Tuple[float, float]:
    """
    Return (breadth_fraction, mean_depth).

    breadth_fraction: fraction of sites with depth >= min_depth across all reported sites.
    mean_depth: arithmetic mean depth across all sites reported by 'samtools depth -a'.
    """
    total_sites = 0
    covered_sites = 0
    sum_depth = 0.0
    with open(depth_txt_path) as f:
        for line in f:
            parts = line.rstrip("\n").split()
            if len(parts) != 3:
                continue
            d = int(parts[2])
            total_sites += 1
            sum_depth += d
            if d >= min_depth:
                covered_sites += 1
    breadth = (covered_sites / total_sites) if total_sites > 0 else 0.0
    mean_depth = (sum_depth / total_sites) if total_sites > 0 else 0.0
    return breadth, mean_depth


# ---------------- intervals / strands ----------------
def norm_interval(a: int, b: int) -> Tuple[int, int]:
    """Return (start,end) 1-based inclusive with start <= end."""
    return (a, b) if a <= b else (b, a)


def overlap_len(a1: int, a2: int, b1: int, b2: int) -> int:
    """Length of overlap (in nt) between [a1,a2] and [b1,b2] (1-based inclusive)."""
    s = max(min(a1, a2), min(b1, b2))
    e = min(max(a1, a2), max(b1, b2))
    return max(0, e - s + 1)
