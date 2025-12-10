#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single-FASTA BoxB classifier (rcore version):
- Input: nucleotide FASTA (one record; pspA coding sequence, NO upstream).
- BLASTX vs BoxB protein set (DB preferred, fallback to subject FASTA).
- Output: one-row CSV and JSON with the best AA BoxB HSP and cluster attribution.
Status labels (harmonized):
 - "Matched"
 - "No_BoxB_match"
This module can be imported (classify_boxb()) or used as CLI.
"""
import os
import argparse
import tempfile
from typing import Dict, Any, Optional

import pandas as pd
from pandas.errors import EmptyDataError
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from pspA_utils import (
    run_cmd,
    blast_db_exists,
    parse_cluster_from_protein_id,
    safe_id,
    write_one_row_csv,
    write_one_row_json,
    normalize_sample_id,
)

HEADER = [
    "Sample", "CodingID", "Status",
    "Matching_pspA_reference", "Hollingshead_family",
    "pident", "bitscore", "AA_Length",
    "AA_BoxB", "evalue", "qstart_nt", "qend_nt", "qframe"
]


def _cap_threads(t: int) -> int:
    return max(1, min(12, int(t)))


def classify_boxb(
    nt_fasta: str,
    sample_id: Optional[str],
    boxb_db: str,
    boxb_subject_fasta: str,
    out_dir: str,
    min_aa_len: int = 60,
    min_pident: float = 80.0,
    sort_field: str = "bitscore",
    threads: int = 12
) -> Dict[str, Any]:
    assert os.path.exists(nt_fasta), f"Input FASTA not found: {nt_fasta}"
    os.makedirs(out_dir, exist_ok=True)

    # Load the (single) record
    recs = list(SeqIO.parse(nt_fasta, "fasta"))
    if len(recs) == 0:
        raise ValueError(f"No records in {nt_fasta}")
    if len(recs) > 1:
        # Proceed with the first record; orchestrator normally supplies single-record FASTAs
        recs = recs[:1]
    rec: SeqRecord = recs[0]
    coding_id = rec.id
    sid = sample_id if sample_id else normalize_sample_id(coding_id)

    # BLASTX command
    outfmt = "6 qseqid sseqid pident length evalue bitscore qstart qend sstart send qframe qseq"
    threads = _cap_threads(threads)
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".fa") as tmpq:
        SeqIO.write(rec, tmpq.name, "fasta")
        tmp_out = os.path.join(out_dir, f"{safe_id(coding_id)}.blastx.txt")
        if blast_db_exists(boxb_db):
            cmd = [
                "blastx", "-query", tmpq.name, "-db", boxb_db,
                "-out", tmp_out, "-outfmt", outfmt, "-num_threads", str(threads)
            ]
        else:
            cmd = [
                "blastx", "-query", tmpq.name, "-subject", boxb_subject_fasta,
                "-out", tmp_out, "-outfmt", outfmt, "-num_threads", str(threads)
            ]
        run_cmd(cmd, check=True)

    # Initialize record (harmonized labels)
    record: Dict[str, Any] = {
        "Sample": sid,
        "CodingID": coding_id,
        "Status": "No_BoxB_match",
        "Matching_pspA_reference": "",
        "Hollingshead_family": "",
        "pident": "",
        "bitscore": "",
        "AA_Length": "",
        "AA_BoxB": "",
        "evalue": "",
        "qstart_nt": "",
        "qend_nt": "",
        "qframe": ""
    }

    # Parse BLAST output (handle empty/no-HSP case explicitly)
    try:
        if (not os.path.exists(tmp_out)) or (os.path.getsize(tmp_out) == 0):
            df = None  # no lines -> no HSPs -> keep No_BoxB_match
        else:
            df = pd.read_csv(tmp_out, sep="\t", header=None)
    except EmptyDataError:
        df = None  # treat empty file as no hits
    except Exception as e:
        record["Status"] = f"ParseError: {e}"
        df = None

    if isinstance(df, pd.DataFrame) and not df.empty:
        df.columns = [
            "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
            "qstart", "qend", "sstart", "send", "qframe", "qseq"
        ]
        df = df[df["length"] >= int(min_aa_len)]
        df = df[df["pident"] >= float(min_pident)]
        if not df.empty:
            if sort_field not in df.columns:
                sort_field_local = "bitscore"
            else:
                sort_field_local = sort_field
            hit = df.sort_values(by=sort_field_local, ascending=False).iloc[0].to_dict()
            aa_seq = str(hit["qseq"]).replace("-", "").replace(".", "")
            cluster = parse_cluster_from_protein_id(str(hit["sseqid"]))
            record.update({
                "Status": "Matched",
                "Matching_pspA_reference": str(hit["sseqid"]),
                "Hollingshead_family": cluster,
                "pident": float(hit["pident"]),
                "bitscore": float(hit["bitscore"]),
                "AA_Length": int(hit["length"]),
                "AA_BoxB": aa_seq,
                "evalue": hit["evalue"],
                "qstart_nt": int(hit["qstart"]),
                "qend_nt": int(hit["qend"]),
                "qframe": int(hit["qframe"])
            })

    # Write outputs
    base = os.path.join(out_dir, f"{safe_id(sid)}__{safe_id(coding_id)}")
    write_one_row_json(record, base + ".boxb.json")
    write_one_row_csv(record, base + ".boxb.csv", HEADER)
    return record


def main():
    ap = argparse.ArgumentParser(description="Single-FASTA pspA BoxB classifier (rcore)")
    ap.add_argument("--input", required=True, help="Nucleotide FASTA (pspA; no upstream)")
    ap.add_argument("--sample-id", default=None, help="Sample ID (default: derive from record id)")
    ap.add_argument("--boxb-db", required=True, help="BLAST protein DB basename for BoxB references")
    ap.add_argument("--boxb-subject-fasta", required=True, help="Fallback FASTA for BoxB references")
    ap.add_argument("--out-dir", required=True, help="Output directory (per-record row CSV+JSON)")
    ap.add_argument("--sort-field", default="bitscore")
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--min-aa-len", type=int, default=60)
    ap.add_argument("--min-pident", type=float, default=80.0)
    args = ap.parse_args()

    classify_boxb(
        nt_fasta=args.input,
        sample_id=args.sample_id,
        boxb_db=args.boxb_db,
        boxb_subject_fasta=args.boxb_subject_fasta,
        out_dir=args.out_dir,
        min_aa_len=args.min_aa_len,
        min_pident=args.min_pident,
        sort_field=args.sort_field,
        threads=args.threads
    )


if __name__ == "__main__":
    main()
