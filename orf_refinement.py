#!/usr/bin/env python3
"""Refine BLASTX pspA HSP coordinates to complete Prodigal ORF coordinates."""
import os
import re
import subprocess
import tempfile
from typing import Dict, Optional, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def _parse_prodigal_header(description: str) -> Optional[Tuple[int, int, str, bool]]:
    match = re.search(r"#\s*(\d+)\s*#\s*(\d+)\s*#\s*(-?\d+)\s*#", description)
    if not match:
        return None
    start, end, strand = int(match.group(1)), int(match.group(2)), match.group(3)
    partial_match = re.search(r"partial=(\d{2})", description)
    complete = not partial_match or partial_match.group(1) == "00"
    return start, end, "+" if strand == "1" else "-", complete


def _map_region_coordinates(
    window_start: int,
    local_start: int,
    local_end: int,
) -> Tuple[int, int]:
    """Map Prodigal's 1-based inclusive region coordinates to assembly coordinates."""
    return window_start + local_start - 1, window_start + local_end - 1


def _blastp_best(
    orf_fasta: str,
    protein_db: Optional[str],
    protein_subject_fasta: Optional[str],
    output_dir: str,
    threads: int,
    eligible_ids: Optional[set] = None,
) -> Optional[str]:
    output = os.path.join(output_dir, "orf_blastp.tsv")
    target = ["-db", protein_db] if protein_db else ["-subject", protein_subject_fasta]
    subprocess.run(
        [
            "blastp", "-query", orf_fasta, *target,
            "-out", output,
            "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
            "-num_threads", str(max(1, threads)),
        ],
        check=True,
    )
    if not os.path.exists(output) or os.path.getsize(output) == 0:
        return None
    columns = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]
    hits = pd.read_csv(output, sep="\t", names=columns)
    if eligible_ids is not None:
        hits = hits[hits["qseqid"].isin(eligible_ids)]
    if hits.empty:
        return None
    hit = hits.sort_values(
        ["bitscore", "length", "pident"],
        ascending=False,
    ).iloc[0]
    return str(hit["qseqid"])


def refine_pspa_orf(
    assembly_fasta: str,
    contig_id: str,
    strand: str,
    blastx_start: int,
    blastx_end: int,
    protein_db: Optional[str],
    protein_subject_fasta: Optional[str] = None,
    threads: int = 1,
    window: int = 2000,
) -> Dict[str, object]:
    """Return refined genomic coordinates and confidence for one confirmed hit."""
    contigs = SeqIO.to_dict(SeqIO.parse(assembly_fasta, "fasta"))
    if contig_id not in contigs:
        raise KeyError(f"Contig {contig_id} is absent in {assembly_fasta}")
    contig: SeqRecord = contigs[contig_id]
    contig_length = len(contig.seq)
    window_start = max(1, blastx_start - window)
    window_end = min(contig_length, blastx_end + window)
    region = contig.seq[window_start - 1:window_end]
    prodigal_mode = "single" if len(region) >= 20000 else "meta"

    with tempfile.TemporaryDirectory(prefix="pspa_orf_") as temp_dir:
        region_fasta = os.path.join(temp_dir, "region.fna")
        orf_fasta = os.path.join(temp_dir, "region.ffn")
        orf_protein = os.path.join(temp_dir, "region.faa")
        record = SeqRecord(region, id=contig_id, description="")
        with open(region_fasta, "w") as handle:
            SeqIO.write(record, handle, "fasta")
        try:
            subprocess.run(
                [
                    "prodigal", "-i", region_fasta, "-p", prodigal_mode,
                    "-d", orf_fasta, "-a", orf_protein, "-q",
                ],
                check=True,
            )
        except subprocess.CalledProcessError as exc:
            print(
                f"[WARN] Prodigal ORF refinement failed for {contig_id}: "
                f"mode={prodigal_mode}, exit={exc.returncode}; using BLASTX coordinates"
            )
            return _fallback(blastx_start, blastx_end, contig_length, strand)
        orfs = list(SeqIO.parse(orf_protein, "fasta"))
        if not orfs:
            return _fallback(blastx_start, blastx_end, contig_length, strand)
        parsed_orfs = {}
        for orf in orfs:
            parsed = _parse_prodigal_header(orf.description)
            if parsed is None:
                continue
            local_start, local_end, predicted_strand, complete_orf = parsed
            candidate_start, candidate_end = _map_region_coordinates(
                window_start, local_start, local_end
            )
            if predicted_strand == strand and (
                candidate_start <= blastx_start <= candidate_end
                or candidate_start <= blastx_end <= candidate_end
            ):
                parsed_orfs[orf.id] = (
                    candidate_start, candidate_end, predicted_strand, complete_orf
                )
        if not parsed_orfs:
            return _fallback(blastx_start, blastx_end, contig_length, strand)
        best_id = _blastp_best(
            orf_protein, protein_db, protein_subject_fasta,
            temp_dir, threads, set(parsed_orfs)
        )
        if best_id is None:
            return _fallback(blastx_start, blastx_end, contig_length, strand)
        parsed = parsed_orfs.get(best_id)
        if parsed is None:
            return _fallback(blastx_start, blastx_end, contig_length, strand)
    orf_start, orf_end, predicted_strand, complete_orf = parsed
    confidence = "High" if prodigal_mode == "single" else "Medium"
    return _result(orf_start, orf_end, contig_length, confidence, complete_orf, False)


def _result(
    start: int,
    end: int,
    contig_length: int,
    confidence: str,
    complete_orf: bool,
    ambiguous: bool,
) -> Dict[str, object]:
    at_start = start < 200
    at_end = contig_length - end < 200
    upstream_break = at_start or ambiguous
    if at_start or at_end or not complete_orf:
        confidence = "Low" if at_start or at_end or not complete_orf else "Medium"
    return {
        "ORF_start": start,
        "ORF_end": end,
        "ORF_confidence": confidence,
        "PSPA_CONTIG_START": at_start,
        "PSPA_CONTIG_END": at_end,
        "UPSTREAM_CONTIG_BREAK": upstream_break,
    }


def _fallback(
    blastx_start: int,
    blastx_end: int,
    contig_length: int,
    strand: str,
) -> Dict[str, object]:
    return _result(
        blastx_start, blastx_end, contig_length,
        "Low", False, True,
    )
