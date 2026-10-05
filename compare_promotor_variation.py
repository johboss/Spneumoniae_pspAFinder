#!/usr/bin/env python3
"""Compare concatenated upstream + pspA sequences within selected strain groups."""

import argparse
import csv
import itertools
import os
import re
import shutil
import subprocess
from typing import Dict, Iterable, List, Sequence, Tuple

from Bio import SeqIO


DEFAULT_UPSTREAM = "/Volumes/T7/pspA/02_JPIAMR_B/08_upstream_phylogeny/upstream_sequences.fasta"
DEFAULT_GENE = "/Volumes/T7/pspA/02_JPIAMR_B/06_boxb_reports/confirmed_pspA_sequences.fasta"
DEFAULT_OUTPUT = "/Volumes/T7/pspA/02_JPIAMR_B/12_promotor_variation_analysis"
STRAIN_PATTERN = re.compile(r"^(j\d+)", re.IGNORECASE)
GROUPS = {
    "group1": ("j484", "j539"),
    "group2": ("j319", "j283"),
    "group3": ("j320", "j242", "j327"),
    "group4": ("j355", "j571", "j528"),
}


def strain_key(identifier: str) -> str:
    match = STRAIN_PATTERN.match(identifier.split()[0])
    if not match:
        raise ValueError(f"Could not determine j### strain prefix from FASTA ID {identifier!r}")
    return match.group(1).lower()


def read_unique_sequences(path: str) -> Dict[str, Tuple[str, str]]:
    records: Dict[str, Tuple[str, str]] = {}
    for record in SeqIO.parse(path, "fasta"):
        key = strain_key(record.id)
        if key in records:
            raise ValueError(f"Multiple FASTA records found for strain {key} in {path}")
        records[key] = (record.id, str(record.seq).upper())
    return records


def write_combined_group(
    group_name: str,
    strains: Sequence[str],
    upstream: Dict[str, Tuple[str, str]],
    genes: Dict[str, Tuple[str, str]],
    output_dir: str,
) -> str:
    path = os.path.join(output_dir, f"{group_name}_combined.fasta")
    with open(path, "w") as handle:
        for strain in strains:
            key = strain.lower()
            if key not in upstream:
                raise KeyError(f"Strain {strain} is missing from upstream FASTA")
            if key not in genes:
                raise KeyError(f"Strain {strain} is missing from confirmed pspA FASTA")
            handle.write(f">{strain}\n{upstream[key][1]}{genes[key][1]}\n")
    return path


def align_group(combined_path: str, output_path: str, threads: int) -> None:
    mafft = shutil.which("mafft")
    if not mafft:
        raise FileNotFoundError("MAFFT executable not found on PATH")
    command = [mafft, "--auto", "--thread", str(max(1, threads)), combined_path]
    with open(output_path, "w") as handle:
        subprocess.run(command, stdout=handle, check=True)


def pairwise_difference(left: str, right: str) -> Tuple[int, int, int]:
    if len(left) != len(right):
        raise ValueError("Aligned sequences must have equal lengths")
    substitutions = 0
    indel_columns = 0
    for left_base, right_base in zip(left.upper(), right.upper()):
        if left_base == right_base:
            continue
        if left_base == "-" or right_base == "-":
            indel_columns += 1
        else:
            substitutions += 1
    return substitutions + indel_columns, substitutions, indel_columns


def summarize_group(group_name: str, alignment_path: str) -> List[Dict[str, object]]:
    records = list(SeqIO.parse(alignment_path, "fasta"))
    rows = []
    for left, right in itertools.combinations(records, 2):
        total, substitutions, indels = pairwise_difference(str(left.seq), str(right.seq))
        rows.append({
            "Group": group_name,
            "Strain_1": left.id,
            "Strain_2": right.id,
            "Aligned_length": len(left.seq),
            "Total_differences": total,
            "SNP_differences": substitutions,
            "Indel_differences": indels,
        })
    return rows


def compare_groups(
    upstream_fasta: str,
    gene_fasta: str,
    output_dir: str,
    threads: int = 1,
    groups: Dict[str, Sequence[str]] = GROUPS,
) -> str:
    os.makedirs(output_dir, exist_ok=True)
    upstream = read_unique_sequences(upstream_fasta)
    genes = read_unique_sequences(gene_fasta)
    summary: List[Dict[str, object]] = []
    for group_name, strains in groups.items():
        combined = write_combined_group(group_name, strains, upstream, genes, output_dir)
        aligned = os.path.join(output_dir, f"{group_name}_aligned.fasta")
        align_group(combined, aligned, threads)
        summary.extend(summarize_group(group_name, aligned))

    summary_path = os.path.join(output_dir, "pairwise_difference_summary.csv")
    with open(summary_path, "w", newline="") as handle:
        fields = [
            "Group", "Strain_1", "Strain_2", "Aligned_length",
            "Total_differences", "SNP_differences", "Indel_differences",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary)
    with open(os.path.join(output_dir, "comparison_parameters.txt"), "w") as handle:
        handle.write(f"Upstream FASTA: {upstream_fasta}\n")
        handle.write(f"Confirmed pspA FASTA: {gene_fasta}\n")
        handle.write("Difference count: aligned columns differing; substitutions and gap columns included.\n")
        for group_name, strains in groups.items():
            handle.write(f"{group_name}: {', '.join(strains)}\n")
    return summary_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-fasta", default=DEFAULT_UPSTREAM)
    parser.add_argument("--gene-fasta", default=DEFAULT_GENE)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    summary_path = compare_groups(
        args.upstream_fasta,
        args.gene_fasta,
        args.output_dir,
        args.threads,
    )
    print(f"[INFO] Wrote pairwise summary to {summary_path}")


if __name__ == "__main__":
    main()
