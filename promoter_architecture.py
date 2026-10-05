#!/usr/bin/env python3
"""Promoter architecture and metadata association analysis for pspA upstream alignments."""
import argparse
import csv
import os
from collections import Counter
from typing import Dict, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO
from scipy.stats import chi2_contingency, kruskal


MOTIFS = {"minus35": "TTGACA", "minus10": "TATAAT"}


def _motif_candidates(sequence: str, motif: str) -> List[Tuple[int, str, int]]:
    candidates = []
    for position in range(len(sequence) - len(motif) + 1):
        observed = sequence[position:position + len(motif)].upper()
        if "-" in observed or "N" in observed:
            continue
        mismatches = sum(a != b for a, b in zip(observed, motif))
        if mismatches <= 2:
            candidates.append((position, observed, mismatches))
    return candidates


def find_best_promoter(sequence: str) -> Dict[str, object]:
    minus35 = _motif_candidates(sequence, MOTIFS["minus35"])
    minus10 = _motif_candidates(sequence, MOTIFS["minus10"])
    pairs = []
    for pos35, seq35, mm35 in minus35:
        for pos10, seq10, mm10 in minus10:
            spacer = pos10 - (pos35 + len(MOTIFS["minus35"]))
            if 14 <= spacer <= 20:
                pairs.append((mm35 + mm10, abs(spacer - 17), mm35, mm10,
                              pos35, seq35, pos10, seq10, spacer))
    if not pairs:
        return {
            "promoter_type": "NoPromoterDetected",
            "minus35_seq": "", "minus35_pos": "", "minus35_mismatches": "",
            "minus35_start": "", "minus35_end": "",
            "minus10_seq": "", "minus10_pos": "", "minus10_mismatches": "",
            "minus10_start": "", "minus10_end": "",
            "spacer_length": "", "spacer_start": "", "spacer_end": "",
        }
    best = min(pairs)
    total_mm, _, mm35, mm10, pos35, seq35, pos10, seq10, spacer = best
    if total_mm == 0 and 16 <= spacer <= 18:
        promoter_type = "Canonical"
    elif total_mm <= 2 and 16 <= spacer <= 18:
        promoter_type = "NearCanonical"
    else:
        promoter_type = "Weak"
    return {
        "promoter_type": promoter_type,
        "minus35_seq": seq35,
        "minus35_pos": pos35 + 1,
        "minus35_start": pos35 + 1,
        "minus35_end": pos35 + len(MOTIFS["minus35"]),
        "minus35_mismatches": mm35,
        "minus10_seq": seq10,
        "minus10_pos": pos10 + 1,
        "minus10_start": pos10 + 1,
        "minus10_end": pos10 + len(MOTIFS["minus10"]),
        "minus10_mismatches": mm10,
        "spacer_length": spacer,
        "spacer_start": pos35 + len(MOTIFS["minus35"]) + 1,
        "spacer_end": pos10,
    }


def _gff_attributes(**values: object) -> str:
    return ";".join(
        f"{key}={str(value).replace(';', '%3B').replace(' ', '_')}"
        for key, value in values.items()
        if value is not None and str(value) != ""
    )


def write_promoter_annotation(
    alignment_fasta: str,
    output_path: str,
    reference_id: Optional[str] = None,
    pspa_length: Optional[int] = None,
) -> Dict[str, object]:
    """Write promoter features in the alignment's 1-based, gap-preserving coordinates."""
    records = list(SeqIO.parse(alignment_fasta, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in alignment: {alignment_fasta}")
    reference = next(
        (record for record in records if reference_id and record.id == reference_id),
        records[0],
    )
    if reference_id and reference.id != reference_id:
        raise ValueError(f"Reference ID {reference_id!r} is not present in the alignment")
    ungapped_reference = str(reference.seq).replace("-", "")
    result = find_best_promoter(ungapped_reference)
    alignment_length = len(reference.seq)
    pspa_start = alignment_length + 1
    with open(output_path, "w", newline="") as handle:
        handle.write("##gff-version 3\n")
        handle.write("##sequence-region {} 1 {}\n".format(
            reference.id,
            alignment_length + (pspa_length or 0)
        ))

        def feature(start: int, end: int, feature_type: str, name: str, **extra: object) -> None:
            handle.write(
                "{}\tpspA_promoter\t{}\t{}\t{}\t.\t+\t.\t{}\n".format(
                    reference.id,
                    feature_type, start, end,
                    _gff_attributes(ID=name, Name=name, **extra),
                )
            )

        if result["promoter_type"] != "NoPromoterDetected":
            alignment_positions = [
                index + 1 for index, base in enumerate(str(reference.seq))
                if base != "-"
            ]

            def map_position(ungapped_position: int) -> int:
                return alignment_positions[ungapped_position - 1]

            minus35_start = map_position(int(result["minus35_pos"]))
            minus35_end = map_position(int(result["minus35_pos"]) + len(MOTIFS["minus35"]) - 1)
            minus10_start = map_position(int(result["minus10_pos"]))
            minus10_end = map_position(int(result["minus10_pos"]) + len(MOTIFS["minus10"]) - 1)
            spacer_start = minus35_end + 1
            spacer_end = minus10_start - 1
            feature(minus35_start, minus35_end, "promoter_minus35", "minus35")
            if spacer_start <= spacer_end:
                feature(spacer_start, spacer_end, "promoter_spacer", "spacer")
            feature(minus10_start, minus10_end, "promoter_minus10", "minus10")

        feature(pspa_start, pspa_start, "gene_start", "pspA_start")
        if pspa_length:
            if pspa_length < 1:
                raise ValueError("pspa_length must be at least 1")
            feature(pspa_start, pspa_start + pspa_length - 1, "CDS", "pspA_ORF")
    return result


def _consensus_column(column: str) -> str:
    symbols = [base.upper() for base in column if base.upper() in "ACGT-"]
    return Counter(symbols).most_common(1)[0][0] if symbols else "-"


def write_alignment_variants(
    alignment_fasta: str,
    variants_tsv: str,
    variants_gff: str,
    reference_id: Optional[str] = None,
) -> int:
    """Write all alignment-column differences from the majority consensus."""
    records = list(SeqIO.parse(alignment_fasta, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in alignment: {alignment_fasta}")
    lengths = {len(record.seq) for record in records}
    if len(lengths) != 1:
        raise ValueError("Alignment sequences must all have the same length")
    reference = next(
        (record for record in records if reference_id and record.id == reference_id),
        records[0],
    )
    if reference_id and reference.id != reference_id:
        raise ValueError(f"Reference ID {reference_id!r} is not present in the alignment")
    consensus = "".join(
        _consensus_column("".join(str(record.seq)[index] for record in records))
        for index in range(len(reference.seq))
    )
    events = []
    for index in range(len(consensus)):
        ref = consensus[index]
        if ref == "-":
            # Alignment columns with a consensus gap are insertion columns.
            for record in records:
                alt = str(record.seq)[index].upper()
                if alt in "ACGT":
                    events.append((record.id, index + 1, ref, alt, "Insertion"))
        else:
            for record in records:
                alt = str(record.seq)[index].upper()
                if alt == ref or alt not in "-ACGT":
                    continue
                event_type = "Deletion" if alt == "-" else "SNP"
                events.append((record.id, index + 1, ref, alt, event_type))
    with open(variants_tsv, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["ID", "Position", "Reference", "Alternate", "Type"])
        writer.writerows(events)
    with open(variants_gff, "w", newline="") as handle:
        handle.write("##gff-version 3\n")
        handle.write(f"##sequence-region {reference.id} 1 {len(reference.seq)}\n")
        for number, (sample_id, position, ref, alt, event_type) in enumerate(events, 1):
            attributes = _gff_attributes(
                ID=f"variant_{number}",
                Name=f"{event_type}_{position}",
                Sample=sample_id,
                Reference=ref,
                Alternate=alt,
            )
            handle.write(
                f"{reference.id}\tpspA_variation\t{event_type}\t{position}\t{position}"
                f"\t.\t+\t.\t{attributes}\n"
            )
    with open(os.path.join(os.path.dirname(variants_tsv), "consensus.fasta"), "w") as handle:
        handle.write(f">{reference.id}_consensus\n{consensus}\n")
    return len(events)


def _read_metadata(path: str) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    raw = open(path, "rb").read()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    delimiter = ";" if text.splitlines()[0].count(";") > text.splitlines()[0].count(",") else ","
    return pd.read_csv(pd.io.common.StringIO(text), sep=delimiter, dtype=str).fillna("")


def _write_associations(path: str, table: pd.DataFrame) -> None:
    rows = []
    def categorical(left: str, right: str) -> None:
        if left not in table.columns or right not in table.columns:
            rows.append((f"{left}_vs_{right}", left, right, 0, "", "insufficient data"))
            return
        subset = table[[left, right]].replace("", pd.NA).dropna()
        if subset.empty or subset[left].nunique() < 2 or subset[right].nunique() < 2:
            rows.append((f"{left}_vs_{right}", left, right, len(subset), "", "insufficient data"))
            return
        contingency = pd.crosstab(subset[left], subset[right])
        statistic, pvalue, _, _ = chi2_contingency(contingency)
        rows.append((f"{left}_vs_{right}", left, right, len(subset),
                     float(statistic), float(pvalue)))
    categorical("Phenotype", "promoter_type")
    categorical("pspA_family", "promoter_type")
    categorical("Phenotype", "pspA_family")

    if "Phenotype" not in table.columns or "spacer_length" not in table.columns:
        subset = pd.DataFrame()
    else:
        subset = table[["Phenotype", "spacer_length"]].copy()
    if not subset.empty:
        subset["spacer_length"] = pd.to_numeric(subset["spacer_length"], errors="coerce")
        subset = subset.dropna()
    groups = (
        [values["spacer_length"].to_numpy() for _, values in subset.groupby("Phenotype")
         if len(values) >= 2]
        if not subset.empty else []
    )
    if len(groups) >= 2:
        if subset["spacer_length"].nunique() == 1:
            rows.append(("Phenotype_vs_spacer_length", "Phenotype", "spacer_length",
                         len(subset), "", "not testable: identical spacer lengths"))
        else:
            try:
                statistic, pvalue = kruskal(*groups)
            except ValueError as error:
                if "identical" not in str(error).lower():
                    raise
                rows.append(("Phenotype_vs_spacer_length", "Phenotype", "spacer_length",
                             len(subset), "", "not testable: identical spacer lengths"))
            else:
                rows.append(("Phenotype_vs_spacer_length", "Phenotype", "spacer_length",
                             len(subset), float(statistic), float(pvalue)))
    else:
        rows.append(("Phenotype_vs_spacer_length", "Phenotype", "spacer_length",
                     len(subset), "", "insufficient data"))
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["analysis", "variable_a", "variable_b", "n", "statistic", "pvalue"])
        writer.writerows(rows)


def build_promoter_architecture(
    alignment_fasta: str,
    output_dir: str,
    phandango_metadata: Optional[str] = None,
    reference_id: Optional[str] = None,
    pspa_length: Optional[int] = None,
) -> str:
    os.makedirs(output_dir, exist_ok=True)
    promoter_rows = []
    for record in SeqIO.parse(alignment_fasta, "fasta"):
        row = {"ID": record.id}
        row.update(find_best_promoter(str(record.seq)))
        promoter_rows.append(row)
    promoter = pd.DataFrame(promoter_rows)
    promoter.to_csv(os.path.join(output_dir, "promoter_metadata.csv"), index=False)
    write_promoter_annotation(
        alignment_fasta,
        os.path.join(output_dir, "promoter_annotation.gff"),
        reference_id=reference_id,
        pspa_length=pspa_length,
    )
    write_alignment_variants(
        alignment_fasta,
        os.path.join(output_dir, "upstream_variants.tsv"),
        os.path.join(output_dir, "variation_annotation.gff"),
        reference_id=reference_id,
    )
    metadata = _read_metadata(phandango_metadata) if phandango_metadata else pd.DataFrame()
    if not metadata.empty:
        metadata.columns = [str(column).strip() for column in metadata.columns]
        if "ID" not in metadata.columns:
            raise ValueError("Phandango metadata must contain an ID column matching FASTA headers")
        metadata = metadata[metadata["ID"].isin(promoter["ID"])]
        merged = metadata.merge(promoter, on="ID", how="inner")
        merged.to_csv(os.path.join(output_dir, "phandango_promoter_metadata.csv"), index=False)
        _write_associations(os.path.join(output_dir, "promoter_associations.csv"), merged)
    else:
        promoter.to_csv(os.path.join(output_dir, "phandango_promoter_metadata.csv"), index=False)
        _write_associations(os.path.join(output_dir, "promoter_associations.csv"), promoter)
    promoter.to_csv(os.path.join(output_dir, "promoter_summary.csv"), index=False)
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyse pspA upstream promoter architecture")
    parser.add_argument("--alignment", required=True, help="upstream_alignment.fasta")
    parser.add_argument("--metadata", default=None, help="phandango_metadata.csv")
    parser.add_argument(
        "--reference-id",
        default=None,
        help="FASTA ID whose alignment coordinates define promoter annotations (default: first record)",
    )
    parser.add_argument(
        "--pspa-length",
        type=int,
        default=None,
        help="optional pspA ORF length in alignment-coordinate bases",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="output directory (default: sibling 11_promoter_analysis beside the upstream folder)",
    )
    args = parser.parse_args()
    output_dir = args.output_dir
    if output_dir is None:
        upstream_dir = os.path.dirname(os.path.abspath(args.alignment))
        output_dir = os.path.join(os.path.dirname(upstream_dir), "11_promoter_analysis")
    build_promoter_architecture(
        args.alignment,
        output_dir,
        args.metadata,
        reference_id=args.reference_id,
        pspa_length=args.pspa_length,
    )
    print(f"[INFO] Wrote promoter analysis outputs to {output_dir}")


if __name__ == "__main__":
    main()
