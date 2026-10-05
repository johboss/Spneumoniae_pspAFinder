#!/usr/bin/env python3
"""Extract and analyse upstream regions of BoxB-confirmed pspA loci."""
import csv
import io
import os
import re
import shutil
import subprocess
from typing import Dict, Iterable, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ConvertItolLabelColor import convert_colorstrips


def _find_by_stem(root: str, stem: str, suffixes: Iterable[str]) -> str:
    for directory, _, filenames in os.walk(root):
        for filename in filenames:
            if filename.startswith("._"):
                continue
            base, ext = os.path.splitext(filename)
            if base == stem and ext.lower() in suffixes:
                return os.path.join(directory, filename)
    raise FileNotFoundError(f"Could not find {stem} ({', '.join(sorted(suffixes))}) below {root}")


def _parse_attributes(value: str) -> Dict[str, str]:
    attrs = {}
    for item in value.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            attrs[key] = val
    return attrs


def _read_cds(gff_path: str) -> Dict[str, List[Tuple[int, int]]]:
    cds_by_contig: Dict[str, List[Tuple[int, int]]] = {}
    gff_bytes = open(gff_path, "rb").read()
    try:
        gff_text = gff_bytes.decode("utf-8-sig")
    except UnicodeDecodeError:
        gff_text = gff_bytes.decode("latin-1")
    for line in gff_text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 5 or fields[2].lower() != "cds":
            continue
        try:
            start, end = int(fields[3]), int(fields[4])
        except ValueError:
            continue
        cds_by_contig.setdefault(fields[0], []).append((start, end))
    for intervals in cds_by_contig.values():
        intervals.sort()
    return cds_by_contig


def _extract_upstream(
    contig: SeqRecord,
    cds_by_contig: Dict[str, List[Tuple[int, int]]],
    strand: str,
    start: int,
    end: int,
    max_length: int = 500,
) -> Tuple[Seq, int, Optional[Tuple[int, int]]]:
    intervals = cds_by_contig.get(contig.id, [])
    if strand == "+":
        previous = [interval for interval in intervals if interval[1] < start]
        boundary = max(previous, key=lambda interval: interval[1]) if previous else None
        region_start = (boundary[1] + 1) if boundary else 1
        region_end = start - 1
        sequence = contig.seq[region_start - 1:region_end]
        if len(sequence) > max_length:
            sequence = sequence[-max_length:]
        return sequence, len(sequence), boundary

    following = [interval for interval in intervals if interval[0] > end]
    boundary = min(following, key=lambda interval: interval[0]) if following else None
    region_start = end + 1
    region_end = (boundary[0] - 1) if boundary else len(contig.seq)
    sequence = contig.seq[region_start - 1:region_end].reverse_complement()
    if len(sequence) > max_length:
        sequence = sequence[:max_length]
    return sequence, len(sequence), boundary


def _family_clade(value: str) -> Tuple[str, str]:
    match = re.search(r"(F\d+)_?(C\d+)", value or "")
    return (match.group(1), match.group(2)) if match else ("", "")


def _run_alignment_and_tree(output_dir: str, threads: int) -> None:
    sequences = os.path.join(output_dir, "upstream_sequences.fasta")
    alignment = os.path.join(output_dir, "upstream_alignment.fasta")
    tree = os.path.join(output_dir, "upstream_tree.nwk")
    if not os.path.exists(sequences) or os.path.getsize(sequences) == 0:
        return
    with open(alignment, "w") as alignment_handle:
        subprocess.run(
            ["mafft", "--thread", str(max(1, threads)), "--auto", sequences],
            stdout=alignment_handle,
            check=True,
        )
    fasttree = shutil.which("FastTree") or shutil.which("fasttree")
    if not fasttree:
        raise FileNotFoundError("FastTree executable not found on PATH")
    with open(tree, "w") as tree_handle:
        subprocess.run([fasttree, "-nt", "-gtr", alignment], stdout=tree_handle, check=True)


def _write_itol(output_dir: str, rows: List[Dict[str, str]]) -> None:
    itol_dir = os.path.join(output_dir, "itol")
    os.makedirs(itol_dir, exist_ok=True)
    # Bright, high-contrast colors remain legible on iTOL's light background.
    palette = [
        "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728",
        "#9467BD", "#8C564B", "#E377C2", "#7F7F7F",
        "#BCBD22", "#17BECF",
    ]
    colorstrip_fields = {
        "pspA_clade": "pspA clade",
        "pspA_family": "pspA family",
        "Serotype": "Serotype",
        "GPSC": "GPSC",
        "ST": "ST",
        "Year": "Year",
        "Phenotype": "Phenotype",
    }
    for column, label in colorstrip_fields.items():
        values = sorted({str(row.get(column, "")).strip() for row in rows
                         if str(row.get(column, "")).strip()})
        if not values:
            continue
        colors = {value: palette[index % len(palette)] for index, value in enumerate(values)}
        path = os.path.join(itol_dir, f"{column}.txt")
        with open(path, "w", newline="") as handle:
            handle.write(
                "DATASET_COLORSTRIP\n"
                "SEPARATOR TAB\n"
                f"DATASET_LABEL\t{label}\n"
                "COLOR\t#333333\n"
                f"LEGEND_TITLE\t{label}\n"
                f"LEGEND_SHAPES\t{' '.join(['1'] * len(values))}\n"
                f"LEGEND_COLORS\t{' '.join(colors[value] for value in values)}\n"
                f"LEGEND_LABELS\t{' '.join(values)}\n"
                "DATA\n"
            )
            for row in rows:
                value = str(row.get(column, "")).strip()
                if value:
                    handle.write(f"{row['CodingID']}\t{colors[value]}\t{value}\n")

    lengths = [
        float(row["upstream_length"])
        for row in rows
        if str(row.get("upstream_length", "")).strip()
    ]
    if lengths:
        minimum, maximum = min(lengths), max(lengths)
        if minimum == maximum:
            maximum = minimum + 1
        with open(os.path.join(itol_dir, "upstream_length.txt"), "w", newline="") as handle:
            handle.write(
                "DATASET_GRADIENT\n"
                "SEPARATOR TAB\n"
                "DATASET_LABEL\tUpstream length\n"
                "COLOR_MIN\t#56B4E9\n"
                "COLOR_MAX\t#D55E00\n"
                "LEGEND_TITLE\tUpstream length\n"
                f"LEGEND_LABELS\t{minimum:g} {maximum:g}\n"
                "LEGEND_COLORS\t#56B4E9 #D55E00\n"
                "DATA\n"
            )
            for row in rows:
                value = str(row.get("upstream_length", "")).strip()
                if value:
                    handle.write(f"{row['CodingID']}\t{float(value):g}\n")
    _write_node_popup(itol_dir, rows)


def _write_node_popup(itol_dir: str, rows: List[Dict[str, str]]) -> None:
    """Write an iTOL text dataset shown when a node is selected."""
    metadata_fields = [
        ("Assembly", "Assembly"),
        ("CodingID", "CodingID"),
        ("GPSC", "GPSC"),
        ("ST", "ST"),
        ("Serotype", "Serotype"),
        ("Year", "Year"),
        ("Phenotype", "Phenotype"),
        ("pspA_family", "pspA family"),
        ("pspA_clade", "pspA clade"),
        ("pspA_hollingshead", "Hollingshead family"),
        ("upstream_length", "Upstream length"),
        ("ORF_confidence", "ORF confidence"),
    ]
    path = os.path.join(itol_dir, "node_popup.txt")
    with open(path, "w", newline="") as handle:
        handle.write(
            "DATASET_TEXT\nSEPARATOR TAB\n"
            "DATASET_LABEL\tpspA metadata popup\n"
            "COLOR\t#333333\nSHOW_LABELS\t0\nDATA\n"
        )
        for row in rows:
            values = [
                f"{label}: {row.get(column, '')}"
                for column, label in metadata_fields
                if str(row.get(column, "")).strip()
            ]
            handle.write(f"{row['CodingID']}\t{' | '.join(values)}\n")


def _write_phandango_metadata(output_dir: str, rows: List[Dict[str, str]]) -> None:
    """Write metadata with IDs matching the upstream tree tip names exactly."""
    fields = [
        "ID", "Filename", "Phenotype", "GPSC", "ST", "Serotype", "Year",
        "Project", "pspA_family", "pspA_clade", "pspA_hollingshead",
        "upstream_length",
    ]
    with open(os.path.join(output_dir, "phandango_metadata.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                field: row.get("CodingID", "") if field == "ID" else row.get(field, "")
                for field in fields
            })


def build_upstream_phylogeny(
    confirmed_hits_csv: str,
    assemblies_root: str,
    gff_root: str,
    metadata_csv: str,
    output_dir: str,
    threads: int = 12,
    max_upstream_length: int = 500,
) -> str:
    """Extract confirmed upstream regions, write metadata, alignment, tree and iTOL files."""
    os.makedirs(output_dir, exist_ok=True)
    if max_upstream_length < 1:
        raise ValueError("max_upstream_length must be at least 1")
    with open(confirmed_hits_csv, newline="") as handle:
        hits = list(csv.DictReader(handle))
    confirmed = [row for row in hits if row.get("Status") in ("Matched", "Rescue_Matched")]
    copy_numbers: Dict[str, int] = {}
    for row in confirmed:
        copy_numbers[row["Assembly"]] = copy_numbers.get(row["Assembly"], 0) + 1

    metadata_by_filename: Dict[str, Dict[str, str]] = {}
    metadata_bytes = open(metadata_csv, "rb").read()
    try:
        metadata_text = metadata_bytes.decode("utf-8-sig")
    except UnicodeDecodeError:
        metadata_text = metadata_bytes.decode("latin-1")
    first_line = metadata_text.splitlines()[0] if metadata_text.splitlines() else ""
    delimiter = ";" if first_line.count(";") > first_line.count(",") else ","
    metadata_reader = csv.DictReader(io.StringIO(metadata_text, newline=""), delimiter=delimiter)
    metadata_reader.fieldnames = [
        field.strip().lstrip("\ufeff") if field else field
        for field in (metadata_reader.fieldnames or [])
    ]
    for raw_row in metadata_reader:
        row = {
            key.strip().lstrip("\ufeff"): value.strip()
            for key, value in raw_row.items()
            if key is not None
        }
        aliases = {key.lower().replace("_", "").replace(" ", ""): key for key in row}
        for canonical in ("Filename", "Phenotype", "GPSC", "ST", "Serotype", "Year", "Project"):
            source = aliases.get(canonical.lower())
            if source and canonical not in row:
                row[canonical] = row[source]
        filename = row.get("Filename", "")
        metadata_by_filename[os.path.splitext(filename)[0]] = row

    upstream_rows: List[Dict[str, str]] = []
    fasta_path = os.path.join(output_dir, "upstream_sequences.fasta")
    with open(fasta_path, "w") as fasta:
        for hit in confirmed:
            sample = hit["Assembly"]
            try:
                assembly_path = _find_by_stem(assemblies_root, sample, {".fasta", ".fa", ".fna"})
                gff_path = _find_by_stem(gff_root, sample, {".gff", ".gff3"})
            except FileNotFoundError as exc:
                raise FileNotFoundError(f"{exc} (required for confirmed pspA {sample})") from exc
            contigs = SeqIO.to_dict(SeqIO.parse(assembly_path, "fasta"))
            if hit["Contig"] not in contigs:
                raise KeyError(f"Contig {hit['Contig']} from {sample} is absent in {assembly_path}")
            orf_start = int(hit.get("ORF_start") or hit["qstart"])
            orf_end = int(hit.get("ORF_end") or hit["qend"])
            sequence, length, _ = _extract_upstream(
                contigs[hit["Contig"]], _read_cds(gff_path), hit["Strand"],
                orf_start, orf_end, max_upstream_length,
            )
            coding_id = hit["CodingID"]
            fasta.write(f">{coding_id}\n{sequence}\n")
            family, clade = _family_clade(hit.get("Hollingshead_family", ""))
            row = dict(metadata_by_filename.get(sample, {}))
            row.update({
                "Assembly": sample,
                "CodingID": coding_id,
                "Contig": hit["Contig"],
                "Strand": hit["Strand"],
                "qstart": hit["qstart"],
                "qend": hit["qend"],
                "BLASTX_start": hit.get("BLASTX_start", hit["qstart"]),
                "BLASTX_end": hit.get("BLASTX_end", hit["qend"]),
                "ORF_start": str(orf_start),
                "ORF_end": str(orf_end),
                "ORF_confidence": hit.get("ORF_confidence", ""),
                "PSPA_CONTIG_START": hit.get("PSPA_CONTIG_START", ""),
                "PSPA_CONTIG_END": hit.get("PSPA_CONTIG_END", ""),
                "UPSTREAM_CONTIG_BREAK": hit.get("UPSTREAM_CONTIG_BREAK", ""),
                "GPSC": row.get("GPSC", ""),
                "ST": row.get("ST", ""),
                "Serotype": row.get("Serotype", ""),
                "Phenotype": row.get("Phenotype", ""),
                "pspA_detected": "yes",
                "pspA_hollingshead": hit.get("Hollingshead_family", ""),
                "pspA_family": hit.get("pspA_family", "") or family,
                "pspA_clade": hit.get("pspA_clade", "") or clade,
                "pspA_reference": hit.get("Matching_pspA_reference", ""),
                "pspA_copy_number": str(copy_numbers[sample]),
                "upstream_length": str(length),
            })
            upstream_rows.append(row)

    upstream_metadata = os.path.join(output_dir, "upstream_metadata.csv")
    base_metadata_columns = (
        list(metadata_by_filename[next(iter(metadata_by_filename))].keys())
        if metadata_by_filename else []
    )
    columns = list(dict.fromkeys(base_metadata_columns + [
        "Assembly", "CodingID", "Contig", "Strand", "qstart", "qend",
        "BLASTX_start", "BLASTX_end", "ORF_start", "ORF_end",
        "ORF_confidence", "PSPA_CONTIG_START", "PSPA_CONTIG_END",
        "UPSTREAM_CONTIG_BREAK",
        "pspA_detected", "pspA_hollingshead", "pspA_family", "pspA_clade",
        "pspA_reference", "pspA_copy_number", "upstream_length"
    ]))
    if upstream_rows:
        columns = list(dict.fromkeys(columns + list(upstream_rows[0].keys())))
    with open(upstream_metadata, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows({column: row.get(column, "") for column in columns} for row in upstream_rows)

    metadata_output = os.path.join(output_dir, "metadata_ALL_with_pspA.csv")
    all_columns = list(dict.fromkeys(
        (list(metadata_by_filename[next(iter(metadata_by_filename))].keys())
         if metadata_by_filename else [])
        + ["pspA_detected", "pspA_hollingshead", "pspA_family", "pspA_clade",
           "pspA_reference", "pspA_copy_number", "upstream_length"]
    ))
    by_sample = {row["Assembly"]: row for row in upstream_rows}
    with open(metadata_output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_columns)
        writer.writeheader()
        for sample, original in metadata_by_filename.items():
            update = by_sample.get(sample, {})
            result = dict(original)
            result.update({column: update.get(column, "") for column in all_columns if column not in original})
            if update:
                result.update({column: update.get(column, "") for column in
                               ("pspA_detected", "pspA_hollingshead", "pspA_family",
                                "pspA_clade", "pspA_reference", "pspA_copy_number",
                                "upstream_length")})
            writer.writerow({column: result.get(column, "") for column in all_columns})

    _write_itol(output_dir, upstream_rows)
    _write_phandango_metadata(output_dir, upstream_rows)
    convert_colorstrips(os.path.join(output_dir, "itol"))
    _run_alignment_and_tree(output_dir, threads)
    return output_dir
