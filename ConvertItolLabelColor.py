#!/usr/bin/env python3
"""Convert iTOL COLORSTRIP datasets to per-label DATASET_RANGE files."""
import argparse
from pathlib import Path


def convert_colorstrips(input_dir: str) -> int:
    directory = Path(input_dir)
    created = 0
    for infile in sorted(directory.glob("*.txt")):
        lines = infile.read_text(encoding="utf-8", errors="replace").splitlines()
        if "DATASET_COLORSTRIP" not in lines:
            continue
        dataset_label = infile.stem
        legend_title = infile.stem
        dataset_color = "#333333"
        legend_colors = []
        legend_labels = []
        data_start = None
        for index, line in enumerate(lines):
            fields = line.split("\t", 1)
            if line.startswith("DATASET_LABEL") and len(fields) > 1:
                dataset_label = fields[1]
            elif line.startswith("COLOR") and len(fields) > 1:
                dataset_color = fields[1]
            elif line.startswith("LEGEND_TITLE") and len(fields) > 1:
                legend_title = fields[1]
            elif line.startswith("LEGEND_COLORS"):
                legend_colors = line.split()[1:]
            elif line.startswith("LEGEND_LABELS"):
                legend_labels = line.split()[1:]
            elif line.strip() == "DATA":
                data_start = index + 1
                break
        if data_start is None:
            continue
        output = infile.with_name(f"{infile.stem}_range.txt")
        with output.open("w", encoding="utf-8", newline="") as handle:
            handle.write("DATASET_RANGE\n")
            handle.write("SEPARATOR COMMA\n")
            handle.write(f"DATASET_LABEL,{dataset_label}\n")
            handle.write(f"COLOR,{dataset_color}\n")
            handle.write("RANGE_TYPE,box\nRANGE_COVER,label\n")
            handle.write("UNROOTED_SMOOTH,simplify\nCOVER_LABELS,0\n")
            handle.write("COVER_DATASETS,0\nFIT_LABELS,0\nSHOW_LABELS,1\n")
            if legend_colors:
                handle.write(f"LEGEND_TITLE,{legend_title}\n")
                handle.write("LEGEND_HORIZONTAL,1\n")
                handle.write(f"LEGEND_SHAPES,{','.join(['1'] * len(legend_colors))}\n")
                handle.write(f"LEGEND_COLORS,{','.join(legend_colors)}\n")
                if legend_labels:
                    handle.write(f"LEGEND_LABELS,{','.join(legend_labels)}\n")
            handle.write("DATA\n")
            for line in lines[data_start:]:
                if not line.strip():
                    continue
                fields = line.split("\t")
                if len(fields) < 3:
                    fields = line.split()
                if len(fields) >= 3:
                    handle.write(f"{fields[0]},{fields[0]},{fields[1]}\n")
        created += 1
    return created


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir")
    args = parser.parse_args()
    print(f"Created {convert_colorstrips(args.input_dir)} iTOL range files")


if __name__ == "__main__":
    main()
