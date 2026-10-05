#!/usr/bin/env python3

import argparse
import pandas as pd


def main():

    ap = argparse.ArgumentParser(
        description="Generate Phandango metadata table for pspA upstream analysis"
    )

    ap.add_argument(
        "--metadata",
        required=True,
        help="Metadata_ALL.csv"
    )

    ap.add_argument(
        "--confirmed",
        required=True,
        help="confirmed_pspA_hits.csv"
    )

    ap.add_argument(
        "--upstream-metadata",
        required=True,
        help="upstream_metadata.csv"
    )

    ap.add_argument(
        "--out",
        required=True,
        help="phandango_metadata.csv"
    )

    args = ap.parse_args()

    def read_table(path):
        raw = open(path, "rb").read()
        try:
            text = raw.decode("utf-8-sig")
        except UnicodeDecodeError:
            text = raw.decode("latin-1")
        delimiter = ";" if text.splitlines()[0].count(";") > text.splitlines()[0].count(",") else ","
        return pd.read_csv(pd.io.common.StringIO(text), sep=delimiter, dtype=str)

    md = read_table(args.metadata).fillna("")
    confirmed = read_table(args.confirmed).fillna("")
    upstream = read_table(args.upstream_metadata).fillna("")

    md.columns = md.columns.str.strip()
    confirmed.columns = confirmed.columns.str.strip()
    upstream.columns = upstream.columns.str.strip()
    for frame in (md, confirmed, upstream):
        aliases = {
            str(column).lower().replace("_", "").replace(" ", ""): column
            for column in frame.columns
        }
        for canonical in ("Filename", "Phenotype", "GPSC", "ST", "Serotype", "Year", "Project"):
            source = aliases.get(canonical.lower())
            if source and canonical not in frame.columns:
                frame[canonical] = frame[source]

    if "CodingID" not in confirmed.columns:
        raise ValueError("confirmed hits must contain CodingID")
    confirmed_subset = confirmed.copy()
    confirmed_subset["ID"] = confirmed_subset["CodingID"]
    if "Assembly" in confirmed_subset.columns:
        confirmed_subset["Filename"] = confirmed_subset["Assembly"]
    if "CodingID" in upstream.columns:
        upstream_subset = upstream.drop_duplicates("CodingID").rename(columns={"CodingID": "ID"})
        out = confirmed_subset.merge(upstream_subset, on="ID", how="left", suffixes=("", "_upstream"))
    else:
        out = confirmed_subset
    if "Filename" not in out.columns:
        out["Filename"] = ""
    if "Filename" in md.columns:
        md_subset = md.drop_duplicates("Filename")
        metadata_columns = [
            c for c in ["Phenotype", "GPSC", "ST", "Serotype", "Year", "Project"]
            if c in md_subset.columns
        ]
        if metadata_columns:
            out = out.merge(
                md_subset[["Filename"] + metadata_columns],
                on="Filename",
                how="left",
                suffixes=("", "_metadata"),
            )
            for column in metadata_columns:
                metadata_column = f"{column}_metadata"
                if metadata_column in out.columns:
                    out[column] = out[column].where(out[column].ne(""), out[metadata_column])
                    out = out.drop(columns=[metadata_column])
    if "Upstream_Length" in out.columns and "upstream_length" not in out.columns:
        out["upstream_length"] = out["Upstream_Length"]
    out["pspA_detected"] = "Yes"

    wanted = [
        "ID",
        "Filename",
        "Phenotype",
        "GPSC",
        "ST",
        "Serotype",
        "Year",
        "Project",
        "pspA_detected",
        "Hollingshead_family",
        "pspA_family",
        "pspA_clade",
        "Matching_pspA_reference",
        "upstream_length"
    ]

    existing = [
        c for c in wanted
        if c in out.columns
    ]

    out = out[existing]

    out.to_csv(
        args.out,
        index=False
    )

    print(
        f"[INFO] Wrote {len(out)} rows to {args.out}"
    )


if __name__ == "__main__":
    main()