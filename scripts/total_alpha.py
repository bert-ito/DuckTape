#!/usr/bin/env python3
"""
Calculate alpha MLE by expression bin.
"""

import argparse

import numpy as np
import pandas as pd


def featurecounts_columns(path):
    """Return the header fields of a featureCounts table."""
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            return line.rstrip("\n").split("\t")
    raise ValueError(f"No header row found in {path}")


def read_counts(path, want_length):
    """
    Read a featureCounts table, keeping only the columns we need.

    The count column is always the last one; Length is column 5 (0-based).
    """
    header = featurecounts_columns(path)
    n_cols = len(header)

    if n_cols < 7:
        raise ValueError(
            f"{path}: expected at least 7 columns, found {n_cols}"
        )

    count_idx = n_cols - 1

    if want_length:
        usecols = [0, 5, count_idx]
        names = ["Geneid", "Length", "count"]
    else:
        usecols = [0, count_idx]
        names = ["Geneid", "count"]

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=0,
        usecols=usecols,
        names=names,
        dtype={"Geneid": "string"},
        engine="c",
        low_memory=False,
    )

    df["count"] = pd.to_numeric(df["count"], errors="coerce")
    if want_length:
        df["Length"] = pd.to_numeric(df["Length"], errors="coerce")

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Calculate median alpha by expression bin."
    )
    parser.add_argument("--sense", required=True,
                        help="featureCounts sense file")
    parser.add_argument("--antisense", required=True,
                        help="featureCounts antisense file")
    parser.add_argument("--output", required=True,
                        help="Output TSV")
    parser.add_argument("--min-expression", type=int, default=50,
                        help="Minimum sense counts (default: 50)")
    parser.add_argument("--min-genes-per-bin", type=int, default=1,
                        help="Minimum genes per bin (default: 1)")

    args = parser.parse_args()

    sense = read_counts(args.sense, want_length=True)
    anti = read_counts(args.antisense, want_length=False)

    sense = sense.rename(columns={"count": "sense"})
    anti = anti.rename(columns={"count": "antisense"})

    df = pd.merge(sense, anti, on="Geneid", how="inner")

    df = df.dropna(subset=["sense", "antisense"])
    df = df[df["sense"] >= args.min_expression]
    df = df[df["antisense"] >= 0].copy()

    if df.empty:
        raise SystemExit(
            "No genes passed the expression filter "
            f"(--min-expression {args.min_expression})"
        )

    df["ratio"] = df["antisense"] / (df["sense"] + df["antisense"])
    df = df[df["ratio"] <= 0.5].copy()

    if df.empty:
        raise SystemExit("No genes remained after the ratio <= 0.5 filter")

    # Expression bins of width 100, identical to the original
    bins = np.arange(0, df["sense"].max() + 100, 100)
    df["bin"] = pd.cut(df["sense"], bins=bins, include_lowest=True)

    grouped = df.groupby("bin", observed=False)

    count_by_bin = grouped.size()

    antisense_sum = grouped["antisense"].sum()
    sense_sum = grouped["sense"].sum()

    total_reads = antisense_sum + sense_sum

    alpha_by_bin = antisense_sum / total_reads

    summary = pd.DataFrame({
        "bin": alpha_by_bin.index.astype(str),
        "n_genes": count_by_bin.to_numpy(),
        "antisense_reads": antisense_sum.to_numpy(),
        "total_reads": total_reads.to_numpy(),
        "alpha_mle": alpha_by_bin.to_numpy(),
    })

    summary = summary[
        summary["n_genes"] >= args.min_genes_per_bin
    ].reset_index(drop=True)

    summary.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
