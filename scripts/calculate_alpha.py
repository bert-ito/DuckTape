#!/usr/bin/env python3
"""
Estimate global antisense leakage (p_e) for strand-specific RNA-seq.

Optimisations over the original:
  * Only Geneid and the count column are parsed. featureCounts stores
    Chr/Start/End/Strand as semicolon-delimited per-exon lists, which for a
    multi-exon gene is ~97% of the row's bytes and is never used here.
  * Explicit dtypes, avoiding object-dtype fallback on large tables.

Numerical output is identical to the original implementation.
"""

import argparse

import pandas as pd


def featurecounts_columns(path):
    """Return the header fields of a featureCounts table."""
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            return line.rstrip("\n").split("\t")
    raise ValueError(f"No header row found in {path}")


def read_counts(path, label):
    """Read a featureCounts table keeping only Geneid and the count column."""
    header = featurecounts_columns(path)
    n_cols = len(header)

    if n_cols < 2:
        raise ValueError(f"{path}: expected at least 2 columns, found {n_cols}")

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=0,
        usecols=[0, n_cols - 1],
        names=["Geneid", label],
        dtype={"Geneid": "string"},
        engine="c",
        low_memory=False,
    )

    df[label] = pd.to_numeric(df[label], errors="coerce")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Estimate global antisense leakage (p_e) for "
                    "strand-specific RNA-seq."
    )
    parser.add_argument("-s", "--sense", required=True,
                        help="Sense featureCounts table")
    parser.add_argument("-a", "--antisense", required=True,
                        help="Antisense featureCounts table")
    parser.add_argument("-m", "--min-expression", type=float, default=50,
                        help="Minimum sense counts (default: 50)")
    parser.add_argument("-t", "--min-total", type=float, default=20,
                        help="Minimum total counts (default: 20)")
    parser.add_argument("-r", "--max-ratio", type=float, default=1.0,
                        help="Maximum antisense/sense ratio to keep "
                             "(default: 1.0)")
    parser.add_argument("--median-check", action="store_true",
                        help="Also report median per-gene antisense fraction")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file to save p_e")

    args = parser.parse_args()

    sense = read_counts(args.sense, "sense")
    antisense = read_counts(args.antisense, "antisense")

    df = pd.merge(sense, antisense, on="Geneid", how="inner")
    df = df.dropna(subset=["sense", "antisense"])

    # 1. Minimum sense expression
    df = df[df["sense"] >= args.min_expression].copy()

    # 2. Antisense/sense ratio; 3. drop likely NATs or misannotation
    df["ratio"] = df["antisense"] / df["sense"]
    df = df[df["ratio"] <= args.max_ratio].copy()

    # 4. Total reads; 5. minimum total coverage
    df["total"] = df["sense"] + df["antisense"]
    df = df[df["total"] >= args.min_total].copy()

    total_antisense = df["antisense"].sum()
    total_reads = df["total"].sum()

    if total_reads == 0:
        raise ValueError("No reads remaining after filtering.")

    p_e = total_antisense / total_reads

    median_pe = None
    if args.median_check:
        median_pe = (df["antisense"] / df["total"]).median()

    with open(args.output, "w") as handle:
        handle.write(f"{p_e:.10f}\n")

    print("===== Alpha (p_e) estimation =====")
    print(f"Genes used: {len(df)}")
    print(f"Total reads: {int(total_reads)}")
    print(f"Antisense reads: {int(total_antisense)}")
    print(f"\np_e (global MLE) = {p_e:.6f}")

    if median_pe is not None:
        print(f"p_e (median per-gene) = {median_pe:.6f}")

    print(f"\nSaved to: {args.output}")


if __name__ == "__main__":
    main()
