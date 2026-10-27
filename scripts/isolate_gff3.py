#!/usr/bin/env python3
"""
Identify isolated genes from a GFF3 file.

A gene is "isolated" when the gap to its immediate neighbour on the same
chromosome, on either side, exceeds the buffer.

Optimisations over the original:
  * GFF parsing via pandas read_csv + vectorised regex instead of a
    per-line Python loop.
  * Neighbour distances via groupby().shift() instead of df.iloc lookups
    inside a Python loop, which was O(n) pandas scalar indexing.
  * Transparent .gz support.

Selection logic and output are identical to the original.
"""

import argparse

import pandas as pd

GFF_COLUMNS = [
    "chr", "source", "feature", "start", "end",
    "score", "strand", "frame", "attributes",
]


def load_genes(path):
    """Read gene features from a GFF3 into a tidy DataFrame."""
    # comment= is deliberately NOT used: it truncates at '#' anywhere in the
    # line, which would corrupt attribute fields containing '#'.
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=GFF_COLUMNS,
        dtype=str,
        engine="c",
        on_bad_lines="skip",
        low_memory=False,
    )

    # Header lines (##gff-version) and any FASTA section have no feature field
    df = df[df["feature"] == "gene"]

    if df.empty:
        raise ValueError(
            "No gene features were found. Check that the GFF3 contains "
            "entries with feature type 'gene'."
        )

    attrs = df["attributes"]

    # Same precedence as the original: gene_id, then ID, then Name
    gene_id = attrs.str.extract(r"(?:^|;)gene_id=([^;]+)", expand=False)
    gene_id = gene_id.fillna(attrs.str.extract(r"(?:^|;)ID=([^;]+)", expand=False))
    gene_id = gene_id.fillna(attrs.str.extract(r"(?:^|;)Name=([^;]+)", expand=False))

    out = pd.DataFrame({
        "gene_id": gene_id,
        "chr": df["chr"],
        "start": pd.to_numeric(df["start"], errors="coerce"),
        "end": pd.to_numeric(df["end"], errors="coerce"),
        "strand": df["strand"],
    })

    out = out.dropna(subset=["gene_id", "start", "end"])
    out["start"] = out["start"].astype("int64")
    out["end"] = out["end"].astype("int64")

    # Matches the original's gene_id.replace("gene:", "")
    out["gene_id"] = out["gene_id"].str.replace("gene:", "", regex=False)

    return out.reset_index(drop=True)


def find_isolated(df, buffer):
    """
    Return gene IDs whose neighbours on both sides are further than `buffer`.

    Equivalent to the original row-by-row comparison against the previous and
    next row in (chr, start) order, restricted to the same chromosome.
    """
    df = df.sort_values(["chr", "start"], kind="mergesort").reset_index(drop=True)

    grouped = df.groupby("chr", sort=False)
    prev_end = grouped["end"].shift(1)
    next_start = grouped["start"].shift(-1)

    # NaN at chromosome boundaries means "no neighbour" -> not too close
    close_left = (df["start"] - prev_end) <= buffer
    close_right = (next_start - df["end"]) <= buffer

    close_left = close_left.fillna(False).astype(bool)
    close_right = close_right.fillna(False).astype(bool)

    isolated = df.loc[~(close_left | close_right), "gene_id"]

    return df, isolated


def main():
    parser = argparse.ArgumentParser(
        description="Identify isolated genes from a GFF3 file."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input GFF3 file (may be gzipped)")
    parser.add_argument("-b", "--buffer", type=int, required=True,
                        help="Minimum distance (bp) from neighbouring genes")
    parser.add_argument("-o", "--output", required=True,
                        help="Output text file containing isolated gene IDs")

    args = parser.parse_args()

    df = load_genes(args.input)
    df, isolated = find_isolated(df, args.buffer)

    with open(args.output, "w") as out:
        out.write("\n".join(isolated.tolist()))
        if len(isolated):
            out.write("\n")

    print(f"Total genes: {len(df)}")
    print(f"Isolated genes: {len(isolated)}")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
