#!/usr/bin/env python3
"""
Weighted correlation between expression bin midpoint and median alpha.

Changes over the original:
  * Guards against division by zero when |r| == 1 or when the effective
    number of bins drops to 2 or fewer, which previously produced inf/nan
    and a confusing traceback.
  * Guards against zero variance in either axis.
  * scipy is imported at the top rather than mid-script.

Numerical output is identical to the original on well-behaved input.
"""

import argparse
import re

import numpy as np
import pandas as pd
from scipy.stats import t as student_t


def midpoint(bin_str):
    """Midpoint of a pandas Interval rendered as a string, e.g. '(100, 200]'."""
    nums = re.findall(r"[0-9.]+", str(bin_str))
    if len(nums) != 2:
        return np.nan
    return (float(nums[0]) + float(nums[1])) / 2


def main():
    parser = argparse.ArgumentParser(
        description="Weighted correlation of alpha against expression."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    if df.empty:
        raise SystemExit(f"{args.input} contains no rows")

    df["bin_mid"] = df["bin"].apply(midpoint)
    df = df.dropna()

    if len(df) < 3:
        raise SystemExit(
            f"Need at least 3 usable expression bins, found {len(df)}. "
            "Check that the alpha table was generated correctly."
        )

    x = df["bin_mid"].to_numpy(dtype=float)
    y = df["alpha_mle"].to_numpy(dtype=float)
    w = df["n_genes"].to_numpy(dtype=float)

    w_sum = w.sum()
    if w_sum <= 0:
        raise SystemExit("All bins have zero weight (n_genes)")

    mx = np.sum(w * x) / w_sum
    my = np.sum(w * y) / w_sum

    cov_xy = np.sum(w * (x - mx) * (y - my)) / w_sum
    var_x = np.sum(w * (x - mx) ** 2) / w_sum
    var_y = np.sum(w * (y - my) ** 2) / w_sum

    if var_x <= 0 or var_y <= 0:
        raise SystemExit(
            "Zero variance in expression or alpha; correlation is undefined"
        )

    r = cov_xy / np.sqrt(var_x * var_y)
    r = float(np.clip(r, -1.0, 1.0))
    r2 = r ** 2

    # Effective number of bins (Kish), accounting for unequal weights
    n_eff = (w_sum ** 2) / np.sum(w ** 2)

    if n_eff <= 2:
        p = float("nan")
        note = "  (too few effective bins for a p-value)"
    elif r2 >= 1.0:
        p = 0.0
        note = "  (perfect correlation)"
    else:
        t_stat = r * np.sqrt((n_eff - 2) / (1 - r2))
        p = float(2 * student_t.sf(abs(t_stat), df=n_eff - 2))
        note = ""

    with open(args.output, "w") as out:
        out.write(f"weighted_R_squared\t{r2:.6f}\n")
        out.write(f"p_value\t{p:.6e}\n")
        out.write(f"effective_bins\t{n_eff:.2f}\n")

    print(f"Weighted R squared = {r2:.6f}")
    print(f"p-value = {p:.6e}{note}")


if __name__ == "__main__":
    main()
