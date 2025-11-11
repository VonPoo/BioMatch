#!/usr/bin/env python3
"""
BioMatch Analysis Script: special_sequence_method_snp_filter

Purpose
- Filter PLINK `.bim` variants by an allowed base set (keep-base), keeping only
  SNPs whose alleles (A1/A2) are composed entirely of specified bases.
- Normalize missing/placeholder BIM IDs in-place to `CHR_BP` format to ensure
  PLINK `--extract` behaves consistently.

Inputs
- `bim`: Path to PLINK `.bim` file to inspect and normalize.
- `--output`: Path to write kept SNP IDs, one per line.
- `--keep-base`: Allowed bases, compact (e.g., `ATC`) or comma-separated
  (e.g., `A,T,C`). Case-insensitive; only `A/C/G/T` are accepted.

Outputs
- Text file containing kept SNP IDs, one per line. The input `.bim` may be
  updated in-place to normalize missing IDs prior to filtering.

Usage Example
    python special_sequence_method_snp_filter.py \
        path/to/combined_data.bim \
        --output path/to/combined_data_filtered_snps.txt \
        --keep-base AT

Notes
- Treats `'.'`, `'0'`, and `'NA'` as disallowed or missing alleles/IDs.
- Multi-letter alleles (e.g., short indels) are evaluated per character and
  must consist entirely of allowed bases to pass.
- Requires Python 3 and `pandas`.
"""
import argparse
import sys
import os
import pandas as pd


def parse_args():
    """Parse CLI arguments for BIM filtering.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with `bim`, `output`, and `keep_base` attributes.
    """
    p = argparse.ArgumentParser(
        description=(
            "Filter PLINK BIM variants by allowed base-set (keep-base). "
            "Keeps only rows whose alleles (A1/A2) consist entirely of allowed bases."
        )
    )
    p.add_argument("bim", help="Path to PLINK .bim file")
    p.add_argument("--output", required=True, help="Output file path for kept SNP IDs (one per line)")
    p.add_argument(
        "--keep-base",
        required=True,
        help=(
            "Allowed bases set. Accepts compact form (e.g., ATC) or comma-separated (e.g., A,T,C). "
            "Case-insensitive; any allele containing disallowed bases is excluded."
        ),
    )
    return p.parse_args()


def read_bim(path: str) -> pd.DataFrame:
    """Read a whitespace-separated PLINK BIM file into a DataFrame.

    Parameters
    ----------
    path : str
        Path to the `.bim` file.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: `CHR`, `SNP`, `CM`, `BP`, `A1`, `A2`.
    """
    # Robust whitespace-separated BIM reading
    # Keep CM as string to avoid formatting changes like 0 -> 0.0
    cols = ["CHR", "SNP", "CM", "BP", "A1", "A2"]
    dtypes = {"CHR": str, "SNP": str, "CM": object, "BP": int, "A1": str, "A2": str}
    return pd.read_csv(path, sep=r"\s+", names=cols, dtype=dtypes, engine="python")


def normalize_bim_ids_inplace(df: pd.DataFrame, bim_path: str) -> pd.DataFrame:
    """Normalize missing SNP IDs in a BIM DataFrame and write back in-place.

    Missing or placeholder IDs (e.g., `.`, `0`, `NA`, empty) are replaced with
    `CHR_BP` format to ensure PLINK ID-based operations work reliably.

    Parameters
    ----------
    df : pandas.DataFrame
        BIM DataFrame.
    bim_path : str
        Path to the source `.bim` file to rewrite.

    Returns
    -------
    pandas.DataFrame
        Updated DataFrame with normalized `SNP` column.
    """
    # Normalize missing/placeholder SNP IDs to CHR_BP (e.g., '1_28042')
    snp = df["SNP"].astype(str)
    # Treat '.', empty, '0', 'NA', and 'nan' as missing IDs
    missing_vals = {".", "", "0", "NA", "nan", "None"}
    mask = snp.str.lower().isin(missing_vals) | df["SNP"].isna()
    # Precompute normalized IDs
    chr_str = df["CHR"].astype(str)
    bp_str = df["BP"].astype(int).astype(str)
    new_ids = chr_str + "_" + bp_str
    replaced_count = int(mask.sum())
    if replaced_count > 0:
        df.loc[mask, "SNP"] = new_ids[mask]
        # Write back to the BIM file to ensure PLINK --extract works with normalized IDs
        try:
            df.to_csv(bim_path, sep="\t", header=False, index=False)
        except Exception as e:
            print(f"Error writing normalized BIM: {e}", file=sys.stderr)
            raise
        print(f"Normalized BIM IDs in-place: replaced={replaced_count}", file=sys.stderr)
    else:
        print("No BIM IDs required normalization", file=sys.stderr)
    return df


def normalize_keep_bases(raw: str) -> set:
    """Normalize `--keep-base` input into a set of allowed bases.

    Parameters
    ----------
    raw : str
        Raw input string (compact or comma-separated).

    Returns
    -------
    set
        Uppercase allowed bases subset of `{A,C,G,T}`.
    """
    s = (raw or "").upper().replace(" ", "").replace(",", "")
    allowed = set(ch for ch in s if ch in {"A", "C", "G", "T"})
    if not allowed:
        raise ValueError("--keep-base must specify at least one of A,C,G,T")
    return allowed


def allele_allowed(allele: str, allowed: set) -> bool:
    """Return True if all characters in the allele are allowed bases.

    Treats `None`, `'.'`, `'0'`, and `'NA'` as disallowed.
    """
    if allele is None:
        return False
    a = allele.upper()
    # Treat '.' or missing as disallowed
    if a == "." or a == "0" or a == "NA":
        return False
    # Each character must be in allowed set; multi-letter alleles (e.g., indels) are handled char-wise
    return all(ch in allowed for ch in a)


def main():
    """Entry point: normalize BIM IDs and filter SNPs by allowed bases."""
    args = parse_args()
    bim_path = args.bim
    out_path = args.output
    try:
        allowed = normalize_keep_bases(args.keep_base)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    if not os.path.exists(bim_path):
        print(f"Error: BIM file not found: {bim_path}", file=sys.stderr)
        return 1

    try:
        df = read_bim(bim_path)
    except Exception as e:
        print(f"Error reading BIM: {e}", file=sys.stderr)
        return 1

    # Normalize BIM IDs in-place before filtering so that PLINK --extract matches
    try:
        df = normalize_bim_ids_inplace(df, bim_path)
    except Exception:
        return 1

    # Apply filtering mask
    mask = df["A1"].apply(lambda x: allele_allowed(x, allowed)) & df["A2"].apply(lambda x: allele_allowed(x, allowed))
    kept = df.loc[mask, "SNP"].dropna()

    # Write output as one ID per line
    try:
        with open(out_path, "w") as f:
            for snp in kept:
                f.write(str(snp) + "\n")
    except Exception as e:
        print(f"Error writing output: {e}", file=sys.stderr)
        return 1

    # Summary to stderr to avoid mixing with output
    total = len(df)
    n_kept = len(kept)
    n_removed = total - n_kept
    print(
        f"Filtered BIM: total={total}, kept={n_kept}, removed={n_removed}; keep-bases={''.join(sorted(allowed))}",
        file=sys.stderr,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
