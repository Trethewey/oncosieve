#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — annotate_primateai.py
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
#
# Annotates the pan-cancer whitelist with PrimateAI-3D pathogenicity scores.
# PrimateAI-3D provides deep-learning-based pathogenicity predictions for
# all possible missense variants in the human genome (70.6M variants).
#
# Adds three columns to the whitelist:
#   - primateai_score       : Raw pathogenicity score (0.0-1.0)
#   - primateai_percentile  : Genome-wide percentile ranking
#   - primateai_prediction  : Binary classification (benign/pathogenic)
#
# Usage:
#   # Standalone:
#   python3 tools/annotate_primateai.py \
#       --whitelist output/pan_cancer_whitelist_GRCh38_full.tsv.gz \
#       --primateai data/PRIMATE_AI/PrimateAI-3D.hg38.txt.gz \
#       --output output/pan_cancer_whitelist_GRCh38_full_pai3d.tsv.gz
#
#   # Integrated into post_pipeline.py (see integration notes at bottom)
# =============================================================================

import argparse
import os
import sys

import pandas as pd

# Try Polars first (much faster for 70M row file); fall back to pandas
try:
    import polars as pl
    HAS_POLARS = True
except ImportError:
    HAS_POLARS = False


# =============================================================================
# PrimateAI-3D column definitions
# =============================================================================

# Columns in PrimateAI-3D.hg38.txt.gz (no header row in file)
PAI3D_COLUMNS = [
    'chr', 'pos', 'non_flipped_ref', 'non_flipped_alt',
    'gene_name', 'change_position_1based', 'ref_aa', 'alt_aa',
    'score_PAI3D', 'percentile_PAI3D', 'refseq', 'prediction',
]

# Columns we add to the whitelist
OUTPUT_COLS = ['primateai_score', 'primateai_percentile', 'primateai_prediction']


# =============================================================================
# Polars-based annotation (recommended: handles 70M rows efficiently)
# =============================================================================

def annotate_primateai_polars(
    wl: pd.DataFrame,
    primateai_path: str,
) -> pd.DataFrame:
    """
    Annotate whitelist DataFrame with PrimateAI-3D scores using Polars
    lazy scanning for memory efficiency.

    Parameters
    ----------
    wl : pd.DataFrame
        Whitelist with columns: chrom, pos, ref, alt (at minimum).
    primateai_path : str
        Path to PrimateAI-3D.hg38.txt.gz.

    Returns
    -------
    pd.DataFrame
        Whitelist with primateai_score, primateai_percentile,
        primateai_prediction columns added.
    """
    print(f'  Loading PrimateAI-3D: {primateai_path}')
    print('  This may take 2-4 minutes for 70M variants...')

    # Lazy scan -- Polars reads only the columns we need and streams the data
    pai = pl.scan_csv(
        primateai_path,
        separator='\t',
        has_header=True,
        schema_overrides={
            'chr':              pl.Utf8,
            'pos':              pl.Int64,
            'non_flipped_ref':  pl.Utf8,
            'non_flipped_alt':  pl.Utf8,
            'score_PAI3D':      pl.Float64,
            'percentile_PAI3D': pl.Float64,
            'prediction':       pl.Utf8,
        },
    ).select([
        'chr', 'pos', 'non_flipped_ref', 'non_flipped_alt',
        'score_PAI3D', 'percentile_PAI3D', 'prediction',
    ]).filter(
        pl.col('pos').is_not_null()
    )

    # Strip 'chr' prefix from PrimateAI-3D chromosomes to match whitelist
    pai = pai.with_columns(
        pl.col('chr').str.replace(r'^chr', '').alias('chr')
    )

    # Deduplicate: take max score per (chr, pos, ref, alt) if duplicates exist
    pai = pai.group_by(['chr', 'pos', 'non_flipped_ref', 'non_flipped_alt']).agg([
        pl.col('score_PAI3D').max().alias('primateai_score'),
        pl.col('percentile_PAI3D').max().alias('primateai_percentile'),
        pl.col('prediction').first().alias('primateai_prediction'),
    ])

    # Collect into memory
    pai_df = pai.collect().to_pandas()
    # Rename join keys to avoid collision with whitelist columns
    pai_df = pai_df.rename(columns={
        'chr': '_pai_chr', 'pos': '_pai_pos',
        'non_flipped_ref': '_pai_ref', 'non_flipped_alt': '_pai_alt',
    })
    print(f'  {len(pai_df):,} unique PrimateAI-3D entries loaded')

    # Prepare whitelist join keys (strip chr prefix, ensure numeric pos)
    wl['_pai_chrom'] = wl['chrom'].astype(str).str.replace(r'^chr', '', regex=True)
    wl['_pai_pos_key'] = pd.to_numeric(wl['pos'], errors='coerce')

    # Left merge
    merged = wl.merge(
        pai_df,
        left_on=['_pai_chrom', '_pai_pos_key', 'ref', 'alt'],
        right_on=['_pai_chr', '_pai_pos', '_pai_ref', '_pai_alt'],
        how='left',
    ).drop(columns=[
        '_pai_chrom', '_pai_pos_key', '_pai_chr', '_pai_pos',
        '_pai_ref', '_pai_alt',
    ], errors='ignore')

    # Format scores
    merged['primateai_score'] = merged['primateai_score'].round(4).astype(str)
    merged.loc[merged['primateai_score'] == 'nan', 'primateai_score'] = ''

    merged['primateai_percentile'] = merged['primateai_percentile'].round(4).astype(str)
    merged.loc[merged['primateai_percentile'] == 'nan', 'primateai_percentile'] = ''

    merged['primateai_prediction'] = merged['primateai_prediction'].fillna('')

    n_scored = (merged['primateai_score'] != '').sum()
    print(f'  Variants with PrimateAI-3D score: {n_scored:,} / {len(merged):,}')

    return merged


# =============================================================================
# Pandas-based fallback (slower, higher memory, but works without Polars)
# =============================================================================

def annotate_primateai_pandas(
    wl: pd.DataFrame,
    primateai_path: str,
) -> pd.DataFrame:
    """
    Pandas fallback for PrimateAI-3D annotation.
    WARNING: This loads the full 1.7GB gzipped file into memory (~6GB uncompressed).
    Requires 16+ GB RAM.
    """
    print(f'  Loading PrimateAI-3D (pandas fallback): {primateai_path}')
    print('  WARNING: This requires 16+ GB RAM. Consider installing Polars.')

    pai = pd.read_csv(
        primateai_path,
        sep='\t',
        usecols=['chr', 'pos', 'non_flipped_ref', 'non_flipped_alt',
                 'score_PAI3D', 'percentile_PAI3D', 'prediction'],
        dtype={
            'chr': str, 'pos': 'Int64',
            'non_flipped_ref': str, 'non_flipped_alt': str,
            'score_PAI3D': float, 'percentile_PAI3D': float,
            'prediction': str,
        },
        compression='gzip',
    )

    pai['chr'] = pai['chr'].str.replace(r'^chr', '', regex=True)
    pai = pai.groupby(['chr', 'pos', 'non_flipped_ref', 'non_flipped_alt']).agg({
        'score_PAI3D': 'max',
        'percentile_PAI3D': 'max',
        'prediction': 'first',
    }).reset_index()

    pai = pai.rename(columns={
        'score_PAI3D': 'primateai_score',
        'percentile_PAI3D': 'primateai_percentile',
        'prediction': 'primateai_prediction',
    })

    wl['_pai_chrom'] = wl['chrom'].astype(str).str.replace(r'^chr', '', regex=True)
    wl['_pai_pos'] = pd.to_numeric(wl['pos'], errors='coerce')

    merged = wl.merge(
        pai,
        left_on=['_pai_chrom', '_pai_pos', 'ref', 'alt'],
        right_on=['chr', 'pos', 'non_flipped_ref', 'non_flipped_alt'],
        how='left',
    ).drop(columns=[
        '_pai_chrom', '_pai_pos', 'chr', 'pos',
        'non_flipped_ref', 'non_flipped_alt',
    ], errors='ignore')

    merged['primateai_score'] = merged['primateai_score'].round(4).astype(str)
    merged.loc[merged['primateai_score'] == 'nan', 'primateai_score'] = ''
    merged['primateai_percentile'] = merged['primateai_percentile'].round(4).astype(str)
    merged.loc[merged['primateai_percentile'] == 'nan', 'primateai_percentile'] = ''
    merged['primateai_prediction'] = merged['primateai_prediction'].fillna('')

    n_scored = (merged['primateai_score'] != '').sum()
    print(f'  Variants with PrimateAI-3D score: {n_scored:,} / {len(merged):,}')
    return merged


# =============================================================================
# Public API
# =============================================================================

def annotate_primateai(
    wl: pd.DataFrame,
    primateai_path: str,
) -> pd.DataFrame:
    """
    Annotate whitelist with PrimateAI-3D scores.
    Automatically selects Polars (preferred) or pandas backend.
    """
    if not os.path.isfile(primateai_path):
        print(f'  PrimateAI-3D file not found: {primateai_path}')
        print('  Skipping PrimateAI-3D annotation.')
        for col in OUTPUT_COLS:
            wl[col] = ''
        return wl

    if HAS_POLARS:
        return annotate_primateai_polars(wl, primateai_path)
    else:
        return annotate_primateai_pandas(wl, primateai_path)


# =============================================================================
# Standalone CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Annotate ONCOSIEVE whitelist with PrimateAI-3D scores',
    )
    parser.add_argument(
        '--whitelist', required=True,
        help='Input whitelist TSV (.tsv.gz)',
    )
    parser.add_argument(
        '--primateai', required=True,
        help='Path to PrimateAI-3D.hg38.txt.gz',
    )
    parser.add_argument(
        '--output', required=True,
        help='Output annotated TSV (.tsv.gz)',
    )
    args = parser.parse_args()

    print('ONCOSIEVE -- PrimateAI-3D Annotation')
    print(f'  Input:     {args.whitelist}')
    print(f'  PrimateAI: {args.primateai}')
    print(f'  Output:    {args.output}')
    print()

    wl = pd.read_csv(args.whitelist, sep='\t', dtype=str, low_memory=False)
    print(f'  Loaded {len(wl):,} whitelist variants')

    wl = annotate_primateai(wl, args.primateai)

    compression = 'gzip' if args.output.endswith('.gz') else None
    wl.to_csv(args.output, sep='\t', index=False, compression=compression)
    print(f'  Written: {args.output}')
    print('  Done.')


if __name__ == '__main__':
    main()


# =============================================================================
# INTEGRATION NOTES
# =============================================================================
#
# To integrate PrimateAI-3D into the existing post_pipeline.py workflow:
#
# 1. Add to config.yaml:
#
#    primateai:
#      enabled: true
#      scores: data/PRIMATE_AI/PrimateAI-3D.hg38.txt.gz
#
# 2. In tools/post_pipeline.py, after the REVEL annotation step, add:
#
#    from tools.annotate_primateai import annotate_primateai
#
#    # After: wl = annotate_revel(wl, args.revel)
#    pai_path = 'data/PRIMATE_AI/PrimateAI-3D.hg38.txt.gz'
#    if os.path.isfile(pai_path):
#        print('\n--- PrimateAI-3D Annotation ---')
#        wl = annotate_primateai(wl, pai_path)
#
# 3. In tools/generate_report.py, add 'primateai_score' to DISPLAY_COLS
#    and include it in the KPI card section.
#
# 4. In run_oncosieve.sh, add a --primateai flag:
#
#    PRIMATEAI_FILE="${DATA_DIR:-data}/PRIMATE_AI/PrimateAI-3D.hg38.txt.gz"
#    if [[ -f "$PRIMATEAI_FILE" ]]; then
#        "$PYTHON" tools/annotate_primateai.py \
#            --whitelist "$WL_FULL" \
#            --primateai "$PRIMATEAI_FILE" \
#            --output "${WL_FULL%.tsv.gz}_pai3d.tsv.gz"
#    fi
#
# 5. Tier promotion (v2.1: implemented in post_pipeline.py after the
#    annotate_primateai() call — a Tier 3 variant with
#    primateai_prediction == 'pathogenic' and primateai_percentile >= 0.9
#    is rescued to Tier 2. Kept out of build_whitelist.py because tiering
#    already runs before PrimateAI annotation.
#
# =============================================================================
