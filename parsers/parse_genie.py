#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Parser: AACR GENIE MAF
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
parse_genie.py
Parse AACR Project GENIE data_mutations_extended.txt (MAF format, GRCh38).

Required file:
  data_mutations_extended.txt
  Download: https://www.synapse.org/#!Synapse:syn7222066
  Requires a free Synapse account and GENIE data access agreement.
  Select the latest public release (e.g. GENIE 15.0-public).

The MAF file follows the GDC MAF format specification. Key columns used:
  Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2,
  Hugo_Symbol, HGVSc, HGVSp_Short, Variant_Classification,
  Tumor_Sample_Barcode, Center (maps to cancer centre / panel)

GENIE does not annotate cancer type per mutation row directly. Cancer type
is retrieved from the companion clinical file: data_clinical_sample.txt
using SAMPLE_ID -> ONCOTREE_CODE.
"""

import os

import pandas as pd

from parsers.common import (
    CONSEQUENCE_MAP,
    STANDARD_COLS,
    clean_allele,
    empty_standard_df,
    is_valid_allele,
    map_consequence,
    normalise_chrom,
    setup_logger,
)

log = setup_logger('GENIE')

_MAF_COLS = {
    'chrom':       'Chromosome',
    'pos':         'Start_Position',
    'ref':         'Reference_Allele',
    'alt':         'Tumor_Seq_Allele2',
    'gene':        'Hugo_Symbol',
    'hgvsc':       'HGVSc',
    'hgvsp':       'HGVSp_Short',
    'consequence': 'Variant_Classification',
    'sample_id':   'Tumor_Sample_Barcode',
}

# Optional MAF columns — used if present, ignored if absent so the parser
# still works on older MAF releases that lack them.
_MAF_OPT_COL_CSQ    = 'Consequence'      # VEP term (multi-term '&'-joined); richer than Variant_Classification
_MAF_OPT_COL_TID    = 'Transcript_ID'    # ENST — prepended to bare c. HGVSc so MANE Select can resolve
_MAF_OPT_COL_FILTER = 'FILTER'           # only PASS / . / empty is kept

_CLINICAL_SAMPLE_ID_COL    = 'SAMPLE_ID'
_CLINICAL_ONCOTREE_COL     = 'ONCOTREE_CODE'
_CLINICAL_CANCER_TYPE_COL  = 'CANCER_TYPE'


def parse_genie(maf_path: str,
                clinical_sample_path: str | None = None) -> pd.DataFrame:
    """
    Parse GENIE MAF file.

    Parameters
    ----------
    maf_path              : path to data_mutations_extended.txt
    clinical_sample_path  : path to data_clinical_sample.txt (optional;
                            enables cancer type annotation per sample)

    Returns
    -------
    DataFrame with STANDARD_COLS schema.
    One row per mutation per sample. Aggregation to unique samples per variant
    is performed at the merge step.
    """
    if not os.path.exists(maf_path):
        log.warning('GENIE MAF not found: %s  — skipping', maf_path)
        return empty_standard_df()

    log.info('Parsing GENIE MAF: %s', maf_path)

    # Build cancer type lookup from clinical file if available
    cancer_type_map: dict[str, str] = {}
    if clinical_sample_path and os.path.exists(clinical_sample_path):
        cancer_type_map = _load_clinical_sample(clinical_sample_path)
        log.info('GENIE clinical lookup: %d samples', len(cancer_type_map))

    # Read MAF; skip metadata lines starting with '#'
    try:
        df_maf = pd.read_csv(
            maf_path, sep='\t', comment='#', dtype=str,
            low_memory=False, na_filter=False,
        )
    except Exception as e:
        log.error('Failed to read GENIE MAF: %s', e)
        return empty_standard_df()

    # Validate required columns
    required = list(_MAF_COLS.values())
    missing  = [c for c in required if c not in df_maf.columns]
    if missing:
        log.error('GENIE MAF missing columns: %s', missing)
        return empty_standard_df()

    # FILTER: only PASS calls contribute to the whitelist. Mirrors parse_tcga.
    if _MAF_OPT_COL_FILTER in df_maf.columns:
        n_before = len(df_maf)
        flt = df_maf[_MAF_OPT_COL_FILTER].astype(str).str.strip().str.upper()
        df_maf = df_maf[flt.isin({'PASS', '.', ''})].copy()
        n_after = len(df_maf)
        if n_before != n_after:
            log.info('GENIE FILTER: %d rows excluded (%d kept)',
                     n_before - n_after, n_after)
    else:
        log.warning('GENIE MAF has no FILTER column — assuming all calls PASS')

    # Vectorized processing — avoids iterrows on ~3.4M rows
    working_cols = [
        _MAF_COLS['chrom'], _MAF_COLS['pos'], _MAF_COLS['ref'], _MAF_COLS['alt'],
        _MAF_COLS['gene'], _MAF_COLS['hgvsc'], _MAF_COLS['hgvsp'],
        _MAF_COLS['consequence'], _MAF_COLS['sample_id'],
    ]
    if _MAF_OPT_COL_CSQ in df_maf.columns:
        working_cols.append(_MAF_OPT_COL_CSQ)
    if _MAF_OPT_COL_TID in df_maf.columns:
        working_cols.append(_MAF_OPT_COL_TID)

    df = df_maf[working_cols].copy()
    new_names = ['chrom', 'pos', 'ref', 'alt', 'gene', 'hgvsc', 'hgvsp',
                 'variant_classification', 'sample_id']
    if _MAF_OPT_COL_CSQ in df_maf.columns:
        new_names.append('vep_consequence')
    if _MAF_OPT_COL_TID in df_maf.columns:
        new_names.append('transcript_id')
    df.columns = new_names

    # Clean alleles
    df['ref'] = df['ref'].str.strip().str.upper()
    df['alt'] = df['alt'].str.strip().str.upper()

    # Filter invalid alleles
    valid_re = r'^[ACGTNacgtn*\-]+$'
    mask = (
        df['ref'].str.match(valid_re, na=False) &
        df['alt'].str.match(valid_re, na=False) &
        ~df['alt'].isin(['', '-', '.'])
    )
    df = df[mask].copy()

    # Normalise chromosome: strip 'chr' prefix then re-add, handle MT
    raw_chrom = df['chrom'].str.strip().str.replace(r'^chr', '', regex=True)
    raw_chrom = raw_chrom.replace({'MT': 'M', 'mt': 'M'})
    df['chrom'] = 'chr' + raw_chrom

    # Position
    df['pos'] = pd.to_numeric(df['pos'], errors='coerce')
    df = df.dropna(subset=['pos'])
    df['pos'] = df['pos'].astype(int)

    # Strip string columns
    for col in ('gene', 'hgvsc', 'hgvsp'):
        df[col] = df[col].str.strip()

    # Prepend Transcript_ID to bare c. HGVSc so the MANE Select lookup
    # downstream can find an ENST accession. Without this, GENIE (transcript
    # priority 1) contributes bare c.-strings that never resolve to a MANE
    # transcript and is_mane_select stays False for most whitelist rows.
    if 'transcript_id' in df.columns:
        tid = df['transcript_id'].astype(str).str.strip()
        bare_c = df['hgvsc'].str.startswith('c.', na=False) & tid.ne('') & tid.ne('nan')
        df.loc[bare_c, 'hgvsc'] = tid[bare_c] + ':' + df.loc[bare_c, 'hgvsc']
        df = df.drop(columns=['transcript_id'])

    # Consequence mapping. Prefer VEP Consequence (multi-term '&'-joined) via
    # map_consequence, fall back to MAF Variant_Classification via CONSEQUENCE_MAP.
    if 'vep_consequence' in df.columns:
        df['consequence'] = df['vep_consequence'].fillna('').astype(str).str.strip().apply(map_consequence)
        # Backfill any 'unknown' from the coarser MAF column
        mask = df['consequence'].eq('unknown')
        if mask.any():
            df.loc[mask, 'consequence'] = (
                df.loc[mask, 'variant_classification']
                  .str.strip().map(CONSEQUENCE_MAP).fillna('unknown')
            )
        df = df.drop(columns=['vep_consequence'])
    else:
        df['consequence'] = df['variant_classification'].str.strip().map(CONSEQUENCE_MAP).fillna('unknown')
    df = df.drop(columns=['variant_classification'])

    # Cancer type from clinical lookup
    df['sample_id'] = df['sample_id'].str.strip()
    if cancer_type_map:
        df['cancer_type'] = df['sample_id'].map(cancer_type_map).fillna('unspecified')
    else:
        df['cancer_type'] = 'unspecified'

    df['n_samples'] = 1
    df['source'] = 'GENIE'

    # Drop sample_id and align to STANDARD_COLS
    df = df.drop(columns=['sample_id'])
    for col in STANDARD_COLS:
        if col not in df.columns:
            df[col] = ''
    df = df[STANDARD_COLS]

    if df.empty:
        log.warning('GENIE: no rows produced')
        return empty_standard_df()

    log.info('GENIE: %d rows', len(df))
    return df


def _load_clinical_sample(clinical_path: str) -> dict[str, str]:
    """
    Load GENIE data_clinical_sample.txt and return a dict of
    SAMPLE_ID -> cancer type label (CANCER_TYPE preferred, else ONCOTREE_CODE).
    """
    try:
        df = pd.read_csv(
            clinical_path, sep='\t', comment='#',
            dtype=str, low_memory=False, na_filter=False,
        )
    except Exception as e:
        log.warning('Could not load clinical sample file: %s', e)
        return {}

    if _CLINICAL_SAMPLE_ID_COL not in df.columns:
        log.warning('Clinical file missing %s column', _CLINICAL_SAMPLE_ID_COL)
        return {}

    # Prefer CANCER_TYPE, fall back to ONCOTREE_CODE
    if _CLINICAL_CANCER_TYPE_COL in df.columns:
        ctype_col = _CLINICAL_CANCER_TYPE_COL
    elif _CLINICAL_ONCOTREE_COL in df.columns:
        ctype_col = _CLINICAL_ONCOTREE_COL
    else:
        log.warning('Clinical file has no cancer type column')
        return {}

    df = df[[_CLINICAL_SAMPLE_ID_COL, ctype_col]].copy()
    df[_CLINICAL_SAMPLE_ID_COL] = df[_CLINICAL_SAMPLE_ID_COL].str.strip()
    df[ctype_col] = df[ctype_col].str.strip().str.lower()
    df = df[df[_CLINICAL_SAMPLE_ID_COL] != '']
    return df.set_index(_CLINICAL_SAMPLE_ID_COL)[ctype_col].to_dict()
