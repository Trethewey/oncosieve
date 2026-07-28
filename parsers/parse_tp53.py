#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Parser: TP53 somatic variants database
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
parse_tp53.py
Parse the NCI TP53 database (GRCh38) somatic and germline variant files.

Required files:
  SomaticVariants_GRCh38.csv   (primary)
  GermlineVariants_GRCh38.csv  (optional, if include_germline=True)
  reference_fasta               (optional but strongly recommended) — used to
                                 build proper VCF-style anchor-base indels.
                                 Without it, TP53 indels are DROPPED rather
                                 than fabricated with placeholder anchors.

Download: https://tp53.cancer.gov

Key column usage:
  - Position  : g_description_GRCh38 (e.g. g.7675184A>G)
  - REF/ALT   : parsed from g_description_GRCh38 (SNVs directly; indels via
                the reference FASTA when available, else dropped)
  - Consequence: Effect (case-insensitive)
  - Cancer type: Morphology (histological diagnosis)
  - HGVSc     : c_description
  - HGVSp     : ProtDescription
  - Count     : TCGA_ICGC_GENIE_count if available, else 1 per row
"""

import os
import re

import pandas as pd

from parsers.common import (
    STANDARD_COLS,
    clean_allele,
    empty_standard_df,
    is_valid_allele,
    map_consequence,
    normalise_chrom,
    setup_logger,
)

log = setup_logger('TP53')

# Reference-FASTA singleton. Set by parse_tp53() when a fasta path is provided
# so _parse_gdesc can look up the anchor base for indel normalisation. A None
# value means indels get dropped rather than fabricated with placeholder bases.
_FASTA = None


def _open_fasta(fasta_path: str):
    """Return a pyfaidx.Fasta handle or None."""
    if not fasta_path or not os.path.exists(fasta_path):
        return None
    try:
        import pyfaidx
    except ImportError:
        log.warning('TP53: pyfaidx not installed — indels will be dropped. '
                    'pip install pyfaidx')
        return None
    try:
        return pyfaidx.Fasta(fasta_path)
    except Exception as e:
        log.warning('TP53: could not open reference FASTA %s: %s — indels will be dropped',
                    fasta_path, e)
        return None


def _fasta_base(chrom: str, pos: int) -> str | None:
    """Return the 1-based base at (chrom, pos) or None."""
    if _FASTA is None:
        return None
    for key in (chrom, chrom.lstrip('chr'), f'chr{chrom.lstrip("chr")}'):
        try:
            return str(_FASTA[key][pos - 1]).upper()
        except (KeyError, IndexError):
            continue
    return None

# Regex to parse g_description_GRCh38
# SNV:    g.7674220C>T
_GDESC_SNV_RE = re.compile(r'g\.(\d+)([ACGTacgt])>([ACGTacgt])')
# Deletion with bases: g.7674220_7674225delACGTAC  or  g.7674220delC
_GDESC_DEL_RE = re.compile(r'g\.(\d+)(?:_\d+)?del([ACGTacgt]+)', re.IGNORECASE)
# Insertion: g.7674220_7674221insACG
_GDESC_INS_RE = re.compile(r'g\.(\d+)_\d+ins([ACGTacgt]+)', re.IGNORECASE)
# Delins with bases: g.7674220_7674225delACGTACinsT  or  g.7674220delCinsT
_GDESC_DELINS_RE = re.compile(r'g\.(\d+)(?:_\d+)?del([ACGTacgt]+)ins([ACGTacgt]+)', re.IGNORECASE)
# Duplication: g.7674220_7674222dup  or  g.7674220_7674222dupACG
_GDESC_DUP_RE = re.compile(r'g\.(\d+)(?:_\d+)?dup([ACGTacgt]*)', re.IGNORECASE)


def _parse_gdesc(gdesc: str):
    """
    Parse g_description_GRCh38 into (pos, ref, alt) in VCF-style coordinates.
    Returns (pos, ref, alt) or None if unparseable.

    SNVs are resolved directly. Indels (deletion, insertion, delins,
    duplication) require an anchor base from the reference FASTA — when
    the module-level `_FASTA` singleton is None (no reference configured
    or pyfaidx missing), indels are dropped rather than fabricated with
    placeholder bases. The v2 fabrication produced REF==ALT deletions
    and REF='N' insertions that clinicians could not trust.
    """
    gdesc = gdesc.strip()

    # SNV: g.NNNNref>alt  — no FASTA lookup needed.
    m = _GDESC_SNV_RE.search(gdesc)
    if m:
        return int(m.group(1)), m.group(2).upper(), m.group(3).upper()

    # Try delins first (before plain del, since delins contains 'del').
    # delins is representable directly (both alleles given), but VCF convention
    # is that ref and alt should start with the same anchor base. Some
    # downstream tools also accept the direct form; we emit it as-given.
    m = _GDESC_DELINS_RE.search(gdesc)
    if m:
        pos = int(m.group(1))
        ref = m.group(2).upper()
        alt = m.group(3).upper()
        return pos, ref, alt

    # For everything below, we need a reference base at pos-1 (VCF anchor).
    # If FASTA is not available, drop the row.
    if _FASTA is None:
        return None

    # Deletion with explicit bases: g.NNNNdelBASES
    # VCF: pos=pos-1, ref=anchor+deleted, alt=anchor
    m = _GDESC_DEL_RE.search(gdesc)
    if m:
        pos = int(m.group(1))
        deleted = m.group(2).upper()
        anchor = _fasta_base('chr17', pos - 1)
        if anchor is None:
            return None
        return pos - 1, anchor + deleted, anchor

    # Insertion: g.NNNN_NNNN+1insBASES
    # VCF: pos=pos, ref=anchor(pos), alt=anchor+inserted
    m = _GDESC_INS_RE.search(gdesc)
    if m:
        pos = int(m.group(1))
        inserted = m.group(2).upper()
        anchor = _fasta_base('chr17', pos)
        if anchor is None:
            return None
        return pos, anchor, anchor + inserted

    # Duplication with explicit bases: g.NNNNdupBASES
    # A dup at [pos..pos+k-1] is equivalent to an insertion of the same bases
    # immediately after pos+k-1. VCF: pos=pos+k-1, ref=anchor(pos+k-1),
    # alt=anchor+duped.
    m = _GDESC_DUP_RE.search(gdesc)
    if m and m.group(2):
        pos = int(m.group(1))
        duped = m.group(2).upper()
        last_dup_pos = pos + len(duped) - 1
        anchor = _fasta_base('chr17', last_dup_pos)
        if anchor is None:
            return None
        return last_dup_pos, anchor, anchor + duped

    return None

# Consequence mapping — case-insensitive lookup applied at parse time.
# 'deletion' and 'insertion' resolve at call-site by frame-checking (len(ref) - len(alt)) % 3,
# not by unconditionally labelling every indel as frameshift.
_TP53_CONSEQUENCE_MAP = {
    'missense':   'missense',
    'nonsense':   'nonsense',
    'silent':     'synonymous',
    'splice':     'splice_site',
    'frameshift': 'frameshift',
    'in-frame':   'inframe_indel',
    'in frame':   'inframe_indel',
    'complex':    'frameshift',
    'intronic':   'intronic',
    'other':      'unknown',
    'unknown':    'unknown',
    '':           'unknown',
}


def _indel_consequence(ref: str, alt: str) -> str:
    """Return 'frameshift' or 'inframe_indel' based on length difference."""
    delta = len(ref) - len(alt)
    return 'inframe_indel' if delta % 3 == 0 else 'frameshift'

# DNE_LOFclass — retained as metadata but no longer used as a filter.
# The TP53 database is comprehensive; gating on functional class dropped 45% of
# variants and caused TP53-gene entries from other sources to lack TP53 attribution.


def parse_tp53(somatic_tsv: str,
               germline_tsv: str | None = None,
               include_germline: bool = False,
               reference_fasta: str | None = None) -> pd.DataFrame:
    """
    Parse TP53 somatic (and optionally germline) variant files.

    If `reference_fasta` is provided and pyfaidx is available, indels are
    normalised with a real anchor base from the FASTA. Otherwise indels
    are dropped and only SNVs and delins entries reach the output.
    """
    global _FASTA
    _FASTA = _open_fasta(reference_fasta) if reference_fasta else None
    if reference_fasta and _FASTA is None:
        log.warning('TP53: reference FASTA unavailable — indels will be dropped')
    elif _FASTA is None:
        log.info('TP53: no reference_fasta configured — indels will be dropped '
                 '(SNVs and delins still processed)')
    else:
        log.info('TP53: reference FASTA opened for indel anchor normalisation')

    frames = []

    if os.path.exists(somatic_tsv):
        log.info('Parsing TP53 somatic file: %s', somatic_tsv)
        df = _parse_tp53_tsv(somatic_tsv, source_label='TP53_somatic')
        if df is not None and not df.empty:
            frames.append(df)
    else:
        log.warning('TP53 somatic file not found: %s', somatic_tsv)

    if include_germline and germline_tsv and os.path.exists(germline_tsv):
        log.info('Parsing TP53 germline file: %s', germline_tsv)
        df = _parse_tp53_tsv(germline_tsv, source_label='TP53_germline')
        if df is not None and not df.empty:
            frames.append(df)

    if not frames:
        return empty_standard_df()

    result = pd.concat(frames, ignore_index=True)
    log.info('TP53: %d total rows', len(result))
    return result


def _parse_tp53_tsv(path: str, source_label: str) -> pd.DataFrame | None:
    sep = ',' if path.endswith('.csv') else '\t'
    try:
        df_raw = pd.read_csv(path, sep=sep, dtype=str, low_memory=False,
                             na_filter=False)
        if len(df_raw.columns) == 1:
            sep = ',' if sep == '\t' else '\t'
            df_raw = pd.read_csv(path, sep=sep, dtype=str, low_memory=False,
                                 na_filter=False)
    except Exception as e:
        log.error('Failed to read TP53 file %s: %s', path, e)
        return None

    log.info('TP53 columns found: %s', list(df_raw.columns))
    df_raw = df_raw.fillna('')

    # Locate required columns
    pos_col      = _find_col(df_raw, ['g_description_GRCh38'])
    effect_col   = _find_col(df_raw, ['Effect', 'Mutation_type', 'mutation_type'])
    dnelof_col   = _find_col(df_raw, ['DNE_LOFclass'])
    morpho_col   = _find_col(df_raw, ['Morphology', 'Short_topo', 'Topography'])
    count_col    = _find_col(df_raw, ['TCGA_ICGC_GENIE_count', 'Somatic_count', 'Count'])
    hgvsp_col    = _find_col(df_raw, ['ProtDescription', 'HGVSp', 'AAchange'])
    hgvsc_col    = _find_col(df_raw, ['c_description', 'HGVSc'])

    if not pos_col:
        log.error('TP53: g_description_GRCh38 column not found — cannot parse positions')
        return None

    n_skipped_pos    = 0
    rows = []

    for _, row in df_raw.iterrows():
        try:
            # --- Position and alleles from g_description_GRCh38 ---
            gdesc = str(row[pos_col]).strip()
            parsed = _parse_gdesc(gdesc)
            if parsed is None:
                n_skipped_pos += 1
                continue

            pos, ref, alt = parsed

            if not is_valid_allele(ref) or not is_valid_allele(alt):
                continue

            # TP53 is always chr17
            chrom = 'chr17'

            # --- Consequence (case-insensitive) ---
            effect_raw  = str(row[effect_col]).strip().lower() if effect_col else ''
            if effect_raw in ('deletion', 'insertion', 'duplication', 'indel'):
                # Frame-derive from allele-length delta rather than assuming
                # every deletion or insertion is a frameshift.
                consequence = _indel_consequence(ref, alt)
            else:
                consequence = _TP53_CONSEQUENCE_MAP.get(effect_raw,
                              map_consequence(effect_raw))

            # --- Cancer type from Morphology ---
            # Kept lowercase for legacy compatibility with the TP53 pathway,
            # but the aggregator's _clean_cancer_type does canonical case
            # normalisation across all sources.
            cancer_type = str(row[morpho_col]).strip() if morpho_col else ''
            if not cancer_type or cancer_type.lower() in ('', 'na', 'nan'):
                cancer_type = 'unspecified'

            # --- Sample count ---
            n_samples = 1
            if count_col:
                count_val = str(row[count_col]).strip()
                if count_val.isdigit():
                    n_samples = int(count_val)

            # --- HGVSc and HGVSp ---
            hgvsc = str(row[hgvsc_col]).strip() if hgvsc_col else ''
            hgvsp = str(row[hgvsp_col]).strip() if hgvsp_col else ''

            # --- DNE_LOFclass (metadata, not a filter) ---
            dnelof = str(row[dnelof_col]).strip() if dnelof_col else ''

            rows.append({
                'chrom':       chrom,
                'pos':         pos,
                'ref':         ref,
                'alt':         alt,
                'gene':        'TP53',
                'hgvsc':       hgvsc,
                'hgvsp':       hgvsp,
                'consequence': consequence,
                'cancer_type': cancer_type,
                'n_samples':   n_samples,
                'source':      source_label,
                'tp53_class':  dnelof,
            })
        except (ValueError, TypeError, KeyError):
            continue

    log.info(
        'TP53 (%s): %d rows kept | %d excluded by unparseable position',
        source_label, len(rows), n_skipped_pos
    )

    if not rows:
        log.warning('TP53: no rows produced from %s', path)
        return None

    df = pd.DataFrame(rows, columns=STANDARD_COLS + ['tp53_class'])
    return df


def _find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None
