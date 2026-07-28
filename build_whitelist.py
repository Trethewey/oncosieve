#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Main pipeline: aggregate multi-source somatic variant data and build the whitelist
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
build_whitelist.py
Main pipeline script. Runs all parsers, merges outputs, applies filters,
assigns tiers, and writes the final whitelist in TSV and VCF formats.

Usage:
    python build_whitelist.py --config config.yaml [--skip-sources cosmic,genie]
"""

import argparse
import logging
import os
import sys
from datetime import datetime

import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import yaml

from parsers.common import INCLUDED_CONSEQUENCES, STANDARD_COLS, setup_logger

log = setup_logger('build_whitelist')


# ── MANE Select lookup ──────────────────────────────────────────────────────

def load_mane_lookup(mane_path: str = 'data/reference/mane_select.tsv.gz') -> dict:
    """
    Build a lookup from the MANE Select reference file.
    Returns {ENST_base: (refseq_nuc, ensp_full, mane_status, gene)} where
    ENST_base is the version-stripped Ensembl transcript (e.g. 'ENST00000381652').
    """
    import gzip, re
    lookup = {}
    opener = gzip.open if mane_path.endswith('.gz') else open
    with opener(mane_path, 'rt') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        ci = {c.lstrip('#'): i for i, c in enumerate(header)}
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            enst_full = cols[ci['Ensembl_nuc']]     # e.g. ENST00000381652.4
            refseq    = cols[ci['RefSeq_nuc']]       # e.g. NM_004972.4
            ensp_full = cols[ci['Ensembl_prot']]     # e.g. ENSP00000371067.4
            status    = cols[ci['MANE_status']]       # MANE Select or MANE Plus Clinical
            gene      = cols[ci['symbol']]
            m = re.match(r'(ENST\d+)', enst_full)
            if m:
                enst_base = m.group(1)
                lookup[enst_base] = (refseq, ensp_full, status, gene)
    log.info('Loaded MANE lookup: %d transcripts (%s)',
             len(lookup), mane_path)
    return lookup


def load_ensembl_xref(xref_path: str = 'data/reference/ensembl_transcript_xref.tsv') -> dict:
    """
    Load supplementary Ensembl transcript cross-references for non-MANE
    transcripts. Returns {ENST_base: (refseq_nuc, ensp_id)}.
    """
    lookup = {}
    if not os.path.exists(xref_path):
        return lookup
    with open(xref_path, 'rt') as fh:
        fh.readline()  # skip header
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) >= 3:
                enst, ensp, refseq = cols[0], cols[1], cols[2]
                if enst:
                    lookup[enst] = (refseq, ensp)
    log.info('Loaded Ensembl xref: %d transcripts (%s)',
             len(lookup), xref_path)
    return lookup


def apply_mane_lookup(df: pd.DataFrame, mane_lookup: dict,
                      ensembl_xref: dict | None = None) -> pd.DataFrame:
    """
    Using the MANE lookup (and optional Ensembl xref fallback), set:
      - is_mane_select: True if the transcript is MANE Select
      - refseq_id: NM_ accession from MANE or Ensembl xref
      - hgvsp: prepend ENSP accession if bare p. notation
    """
    import re
    _RE_ENST = re.compile(r'(ENST\d+)')
    _RE_ENSP = re.compile(r'ENSP\d+')
    xref = ensembl_xref or {}

    def _enrich_row(hgvsc, hgvsp):
        refseq_id = ''
        is_mane = False
        ensp_prefix = ''

        if isinstance(hgvsc, str) and hgvsc:
            m = _RE_ENST.search(hgvsc)
            if m:
                enst_base = m.group(1)
                info = mane_lookup.get(enst_base)
                if info:
                    refseq_id, ensp_full, status, _ = info
                    is_mane = (status == 'MANE Select')
                    ensp_prefix = ensp_full
                elif enst_base in xref:
                    # Fallback to Ensembl xref for non-MANE transcripts
                    refseq_id, ensp_id = xref[enst_base]
                    ensp_prefix = ensp_id  # no version in xref

        # Prepend ENSP to bare p. hgvsp if we have a mapping
        new_hgvsp = hgvsp
        if (ensp_prefix and isinstance(hgvsp, str) and hgvsp
                and hgvsp not in ('', 'nan')
                and not _RE_ENSP.search(hgvsp)):
            new_hgvsp = f'{ensp_prefix}:{hgvsp}'

        return refseq_id, is_mane, new_hgvsp

    results = df.apply(
        lambda r: _enrich_row(
            r.get('hgvsc', ''),
            r.get('hgvsp', '')),
        axis=1, result_type='expand')
    df['refseq_id'] = results[0]
    df['is_mane_select'] = results[1]
    df['hgvsp'] = results[2]

    n_mane = df['is_mane_select'].sum()
    n_refseq = (df['refseq_id'] != '').sum()
    n_ensp = df['hgvsp'].fillna('').str.contains('ENSP', na=False).sum()
    log.info('MANE lookup: %d MANE Select, %d with RefSeq ID, %d with ENSP in hgvsp',
             n_mane, n_refseq, n_ensp)
    return df

# Sources that contribute to the sample count threshold
COUNT_SOURCES = {'COSMIC', 'TCGA', 'GENIE'}

# Sources that are annotation/curated (not count-based)
CURATED_SOURCES = {'OncoKB', 'ClinVar', 'CancerHotspots', 'TP53_somatic', 'TP53_germline'}

# Priority order for hgvsc/transcript selection during aggregation
# Lower number = higher priority
TRANSCRIPT_SOURCE_PRIORITY = {
    'CancerHotspots': 0,
    'GENIE':          1,
    'ClinVar':        2,
    'COSMIC':         3,
    'TCGA':           4,
    'TP53_somatic':   5,
    'TP53_germline':  5,
}

# VCF header template
_VCF_HEADER = """\
##fileformat=VCFv4.2
##fileDate={date}
##source=pan_cancer_whitelist_pipeline
##reference=GRCh38/hg38
##contig=<ID=1,length=248956422,assembly=GRCh38>
##contig=<ID=2,length=242193529,assembly=GRCh38>
##contig=<ID=3,length=198295559,assembly=GRCh38>
##contig=<ID=4,length=190214555,assembly=GRCh38>
##contig=<ID=5,length=181538259,assembly=GRCh38>
##contig=<ID=6,length=170805979,assembly=GRCh38>
##contig=<ID=7,length=159345973,assembly=GRCh38>
##contig=<ID=8,length=145138636,assembly=GRCh38>
##contig=<ID=9,length=138394717,assembly=GRCh38>
##contig=<ID=10,length=133797422,assembly=GRCh38>
##contig=<ID=11,length=135086622,assembly=GRCh38>
##contig=<ID=12,length=133275309,assembly=GRCh38>
##contig=<ID=13,length=114364328,assembly=GRCh38>
##contig=<ID=14,length=107043718,assembly=GRCh38>
##contig=<ID=15,length=101991189,assembly=GRCh38>
##contig=<ID=16,length=90338345,assembly=GRCh38>
##contig=<ID=17,length=83257441,assembly=GRCh38>
##contig=<ID=18,length=80373285,assembly=GRCh38>
##contig=<ID=19,length=58617616,assembly=GRCh38>
##contig=<ID=20,length=64444167,assembly=GRCh38>
##contig=<ID=21,length=46709983,assembly=GRCh38>
##contig=<ID=22,length=50818468,assembly=GRCh38>
##contig=<ID=X,length=156040895,assembly=GRCh38>
##contig=<ID=Y,length=57227415,assembly=GRCh38>
##contig=<ID=MT,length=16569,assembly=GRCh38>
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=CSQ,Number=1,Type=String,Description="Standardised consequence term">
##INFO=<ID=HGVSc,Number=1,Type=String,Description="HGVSc notation">
##INFO=<ID=HGVSp,Number=1,Type=String,Description="HGVSp notation">
##INFO=<ID=N_SAMPLES,Number=1,Type=Integer,Description="Total samples across count-based sources">
##INFO=<ID=N_CANCER_TYPES,Number=1,Type=Integer,Description="Number of distinct cancer types">
##INFO=<ID=CANCER_TYPES,Number=.,Type=String,Description="Comma-delimited array of contributing cancer types (VCFv4.2 Number=.)">
##INFO=<ID=SOURCES,Number=.,Type=String,Description="Comma-delimited array of contributing source databases (VCFv4.2 Number=.)">
##INFO=<ID=WL_TIER,Number=1,Type=Integer,Description="Whitelist tier: 1=highest, 3=minimum threshold">
##INFO=<ID=ONCOKB,Number=1,Type=String,Description="OncoKB oncogenicity classification (if available)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def _apply_data_dir(cfg: dict, data_dir: str) -> dict:
    """
    Replace the leading path component of all relative file paths in
    data_sources with data_dir. The config.yaml uses paths like
    'data/cosmic/file.tsv' — this replaces 'data/' with the supplied
    data_dir so paths are never doubled.
    """
    if not data_dir:
        return cfg
    data_dir = data_dir.rstrip('/')
    path_keys = {'tsv', 'vcf', 'maf', 'clinical_sample', 'somatic_tsv',
                 'germline_tsv', 'variants_file', 'dir', 'chain'}

    def _repath(v: str) -> str:
        if not v or os.path.isabs(v):
            return v
        # Strip any leading path component (e.g. 'data/') then join with data_dir
        parts = v.replace('\\', '/').split('/', 1)
        remainder = parts[1] if len(parts) > 1 else parts[0]
        return os.path.join(data_dir, remainder)

    for source in cfg.get('data_sources', {}).values():
        if not isinstance(source, dict):
            continue
        for k, v in source.items():
            if k in path_keys and isinstance(v, str):
                source[k] = _repath(v)
    for k, v in cfg.get('reference', {}).items():
        if isinstance(v, str):
            cfg['reference'][k] = _repath(v)
    return cfg


def load_config(path: str, data_dir: str = '') -> dict:
    with open(path) as fh:
        cfg = yaml.safe_load(fh)

    # Load and merge settings.yaml if specified
    settings_file = cfg.pop('settings_file', 'settings.yaml')
    cfg['settings_file'] = settings_file
    if os.path.exists(settings_file):
        with open(settings_file) as fh:
            settings = yaml.safe_load(fh) or {}
        cfg['thresholds']            = settings.get('thresholds', {})
        cfg['tiering']               = settings.get('tiering', {})
        cfg['vaf_rescue']            = settings.get('vaf_rescue', {})
        cfg['included_consequences'] = settings.get('included_consequences', [])
        cfg['log_level']             = settings.get('log_level', 'INFO')
        cfg['performance']           = settings.get('performance', {})
        for src in ('oncokb', 'clinvar', 'cancer_hotspots', 'cbioportal'):
            if src in settings and src in cfg.get('data_sources', {}):
                cfg['data_sources'][src].update(settings[src])
        log.info('Settings loaded from: %s', settings_file)
    else:
        log.warning('settings.yaml not found — using defaults')

    return cfg


def _finalise_config(cfg: dict, data_dir: str = '') -> dict:
    """Apply data_dir prefix after config is fully loaded."""
    if data_dir:
        cfg = _apply_data_dir(cfg, data_dir)
        log.info('Data directory: %s', data_dir)
    return cfg


def _oncokb_prefilter(merged_so_far: pd.DataFrame, thr: dict, logger) -> pd.DataFrame:
    """
    Aggregate raw merged rows by (gene, hgvsp) using Polars and return only
    pairs that pass count thresholds or come from curated sources.
    Falls back to a simple deduplication if Polars is unavailable.
    """
    import logging
    min_samples      = thr.get('min_samples_total', 10)
    min_cancer_types = thr.get('min_cancer_types', 1)
    count_pattern    = '|'.join(COUNT_SOURCES)

    logger.info('OncoKB pre-filter: aggregating %d raw rows via Polars...', len(merged_so_far))

    try:
        import polars as pl

        lf = pl.from_pandas(
            merged_so_far[['gene', 'hgvsp', 'source', 'cancer_type', 'n_samples']]
            .fillna('')
            .assign(n_samples=lambda d: pd.to_numeric(d['n_samples'], errors='coerce').fillna(0))
        ).lazy()

        # Count-source aggregation
        count_agg = (
            lf.filter(pl.col('source').str.contains(count_pattern))
            .group_by(['gene', 'hgvsp'])
            .agg([
                pl.col('n_samples').sum().alias('total_samples'),
                pl.col('cancer_type')
                  .filter(pl.col('cancer_type').str.to_lowercase() != 'unspecified')
                  .n_unique()
                  .alias('n_ct'),
            ])
            .filter(
                (pl.col('total_samples') >= min_samples) &
                (pl.col('n_ct') >= min_cancer_types)
            )
            .select(['gene', 'hgvsp'])
            .collect()
        )

        # Curated source pairs
        curated = (
            lf.filter(pl.col('source').str.contains('ClinVar|CancerHotspots'))
            .select(['gene', 'hgvsp'])
            .unique()
            .collect()
        )

        pair_pl = pl.concat([count_agg, curated]).unique()
        pair_df = pair_pl.to_pandas()

        logger.info('OncoKB pre-filter: %d count-passing + %d curated = %d unique pairs',
                    len(count_agg), len(curated), len(pair_df))

    except ImportError:
        logger.warning('Polars not available — using simple deduplication for OncoKB pre-filter')
        pair_df = (
            merged_so_far[['gene', 'hgvsp']]
            .dropna()
            .drop_duplicates()
            .query("gene != '' and hgvsp != ''")
            .reset_index(drop=True)
        )
        logger.info('OncoKB pre-filter (fallback): %d unique pairs', len(pair_df))

    pair_df = pair_df[
        pair_df['gene'].str.strip().ne('') &
        pair_df['hgvsp'].str.strip().ne('')
    ].reset_index(drop=True)

    logger.info('OncoKB pre-filter: %d pairs to query', len(pair_df))
    return pair_df


def run_parsers(cfg: dict, skip: set, inter_dir: str = 'intermediate') -> dict[str, pd.DataFrame]:
    """Run all enabled parsers and save each intermediate immediately on completion."""
    from parsers.parse_cosmic      import parse_cosmic
    from parsers.parse_oncokb      import parse_oncokb
    from parsers.parse_cbioportal  import parse_cbioportal
    from parsers.parse_clinvar     import parse_clinvar
    from parsers.parse_genie       import parse_genie
    from parsers.parse_tp53        import parse_tp53
    from parsers.parse_hotspots    import parse_hotspots
    from parsers.parse_tcga        import parse_tcga

    frames: dict[str, pd.DataFrame] = {}
    ds = cfg.get('data_sources', {})
    n_threads = cfg.get('performance', {}).get('threads', 4)
    os.makedirs(inter_dir, exist_ok=True)

    def _save(name: str, df: pd.DataFrame) -> None:
        if not df.empty:
            ipath = os.path.join(inter_dir, f'{name.lower()}.tsv.gz')
            df.to_csv(ipath, sep='\t', index=False, compression='gzip')
            log.info('Saved intermediate: %s  (%d rows)', ipath, len(df))

    def _should_run(name: str) -> bool:
        return name not in skip and ds.get(name, {}).get('enabled', True)

    if _should_run('cosmic'):
        log.info('=== Running COSMIC parser ===')
        frames['COSMIC'] = parse_cosmic(
            tsv_path            = ds['cosmic']['tsv'],
            vcf_path            = ds['cosmic'].get('vcf'),
            classification_path = ds['cosmic'].get('classification'),
            chunk_size          = ds['cosmic'].get('chunk_size', 200_000),
            n_threads           = n_threads,
        )
        _save('COSMIC', frames['COSMIC'])

    if _should_run('cbioportal'):
        log.info('=== Running cBioPortal parser ===')
        frames['cBioPortal'] = parse_cbioportal(
            api_base                = ds['cbioportal'].get('api_base', 'https://www.cbioportal.org/api'),
            use_tcga_pancancer_only = ds['cbioportal'].get('use_tcga_pancancer_only', False),
            max_studies             = ds['cbioportal'].get('max_studies', 0),
            request_delay_s         = ds['cbioportal'].get('request_delay_s', 0.5),
        )
        _save('cBioPortal', frames['cBioPortal'])

    if _should_run('tcga'):
        log.info('=== Running TCGA parser ===')
        frames['TCGA'] = parse_tcga(
            maf_path = ds['tcga']['maf'],
        )
        _save('TCGA', frames['TCGA'])

    if _should_run('clinvar'):
        log.info('=== Running ClinVar parser ===')
        frames['ClinVar'] = parse_clinvar(
            vcf_path       = ds['clinvar']['vcf'],
            include_clinsig = ds['clinvar'].get('include_clinsig'),
            somatic_only    = ds['clinvar'].get('somatic_only', True),
        )
        _save('ClinVar', frames['ClinVar'])

    if _should_run('genie'):
        log.info('=== Running GENIE parser ===')
        genie_maf = ds['genie']['maf']
        genie_dir = os.path.dirname(genie_maf)
        clinical  = ds['genie'].get('clinical_sample',
                                    os.path.join(genie_dir, 'data_clinical_sample.txt'))
        frames['GENIE'] = parse_genie(
            maf_path              = genie_maf,
            clinical_sample_path  = clinical if os.path.exists(clinical) else None,
        )
        _save('GENIE', frames['GENIE'])

    # Build task list for parallel execution
    tasks = {}

    if _should_run('tp53'):
        # Pass the reference FASTA so parse_tp53 can properly anchor indels
        # (deletion/insertion/duplication) at pos-1 with the real reference
        # base. Without this, TP53 indels are dropped rather than fabricated.
        tp53_ref_fasta = (cfg.get('reference') or {}).get('fasta')
        tasks['TP53'] = lambda: parse_tp53(
            somatic_tsv      = ds['tp53']['somatic_tsv'],
            germline_tsv     = ds['tp53'].get('germline_tsv'),
            include_germline = ds['tp53'].get('include_germline', False),
            reference_fasta  = tp53_ref_fasta,
        )

    if _should_run('cancer_hotspots'):
        tasks['CancerHotspots'] = lambda: parse_hotspots(
            tsv_path  = ds['cancer_hotspots'].get('tsv'),
            max_qvalue= ds['cancer_hotspots'].get('max_qvalue', 0.05),
        )


    if tasks:
        log.info('Running %d parsers in parallel (threads=%d)...', len(tasks), n_threads)
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            future_to_name = {executor.submit(fn): name for name, fn in tasks.items()}
            for future in as_completed(future_to_name):
                name = future_to_name[future]
                try:
                    frames[name] = future.result()
                    log.info('=== %s parser complete ===', name)
                    _save(name, frames[name])
                except Exception as exc:
                    log.error('%s parser failed: %s', name, exc)

    if _should_run('oncokb'):
        log.info('=== Running OncoKB parser ===')
        all_so_far = [df for df in frames.values() if not df.empty]
        merged_so_far = pd.concat(all_so_far, ignore_index=True) if all_so_far else pd.DataFrame()
        thr = cfg.get('thresholds', {})
        pair_df = _oncokb_prefilter(merged_so_far, thr, log)
        frames['OncoKB'] = parse_oncokb(
            variants_file        = ds['oncokb']['variants_file'],
            merged_df            = pair_df,
            include_oncogenicity = ds['oncokb'].get('include_oncogenicity'),
            api_token            = ds['oncokb'].get('api_token'),
        )
        _save('OncoKB', frames['OncoKB'])
    return frames



def merge_and_aggregate(frames: dict[str, pd.DataFrame],
                         included_consequences: set[str]) -> pd.DataFrame:
    """
    Merge all parser outputs and aggregate by variant key.

    Uses Polars for multi-threaded groupby, falling back to pandas if Polars
    is not installed. On a 40-core server, Polars completes in 2-3 minutes
    vs 6+ hours for the pandas single-threaded version.
    """
    all_dfs = [df for df in frames.values() if not df.empty]
    if not all_dfs:
        log.error('No parser produced any output')
        sys.exit(1)

    raw = pd.concat(all_dfs, ignore_index=True)
    log.info('Total raw rows across all sources: %d', len(raw))

    # Filter consequence types
    raw = raw[raw['consequence'].isin(included_consequences)].copy()
    log.info('After consequence filter: %d rows', len(raw))

    # Separate coordinate-resolved and coordinate-missing rows
    has_coords = raw['chrom'].notna() & (raw['chrom'] != '') & \
                 raw['pos'].notna() & (raw['ref'] != '') & (raw['alt'] != '')
    df_coord   = raw[has_coords].copy()
    df_nocoord = raw[~has_coords].copy()

    log.info('Rows with coordinates: %d', len(df_coord))
    log.info('Rows without coordinates (annotation-only): %d', len(df_nocoord))

    # Normalise chrom — strip chr prefix for consistency
    df_coord['chrom'] = df_coord['chrom'].astype(str).str.replace(r'^chr', '', regex=True)
    df_coord['pos']   = pd.to_numeric(df_coord['pos'], errors='coerce')
    df_coord['n_samples'] = pd.to_numeric(df_coord['n_samples'], errors='coerce').fillna(0)

    # Tag count-source rows
    count_pattern = '|'.join(COUNT_SOURCES)
    df_coord['_is_count_source'] = df_coord['source'].str.contains(count_pattern, na=False)

    log.info('Aggregating by variant key...')

    try:
        import polars as pl
        log.info('Using Polars for multi-threaded aggregation')
        df_agg = _aggregate_polars(df_coord)
    except ImportError:
        log.warning('Polars not installed — falling back to pandas (slow). Run: pip install polars')
        df_agg = _aggregate_pandas(df_coord)

    log.info('Aggregated: %d unique variants', len(df_agg))

    # -------------------------------------------------------------------------
    # OncoKB join — coordinate-based primary, (gene, hgvsp) fallback.
    # COSMIC and TCGA use different reference transcripts, causing Polars
    # first_nonempty to pick an incompatible hgvsp. Fix: resolve via trusted
    # sources (GENIE, ClinVar, CancerHotspots) before aggregation scrambles
    # the hgvsp values, then join by coordinate key post-aggregation.
    # -------------------------------------------------------------------------

    ONCOKB_TRUSTED_SOURCES = {'GENIE', 'ClinVar', 'CancerHotspots'}

    if not df_nocoord.empty and 'source' in df_nocoord.columns:
        oncokb_rows = df_nocoord[
            df_nocoord['source'].str.startswith('OncoKB:', na=False)
        ][['gene', 'hgvsp', 'oncokb_oncogenicity']].copy()

        if not oncokb_rows.empty:
            oncokb_rows = oncokb_rows.dropna(subset=['gene', 'hgvsp'])
            oncokb_rows = oncokb_rows[
                oncokb_rows['gene'].str.strip().ne('') &
                oncokb_rows['hgvsp'].str.strip().ne('')
            ]

            onco_priority = {'Oncogenic': 0, 'Likely Oncogenic': 1, 'Predicted Oncogenic': 2}
            oncokb_rows['_priority'] = oncokb_rows['oncokb_oncogenicity'].map(onco_priority).fillna(99)

            oncokb_by_hgvsp = (
                oncokb_rows.sort_values('_priority')
                .drop_duplicates(subset=['gene', 'hgvsp'], keep='first')
                .set_index(['gene', 'hgvsp'])['oncokb_oncogenicity']
                .to_dict()
            )
            log.info('OncoKB: %d (gene, hgvsp) entries in lookup', len(oncokb_by_hgvsp))

            trusted_coord = df_coord.loc[
                df_coord['source'].isin(ONCOKB_TRUSTED_SOURCES),
                ['chrom', 'pos', 'ref', 'alt', 'gene', 'hgvsp']
            ].copy()
            trusted_coord = trusted_coord[
                trusted_coord['hgvsp'].notna() &
                trusted_coord['hgvsp'].str.strip().ne('')
            ]
            trusted_coord['pos'] = trusted_coord['pos'].astype(int)
            # Vectorized OncoKB lookup via merge instead of apply
            _okb_df = pd.DataFrame(
                [(g, h, v) for (g, h), v in oncokb_by_hgvsp.items()],
                columns=['gene', 'hgvsp', '_oncokb']
            )
            trusted_coord = trusted_coord.merge(_okb_df, on=['gene', 'hgvsp'], how='left')
            trusted_coord['_oncokb'] = trusted_coord['_oncokb'].fillna('')
            trusted_coord = trusted_coord[trusted_coord['_oncokb'] != '']
            trusted_coord['_priority'] = trusted_coord['_oncokb'].map(onco_priority).fillna(99)
            oncokb_by_coord = (
                trusted_coord.sort_values('_priority')
                .drop_duplicates(subset=['chrom', 'pos', 'ref', 'alt'], keep='first')
                .set_index(['chrom', 'pos', 'ref', 'alt'])['_oncokb']
                .to_dict()
            )
            log.info('OncoKB: %d variants resolved to coordinates via trusted sources',
                     len(oncokb_by_coord))

            # Vectorized OncoKB join: coordinate-based, then gene+hgvsp fallback
            _coord_df = pd.DataFrame(
                [(c, p, r, a, v) for (c, p, r, a), v in oncokb_by_coord.items()],
                columns=['chrom', 'pos', 'ref', 'alt', '_okb_coord']
            )
            _coord_df['pos'] = _coord_df['pos'].astype(int)
            df_agg['pos'] = df_agg['pos'].astype(int)
            df_agg = df_agg.merge(_coord_df, on=['chrom', 'pos', 'ref', 'alt'], how='left')
            df_agg['_okb_coord'] = df_agg['_okb_coord'].fillna('')

            df_agg = df_agg.merge(_okb_df.rename(columns={'_oncokb': '_okb_hgvsp'}),
                                  on=['gene', 'hgvsp'], how='left')
            df_agg['_okb_hgvsp'] = df_agg['_okb_hgvsp'].fillna('')

            # Priority: existing > coordinate > hgvsp fallback
            existing = df_agg['oncokb_oncogenicity'].fillna('').astype(str)
            has_existing = existing.str.strip().ne('') & existing.ne('nan') & existing.ne('None')
            coord_hit = df_agg['_okb_coord'].ne('')
            df_agg['oncokb_oncogenicity'] = existing
            df_agg.loc[~has_existing & coord_hit, 'oncokb_oncogenicity'] = df_agg['_okb_coord']
            df_agg.loc[~has_existing & ~coord_hit, 'oncokb_oncogenicity'] = df_agg['_okb_hgvsp']

            n_annotated = (df_agg['oncokb_oncogenicity'] != '').sum()
            n_coord = coord_hit.sum()
            df_agg = df_agg.drop(columns=['_okb_coord', '_okb_hgvsp'])
            log.info('OncoKB: %d variants annotated after join', n_annotated)
            log.info('OncoKB: %d via coordinates, %d via (gene, hgvsp) fallback',
                     n_coord, n_annotated - n_coord)

    # -------------------------------------------------------------------------
    # ClinVar significance join — post-aggregation, coordinate-based.
    # Polars aggregation does not reliably carry this column through.
    # -------------------------------------------------------------------------
    if 'clinvar_clinical_significance' in df_coord.columns:
        clinvar_rows = df_coord[
            df_coord['source'] == 'ClinVar'
        ][['chrom', 'pos', 'ref', 'alt', 'clinvar_clinical_significance']].copy()
        clinvar_rows = clinvar_rows[
            clinvar_rows['clinvar_clinical_significance'].str.strip().ne('')
        ]
        if not clinvar_rows.empty:
            clinvar_rows['pos'] = clinvar_rows['pos'].astype(int)
            clinvar_dedup = clinvar_rows.drop_duplicates(
                subset=['chrom', 'pos', 'ref', 'alt'], keep='first'
            ).rename(columns={'clinvar_clinical_significance': '_cv_sig'})
            log.info('ClinVar significance lookup: %d entries', len(clinvar_dedup))

            # Vectorized merge instead of row-by-row apply
            df_agg['pos'] = df_agg['pos'].astype(int)
            df_agg = df_agg.merge(
                clinvar_dedup[['chrom', 'pos', 'ref', 'alt', '_cv_sig']],
                on=['chrom', 'pos', 'ref', 'alt'], how='left'
            )
            # Use merged value where existing is empty
            existing_cv = df_agg['clinvar_clinical_significance'].fillna('')
            merged_cv = df_agg['_cv_sig'].fillna('')
            df_agg['clinvar_clinical_significance'] = existing_cv.where(
                existing_cv != '', merged_cv
            )
            df_agg = df_agg.drop(columns=['_cv_sig'])
            n_cv = (df_agg['clinvar_clinical_significance'] != '').sum()
            log.info('ClinVar significance: %d variants annotated after join', n_cv)

    return df_agg


def _clean_cancer_type(val: str) -> str:
    """Sanitise a single cancer type string."""
    return str(val).replace('\t', ' ').replace('\n', ' ').replace('|', '-').strip()


def _aggregate_polars(df_coord: pd.DataFrame) -> pd.DataFrame:
    """Multi-threaded aggregation using Polars."""
    import polars as pl

    # Ensure annotation columns exist
    if 'oncokb_oncogenicity' not in df_coord.columns:
        df_coord['oncokb_oncogenicity'] = ''
    if 'clinvar_clinical_significance' not in df_coord.columns:
        df_coord['clinvar_clinical_significance'] = ''
    if 'tp53_class' not in df_coord.columns:
        df_coord['tp53_class'] = ''

    # Fill nulls for string columns before converting
    str_cols = ['chrom', 'ref', 'alt', 'source', 'cancer_type', 'gene',
                'consequence', 'hgvsc', 'hgvsp', 'oncokb_oncogenicity',
                'clinvar_clinical_significance', 'tp53_class']
    for c in str_cols:
        if c in df_coord.columns:
            df_coord[c] = df_coord[c].fillna('').astype(str)

    # Add source priority for deterministic hgvsc selection
    df_coord['_src_priority'] = df_coord['source'].map(
        TRANSCRIPT_SOURCE_PRIORITY).fillna(99).astype(int)

    lf = pl.from_pandas(df_coord).lazy()

    # Sort by priority so first_nonempty picks the best source's transcript
    lf = lf.sort(['chrom', 'pos', 'ref', 'alt', '_src_priority'])

    key_cols = ['chrom', 'pos', 'ref', 'alt']

    def first_nonempty(col):
        """Return first non-empty, non-null, non-'nan' value."""
        return (
            pl.col(col)
            .filter(
                pl.col(col).is_not_null() &
                (pl.col(col).str.strip_chars() != '') &
                (pl.col(col).str.to_lowercase() != 'nan') &
                (pl.col(col).str.to_lowercase() != 'none') &
                (pl.col(col).str.to_lowercase() != 'na')
            )
            .first()
            .fill_null('')
            .alias(col)
        )

    # Sources: sorted unique pipe-joined
    sources_expr = (
        pl.col('source')
        .filter(pl.col('source').is_not_null() & (pl.col('source') != ''))
        .unique()
        .sort()
        .str.join('|')
        .alias('sources')
    )

    # Cancer types: sanitise, sorted unique pipe-joined
    cancer_expr = (
        pl.col('cancer_type')
        .filter(
            pl.col('cancer_type').is_not_null() &
            (pl.col('cancer_type') != '') &
            (pl.col('cancer_type').str.to_lowercase() != 'unspecified')
        )
        .map_elements(lambda v: _clean_cancer_type(v), return_dtype=pl.String)
        .unique()
        .sort()
        .str.join('|')
        .fill_null('unspecified')
        .alias('cancer_types')
    )

    n_cancer_expr = (
        pl.col('cancer_type')
        .filter(
            pl.col('cancer_type').is_not_null() &
            (pl.col('cancer_type') != '') &
            (pl.col('cancer_type').str.to_lowercase() != 'unspecified')
        )
        .n_unique()
        .alias('n_cancer_types')
    )

    # n_samples from count sources only
    n_samples_expr = (
        pl.when(pl.col('_is_count_source'))
        .then(pl.col('n_samples'))
        .otherwise(0.0)
        .sum()
        .cast(pl.Int64)
        .alias('n_samples')
    )

    # hgvsc/hgvsp: priority-aware selection (sort guarantees best source first)
    _hgvsc_filter = (
        pl.col('hgvsc').is_not_null() &
        (pl.col('hgvsc').str.strip_chars() != '') &
        (pl.col('hgvsc').str.to_lowercase() != 'nan')
    )
    hgvsc_expr = (
        pl.col('hgvsc')
        .filter(_hgvsc_filter)
        .first()
        .fill_null('')
        .alias('hgvsc')
    )
    hgvsp_expr = (
        pl.col('hgvsp')
        .filter(
            pl.col('hgvsp').is_not_null() &
            (pl.col('hgvsp').str.strip_chars() != '') &
            (pl.col('hgvsp').str.to_lowercase() != 'nan')
        )
        .first()
        .fill_null('')
        .alias('hgvsp')
    )
    # Source that provided the winning hgvsc
    transcript_source_expr = (
        pl.col('source')
        .filter(_hgvsc_filter)
        .first()
        .fill_null('')
        .alias('transcript_source')
    )
    df_agg = (
        lf.group_by(key_cols)
        .agg([
            sources_expr,
            cancer_expr,
            n_cancer_expr,
            n_samples_expr,
            first_nonempty('gene'),
            first_nonempty('consequence'),
            hgvsc_expr,
            hgvsp_expr,
            transcript_source_expr,
            first_nonempty('oncokb_oncogenicity'),
            first_nonempty('clinvar_clinical_significance'),
            first_nonempty('tp53_class'),
        ])
        .collect(engine="streaming")
    ).to_pandas()

    df_agg['n_cancer_types'] = df_agg['n_cancer_types'].fillna(1).astype(int)
    df_agg['n_samples']      = df_agg['n_samples'].fillna(0).astype(int)
    df_agg['pos']            = df_agg['pos'].astype(int)
    return df_agg


def _aggregate_pandas(df_coord: pd.DataFrame) -> pd.DataFrame:
    """Single-threaded pandas fallback aggregation."""
    key_cols = ['chrom', 'pos', 'ref', 'alt']

    if 'oncokb_oncogenicity' not in df_coord.columns:
        df_coord['oncokb_oncogenicity'] = ''
    if 'tp53_class' not in df_coord.columns:
        df_coord['tp53_class'] = ''
    if 'clinvar_clinical_significance' not in df_coord.columns:
        df_coord['clinvar_clinical_significance'] = ''

    # Add source priority for deterministic hgvsc selection
    df_coord['_src_priority'] = df_coord['source'].map(
        TRANSCRIPT_SOURCE_PRIORITY).fillna(99).astype(int)
    df_coord = df_coord.sort_values(
        ['chrom', 'pos', 'ref', 'alt', '_src_priority'])

    def first_nonempty_agg(s):
        vals = s.dropna()
        vals = vals[~vals.str.strip().str.lower().isin(['', 'nan', 'none', 'na'])]
        return vals.iloc[0] if len(vals) else ''

    def clean_ct_agg(x):
        cleaned = []
        for t in x.dropna().unique():
            t = _clean_cancer_type(str(t))
            if t and t.lower() != 'unspecified':
                cleaned.append(t)
        return '|'.join(sorted(cleaned)) or 'unspecified'

    sources_agg        = df_coord.groupby(key_cols)['source'].agg(lambda x: '|'.join(sorted(x.dropna().unique())))
    cancer_type_agg    = df_coord.groupby(key_cols)['cancer_type'].agg(clean_ct_agg)
    n_cancer_types_agg = df_coord.groupby(key_cols)['cancer_type'].agg(
        lambda x: len(set(t for t in x.dropna().unique() if t and t != 'unspecified')) or 1)
    n_samples_agg      = df_coord[df_coord['_is_count_source']].groupby(key_cols)['n_samples'].sum()
    gene_agg           = df_coord.groupby(key_cols)['gene'].agg(first_nonempty_agg)
    consequence_agg    = df_coord.groupby(key_cols)['consequence'].agg(first_nonempty_agg)
    hgvsc_agg          = df_coord.groupby(key_cols)['hgvsc'].agg(first_nonempty_agg)
    hgvsp_agg          = df_coord.groupby(key_cols)['hgvsp'].agg(first_nonempty_agg)
    oncokb_agg      = df_coord.groupby(key_cols)['oncokb_oncogenicity'].agg(first_nonempty_agg)
    clinvar_sig_agg = df_coord.groupby(key_cols)['clinvar_clinical_significance'].agg(first_nonempty_agg)
    tp53_class_agg  = df_coord.groupby(key_cols)['tp53_class'].agg(first_nonempty_agg)

    # transcript_source: which source provided the winning hgvsc
    _hgvsc_valid = df_coord[
        df_coord['hgvsc'].notna() &
        (df_coord['hgvsc'].str.strip() != '') &
        (~df_coord['hgvsc'].str.lower().isin(['nan', 'none', 'na']))
    ]
    transcript_source_agg = _hgvsc_valid.groupby(key_cols)['source'].first()

    df_agg = pd.DataFrame({
        'sources':                       sources_agg,
        'cancer_types':                  cancer_type_agg,
        'n_cancer_types':                n_cancer_types_agg,
        'n_samples':                     n_samples_agg,
        'gene':                          gene_agg,
        'consequence':                   consequence_agg,
        'hgvsc':                         hgvsc_agg,
        'hgvsp':                         hgvsp_agg,
        'transcript_source':             transcript_source_agg,
        'oncokb_oncogenicity':           oncokb_agg,
        'clinvar_clinical_significance': clinvar_sig_agg,
        'tp53_class':                    tp53_class_agg,
    }).reset_index()

    df_agg['n_samples']   = df_agg['n_samples'].fillna(0).astype(int)
    df_agg['pos']         = df_agg['pos'].astype(int)
    return df_agg


def apply_filters(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """
    Apply inclusion filters.

    A variant passes if ANY of the following conditions is true:
      1. n_samples >= min_samples AND n_cancer_types >= min_cancer_types
         (count-based threshold)
      2. OncoKB oncogenicity = Oncogenic | Likely Oncogenic | Predicted Oncogenic
         (expert-curated override)
      3. ClinVar is a source (pathogenic somatic evidence)
      4. CancerHotspots is a source (statistically significant hotspot)
    """
    thr = cfg.get('thresholds', {})
    min_samples      = thr.get('min_samples_total', 10)
    min_cancer_types = thr.get('min_cancer_types', 1)

    passes_count   = (df['n_samples'] >= min_samples) & \
                     (df['n_cancer_types'] >= min_cancer_types)
    passes_oncokb  = df['oncokb_oncogenicity'].isin(
                         ['Oncogenic', 'Likely Oncogenic', 'Predicted Oncogenic'])
    passes_clinvar = df['sources'].str.contains('ClinVar', na=False)
    passes_hotspot = df["sources"].str.contains("CancerHotspots", na=False)
    passes_tp53    = df["sources"].str.contains("TP53", na=False)

    mask = passes_count | passes_oncokb | passes_clinvar | passes_hotspot | passes_tp53

    df_pass = df[mask].copy()
    log.info('After filters: %d / %d variants retained', len(df_pass), len(df))
    log.info('  Count-based:   %d', passes_count.sum())
    log.info('  OncoKB:        %d', passes_oncokb.sum())
    log.info('  ClinVar:       %d', passes_clinvar.sum())
    log.info("  CancerHotspot: %d", passes_hotspot.sum())
    log.info("  TP53:          %d", passes_tp53.sum())
    return df_pass


def assign_tiers(df: pd.DataFrame, cfg: dict | None = None) -> pd.DataFrame:
    """
    Assign whitelist tier to each variant.

    Tier 1: OncoKB Oncogenic/Likely Oncogenic
            OR (n_samples >= tier1_min_samples AND n_cancer_types >= tier1_min_cancer_types)
    Tier 2: OncoKB Predicted Oncogenic OR ClinVar/CancerHotspots
            OR (n_samples >= tier2_min_samples AND n_cancer_types >= tier2_min_cancer_types)
    Tier 3: passes base thresholds (minimum confidence)

    Thresholds are read from settings.yaml via cfg['tiering'].
    """
    tiering = (cfg or {}).get('tiering', {})

    # Validate that the settings.yaml keys are as expected before falling back
    # to defaults — a single typo silently substituted the built-in defaults
    # and produced a mistiered whitelist with no evidence in the log.
    _EXPECTED_TIER_KEYS = {
        'tier1_min_samples', 'tier1_min_cancer_types',
        'tier2_min_samples', 'tier2_min_cancer_types',
    }
    if tiering:
        extra = set(tiering) - _EXPECTED_TIER_KEYS
        if extra:
            log.warning('Unrecognised keys in settings.yaml tiering block '
                        '(possible typo — will not affect tiering): %s',
                        sorted(extra))
        missing = _EXPECTED_TIER_KEYS - set(tiering)
        if missing:
            log.warning('Missing keys in settings.yaml tiering block '
                        '(using defaults): %s', sorted(missing))

    t1_s  = tiering.get('tier1_min_samples',       50)
    t1_ct = tiering.get('tier1_min_cancer_types',   3)
    t2_s  = tiering.get('tier2_min_samples',        25)
    t2_ct = tiering.get('tier2_min_cancer_types',    2)

    log.info('Tier thresholds: T1 >= %d samples & >= %d cancer types  |  '
             'T2 >= %d samples & >= %d cancer types', t1_s, t1_ct, t2_s, t2_ct)

    # Vectorised tiering — avoids iterrows on large dataframes
    onco = df['oncokb_oncogenicity'].fillna('').astype(str)
    n_s  = df['n_samples'].astype(int)
    n_ct = df['n_cancer_types'].astype(int)
    srcs = df['sources'].fillna('').astype(str)

    is_tier1 = (
        onco.isin(['Oncogenic', 'Likely Oncogenic'])
        | ((n_s >= t1_s) & (n_ct >= t1_ct))
    )
    is_tier2 = (
        (onco == 'Predicted Oncogenic')
        | srcs.str.contains('ClinVar|CancerHotspots', na=False)
        | ((n_s >= t2_s) & (n_ct >= t2_ct))
    )

    df['wl_tier'] = 3
    df.loc[is_tier2 & ~is_tier1, 'wl_tier'] = 2
    df.loc[is_tier1, 'wl_tier'] = 1
    log.info('Tier distribution: Tier1=%d  Tier2=%d  Tier3=%d',
             (df['wl_tier'] == 1).sum(),
             (df['wl_tier'] == 2).sum(),
             (df['wl_tier'] == 3).sum())
    return df


def write_tsv(df: pd.DataFrame, path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    out_cols = [
        'chrom', 'pos', 'ref', 'alt', 'gene', 'hgvsc', 'hgvsp',
        'consequence', 'n_cancer_types', 'cancer_types', 'n_samples',
        'sources', 'oncokb_oncogenicity', 'clinvar_clinical_significance',
        'transcript_source', 'is_mane_select', 'refseq_id',
        'tp53_class', 'wl_tier',
    ]
    # Karyotypic sort so TSV row order matches the VCF (1,2,...,22,X,Y,MT).
    # Lexicographic sort produced 1,10,11,...,2 which surprised users diffing
    # the TSV against the VCF.
    chrom_order = {str(i): i for i in range(1, 23)}
    chrom_order.update({'X': 23, 'Y': 24, 'MT': 25, 'M': 25})
    df_out = df[out_cols].copy()
    chrom_clean = df_out['chrom'].astype(str).str.replace(r'^chr', '', regex=True)
    df_out['_chrom_sort'] = chrom_clean.map(lambda c: chrom_order.get(c, 99))
    df_out = df_out.sort_values(['_chrom_sort', 'pos']).drop(columns='_chrom_sort')
    df_out.to_csv(path, sep='\t', index=False,
                  compression='gzip' if path.endswith('.gz') else None)
    log.info('Wrote TSV: %s  (%d rows)', path, len(df_out))


def write_vcf(df: pd.DataFrame, path: str) -> None:
    """Write whitelist as a sites-only VCF (no sample columns)."""

    def _vcf_escape(val):
        """Percent-encode VCF INFO special characters.

        Encodes the four reserved-in-INFO chars (%, ;, =, space) plus comma
        (which is the VCFv4.2 array separator) and CR. Tabs and newlines are
        percent-encoded rather than silently deleted so their presence
        remains visible for debugging.
        """
        if not isinstance(val, str):
            val = str(val)
        return (val
                .replace('%',  '%25')
                .replace(';',  '%3B')
                .replace('=',  '%3D')
                .replace(',',  '%2C')
                .replace(' ',  '%20')
                .replace('\r', '%0D')
                .replace('\t', '%09')
                .replace('\n', '%0A'))

    def _vcf_array(pipe_joined, sep_out=','):
        """Convert a pipe-delimited aggregator string to a VCF INFO array.

        The aggregator emits `A|B|C`; VCFv4.2 arrays expect `A,B,C`. Each
        element is individually escaped so embedded commas / spaces /
        semicolons in a single cancer_type or source label do not break
        array parsing downstream.
        """
        if pipe_joined is None:
            return ''
        s = str(pipe_joined)
        if not s:
            return ''
        return sep_out.join(_vcf_escape(part) for part in s.split('|'))

    os.makedirs(os.path.dirname(path), exist_ok=True)

    # Sort by chromosome order then position
    chrom_order = {str(i): i for i in range(1, 23)}
    chrom_order.update({'X': 23, 'Y': 24, 'MT': 25, 'M': 25})
    df = df.copy()
    df['chrom'] = df['chrom'].astype(str).str.replace(r'^chr', '', regex=True)
    df.loc[df['chrom'] == 'M', 'chrom'] = 'MT'
    df['_chrom_sort'] = df['chrom'].map(lambda c: chrom_order.get(c, 99))
    df = df.sort_values(['_chrom_sort', 'pos']).drop(columns='_chrom_sort')

    header = _VCF_HEADER.format(date=datetime.today().strftime('%Y%m%d'))

    import subprocess

    # Write to temp uncompressed file then bgzip for tabix compatibility
    tmp_path = path.replace('.gz', '') if path.endswith('.gz') else path + '.tmp'
    with open(tmp_path, 'wt') as fh:
        fh.write(header)
        for row in df.itertuples(index=False):
            info_parts = [
                f'GENE={_vcf_escape(row.gene)}',
                f'CSQ={_vcf_escape(row.consequence)}',
                f'HGVSc={_vcf_escape(row.hgvsc)}',
                f'HGVSp={_vcf_escape(row.hgvsp)}',
                f'N_SAMPLES={row.n_samples}',
                f'N_CANCER_TYPES={row.n_cancer_types}',
                # CANCER_TYPES and SOURCES declared as Number=. arrays; the
                # aggregator stores them pipe-joined so the TSV stays human
                # readable — convert to a comma-separated VCF array here.
                f'CANCER_TYPES={_vcf_array(row.cancer_types)}',
                f'SOURCES={_vcf_array(row.sources)}',
                f'WL_TIER={row.wl_tier}',
            ]
            onco = getattr(row, 'oncokb_oncogenicity', None)
            if isinstance(onco, str) and onco and onco.lower() not in ('nan', 'none'):
                info_parts.append(f'ONCOKB={_vcf_escape(onco)}')

            info = ';'.join(p for p in info_parts if '=' in p and p.split('=', 1)[1]) or '.'
            fh.write('\t'.join([
                str(row.chrom),
                str(row.pos),
                '.',
                str(row.ref),
                str(row.alt),
                '.',
                'PASS',
                info,
            ]) + '\n')

    if path.endswith('.gz'):
        import shutil, gzip as _gzip
        if shutil.which('bgzip'):
            subprocess.run(['bgzip', '-f', tmp_path], check=True)
            # bgzip writes to tmp_path + '.gz' by default, rename if needed
            bgz_path = tmp_path + '.gz'
            if bgz_path != path:
                os.rename(bgz_path, path)
        else:
            log.warning('bgzip not found; using Python gzip (not BGZF — tabix indexing may fail)')
            with open(tmp_path, 'rb') as f_in, _gzip.open(path, 'wb') as f_out:
                f_out.writelines(f_in)
            os.remove(tmp_path)

    log.info('Wrote VCF: %s  (%d variants)', path, len(df))


# ── Entry point ───────────────────────────────────────────────────────────────


def log_database_versions(cfg: dict, settings_file: str = 'settings.yaml') -> None:
    """
    Log the version/date of each database used in the pipeline.
    Writes a human-readable summary to the log and to output/database_versions.txt.
    """
    import re
    import requests
    from datetime import datetime, timezone

    ds      = cfg.get('data_sources', {})
    out_dir = cfg.get('output', {}).get('dir', 'output')
    os.makedirs(out_dir, exist_ok=True)
    lines   = []

    def _add(source, version, path_or_url=''):
        msg = f'  {source:<25} {version}'
        if path_or_url:
            msg += f'  ({path_or_url})'
        log.info(msg)
        lines.append(msg)

    log.info('=' * 60)
    log.info('ONCOSIEVE — pan-cancer variant curation and rescue tool')
    log.info('Author   : Dr Christopher Trethewey')
    log.info('Email    : christopher.trethewey@nhs.net')
    log.info('=' * 60)
    log.info('DATABASE VERSIONS')
    log.info('=' * 60)

    # COSMIC — extract version from filename
    if ds.get('cosmic', {}).get('enabled'):
        tsv = ds['cosmic'].get('tsv', '')
        m = re.search(r'v(\d+)', tsv)
        version = f'v{m.group(1)}' if m else 'unknown'
        _add('COSMIC', version, tsv)

    # GENIE — hardcoded v19
    if ds.get('genie', {}).get('enabled'):
        maf = ds['genie'].get('maf', '')
        _add('GENIE', 'v19', maf)

    # ClinVar — read date from VCF header
    if ds.get('clinvar', {}).get('enabled'):
        vcf = ds['clinvar'].get('vcf', '')
        version = 'unknown'
        if os.path.exists(vcf):
            try:
                import pysam
                v = pysam.VariantFile(vcf)
                for rec in v.header.records:
                    if 'fileDate' in str(rec) or 'dbSNP' in str(rec):
                        version = str(rec).strip()
                        break
                v.close()
            except Exception:
                pass
        _add('ClinVar', version, vcf)

    # OncoKB — query API for data version
    if ds.get('oncokb', {}).get('enabled'):
        version = 'unknown'
        token = cfg.get('data_sources', {}).get('oncokb', {}).get('api_token', '')
        try:
            headers = {'Accept': 'application/json'}
            if token:
                headers['Authorization'] = f'Bearer {token}'
            resp = requests.get(
                'https://www.oncokb.org/api/v1/info',
                headers=headers, timeout=10
            )
            if resp.status_code == 200:
                info = resp.json()
                dv = info.get('dataVersion', {})
                version = f"{dv.get('version','?')} ({dv.get('date','?')})"
        except Exception:
            pass
        _add('OncoKB', version, 'https://www.oncokb.org/api/v1')

    # TP53
    if ds.get('tp53', {}).get('enabled'):
        tsv = ds['tp53'].get('somatic_tsv', '')
        version = ds['tp53'].get('version')
        if not version:
            version = 'unknown'
            if os.path.exists(tsv):
                mtime = os.path.getmtime(tsv)
                version = f'downloaded {datetime.fromtimestamp(mtime).strftime("%Y-%m-%d")}'
        _add('TP53 database', version, tsv)

    # TCGA — version from filename
    if ds.get('tcga', {}).get('enabled'):
        maf = ds['tcga'].get('maf', '')
        _add('TCGA mc3', 'v0.2.8', maf)

    # cBioPortal — live API, log URL
    if ds.get('cbioportal', {}).get('enabled'):
        version = 'live API'
        try:
            resp = requests.get(
                f"{ds['cbioportal'].get('api_base','https://www.cbioportal.org/api')}/info",
                timeout=10
            )
            if resp.status_code == 200:
                info = resp.json()
                version = f"live API ({info.get('portalVersion', '?')})"
        except Exception:
            pass
        _add('cBioPortal', version, ds['cbioportal'].get('api_base', ''))

    # CancerHotspots — live API
    if ds.get('cancer_hotspots', {}).get('enabled'):
        _add('CancerHotspots', 'live API v2', 'https://www.cancerhotspots.org/api')

    log.info('=' * 60)

    # Write to file
    run_date = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
    out_path = os.path.join(out_dir, 'database_versions.txt')
    with open(out_path, 'w') as fh:
        fh.write('ONCOSIEVE — pan-cancer variant curation and rescue tool\n')
        fh.write('Author   : Dr Christopher Trethewey\n')
        fh.write('Email    : christopher.trethewey@nhs.net\n')
        fh.write('=' * 60 + '\n')
        fh.write(f'Run date : {run_date}\n')
        fh.write(f'Settings : {settings_file}\n')
        fh.write('=' * 60 + '\n')
        for line in lines:
            fh.write(line.strip() + '\n')
    log.info('Database version log written to: %s', out_path)

def main():
    parser = argparse.ArgumentParser(
        description='Build pan-cancer whitelist from multiple somatic variant sources.'
    )
    parser.add_argument('--config', default='config.yaml',
                        help='Path to config.yaml (default: config.yaml)')
    parser.add_argument('--skip-sources', default='',
                        help='Comma-separated source names to skip '
                             '(e.g. cosmic,genie)')
    parser.add_argument('--from-intermediates', action='store_true',
                        help='Skip parsers and load from saved intermediate files')
    parser.add_argument('--rerun-oncokb', action='store_true',
                        help='Load all intermediates then re-run OncoKB against merged variants')
    parser.add_argument('--intermediate-only', action='store_true',
                        help='Run parsers and save intermediates; do not merge')
    parser.add_argument('--data-dir', default='',
                        help='Path to reference data directory. Overrides relative '
                             'paths in config.yaml (e.g. /path/to/reference/)')
    args = parser.parse_args()

    cfg  = load_config(args.config)
    cfg  = _finalise_config(cfg, args.data_dir)
    skip = set(s.strip().lower() for s in args.skip_sources.split(',') if s.strip())

    log.setLevel(getattr(logging, cfg.get('log_level', 'INFO')))

    # Log database versions
    log_database_versions(cfg, settings_file=cfg.get('settings_file', 'settings.yaml'))

    inter_dir = cfg.get('intermediate_dir', 'intermediate')
    out_dir   = cfg.get('output', {}).get('dir', 'output')
    prefix    = cfg.get('output', {}).get('prefix', 'pan_cancer_whitelist_GRCh38')
    os.makedirs(inter_dir, exist_ok=True)
    os.makedirs(out_dir,   exist_ok=True)

    # Run all parsers (or load from intermediates)
    if args.rerun_oncokb:
        log.info('--rerun-oncokb: loading all intermediates except OncoKB...')
        frames = {}
        for fname in os.listdir(inter_dir):
            if fname.endswith('.tsv.gz') and 'oncokb' not in fname.lower():
                name = fname.replace('.tsv.gz', '')
                path = os.path.join(inter_dir, fname)
                df = pd.read_csv(path, sep='\t', dtype=str, low_memory=False)
                for col in ('n_samples', 'n_cancer_types', 'pos'):
                    if col in df.columns:
                        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
                frames[name] = df
                log.info('Loaded intermediate: %s  (%d rows)', path, len(df))
        merged_so_far = pd.concat([df for df in frames.values() if not df.empty], ignore_index=True)
        log.info('OncoKB re-run: %d raw rows loaded from intermediates', len(merged_so_far))
        from parsers.parse_oncokb import parse_oncokb
        ds      = cfg['data_sources']
        thr     = cfg.get('thresholds', {})
        pair_df = _oncokb_prefilter(merged_so_far, thr, log)
        frames['OncoKB'] = parse_oncokb(
            variants_file        = ds['oncokb']['variants_file'],
            merged_df            = pair_df,
            include_oncogenicity = ds['oncokb'].get('include_oncogenicity'),
            api_token            = ds['oncokb'].get('api_token'),
        )
        ipath = os.path.join(inter_dir, 'oncokb.tsv.gz')
        frames['OncoKB'].to_csv(ipath, sep='\t', index=False, compression='gzip')
        log.info('Saved intermediate: %s  (%d rows)', ipath, len(frames['OncoKB']))
    elif args.from_intermediates:
        log.info('Loading from intermediate files in: %s', inter_dir)
        frames = {}
        for fname in os.listdir(inter_dir):
            if fname.endswith('.tsv.gz'):
                name = fname.replace('.tsv.gz', '')
                path = os.path.join(inter_dir, fname)
                df = pd.read_csv(path, sep='\t', dtype=str, low_memory=False)
                # Restore numeric dtypes lost when saving/loading as TSV
                for col in ('n_samples', 'n_cancer_types', 'pos'):
                    if col in df.columns:
                        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
                frames[name] = df
                log.info('Loaded intermediate: %s  (%d rows)', path, len(df))

        # Warn loudly if any enabled source has no intermediate.
        # A missing intermediate means that source is silently absent from
        # the whitelist, which can cause severe output degradation.
        ds_check = cfg.get('data_sources', {})
        source_to_frame = {
            'cosmic': 'cosmic', 'genie': 'genie', 'tcga': 'tcga',
            'clinvar': 'clinvar', 'tp53': 'tp53',
            'cancer_hotspots': 'cancerhotspots', 'cbioportal': 'cbioportal',
        }
        for src_key, frame_key in source_to_frame.items():
            if ds_check.get(src_key, {}).get('enabled', False):
                if frame_key not in frames and src_key not in skip:
                    log.warning(
                        '='*60
                    )
                    log.warning(
                        'MISSING intermediate for enabled source: %s '
                        '(expected: %s/%s.tsv.gz)',
                        src_key.upper(), inter_dir, frame_key
                    )
                    log.warning(
                        'This source will be ABSENT from the whitelist. '
                        'Run a full pipeline (without --from-intermediates) '
                        'to reparse this source.'
                    )
                    log.warning('='*60)
    else:
        frames = run_parsers(cfg, skip, inter_dir=inter_dir)

    if args.intermediate_only:
        log.info('--intermediate-only flag set; stopping before merge')
        return

    # Merge and aggregate
    included_csq = set(cfg.get('included_consequences', list(INCLUDED_CONSEQUENCES)))
    df_merged = merge_and_aggregate(frames, included_csq)

    # MANE Select lookup — set is_mane_select and refseq_id from reference
    mane_path = os.path.join(
        cfg.get('data_dir', ''), 'data', 'reference', 'mane_select.tsv.gz')
    if not os.path.exists(mane_path):
        mane_path = 'data/reference/mane_select.tsv.gz'
    if os.path.exists(mane_path):
        mane_lookup = load_mane_lookup(mane_path)
        xref_path = os.path.join(os.path.dirname(mane_path),
                                 'ensembl_transcript_xref.tsv')
        ensembl_xref = load_ensembl_xref(xref_path)
        df_merged = apply_mane_lookup(df_merged, mane_lookup, ensembl_xref)
    else:
        log.warning('MANE reference not found at %s — skipping MANE lookup', mane_path)
        df_merged['is_mane_select'] = False
        df_merged['refseq_id'] = ''

    # Apply filters
    df_filtered = apply_filters(df_merged, cfg)

    # Assign tiers
    df_tiered = assign_tiers(df_filtered, cfg=cfg)

    # Write outputs
    tsv_path = os.path.join(out_dir, f'{prefix}.tsv.gz')
    vcf_path = os.path.join(out_dir, f'{prefix}.vcf.gz')
    write_tsv(df_tiered, tsv_path)
    write_vcf(df_tiered, vcf_path)

    log.info('Pipeline complete. Output in: %s', out_dir)
    log.info('  TSV: %s', tsv_path)
    log.info('  VCF: %s', vcf_path)
    log.info('  Index VCF with: tabix -p vcf %s', vcf_path)


if __name__ == '__main__':
    main()
