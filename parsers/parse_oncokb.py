#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Parser: OncoKB annotation API
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
parse_oncokb.py
Query the OncoKB public API to annotate variants with oncogenicity.

No API token is required for the public endpoint. The endpoint returns
the 'oncogenic' field (Oncogenic, Likely Oncogenic, Predicted Oncogenic,
Likely Neutral, Unknown) for all submitted variants.

If a local allAnnotatedVariants.txt file exists (from a token-based download),
it is used instead to avoid API calls.
"""

import os
import re
import time

import pandas as pd
import requests

from parsers.common import (
    STANDARD_COLS,
    empty_standard_df,
    setup_logger,
)

log = setup_logger('OncoKB')

_API_BASE       = 'https://www.oncokb.org/api/v1'
_API_TOKEN      = None  # set via api_token argument — never hardcode here
_BATCH_ENDPOINT = f'{_API_BASE}/annotate/mutations/byProteinChange'
_BATCH_SIZE     = 50
_RETRY_DELAY    = 2.0
_REQUEST_DELAY  = 0.5
_MAX_RETRIES    = 3

_ONCOGENICITY_INCLUDE = {
    'Oncogenic',
    'Likely Oncogenic',
    'Predicted Oncogenic',
}

_COL_GENE   = 'Hugo Symbol'
_COL_ALT    = 'Alteration'
_COL_EFFECT = 'Mutation Effect'
_COL_ONCO   = 'Oncogenicity'
_COL_CTYPE  = 'Cancer Type'


_AA3TO1 = {
    'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe': 'F',
    'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Lys': 'K', 'Leu': 'L',
    'Met': 'M', 'Asn': 'N', 'Pro': 'P', 'Gln': 'Q', 'Arg': 'R',
    'Ser': 'S', 'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y',
    'Ter': '*', 'Sec': 'U',
}


def _hgvsp_to_oncokb(hgvsp: str) -> str | None:
    """
    Convert an hgvsp string to the single-letter substitution format expected
    by OncoKB (e.g. G12D, W235R).

    Handles:
      ENSP00000484643.1:p.Trp235Arg  -> W235R
      p.Trp235Arg                    -> W235R
      p.Gly12Asp                     -> G12D
      W235R                          -> W235R  (already correct)

    Rejects (returns None) — the OncoKB byProteinChange endpoint interprets
    the query strictly as a substitution, so anything with a suffix (fs, ext,
    del, ins, dup) or a non-substitution class would otherwise silently
    inherit an unrelated variant's oncogenicity:
      p.Ala100Argfs*7                -> None
      p.Trp235TerextTer49            -> None
      p.Lys100_Lys101del             -> None
    """
    s = hgvsp.strip()

    # Strip ENSP prefix
    if ':' in s:
        s = s.split(':', 1)[1].strip()

    # Strip p. prefix
    if s.startswith('p.'):
        s = s[2:]

    # Already single-letter substitution form
    if re.match(r'^[A-Z*]\d+[A-Z*]$', s):
        return s

    # 3-letter form: only accept a clean substitution with an empty suffix.
    m = re.match(r'^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)(.*)$', s)
    if m:
        ref3, pos, alt3, suffix = m.groups()
        if suffix:
            # e.g. fs, ext, del, ins, dup — not a substitution; reject.
            return None
        ref1 = _AA3TO1.get(ref3)
        alt1 = _AA3TO1.get(alt3, alt3) if alt3 != '*' else '*'
        if ref1:
            return f'{ref1}{pos}{alt1}'

    return None


def parse_oncokb(variants_file: str,
                 merged_df=None,
                 include_oncogenicity: list | None = None,
                 api_token: str | None = None) -> pd.DataFrame:
    """
    Annotate variants with OncoKB oncogenicity.

    If variants_file exists on disk, it is parsed directly.
    Otherwise, the public API is queried using (gene, hgvsp) pairs from merged_df.
    """
    if include_oncogenicity is None:
        include_oncogenicity = list(_ONCOGENICITY_INCLUDE)

    if os.path.exists(variants_file):
        log.info('OncoKB: using local file %s', variants_file)
        return _parse_file(variants_file, include_oncogenicity)

    if merged_df is None or merged_df.empty:
        log.warning('OncoKB: no local file and no merged_df provided — skipping')
        return empty_standard_df()

    if not api_token:
        log.warning('OncoKB: no api_token provided — skipping API query')
        return empty_standard_df()
    log.info('OncoKB: querying public API for %d unique variants', len(merged_df))
    return _query_api(merged_df, include_oncogenicity, api_token)


def _parse_file(variants_file: str, include_oncogenicity: list) -> pd.DataFrame:
    try:
        df_raw = pd.read_csv(variants_file, sep='\t', dtype=str, low_memory=False)
    except Exception as e:
        log.error('Failed to read OncoKB file: %s', e)
        return empty_standard_df()

    required = [_COL_GENE, _COL_ALT, _COL_EFFECT, _COL_ONCO]
    missing = [c for c in required if c not in df_raw.columns]
    if missing:
        log.error('OncoKB file missing columns: %s', missing)
        return empty_standard_df()

    df_raw = df_raw.fillna('')
    df_f = df_raw[df_raw[_COL_ONCO].isin(include_oncogenicity)].copy()
    log.info('OncoKB file: %d / %d entries pass filter', len(df_f), len(df_raw))

    return _build_output(
        genes=df_f[_COL_GENE].str.strip().tolist(),
        alts=df_f[_COL_ALT].str.strip().tolist(),
        effects=df_f[_COL_EFFECT].str.strip().tolist(),
        oncogenicities=df_f[_COL_ONCO].tolist(),
        cancer_types=df_f[_COL_CTYPE].str.strip().tolist() if _COL_CTYPE in df_f.columns else [''] * len(df_f),
    )


def _query_api(merged_df: pd.DataFrame, include_oncogenicity: list, api_token: str) -> pd.DataFrame:
    pairs = (
        merged_df[['gene', 'hgvsp']]
        .dropna()
        .drop_duplicates()
        .query("gene != '' and hgvsp != ''")
    )

    if pairs.empty:
        log.warning('OncoKB: no valid (gene, hgvsp) pairs to query')
        return empty_standard_df()

    all_results = []
    total = len(pairs)
    n_batches = (total + _BATCH_SIZE - 1) // _BATCH_SIZE
    log.info('OncoKB API: %d variants in %d batches', total, n_batches)

    for batch_idx in range(n_batches):
        batch = pairs.iloc[batch_idx * _BATCH_SIZE:(batch_idx + 1) * _BATCH_SIZE]
        payload = []
        for i, (_, row) in enumerate(batch.iterrows()):
            alteration = _hgvsp_to_oncokb(str(row['hgvsp']))
            if not alteration:
                continue
            payload.append({
                'id':          str(i),
                'gene':        {'hugoSymbol': str(row['gene']).strip()},
                'alteration':  alteration,
            })
        if not payload:
            continue

        for attempt in range(_MAX_RETRIES):
            try:
                resp = requests.post(
                    _BATCH_ENDPOINT,
                    json=payload,
                    headers={'Content-Type': 'application/json',
                             'Authorization': f'Bearer {api_token}',
                             'Accept': 'application/json'},
                    timeout=30,
                )
                if resp.status_code == 200:
                    all_results.extend(resp.json())
                    break
                elif resp.status_code == 429:
                    wait = _RETRY_DELAY * (attempt + 2)
                    log.warning('OncoKB rate limited — waiting %.1fs', wait)
                    time.sleep(wait)
                elif resp.status_code >= 500:
                    wait = _RETRY_DELAY * (attempt + 1)
                    log.warning('OncoKB batch %d: server error HTTP %d — retrying in %.1fs',
                                batch_idx, resp.status_code, wait)
                    time.sleep(wait)
                else:
                    log.warning('OncoKB batch %d: HTTP %d — %s',
                                batch_idx, resp.status_code, resp.text[:200])
                    break
            except requests.RequestException as e:
                log.warning('OncoKB batch %d attempt %d: %s', batch_idx, attempt + 1, e)
                time.sleep(_RETRY_DELAY)

        if batch_idx % 10 == 0:
            log.info('OncoKB API: processed %d / %d',
                     min((batch_idx + 1) * _BATCH_SIZE, total), total)
        time.sleep(_REQUEST_DELAY)

    if not all_results:
        log.warning('OncoKB API: no results returned')
        return empty_standard_df()

    return _parse_api_results(all_results, include_oncogenicity)


def _parse_api_results(results: list, include_oncogenicity: list) -> pd.DataFrame:
    genes, alts, effects, oncogenicities, cancer_types = [], [], [], [], []

    for item in results:
        onco = item.get('oncogenic', '') or ''
        if onco not in include_oncogenicity:
            continue
        query  = item.get('query', {})
        gene   = query.get('hugoSymbol', '') or ''
        alt    = query.get('alteration', '') or ''
        effect = (item.get('mutationEffect') or {}).get('knownEffect', '') or ''
        genes.append(gene)
        alts.append(alt)
        effects.append(effect)
        oncogenicities.append(onco)
        cancer_types.append('unspecified')

    log.info('OncoKB API: %d / %d pass oncogenicity filter', len(genes), len(results))
    if not genes:
        return empty_standard_df()

    return _build_output(genes, alts, effects, oncogenicities, cancer_types)


def _build_output(genes, alts, effects, oncogenicities, cancer_types) -> pd.DataFrame:
    rows = []
    for gene, alt, effect, onco, ctype in zip(
            genes, alts, effects, oncogenicities, cancer_types):
        rows.append({
            'chrom':       '',
            'pos':         pd.NA,
            'ref':         '',
            'alt':         '',
            'gene':        gene,
            'hgvsc':       '',
            'hgvsp':       _normalise_hgvsp(alt),
            'consequence': _map_effect(effect, alt),
            # Use 'unspecified' (aggregator sentinel) rather than 'pan_cancer'
            # when the row has no cancer-type label; the pan-cancer semantics
            # of OncoKB entries are already carried by the OncoKB source label.
            'cancer_type': (ctype or 'unspecified').lower(),
            'n_samples':   0,
            'source':      f'OncoKB:{onco}',
        })
    df_out = pd.DataFrame(rows, columns=STANDARD_COLS)
    df_out['oncokb_oncogenicity'] = oncogenicities
    log.info('OncoKB: %d entries retained', len(df_out))
    return df_out


def _map_effect(effect: str, alteration: str = '') -> str:
    """
    Derive a standardised consequence term from OncoKB's Mutation Effect plus
    the raw Alteration string.

    OncoKB's Mutation Effect column carries functional effect (Gain-of-function,
    Loss-of-function, Neutral, Inconclusive) — NOT variant class — so labelling
    every gain-/loss-of-function as 'missense' was wrong for truncating,
    fusion, amplification, splice and indel alterations. We inspect the
    Alteration string first for its class, and only fall back to 'missense'
    when the alteration is a clean amino-acid substitution.
    """
    a = alteration.strip()
    e = effect.lower()

    # Class inferred from Alteration string
    if a:
        low = a.lower()
        if 'truncating' in low:
            return 'nonsense'
        if 'amplification' in low:
            return 'structural'
        if 'deletion' in low or 'del' in low.split():
            return 'inframe_indel'
        if 'fusion' in low:
            return 'structural'
        if 'splice' in low:
            return 'splice_site'
        if re.search(r'fs\*?\d*$', a):
            return 'frameshift'
        if a.endswith('del') or 'del' in a:
            return 'inframe_indel'
        if a.endswith('ins') or 'ins' in a:
            return 'inframe_indel'
        if a.endswith('dup') or 'dup' in a:
            return 'inframe_indel'
        # Clean single-substitution: G12D, R175H, *510S
        if re.match(r'^[A-Z*]\d+[A-Z*]$', a):
            return 'missense'
        # 3-letter substitution: Gly12Asp
        if re.match(r'^[A-Z][a-z]{2}\d+([A-Z][a-z]{2}|\*)$', a):
            return 'missense'

    if 'neutral' in e:
        return 'synonymous'
    return 'unknown'


def _normalise_hgvsp(alteration: str) -> str:
    alt = alteration.strip()
    if alt.startswith('p.'):
        return alt
    if re.match(r'^[A-Z]\d+[A-Z*]$', alt):
        return f'p.{alt}'
    return alt
