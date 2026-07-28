#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Shared utilities and constants
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""Shared constants and helper functions used by all parsers."""

import gzip
import logging
import re
import sys
from typing import Optional

import pandas as pd


# ── Standard output schema for every parser ───────────────────────────────────
# All parsers must return a DataFrame with exactly these columns.
STANDARD_COLS = [
    'chrom',        # str  e.g. 'chr1'
    'pos',          # int  1-based
    'ref',          # str  REF allele (normalised, upper-case)
    'alt',          # str  ALT allele (normalised, upper-case)
    'gene',         # str  HGNC gene symbol
    'hgvsc',        # str  HGVSc notation or ''
    'hgvsp',        # str  HGVSp notation or ''
    'consequence',  # str  standardised consequence term (see CONSEQUENCE_MAP)
    'cancer_type',  # str  primary cancer type label
    'n_samples',    # int  number of samples contributing this entry
    'source',       # str  source database label
]

# ── Consequence term mapping ───────────────────────────────────────────────────
# Maps raw strings from each source to standardised internal terms.
CONSEQUENCE_MAP: dict[str, str] = {
    # VEP / Ensembl terms
    'missense_variant':                  'missense',
    'stop_gained':                       'nonsense',
    'splice_donor_variant':              'splice_site',
    'splice_acceptor_variant':           'splice_site',
    'splice_region_variant':             'splice_region',
    'frameshift_variant':                'frameshift',
    'inframe_insertion':                 'inframe_indel',
    'inframe_deletion':                  'inframe_indel',
    'synonymous_variant':                'synonymous',
    'stop_lost':                         'stop_lost',
    'start_lost':                        'start_lost',
    'coding_sequence_variant':           'coding_other',
    'protein_altering_variant':          'missense',
    # COSMIC mutation description terms
    'Substitution - Missense':           'missense',
    'Substitution - Nonsense':           'nonsense',
    'Substitution - synonymous':         'synonymous',
    'Substitution - Silent':             'synonymous',
    'Substitution - coding silent':      'synonymous',
    'Substitution - Coding silent':      'synonymous',
    'Insertion - Frameshift':            'frameshift',
    'Deletion - Frameshift':             'frameshift',
    'Insertion - In frame':              'inframe_indel',
    'Deletion - In frame':               'inframe_indel',
    'Complex - deletion inframe':        'inframe_indel',
    'Complex - insertion inframe':       'inframe_indel',
    'Complex - frameshift':              'frameshift',
    'Complex':                           'unknown',
    'Complex - compound substitution':   'missense',
    'Nonstop extension':                 'stop_lost',
    'Whole gene deletion':               'structural',
    'Splice site':                       'splice_site',
    'Splice - donor':                    'splice_site',
    'Splice - acceptor':                 'splice_site',
    'Unknown':                           'unknown',
    'unknown':                           'unknown',
    # MAF variant classification terms
    'Missense_Mutation':                 'missense',
    'Nonsense_Mutation':                 'nonsense',
    'Splice_Site':                       'splice_site',
    'Splice_Region':                     'splice_region',
    'Frame_Shift_Del':                   'frameshift',
    'Frame_Shift_Ins':                   'frameshift',
    'In_Frame_Del':                      'inframe_indel',
    'In_Frame_Ins':                      'inframe_indel',
    'Silent':                            'synonymous',
    'Stop_Codon_Del':                    'stop_lost',
    'Stop_Codon_Ins':                    'stop_lost',
    'Nonstop_Mutation':                  'stop_lost',
    'Translation_Start_Site':            'start_lost',
    'De_novo_Start_InFrame':             'inframe_indel',
    'De_novo_Start_OutOfFrame':          'frameshift',
    'RNA':                               'noncoding',
    "3'UTR":                             'noncoding',
    "5'UTR":                             'noncoding',
    'Intron':                            'intronic',
    'IGR':                               'intergenic',
    'Targeted_Region':                   'unknown',
    '':                                  'unknown',
}

# Consequences to retain in the whitelist
INCLUDED_CONSEQUENCES: set[str] = {
    'missense',
    'nonsense',
    'splice_site',
    'splice_region',
    'frameshift',
    'inframe_indel',
    'synonymous',
    'stop_lost',
    'start_lost',
    'coding_other',
    'unknown',    # retain unknowns; review manually downstream
}

# Valid nucleotide characters for REF/ALT validation.
# NOTE: '-' is deliberately excluded. MAF-style '-' alleles for indels are not
# valid VCF and must be normalised to an anchor-base representation upstream
# in the parser (or the row must be dropped). Allowing '-' here caused
# malformed VCFs and silent data loss across multiple parsers.
_VALID_ALLELE_CHARS = re.compile(r'^[ACGTNacgtn*]+$')


# ── Chromosome normalisation ───────────────────────────────────────────────────

def normalise_chrom(chrom: str) -> str:
    """
    Normalise a chromosome label to 'chr'-prefixed GRCh38 form.
    Converts '1' -> 'chr1', 'X' -> 'chrX', 'MT' -> 'chrM', etc.
    """
    chrom = str(chrom).strip()
    if chrom.startswith('chr'):
        raw = chrom[3:]
    else:
        raw = chrom
    # Map mitochondrial aliases
    if raw in ('MT', 'mt'):
        raw = 'M'
    return f'chr{raw}'


# ── Allele validation ─────────────────────────────────────────────────────────

def is_valid_allele(allele: str) -> bool:
    return bool(allele) and bool(_VALID_ALLELE_CHARS.match(allele))


def clean_allele(allele: str) -> str:
    return str(allele).strip().upper()


# ── Consequence mapping ────────────────────────────────────────────────────────

def map_consequence(raw: str) -> str:
    """
    Map a raw consequence string to a standardised internal term.
    Returns 'other' for unmapped values; logs a debug message.
    """
    if not raw:
        return 'unknown'
    # Some sources pipe-delimit multiple consequences; take the most severe.
    terms = [t.strip() for t in str(raw).split('&')]
    for t in terms:
        mapped = CONSEQUENCE_MAP.get(t)
        if mapped:
            return mapped
    return 'unknown'


# ── Logging ───────────────────────────────────────────────────────────────────

def setup_logger(name: str, level: str = 'INFO') -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        fmt = logging.Formatter('%(asctime)s [%(name)s] %(levelname)s: %(message)s',
                                datefmt='%Y-%m-%d %H:%M:%S')
        handler.setFormatter(fmt)
        logger.addHandler(handler)
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    return logger


# ── Empty frame factory ────────────────────────────────────────────────────────

def empty_standard_df() -> pd.DataFrame:
    """Return an empty DataFrame with the standard schema."""
    return pd.DataFrame(columns=STANDARD_COLS).astype({
        'chrom':       str,
        'pos':         'Int64',
        'ref':         str,
        'alt':         str,
        'gene':        str,
        'hgvsc':       str,
        'hgvsp':       str,
        'consequence': str,
        'cancer_type': str,
        'n_samples':   'Int64',
        'source':      str,
    })


# ── Simple HGVSc REF/ALT extractor (SNV only) ─────────────────────────────────

_SNV_HGVSC = re.compile(r'c\.\d+([ACGTacgt])>([ACGTacgt])')

def extract_snv_alleles_from_hgvsc(hgvsc: str) -> Optional[tuple[str, str]]:
    """
    Extract (REF, ALT) from a simple HGVSc substitution notation.
    Returns None for indels or unparseable strings.
    Note: alleles are on the coding strand; caller must apply strand correction.
    """
    m = _SNV_HGVSC.search(str(hgvsc))
    if m:
        return m.group(1).upper(), m.group(2).upper()
    return None


# ── Adaptive column resolver ───────────────────────────────────────────────────

def resolve_column(available_cols: list, candidates: list, field_name: str,
                   required: bool = False, log=None) -> str:
    """Return first column from candidates found in available_cols, or '' if missing."""
    for c in candidates:
        if c in available_cols:
            if log:
                log.debug('Field "%s" resolved to column "%s"', field_name, c)
            return c
    msg = f'Field "{field_name}" not found in file. Tried: {candidates}. Available: {available_cols[:20]}'
    if required:
        raise RuntimeError(msg)
    if log:
        log.warning(msg)
    return ''


def detect_separator(path: str) -> str:
    """
    Detect whether a file uses tab or comma as its delimiter.
    Returns '\t' or ','.
    """
    opener = gzip.open if path.endswith('.gz') else open
    try:
        with opener(path, 'rt') as fh:
            header = fh.readline()
        return ',' if header.count(',') > header.count('\t') else '\t'
    except Exception:
        return '\t'


def get_file_columns(path: str, comment: str = '#') -> list:
    """
    Return column names from the first non-comment line of a delimited file.
    Automatically detects separator.
    """
    opener = gzip.open if path.endswith('.gz') else open
    sep = detect_separator(path)
    try:
        with opener(path, 'rt') as fh:
            for line in fh:
                if not line.startswith(comment):
                    return line.rstrip('\n').split(sep)
    except Exception:
        return []
    return []


def detect_chrom_prefix(path: str, chrom_col_idx: int = 0, n: int = 500,
                        comment: str = '#') -> str:
    """
    Sample up to n data rows and return chromosome prefix style.

    Returns one of: 'chr', 'plain', 'mixed', 'unknown'
    """
    opener = gzip.open if path.endswith('.gz') else open
    chr_count   = 0
    plain_count = 0
    seen        = 0
    skip_values = set()

    try:
        with opener(path, 'rt') as fh:
            header = fh.readline()
            # Collect header tokens to skip if they appear in data
            sep = ',' if header.count(',') > header.count('\t') else '\t'
            skip_values = set(header.rstrip('\n').split(sep))
            for line in fh:
                if line.startswith(comment):
                    continue
                parts = line.split(sep)
                if len(parts) <= chrom_col_idx:
                    continue
                chrom = parts[chrom_col_idx].strip()
                if not chrom or chrom in skip_values:
                    continue
                if chrom.startswith('chr'):
                    chr_count += 1
                elif chrom[:1].isdigit() or chrom in ('X', 'Y', 'MT', 'M'):
                    plain_count += 1
                seen += 1
                if seen >= n:
                    break
    except Exception:
        return 'unknown'

    if chr_count > 0 and plain_count == 0:
        return 'chr'
    if plain_count > 0 and chr_count == 0:
        return 'plain'
    if chr_count > 0 and plain_count > 0:
        return 'mixed'
    return 'unknown'
