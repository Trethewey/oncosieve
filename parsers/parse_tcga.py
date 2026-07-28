#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Parser: TCGA PanCancer Atlas mc3 MAF
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
parse_tcga.py
Parse the TCGA PanCancer Atlas mc3 somatic mutation MAF.

Required file:
  mc3.v0.2.8.PUBLIC.GRCh38.maf.gz
  Source: https://gdc.cancer.gov/about-data/publications/pancanatlas
  File:   mc3.v0.2.8.PUBLIC.maf.gz  (md5: 639ad8f8386e98dacc22e439188aa8fa)
  The GRCh37 MAF must be lifted to GRCh38 before use.
  Run: python3 tools/db_fix.py --config config.yaml

The mc3 MAF is the TCGA PanCancer Atlas harmonised somatic mutation call set,
produced by merging calls from 6 variant callers across 10,295 tumour-normal
pairs from 33 cancer types.

Key columns used:
  Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2
  Hugo_Symbol, HGVSc, HGVSp_Short, Transcript_ID
  Consequence      — VEP consequence term (used directly, no remapping needed)
  Tumor_Sample_Barcode — encodes TCGA project code for cancer type
  FILTER           — only PASS variants are included

Cancer type is derived from the Tumor_Sample_Barcode, which follows the
mc3 format TCGA-{TSS}-{Participant}-{Sample}-... where TSS is a two-character
Tissue Source Site code (e.g. '02' -> GBM -> 'glioblastoma multiforme').
Uses GDC's tissueSourceSite table (embedded as _TCGA_TSS_MAP) to go
TSS -> study abbreviation, then _TCGA_PROJECT_MAP to go study -> label.
See _cancer_type_from_barcode() for the full contract.
"""

import gzip
import os

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

log = setup_logger('TCGA')

# Column names in the mc3 MAF
_COL_CHROM    = 'Chromosome'
_COL_POS      = 'Start_Position'
_COL_REF      = 'Reference_Allele'
_COL_ALT      = 'Tumor_Seq_Allele2'
_COL_GENE     = 'Hugo_Symbol'
_COL_HGVSC    = 'HGVSc'
_COL_HGVSP    = 'HGVSp_Short'
_COL_CSQ      = 'Consequence'       # VEP term — preferred over Variant_Classification
_COL_VARCLASS = 'Variant_Classification'  # fallback if Consequence absent
_COL_SAMPLE   = 'Tumor_Sample_Barcode'
_COL_FILTER   = 'FILTER'

# TCGA project code -> cancer type label
# Source: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
_TCGA_PROJECT_MAP: dict[str, str] = {
    'ACC':  'adrenocortical carcinoma',
    'BLCA': 'bladder urothelial carcinoma',
    'BRCA': 'breast invasive carcinoma',
    'CESC': 'cervical squamous cell carcinoma',
    'CHOL': 'cholangiocarcinoma',
    'COAD': 'colon adenocarcinoma',
    'DLBC': 'diffuse large b-cell lymphoma',
    'ESCA': 'esophageal carcinoma',
    'GBM':  'glioblastoma multiforme',
    'HNSC': 'head and neck squamous cell carcinoma',
    'KICH': 'kidney chromophobe',
    'KIRC': 'kidney renal clear cell carcinoma',
    'KIRP': 'kidney renal papillary cell carcinoma',
    'LAML': 'acute myeloid leukemia',
    'LGG':  'brain lower grade glioma',
    'LIHC': 'liver hepatocellular carcinoma',
    'LUAD': 'lung adenocarcinoma',
    'LUSC': 'lung squamous cell carcinoma',
    'MESO': 'mesothelioma',
    'OV':   'ovarian serous cystadenocarcinoma',
    'PAAD': 'pancreatic adenocarcinoma',
    'PCPG': 'pheochromocytoma and paraganglioma',
    'PRAD': 'prostate adenocarcinoma',
    'READ': 'rectum adenocarcinoma',
    'SARC': 'sarcoma',
    'SKCM': 'skin cutaneous melanoma',
    'STAD': 'stomach adenocarcinoma',
    'TGCT': 'testicular germ cell tumors',
    'THCA': 'thyroid carcinoma',
    'THYM': 'thymoma',
    'UCEC': 'uterine corpus endometrial carcinoma',
    'UCS':  'uterine carcinosarcoma',
    'UVM':  'uveal melanoma',
    # Rare / historical / administrative TCGA study codes referenced by
    # some tissue source sites in _TCGA_TSS_MAP.
    'CNTL':     'controls',
    'COADREAD': 'colorectal adenocarcinoma',
    'FPPP':     'formalin fixed paraffin-embedded pilot phase ii',
    'LCML':     'chronic myelogenous leukemia',
    'MISC':     'miscellaneous',
    'STES':     'stomach and esophageal carcinoma',
}


# TSS (Tissue Source Site) code -> TCGA study abbreviation.
# Sourced from GDC's public tissue-source-site table:
#   https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes
# When a TSS contributed to more than one study, the first (primary) study is used.
_TCGA_TSS_MAP: dict[str, str] = {
    '01': 'OV',
    '02': 'GBM',
    '04': 'OV',
    '05': 'LUAD',
    '06': 'GBM',
    '07': 'CNTL',
    '08': 'GBM',
    '09': 'OV',
    '10': 'OV',
    '11': 'LUSC',
    '12': 'GBM',
    '13': 'OV',
    '14': 'GBM',
    '15': 'GBM',
    '16': 'GBM',
    '17': 'LUAD',
    '18': 'LUSC',
    '19': 'GBM',
    '1Z': 'THYM',
    '20': 'OV',
    '21': 'LUSC',
    '22': 'LUSC',
    '23': 'OV',
    '24': 'OV',
    '25': 'OV',
    '26': 'GBM',
    '27': 'GBM',
    '28': 'GBM',
    '29': 'OV',
    '2A': 'PRAD',
    '2E': 'UCEC',
    '2F': 'BLCA',
    '2G': 'TGCT',
    '2H': 'ESCA',
    '2J': 'PAAD',
    '2K': 'KIRP',
    '2L': 'PAAD',
    '2M': 'ESCA',
    '2N': 'STAD',
    '2P': 'PAAD',
    '2V': 'LIHC',
    '2W': 'CESC',
    '2X': 'TGCT',
    '2Y': 'LIHC',
    '2Z': 'KIRP',
    '30': 'OV',
    '31': 'OV',
    '32': 'GBM',
    '33': 'LUSC',
    '34': 'LUSC',
    '35': 'LUAD',
    '36': 'OV',
    '37': 'LUSC',
    '38': 'LUAD',
    '39': 'LUSC',
    '3A': 'PAAD',
    '3B': 'SARC',
    '3C': 'BRCA',
    '3E': 'PAAD',
    '3G': 'THYM',
    '3H': 'MESO',
    '3J': 'BRCA',
    '3K': 'LIHC',
    '3L': 'COAD',
    '3M': 'STAD',
    '3N': 'SKCM',
    '3P': 'OV',
    '3Q': 'THYM',
    '3R': 'SARC',
    '3S': 'THYM',
    '3T': 'THYM',
    '3U': 'MESO',
    '3W': 'SARC',
    '3X': 'CHOL',
    '3Z': 'KIRC',
    '41': 'GBM',
    '42': 'OV',
    '43': 'LUSC',
    '44': 'LUAD',
    '46': 'LUSC',
    '49': 'LUAD',
    '4A': 'KIRP',
    '4B': 'LUAD',
    '4C': 'THCA',
    '4D': 'OV',
    '4E': 'UCEC',
    '4G': 'CHOL',
    '4H': 'BRCA',
    '4J': 'CESC',
    '4K': 'TGCT',
    '4L': 'PRAD',
    '4N': 'COAD',
    '4P': 'HNSC',
    '4Q': 'SARC',
    '4R': 'LIHC',
    '4S': 'PRAD',
    '4T': 'COAD',
    '4V': 'THYM',
    '4W': 'GBM',
    '4X': 'THYM',
    '4Y': 'SARC',
    '4Z': 'BLCA',
    '50': 'LUAD',
    '51': 'LUSC',
    '52': 'LUSC',
    '53': 'LUAD',
    '55': 'LUAD',
    '56': 'LUSC',
    '57': 'OV',
    '58': 'LUSC',
    '59': 'OV',
    '5A': 'CHOL',
    '5B': 'UCEC',
    '5C': 'LIHC',
    '5D': 'SARC',
    '5F': 'THCA',
    '5G': 'THYM',
    '5H': 'UVM',
    '5J': 'LAML',
    '5K': 'THYM',
    '5L': 'BRCA',
    '5M': 'COAD',
    '5N': 'BLCA',
    '5P': 'KIRP',
    '5Q': 'PAAD',
    '5R': 'LIHC',
    '5S': 'UCEC',
    '5T': 'BRCA',
    '5U': 'THYM',
    '5V': 'THYM',
    '5W': 'UCEC',
    '5X': 'OV',
    '60': 'LUSC',
    '61': 'OV',
    '62': 'LUAD',
    '63': 'LUSC',
    '64': 'LUAD',
    '65': 'GBM',
    '66': 'LUSC',
    '67': 'LUAD',
    '68': 'LUSC',
    '69': 'LUAD',
    '6A': 'LUSC',
    '6D': 'KIRC',
    '6G': 'READ',
    '6H': 'LCML',
    '70': 'LUSC',
    '71': 'LUAD',
    '72': 'OV',
    '73': 'LUAD',
    '74': 'GBM',
    '75': 'LUAD',
    '76': 'GBM',
    '77': 'LUSC',
    '78': 'LUAD',
    '79': 'LUSC',
    '80': 'LUAD',
    '81': 'GBM',
    '82': 'LUSC',
    '83': 'LUAD',
    '85': 'LUSC',
    '86': 'LUAD',
    '87': 'GBM',
    '90': 'LUSC',
    '91': 'LUAD',
    '92': 'LUSC',
    '93': 'LUAD',
    '94': 'LUSC',
    '95': 'LUAD',
    '96': 'LUSC',
    '97': 'LUAD',
    '98': 'LUSC',
    '99': 'LUAD',
    'A1': 'BRCA',
    'A2': 'BRCA',
    'A3': 'KIRC',
    'A4': 'KIRP',
    'A5': 'UCEC',
    'A6': 'COAD',
    'A7': 'BRCA',
    'A8': 'BRCA',
    'AA': 'COAD',
    'AB': 'LAML',
    'AC': 'BRCA',
    'AD': 'COAD',
    'AF': 'READ',
    'AG': 'READ',
    'AH': 'READ',
    'AJ': 'UCEC',
    'AK': 'KIRC',
    'AL': 'KIRP',
    'AM': 'COAD',
    'AN': 'BRCA',
    'AO': 'BRCA',
    'AP': 'UCEC',
    'AQ': 'BRCA',
    'AR': 'BRCA',
    'AS': 'KIRC',
    'AT': 'KIRP',
    'AU': 'COAD',
    'AV': 'CNTL',
    'AW': 'UCEC',
    'AX': 'UCEC',
    'AY': 'COAD',
    'AZ': 'COAD',
    'B0': 'KIRC',
    'B1': 'KIRP',
    'B2': 'KIRC',
    'B3': 'KIRP',
    'B4': 'KIRC',
    'B5': 'UCEC',
    'B6': 'BRCA',
    'B7': 'STAD',
    'B8': 'KIRC',
    'B9': 'KIRP',
    'BA': 'HNSC',
    'BB': 'HNSC',
    'BC': 'LIHC',
    'BD': 'LIHC',
    'BF': 'SKCM',
    'BG': 'UCEC',
    'BH': 'BRCA',
    'BI': 'CESC',
    'BJ': 'THCA',
    'BK': 'UCEC',
    'BL': 'BLCA',
    'BM': 'READ',
    'BP': 'KIRC',
    'BQ': 'KIRP',
    'BR': 'STAD',
    'BS': 'UCEC',
    'BT': 'BLCA',
    'BW': 'LIHC',
    'C4': 'BLCA',
    'C5': 'CESC',
    'C8': 'BRCA',
    'C9': 'HNSC',
    'CA': 'COAD',
    'CB': 'KIRC',
    'CC': 'LIHC',
    'CD': 'STAD',
    'CE': 'THCA',
    'CF': 'BLCA',
    'CG': 'STAD',
    'CH': 'PRAD',
    'CI': 'READ',
    'CJ': 'KIRC',
    'CK': 'COAD',
    'CL': 'READ',
    'CM': 'COAD',
    'CN': 'HNSC',
    'CQ': 'HNSC',
    'CR': 'HNSC',
    'CS': 'LGG',
    'CU': 'BLCA',
    'CV': 'HNSC',
    'CW': 'KIRC',
    'CX': 'HNSC',
    'CZ': 'KIRC',
    'D1': 'UCEC',
    'D3': 'SKCM',
    'D5': 'COAD',
    'D6': 'HNSC',
    'D7': 'STAD',
    'D8': 'BRCA',
    'D9': 'SKCM',
    'DA': 'SKCM',
    'DB': 'LGG',
    'DC': 'READ',
    'DD': 'LIHC',
    'DE': 'THCA',
    'DF': 'UCEC',
    'DG': 'CESC',
    'DH': 'LGG',
    'DI': 'UCEC',
    'DJ': 'THCA',
    'DK': 'BLCA',
    'DM': 'COAD',
    'DO': 'THCA',
    'DQ': 'HNSC',
    'DR': 'CESC',
    'DS': 'CESC',
    'DT': 'READ',
    'DU': 'LGG',
    'DV': 'KIRC',
    'DW': 'KIRP',
    'DX': 'SARC',
    'DY': 'READ',
    'DZ': 'KIRP',
    'E1': 'LGG',
    'E2': 'BRCA',
    'E3': 'THCA',
    'E5': 'BLCA',
    'E6': 'UCEC',
    'E7': 'BLCA',
    'E8': 'THCA',
    'E9': 'BRCA',
    'EA': 'CESC',
    'EB': 'SKCM',
    'EC': 'UCEC',
    'ED': 'LIHC',
    'EE': 'SKCM',
    'EF': 'READ',
    'EI': 'READ',
    'EJ': 'PRAD',
    'EK': 'CESC',
    'EL': 'THCA',
    'EM': 'THCA',
    'EO': 'UCEC',
    'EP': 'LIHC',
    'EQ': 'STAD',
    'ER': 'SKCM',
    'ES': 'LIHC',
    'ET': 'THCA',
    'EU': 'KIRC',
    'EV': 'KIRP',
    'EW': 'BRCA',
    'EX': 'CESC',
    'EY': 'UCEC',
    'EZ': 'LGG',
    'F1': 'STAD',
    'F2': 'PAAD',
    'F4': 'COAD',
    'F5': 'READ',
    'F6': 'LGG',
    'F7': 'HNSC',
    'F9': 'KIRP',
    'FA': 'DLBC',
    'FB': 'PAAD',
    'FC': 'PRAD',
    'FD': 'BLCA',
    'FE': 'THCA',
    'FF': 'DLBC',
    'FG': 'LGG',
    'FH': 'THCA',
    'FI': 'UCEC',
    'FJ': 'BLCA',
    'FK': 'THCA',
    'FL': 'UCEC',
    'FM': 'DLBC',
    'FN': 'LGG',
    'FP': 'STAD',
    'FQ': 'PAAD',
    'FR': 'SKCM',
    'FS': 'SKCM',
    'FT': 'BLCA',
    'FU': 'CESC',
    'FV': 'LIHC',
    'FW': 'SKCM',
    'FX': 'SARC',
    'FY': 'THCA',
    'FZ': 'PAAD',
    'G2': 'BLCA',
    'G3': 'LIHC',
    'G4': 'COAD',
    'G5': 'READ',
    'G6': 'KIRC',
    'G7': 'KIRP',
    'G8': 'DLBC',
    'G9': 'PRAD',
    'GC': 'BLCA',
    'GD': 'BLCA',
    'GE': 'THCA',
    'GF': 'SKCM',
    'GG': 'UCEC',
    'GH': 'CESC',
    'GI': 'BRCA',
    'GJ': 'LIHC',
    'GK': 'KIRC',
    'GL': 'KIRP',
    'GM': 'BRCA',
    'GN': 'SKCM',
    'GP': 'LAML',
    'GR': 'DLBC',
    'GS': 'DLBC',
    'GU': 'BLCA',
    'GV': 'BLCA',
    'GZ': 'DLBC',
    'H1': 'STAD',
    'H2': 'THCA',
    'H3': 'DLBC',
    'H4': 'BLCA',
    'H5': 'UCEC',
    'H6': 'PAAD',
    'H7': 'HNSC',
    'H8': 'PAAD',
    'H9': 'PRAD',
    'HA': 'STAD',
    'HB': 'SARC',
    'HC': 'PRAD',
    'HD': 'HNSC',
    'HE': 'KIRP',
    'HF': 'STAD',
    'HG': 'CESC',
    'HH': 'STAD',
    'HI': 'PRAD',
    'HJ': 'STAD',
    'HK': 'LGG',
    'HL': 'HNSC',
    'HM': 'CESC',
    'HN': 'BRCA',
    'HP': 'LIHC',
    'HQ': 'BLCA',
    'HR': 'SKCM',
    'HS': 'SARC',
    'HT': 'LGG',
    'HU': 'STAD',
    'HV': 'PAAD',
    'HW': 'LGG',
    'HZ': 'PAAD',
    'IA': 'KIRP',
    'IB': 'PAAD',
    'IC': 'ESCA',
    'IE': 'SARC',
    'IF': 'SARC',
    'IG': 'ESCA',
    'IH': 'SKCM',
    'IJ': 'LAML',
    'IK': 'LGG',
    'IM': 'THCA',
    'IN': 'STAD',
    'IP': 'STAD',
    'IQ': 'HNSC',
    'IR': 'CESC',
    'IS': 'SARC',
    'IW': 'SARC',
    'IZ': 'KIRP',
    'J1': 'LUSC',
    'J2': 'LUAD',
    'J4': 'PRAD',
    'J7': 'KIRP',
    'J8': 'THCA',
    'J9': 'PRAD',
    'JA': 'HNSC',
    'JL': 'BRCA',
    'JU': 'UCEC',
    'JV': 'SARC',
    'JW': 'CESC',
    'JX': 'CESC',
    'JY': 'ESCA',
    'JZ': 'ESCA',
    'K1': 'SARC',
    'K4': 'BLCA',
    'K6': 'UCEC',
    'K7': 'LIHC',
    'K8': 'SKCM',
    'KA': 'ESCA',
    'KB': 'STAD',
    'KC': 'PRAD',
    'KD': 'SARC',
    'KE': 'UCEC',
    'KF': 'SARC',
    'KG': 'PAAD',
    'KH': 'ESCA',
    'KJ': 'UCEC',
    'KK': 'PRAD',
    'KL': 'KICH',
    'KM': 'KICH',
    'KN': 'KICH',
    'KO': 'KICH',
    'KP': 'UCEC',
    'KQ': 'BLCA',
    'KR': 'LIHC',
    'KS': 'THCA',
    'KT': 'LGG',
    'KU': 'HNSC',
    'KV': 'KIRP',
    'KZ': 'STAD',
    'L1': 'PAAD',
    'L3': 'LUSC',
    'L4': 'LUAD',
    'L5': 'ESCA',
    'L6': 'THCA',
    'L7': 'ESCA',
    'L8': 'KIRP',
    'L9': 'LUAD',
    'LA': 'LUSC',
    'LB': 'PAAD',
    'LC': 'BLCA',
    'LD': 'BRCA',
    'LG': 'LIHC',
    'LH': 'SKCM',
    'LI': 'SARC',
    'LK': 'MESO',
    'LL': 'BRCA',
    'LN': 'ESCA',
    'LP': 'CESC',
    'LQ': 'BRCA',
    'LS': 'CESC',
    'LT': 'BLCA',
    'M7': 'PRAD',
    'M8': 'PAAD',
    'M9': 'ESCA',
    'MA': 'CESC',
    'MB': 'SARC',
    'ME': 'LUAD',
    'MF': 'LUSC',
    'MG': 'PRAD',
    'MH': 'KIRP',
    'MI': 'LIHC',
    'MJ': 'SARC',
    'MK': 'THCA',
    'ML': 'LUSC',
    'MM': 'KIRC',
    'MN': 'LUAD',
    'MO': 'SARC',
    'MP': 'LUAD',
    'MQ': 'MESO',
    'MR': 'LIHC',
    'MS': 'BRCA',
    'MT': 'HNSC',
    'MU': 'CESC',
    'MV': 'BLCA',
    'MW': 'KIRC',
    'MX': 'STAD',
    'MY': 'CESC',
    'MZ': 'HNSC',
    'N1': 'SARC',
    'N5': 'UCS',
    'N6': 'UCS',
    'N7': 'UCS',
    'N8': 'UCS',
    'N9': 'UCS',
    'NA': 'UCS',
    'NB': 'LUAD',
    'NC': 'LUSC',
    'ND': 'UCS',
    'NF': 'UCS',
    'NG': 'UCS',
    'NH': 'COAD',
    'NI': 'LIHC',
    'NJ': 'LUAD',
    'NK': 'LUSC',
    'NM': 'HNSC',
    'NP': 'KICH',
    'NQ': 'MESO',
    'NS': 'SKCM',
    'O1': 'LUAD',
    'O2': 'LUSC',
    'O8': 'LIHC',
    'O9': 'KIRP',
    'OC': 'LUSC',
    'OD': 'SKCM',
    'OE': 'PAAD',
    'OJ': 'THCA',
    'OK': 'BRCA',
    'OL': 'BRCA',
    'OR': 'ACC',
    'OU': 'ACC',
    'OW': 'MISC',
    'OX': 'GBM',
    'OY': 'OV',
    'P3': 'HNSC',
    'P4': 'KIRP',
    'P5': 'LGG',
    'P6': 'ACC',
    'P7': 'PCPG',
    'P8': 'PCPG',
    'P9': 'PAAD',
    'PA': 'ACC',
    'PB': 'DLBC',
    'PC': 'SARC',
    'PD': 'LIHC',
    'PE': 'BRCA',
    'PG': 'UCEC',
    'PH': 'LAML',
    'PJ': 'KIRP',
    'PK': 'ACC',
    'PL': 'BRCA',
    'PN': 'CESC',
    'PQ': 'BLCA',
    'PR': 'PCPG',
    'PT': 'SARC',
    'PZ': 'PAAD',
    'Q1': 'CESC',
    'Q2': 'KIRP',
    'Q3': 'PAAD',
    'Q4': 'LAML',
    'Q9': 'ESCA',
    'QA': 'LIHC',
    'QB': 'SKCM',
    'QC': 'SARC',
    'QD': 'THCA',
    'QF': 'UCEC',
    'QG': 'COAD',
    'QH': 'LGG',
    'QJ': 'OV',
    'QK': 'HNSC',
    'QL': 'COAD',
    'QM': 'UCS',
    'QN': 'UCS',
    'QQ': 'SARC',
    'QR': 'PCPG',
    'QS': 'UCEC',
    'QT': 'PCPG',
    'QU': 'PRAD',
    'QV': 'CESC',
    'QW': 'STAD',
    'R1': 'COAD',
    'R2': 'CESC',
    'R3': 'BLCA',
    'R5': 'STAD',
    'R6': 'ESCA',
    'R7': 'HNSC',
    'R8': 'LGG',
    'R9': 'OV',
    'RA': 'CESC',
    'RB': 'PAAD',
    'RC': 'LIHC',
    'RD': 'STAD',
    'RE': 'ESCA',
    'RG': 'LIHC',
    'RH': 'HNSC',
    'RL': 'PAAD',
    'RM': 'PCPG',
    'RN': 'SARC',
    'RP': 'SKCM',
    'RQ': 'DLBC',
    'RR': 'GBM',
    'RS': 'HNSC',
    'RT': 'PCPG',
    'RU': 'COAD',
    'RV': 'PAAD',
    'RW': 'PCPG',
    'RX': 'PCPG',
    'RY': 'LGG',
    'RZ': 'UVM',
    'S2': 'LUAD',
    'S3': 'BRCA',
    'S4': 'PAAD',
    'S5': 'BLCA',
    'S6': 'TGCT',
    'S7': 'PCPG',
    'S8': 'ESCA',
    'S9': 'LGG',
    'SA': 'PCPG',
    'SB': 'TGCT',
    'SC': 'MESO',
    'SD': 'PAAD',
    'SE': 'PCPG',
    'SG': 'SARC',
    'SH': 'MESO',
    'SI': 'SARC',
    'SJ': 'UCEC',
    'SK': 'COAD',
    'SL': 'UCEC',
    'SN': 'TGCT',
    'SO': 'TGCT',
    'SP': 'PCPG',
    'SQ': 'PCPG',
    'SR': 'PCPG',
    'SS': 'COAD',
    'ST': 'HNSC',
    'SU': 'PRAD',
    'SW': 'STAD',
    'SX': 'KIRP',
    'SY': 'BLCA',
    'T1': 'LIHC',
    'T2': 'HNSC',
    'T3': 'HNSC',
    'T6': 'LUAD',
    'T7': 'KIRC',
    'T9': 'COAD',
    'TE': 'SKCM',
    'TG': 'HNSC',
    'TK': 'PRAD',
    'TL': 'STAD',
    'TM': 'LGG',
    'TN': 'HNSC',
    'TP': 'PRAD',
    'TQ': 'LGG',
    'TR': 'SKCM',
    'TS': 'MESO',
    'TT': 'PCPG',
    'TV': 'BRCA',
    'UB': 'LIHC',
    'UC': 'CESC',
    'UD': 'MESO',
    'UE': 'SARC',
    'UF': 'HNSC',
    'UJ': 'LUSC',
    'UL': 'BRCA',
    'UN': 'KIRP',
    'UP': 'HNSC',
    'UR': 'PRAD',
    'US': 'PAAD',
    'UT': 'MESO',
    'UU': 'BRCA',
    'UV': 'LIHC',
    'UW': 'KICH',
    'UY': 'BLCA',
    'UZ': 'KIRP',
    'V1': 'PRAD',
    'V2': 'PRAD',
    'V3': 'UVM',
    'V4': 'UVM',
    'V5': 'ESCA',
    'V6': 'STAD',
    'V7': 'BRCA',
    'V8': 'KIRC',
    'V9': 'KIRP',
    'VA': 'STAD',
    'VB': 'DLBC',
    'VD': 'UVM',
    'VF': 'TGCT',
    'VG': 'OV',
    'VK': 'COAD',
    'VL': 'READ',
    'VM': 'LGG',
    'VN': 'PRAD',
    'VP': 'PRAD',
    'VQ': 'STAD',
    'VR': 'ESCA',
    'VS': 'CESC',
    'VT': 'SARC',
    'VV': 'LGG',
    'VW': 'LGG',
    'VX': 'STAD',
    'VZ': 'PCPG',
    'W2': 'PCPG',
    'W3': 'SKCM',
    'W4': 'TGCT',
    'W5': 'CHOL',
    'W6': 'CHOL',
    'W7': 'CHOL',
    'W8': 'BRCA',
    'W9': 'LGG',
    'WA': 'HNSC',
    'WB': 'PCPG',
    'WC': 'UVM',
    'WD': 'CHOL',
    'WE': 'SKCM',
    'WF': 'PAAD',
    'WG': 'LUSC',
    'WH': 'LGG',
    'WJ': 'LIHC',
    'WK': 'SARC',
    'WL': 'CESC',
    'WM': 'KIRC',
    'WN': 'KIRP',
    'WP': 'SARC',
    'WQ': 'LIHC',
    'WR': 'OV',
    'WS': 'COAD',
    'WT': 'BRCA',
    'WU': 'COAD',
    'WW': 'PRAD',
    'WX': 'LIHC',
    'WY': 'LGG',
    'WZ': 'TGCT',
    'X2': 'SARC',
    'X3': 'TGCT',
    'X4': 'PRAD',
    'X5': 'BLCA',
    'X6': 'SARC',
    'X7': 'THYM',
    'X8': 'ESCA',
    'X9': 'SARC',
    'XA': 'PRAD',
    'XB': 'ESCA',
    'XC': 'LUSC',
    'XD': 'PAAD',
    'XE': 'TGCT',
    'XF': 'BLCA',
    'XG': 'PCPG',
    'XH': 'THYM',
    'XJ': 'PRAD',
    'XK': 'PRAD',
    'XM': 'THYM',
    'XN': 'PAAD',
    'XP': 'ESCA',
    'XQ': 'PRAD',
    'XR': 'LIHC',
    'XS': 'CESC',
    'XT': 'MESO',
    'XU': 'THYM',
    'XV': 'SKCM',
    'XX': 'BRCA',
    'XY': 'TGCT',
    'Y3': 'LAML',
    'Y5': 'SARC',
    'Y6': 'PRAD',
    'Y8': 'KIRP',
    'YA': 'LIHC',
    'YB': 'PAAD',
    'YC': 'BLCA',
    'YD': 'SKCM',
    'YF': 'BLCA',
    'YG': 'SKCM',
    'YH': 'PAAD',
    'YJ': 'PRAD',
    'YL': 'PRAD',
    'YN': 'SKCM',
    'YR': 'CHOL',
    'YS': 'MESO',
    'YT': 'THYM',
    'YU': 'TGCT',
    'YV': 'UVM',
    'YW': 'SARC',
    'YX': 'STAD',
    'YY': 'PAAD',
    'YZ': 'UVM',
    'Z2': 'SKCM',
    'Z3': 'SARC',
    'Z4': 'SARC',
    'Z5': 'PAAD',
    'Z6': 'ESCA',
    'Z7': 'BRCA',
    'Z8': 'PAAD',
    'ZA': 'STAD',
    'ZB': 'THYM',
    'ZC': 'THYM',
    'ZD': 'CHOL',
    'ZE': 'LUSC',
    'ZF': 'BLCA',
    'ZG': 'PRAD',
    'ZH': 'CHOL',
    'ZJ': 'CESC',
    'ZK': 'CHOL',
    'ZL': 'THYM',
    'ZM': 'TGCT',
    'ZN': 'MESO',
    'ZP': 'LIHC',
    'ZQ': 'STAD',
    'ZR': 'ESCA',
    'ZS': 'LIHC',
    'ZT': 'THYM',
    'ZU': 'CHOL',
    'ZW': 'PAAD',
    'ZX': 'CESC',
}


def _cancer_type_from_barcode(barcode: str) -> str:
    """
    Return the cancer type label for a mc3 tumour sample barcode.

    mc3 barcodes follow the format TCGA-{TSS}-{Participant}-{Sample}-... where
    {TSS} is a two-character Tissue Source Site code (e.g. '02', 'A2', '1Z').
    This is NOT a project abbreviation — the naive `parts[1]` lookup against
    _TCGA_PROJECT_MAP never matched and every TCGA row was silently annotated
    with a junk 'TCGA_{TSS}' label that in turn inflated n_cancer_types and
    over-promoted variants past tier thresholds.

    Lookup chain: TSS -> study abbreviation (via _TCGA_TSS_MAP, sourced from
    GDC's public tissue-source-site table) -> cancer type label (via
    _TCGA_PROJECT_MAP). Returns 'unspecified' (the aggregator sentinel) when
    either step misses so unresolved rows don't count towards n_cancer_types.
    """
    parts = str(barcode).split('-')
    if len(parts) < 2:
        return 'unspecified'
    tss = parts[1].strip().upper()
    study = _TCGA_TSS_MAP.get(tss)
    if not study:
        return 'unspecified'
    return _TCGA_PROJECT_MAP.get(study, 'unspecified')


def parse_tcga(maf_path: str) -> pd.DataFrame:
    """
    Parse TCGA mc3 PanCancer Atlas MAF (GRCh38 lifted).

    Parameters
    ----------
    maf_path : path to mc3.v0.2.8.PUBLIC.GRCh38.maf.gz

    Returns
    -------
    DataFrame with STANDARD_COLS schema.
    One row per mutation per sample. Aggregation to unique samples per variant
    is performed at the merge step.
    """
    if not os.path.exists(maf_path):
        log.warning('TCGA MAF not found: %s  — skipping', maf_path)
        return empty_standard_df()

    log.info('Parsing TCGA mc3 MAF: %s', maf_path)

    opener = gzip.open if maf_path.endswith('.gz') else open

    # Read header to get column indices
    col_idx: dict[str, int] = {}
    try:
        with opener(maf_path, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                col_idx = {c: i for i, c in enumerate(cols)}
                break
    except Exception as e:
        log.error('Failed to read TCGA MAF header: %s', e)
        return empty_standard_df()

    required = [_COL_CHROM, _COL_POS, _COL_REF, _COL_ALT,
                _COL_GENE, _COL_SAMPLE]
    missing = [c for c in required if c not in col_idx]
    if missing:
        log.error('TCGA MAF missing required columns: %s', missing)
        return empty_standard_df()

    chrom_i    = col_idx[_COL_CHROM]
    pos_i      = col_idx[_COL_POS]
    ref_i      = col_idx[_COL_REF]
    alt_i      = col_idx[_COL_ALT]
    gene_i     = col_idx[_COL_GENE]
    hgvsc_i    = col_idx.get(_COL_HGVSC)
    hgvsp_i    = col_idx.get(_COL_HGVSP)
    csq_i      = col_idx.get(_COL_CSQ)
    varclass_i = col_idx.get(_COL_VARCLASS)
    sample_i   = col_idx[_COL_SAMPLE]
    filter_i   = col_idx.get(_COL_FILTER)
    # Transcript_ID is used to prepend a transcript prefix to bare c. HGVSc
    # strings so the downstream MANE Select lookup can resolve them.
    tid_i      = col_idx.get('Transcript_ID')

    rows = []
    n_filtered = 0

    try:
        with opener(maf_path, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')

                # Skip header row
                if parts[chrom_i] == _COL_CHROM:
                    continue

                # FILTER column: keep PASS only
                if filter_i is not None and filter_i < len(parts):
                    if parts[filter_i].strip().upper() not in ('PASS', '.', ''):
                        n_filtered += 1
                        continue

                try:
                    chrom = normalise_chrom(parts[chrom_i])
                    pos   = int(float(parts[pos_i]))
                    ref   = clean_allele(parts[ref_i])
                    alt   = clean_allele(parts[alt_i])

                    if not is_valid_allele(ref) or not is_valid_allele(alt):
                        continue
                    if alt in ('', '-', '.'):
                        continue

                    gene   = str(parts[gene_i]).strip()
                    hgvsc  = str(parts[hgvsc_i]).strip()  if hgvsc_i  is not None and hgvsc_i  < len(parts) else ''
                    hgvsp  = str(parts[hgvsp_i]).strip()  if hgvsp_i  is not None and hgvsp_i  < len(parts) else ''
                    sample = str(parts[sample_i]).strip()

                    # Prepend Transcript_ID to bare c. notation so the MANE
                    # Select lookup downstream can resolve the transcript.
                    if hgvsc.startswith('c.') and tid_i is not None and tid_i < len(parts):
                        tid = str(parts[tid_i]).strip()
                        if tid:
                            hgvsc = f'{tid}:{hgvsc}'

                    # Consequence: prefer VEP Consequence column, fall back to
                    # Variant_Classification mapped through map_consequence()
                    if csq_i is not None and csq_i < len(parts):
                        raw_csq = parts[csq_i].strip()
                        # VEP may return pipe-separated terms; take first
                        consequence = map_consequence(raw_csq.split('&')[0].split('|')[0])
                    elif varclass_i is not None and varclass_i < len(parts):
                        consequence = map_consequence(parts[varclass_i].strip())
                    else:
                        consequence = 'unknown'

                    cancer_type = _cancer_type_from_barcode(sample)

                    rows.append({
                        'chrom':       chrom,
                        'pos':         pos,
                        'ref':         ref,
                        'alt':         alt,
                        'gene':        gene,
                        'hgvsc':       hgvsc,
                        'hgvsp':       hgvsp,
                        'consequence': consequence,
                        'cancer_type': cancer_type,
                        'n_samples':   1,
                        'source':      'TCGA',
                    })

                except (ValueError, TypeError, IndexError):
                    continue

    except Exception as e:
        log.error('Failed to parse TCGA MAF: %s', e)
        return empty_standard_df()

    log.info(
        'TCGA: %d rows parsed, %d excluded by FILTER',
        len(rows), n_filtered
    )

    if not rows:
        log.warning('TCGA: no rows produced — check MAF path and FILTER column')
        return empty_standard_df()

    df = pd.DataFrame(rows, columns=STANDARD_COLS)
    log.info('TCGA: %d rows', len(df))
    return df
