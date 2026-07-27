<div align="center">
  <img src="assets/oncosieve-grid-redblue-text.svg" width="550" alt="OncoSieve logo"/>

  [![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

  <p><strong>Pan-Cancer Somatic Variant Whitelist Curation Tool</strong></p>
  <p><em>7 databases · 46.4 million variants</em></p>
  <p><em>1 curated whitelist</em></p>
  <p>
    <strong>Author:</strong> Dr Christopher Trethewey<br>
    <strong>Email:</strong> christopher.trethewey@nhs.net
  </p>
</div>

---

## Overview

OncoSieve builds a pan-cancer somatic variant whitelist from multiple curated databases and uses it to rescue low-VAF variants from Mutect2 post-filter calls. Whitelist entries are scored for pathogenicity with REVEL and PrimateAI-3D. Designed for use in research NGS pipelines.

> [!WARNING]
> **Research use only.** OncoSieve is an experimental research tool and has not been validated for clinical diagnostic use. It must not be used to inform patient management or clinical decision-making without independent validation in an accredited diagnostic setting.

---

## Data sources

| Source         | Version       | Type | Approx. variants (source total) | Notes                                                                                |
|----------------|---------------|------|---------------------------------|--------------------------------------------------------------------------------------|
| COSMIC         | v103          | File | ~38,000,000                     | GRCh38 TSV + VCF + Classification file; cancer type resolved via PRIMARY_HISTOLOGY  |
| AACR GENIE     | v19.0         | File | ~3,750,000                      | MAF; 271,837 samples / 227,696 patients; 19 cancer centres; GRCh37 lifted to GRCh38 |
| TCGA mc3       | v0.2.8        | File | ~3,600,963                      | PanCancer Atlas MAF; 33 cancer types; 10,295 tumours; GRCh37 lifted to GRCh38       |
| ClinVar        | 2026-03-09    | File | ~1,000,000                      | GRCh38 VCF; pathogenic/likely pathogenic and somatic-flagged filtered                |
| OncoKB         | v7.1 (Apr 2026) | API | ~7,700                         | ~850 genes; 130 cancer types; academic token required                                |
| TP53 database  | R21           | File | ~29,900                         | GRCh38 CSV; functional annotations for >9,000 mutant proteins                        |
| CancerHotspots | v2            | File | ~3,181                          | GRCh38 MANE Select TSV; q-value filtered                                             |
|                |               |      |                                 |                                                                                      |
| **Raw total**  |               |      | **~46,400,000**                 | Pre-deduplication; significant inter-database overlap expected                        |

### Annotation layers

Applied after whitelist construction, in `tools/post_pipeline.py`:

| Layer        | Version           | Coverage          | Output columns                                              |
|--------------|-------------------|-------------------|-------------------------------------------------------------|
| REVEL        | v1.3              | Missense          | `revel_score` (0–1, higher = more pathogenic)               |
| PrimateAI-3D | hg38 (70.6M)      | Missense          | `primateai_score`, `primateai_percentile`, `primateai_prediction` |

### Cancer type annotation

- **COSMIC**: resolved via `Cosmic_Classification_v103_GRCh38.tsv.gz` using `PRIMARY_HISTOLOGY`. This file must be present and configured in `config.yaml`.
- **GENIE**: resolved from `data_clinical_sample.txt` per sample using `CANCER_TYPE`, falling back to `ONCOTREE_CODE`.
- **TCGA**: cancer type is not resolvable from the mc3 MAF barcode (position 1 encodes a Tissue Source Site code, not a project abbreviation). TCGA rows contribute sample counts but leave cancer_type empty. A TSS lookup table will be added in a future release.
- **ClinVar**: primary disease name from `CLNDN` field (first term only).
- **OncoKB / CancerHotspots**: `pan_cancer`.
- **TP53**: morphology from the database `Morphology` field.

### ClinVar-only variants

ClinVar variants with no supporting sample observations from count-based sources have `n_samples = 0`. These represent expert clinical curation rather than observed recurrence. The post-pipeline step automatically produces a high-confidence output (`pan_cancer_whitelist_GRCh38_highconf.tsv.gz`) that excludes these entries. See [Post-pipeline processing](#post-pipeline-processing).

---

## Directory structure

```
oncosieve/
├── build_whitelist.py          # Main pipeline entry point
├── mutect2_rescue.py           # Mutect2 post-filter rescue
├── post_process_whitelist.py   # Adds genome_version, transcript_id, protein_change columns
├── pre_check.py                # Pre-run dependency and data audit
├── test_pipeline.py            # Test harness
├── run_oncosieve.sh            # Orchestrator; run this to start the pipeline
├── config.yaml                 # File paths and source enable/disable flags
├── settings.yaml.example       # Template — copy to settings.yaml and add your OncoKB token (settings.yaml is gitignored)
├── requirements.txt
├── assets/                     # Logos used by the HTML report and README
├── packaging/                  # Docker, conda, pip distribution
│   ├── Dockerfile
│   ├── docker-compose.yml
│   ├── environment.yml         # conda env spec
│   ├── setup.py                # pip-installable package
│   └── Makefile
├── panels/
│   └── lymphoma.bed            # Add panel BED files here
├── parsers/
│   ├── parse_cosmic.py         # COSMIC GenomeScreensMutant parser
│   ├── parse_genie.py          # AACR GENIE MAF parser
│   ├── parse_clinvar.py        # ClinVar VCF parser
│   ├── parse_oncokb.py         # OncoKB API parser
│   ├── parse_tp53.py           # NCI TP53 database parser
│   ├── parse_tcga.py           # TCGA mc3 MAF parser
│   └── parse_hotspots.py       # CancerHotspots TSV parser
├── tools/
│   ├── annotate_panels.py      # Panel BED annotation
│   ├── annotate_revel.py       # Standalone REVEL annotation
│   ├── annotate_primateai.py   # Standalone PrimateAI-3D annotation
│   ├── clinvar_vep_annotate.py
│   ├── db_fix.py               # Liftover GRCh37->GRCh38 for GENIE and TCGA MAFs
│   ├── hotspots_vep_remap.py
│   ├── mane_audit.py
│   ├── mane_remap.py
│   ├── post_pipeline.py        # Post-pipeline: REVEL + PrimateAI-3D, high-confidence filter, xlsx export, HTML report
│   └── generate_report.py      # HTML report with Plotly charts and DataTables
└── logs/                       # Auto-created on first run
```

---

## Requirements

### System dependencies

- Python >= 3.10
- bcftools
- tabix (htslib)

### Python dependencies

```
pandas>=2.0
pyyaml>=6.0
requests>=2.31
pysam>=0.21
polars>=0.20
pyarrow>=14.0
openpyxl>=3.1
xlrd>=2.0
plotly>=5.0
pyliftover>=0.4
```

Install with:

```bash
pip install -r requirements.txt
```

---

## Quick start (new users)

```bash
# 1. Clone and install
git clone https://github.com/Trethewey/oncosieve.git && cd oncosieve
pip install -r requirements.txt

# 2. Copy settings.yaml.example to settings.yaml and add your OncoKB token
cp settings.yaml.example settings.yaml
# Edit settings.yaml: set oncokb.api_token

# 3. Download source data (see Data Preparation below for URLs)
#    Place files in the data/ subdirectories as shown in the directory structure

# 4. Create required directories
mkdir -p data/{cosmic,genie,TCGA,clinvar,tp53,hotspots,oncokb,REVEL,PRIMATE_AI,reference}
mkdir -p output intermediate panels logs

# 5. Liftover GENIE and TCGA from GRCh37 to GRCh38 (one-time)
python3 tools/db_fix.py \
    --genie-maf data/genie/data_mutations_extended.txt \
    --tcga-maf  data/TCGA/mc3.v0.2.8.PUBLIC.maf.gz

# 6. Remap CancerHotspots to GRCh38 MANE Select (one-time, requires internet)
python3 tools/hotspots_vep_remap.py

# 7. Index ClinVar VCF (must be BGZF-compressed with CSI index)
bgzip -f data/clinvar/clinvar.vcf       # if not already BGZF
tabix -p vcf -C data/clinvar/clinvar.vcf.gz

# 8. Run pre-check to verify everything is in place
python3 pre_check.py

# 9. Run the full pipeline
bash run_oncosieve.sh
```

### Expected data directory layout

The `data/` directory is **not included in the repository** — it is gitignored because of size (~30 GB for full COSMIC + GENIE + TCGA + REVEL + PrimateAI-3D). Create it yourself and populate from the [data preparation](#data-preparation) section below.

```
data/
├── reference/
│   ├── GRCh38.fa                    # Reference genome FASTA
│   ├── GRCh38.fa.fai               # samtools index
│   └── hg19ToHg38.over.chain.gz    # UCSC liftover chain
├── cosmic/
│   ├── Cosmic_GenomeScreensMutant_v103_GRCh38.tsv.gz
│   ├── Cosmic_GenomeScreensMutant_v103_GRCh38.vcf.gz
│   └── Cosmic_Classification_v103_GRCh38.tsv.gz
├── genie/
│   ├── data_mutations_extended.txt                    # Original (GRCh37)
│   ├── data_mutations_extended_grch38_lifted.txt      # After db_fix.py
│   └── data_clinical_sample.txt
├── TCGA/
│   ├── mc3.v0.2.8.PUBLIC.maf.gz                      # Original (GRCh37)
│   └── mc3.v0.2.8.PUBLIC.GRCh38.maf.gz               # After db_fix.py
├── clinvar/
│   ├── clinvar.vcf.gz               # Must be BGZF-compressed
│   └── clinvar.vcf.gz.csi           # Must be CSI index (not TBI)
├── tp53/
│   ├── SomaticVariants_GRCh38.csv
│   └── GermlineVariants_GRCh38.csv  # Optional
├── hotspots/
│   ├── hotspots_v1.txt              # Source files
│   ├── hotspots_v2.txt
│   ├── hotspots_v3.txt
│   └── hotspots_grch38_mane.tsv     # After hotspots_vep_remap.py
├── oncokb/
│   └── allAnnotatedVariants.txt     # Optional if using API token
├── REVEL/
│   └── revel_with_transcript_ids    # ~6 GB uncompressed
└── PRIMATE_AI/
    └── PrimateAI-3D.hg38.txt.gz     # ~1.7 GB; 70.6M missense variant scores
```

---

## Data preparation

### COSMIC

Download three files from https://cancer.sanger.ac.uk/cosmic/download:

- `Cosmic_GenomeScreensMutant_v103_GRCh38.tsv.gz`
- `Cosmic_GenomeScreensMutant_v103_GRCh38.vcf.gz`
- `Cosmic_Classification_v103_GRCh38.tsv.gz`

Place all three in `data/cosmic/`. The classification file is required for cancer type resolution. Configure the path in `config.yaml` under `data_sources.cosmic.classification`.

### GENIE and TCGA liftover

GENIE and TCGA MAFs are distributed in GRCh37 and must be lifted to GRCh38 before use. This is a one-time step per data version.

```bash
# 1. Download GENIE MAF from Synapse (requires free account and data agreement)
#    https://www.synapse.org/#!Synapse:syn7222066

# 2. Download TCGA mc3 MAF from GDC PanCancer Atlas
#    https://gdc.cancer.gov/about-data/publications/pancanatlas
#    File: mc3.v0.2.8.PUBLIC.maf.gz  (md5: 639ad8f8386e98dacc22e439188aa8fa)

# 3. Run liftover
python3 tools/db_fix.py \
    --genie-maf data/genie/data_mutations_extended.txt \
    --tcga-maf  data/TCGA/mc3.v0.2.8.PUBLIC.maf.gz
```

### ClinVar

Download the GRCh38 VCF from NCBI FTP:

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz -P data/clinvar/

# IMPORTANT: ClinVar VCF must be BGZF-compressed with a CSI index
# If the downloaded file is plain gzip, recompress with bgzip:
bgzip -d data/clinvar/clinvar.vcf.gz
bgzip data/clinvar/clinvar.vcf

# Create CSI index (not TBI — required for GRCh38)
tabix -p vcf -C data/clinvar/clinvar.vcf.gz
```

### TP53

Download from the NCI TP53 Database (no registration required):

- URL: https://tp53.cancer.gov/download/
- Files: `SomaticVariants_GRCh38.csv` (required), `GermlineVariants_GRCh38.csv` (optional)
- Place in `data/tp53/`

### CancerHotspots

Download residue-level files from https://www.cancerhotspots.org:

- `hotspots_v1.txt`, `hotspots_v2.txt`, `hotspots_v3.txt`
- Place in `data/hotspots/`

Then remap to GRCh38 MANE Select coordinates (one-time, requires internet access for Ensembl VEP API):

```bash
python3 tools/hotspots_vep_remap.py
# Use --resume if interrupted
```

### OncoKB

Register for an academic API token at https://www.oncokb.org/account/register (free, expires every 6 months). Set the token in `settings.yaml` under `oncokb.api_token`.

Alternatively, download `allAnnotatedVariants.txt` from OncoKB and place in `data/oncokb/`.

### REVEL

Download the REVEL v1.3 pre-computed score file:

```bash
wget https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip -P data/REVEL/
cd data/REVEL && unzip revel-v1.3_all_chromosomes.zip
```

The extracted file `revel_with_transcript_ids` (~82 million rows, ~6 GB uncompressed) is required for post-pipeline REVEL annotation. The pipeline picks this up automatically from `data/REVEL/revel_with_transcript_ids`.

### PrimateAI-3D

PrimateAI-3D scores require an academic data licence from Illumina. Request access via the Illumina PrimateAI-3D portal and accept the licence terms before downloading.

Place the GRCh38 score file at `data/PRIMATE_AI/PrimateAI-3D.hg38.txt.gz` (~1.7 GB, 70.6M missense variants). The pipeline picks this up automatically and adds `primateai_score`, `primateai_percentile`, and `primateai_prediction` columns to the post-pipeline output. If the file is missing the step is skipped silently.

---

## Usage

### 1. Build the whitelist

```bash
bash run_oncosieve.sh
```

Optional flags:

```bash
# Skip parsers; re-run merge/filter/output from saved intermediates
bash run_oncosieve.sh --from-intermediates

# Run without OncoKB (no token required)
bash run_oncosieve.sh --skip-sources oncokb

# Skip specific sources
bash run_oncosieve.sh --skip-sources cosmic,tcga
```

### 2. Post-pipeline processing

Post-pipeline processing (REVEL + PrimateAI-3D annotation, high-confidence filtering, Excel export, HTML report) runs automatically as step 7 of `run_oncosieve.sh` if the REVEL file is present at `data/REVEL/revel_with_transcript_ids`. PrimateAI-3D annotation is added if `data/PRIMATE_AI/PrimateAI-3D.hg38.txt.gz` is present.

To run manually:

```bash
python3 tools/post_pipeline.py \
    --whitelist output/pan_cancer_whitelist_GRCh38_annotated.tsv.gz \
    --vcf       output/pan_cancer_whitelist_GRCh38.vcf.gz \
    --revel     data/REVEL/revel_with_transcript_ids \
    --primateai data/PRIMATE_AI/PrimateAI-3D.hg38.txt.gz \
    --out-dir   output/
```

This produces:

| File | Description |
|------|-------------|
| `pan_cancer_whitelist_GRCh38_full.tsv.gz` | Full whitelist with REVEL + PrimateAI-3D scores, restructured columns |
| `pan_cancer_whitelist_GRCh38_full.xlsx` | Full whitelist Excel export with summary and data sheets |
| `pan_cancer_whitelist_GRCh38_highconf.tsv.gz` | High-confidence whitelist (ClinVar-only zero-sample variants removed) |
| `pan_cancer_whitelist_GRCh38_highconf.vcf.gz` | High-confidence VCF |
| `pan_cancer_whitelist_GRCh38_highconf.xlsx` | High-confidence Excel export with summary and data sheets |
| `oncosieve_report_YYYY-MM-DD.html` | Interactive HTML report with sources table, KPI cards, Plotly charts, searchable variant table, CSV download, references |

### 3. Annotate with diagnostic panels

```bash
python3 tools/annotate_panels.py \
    --whitelist  output/pan_cancer_whitelist_GRCh38_full.tsv.gz \
    --panels-dir panels/ \
    --output     output/pan_cancer_whitelist_GRCh38_full_annotated.tsv.gz
```

### 4. Rescue Mutect2 variants

```bash
python3 mutect2_rescue.py \
    --input     sample.vcf.gz \
    --whitelist output/pan_cancer_whitelist_GRCh38_highconf.vcf.gz \
    --output    sample.rescued.vcf.gz \
    --decompose
```

> **Note:** Multi-allelic records are rescued at the record level. Pass `--decompose` to split them into biallelic records first with `bcftools norm -m-` (recommended; requires bcftools on PATH). Without this, if any ALT in a multi-allelic record matches the whitelist, all ALTs in that record are rescued.

---

## Output columns

| # | Column | Type | Description |
|---|--------|------|-------------|
| 1 | chrom | str | Chromosome (1-22, X, Y, M; no chr prefix) |
| 2 | pos | int | 1-based genomic position (GRCh38) |
| 3 | ref | str | Reference allele |
| 4 | alt | str | Alternate allele |
| 5 | gene | str | HGNC gene symbol |
| 6 | transcript_id | str | Ensembl transcript accession (e.g. ENST00000256078) |
| 7 | is_mane_select | bool | True if transcript is MANE Select (source = CancerHotspots) |
| 8 | hgvsc | str | Full HGVSc notation including transcript (e.g. ENST00000256078.4:c.35G>A) |
| 9 | hgvsp | str | Full HGVSp notation (e.g. ENSP00000256078.4:p.Gly12Asp) |
| 10 | protein_change | str | Normalised protein change in 3-letter HGVS format (e.g. p.Gly12Asp) |
| 11 | consequence | str | Standardised consequence term |
| 12 | n_cancer_types | int | Number of distinct cancer types |
| 13 | cancer_types | str | Pipe-delimited cancer type labels |
| 14 | n_samples | int | Total samples across count-based sources (0 for ClinVar/OncoKB) |
| 15 | sources | str | Pipe-delimited contributing source databases |
| 16 | oncokb_oncogenicity | str | OncoKB oncogenicity classification |
| 17 | clinvar_clinical_significance | str | ClinVar clinical significance |
| 18 | transcript_source | str | Database that provided the hgvsc/transcript (CancerHotspots, GENIE, ClinVar, COSMIC, TCGA, or empty) |
| 19 | tp53_class | str | TP53 functional classification (DNE_LOF, notDNE_LOF, notDNE_notLOF, unclass.) |
| 20 | wl_tier | int | Whitelist tier (1=highest confidence, 2=good evidence, 3=count-based below Tier 2) |
| 21 | genome_version | str | Genome build (GRCh38) |
| 22 | revel_score | float | REVEL pathogenicity score (missense only; 0-1) |
| 23 | primateai_score | float | PrimateAI-3D pathogenicity score (missense only; 0-1) |
| 24 | primateai_percentile | float | PrimateAI-3D genome-wide percentile (0-1) |
| 25 | primateai_prediction | str | PrimateAI-3D classification (benign / pathogenic) |

---

## Tiering

| Tier | Criteria | Min VAF for rescue |
|------|----------|--------------------|
| 1 | OncoKB Oncogenic or Likely Oncogenic, OR seen in ≥50 samples across ≥3 cancer types | 0.5% |
| 2 | OncoKB Predicted Oncogenic, OR ClinVar pathogenic/somatic, OR CancerHotspots recurrent mutation, OR seen in ≥25 samples across ≥2 cancer types | 0.5% |
| 3 | Count-based only, below Tier 2 threshold | 1.0% |

Tier 3 contains variants that meet the base inclusion filter but fall below the Tier 2 count threshold. This includes TP53 database variants that pass on functional annotation rather than recurrence alone.

PrimateAI-3D scores are annotated onto each variant but do not feed into tier assignment in v2; tiering remains driven by recurrence and curated evidence. PrimateAI-3D promotion of Tier 3 missense variants is implementable (the logic is sketched in `tools/annotate_primateai.py`) but is held back pending in-house validation against orthogonal pathogenicity calls — once benchmarked it can be enabled without re-running the rest of the pipeline.

VAF floors are configurable in `settings.yaml` under `vaf_rescue`.

---

## Transcript annotation

HGVSc strings are carried through from source databases where available. They are not re-annotated by OncoSieve. Transcript concordance across sources is not guaranteed. Use coordinates (chrom, pos, ref, alt) as the authoritative join key.

---

## Configuration

Edit `settings.yaml` to change thresholds, VAF floors, and OncoKB token.
Edit `config.yaml` to change file paths or disable sources.

---

## OncoKB token

Academic licence required. Token expires every 6 months. Update in `settings.yaml` under `oncokb.api_token`.

Register at: https://www.oncokb.org/account/register

---

## Packaging

The `packaging/` directory ships three distribution paths so the pipeline can be deployed beyond a manual git clone.

### Docker

```bash
docker build -f packaging/Dockerfile -t oncosieve .
docker run --rm -v $PWD/data:/app/data -v $PWD/output:/app/output oncosieve
```

A `docker-compose.yml` is provided in `packaging/` for single-command runs once `data/` is in place.

### Conda

```bash
conda env create -f packaging/environment.yml
conda activate oncosieve
bash run_oncosieve.sh
```

### pip (editable install)

```bash
pip install -e packaging/
oncosieve   # exposed as a CLI entry point
```

All three paths assume the `data/` directory is populated as described in [Data preparation](#data-preparation). Reference data is not bundled with any package.

---

## Licence

OncoSieve is released under the [GNU Affero General Public License v3.0 or later (AGPL-3.0-or-later)](LICENSE). AGPL is a strong copyleft licence: you may copy, modify, and redistribute the code, but any distribution or **network-accessible deployment** of a modified version must offer the corresponding modified source under the same licence. See the [FSF page on the AGPL](https://www.gnu.org/licenses/agpl-3.0.html) for the full terms.

The underlying **data sources** are governed by their own licences (e.g. COSMIC requires a commercial-use licence, OncoKB requires an academic API token, Illumina PrimateAI-3D requires a data licence). You are responsible for complying with each source's terms separately.

---

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for release notes per version.

---

## References

If you use OncoSieve in published work, please cite the underlying databases using the references below.

**COSMIC**
Tate JG, Bamford S, Jubb HC, et al. COSMIC: the Catalogue Of Somatic Mutations In Cancer. Nucleic Acids Res. 2019;47(D1):D941-D947.
[doi:10.1093/nar/gky1015](https://doi.org/10.1093/nar/gky1015) · [PMID:30371878](https://pubmed.ncbi.nlm.nih.gov/30371878/)

**AACR Project GENIE**
The AACR Project GENIE Consortium. AACR Project GENIE: Powering Precision Medicine through an International Consortium. Cancer Discov. 2017;7(8):818-831.
[doi:10.1158/2159-8290.CD-17-0151](https://doi.org/10.1158/2159-8290.CD-17-0151)
Include the version of GENIE data used (e.g. v19.0) in your methods.

**TCGA mc3 PanCancer Atlas**
Ellrott K, Bailey MH, Saksena G, et al. Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines. Cell Syst. 2018;6(3):271-281.e7.
[doi:10.1016/j.cels.2018.03.002](https://doi.org/10.1016/j.cels.2018.03.002) · [PMID:29596782](https://pubmed.ncbi.nlm.nih.gov/29596782/)

**OncoKB**
Suehnholz SP, Nissan MH, Zhang H, et al. Quantifying the Expanding Landscape of Clinical Actionability for Patients with Cancer. Cancer Discov. 2024;14(1):49-65.
[doi:10.1158/2159-8290.CD-23-0467](https://doi.org/10.1158/2159-8290.CD-23-0467) · [PMID:37849038](https://pubmed.ncbi.nlm.nih.gov/37849038/)

Chakravarty D, Gao J, Phillips SM, et al. OncoKB: A Precision Oncology Knowledge Base. JCO Precis Oncol. 2017;1:PO.17.00011.
[doi:10.1200/PO.17.00011](https://doi.org/10.1200/PO.17.00011) · [PMID:28890946](https://pubmed.ncbi.nlm.nih.gov/28890946/)

**ClinVar**
Landrum MJ, Lee JM, Benson M, et al. ClinVar: public archive of interpretations of clinically relevant variants. Nucleic Acids Res. 2016;44(D1):D862-868.
[doi:10.1093/nar/gkv1222](https://doi.org/10.1093/nar/gkv1222) · [PMID:26582918](https://pubmed.ncbi.nlm.nih.gov/26582918/)

**NCI TP53 Database**
de Andrade KC, Lee EE, Tookmanian EM, et al. The TP53 Database: transition from the International Agency for Research on Cancer to the US National Cancer Institute. Cell Death Differ. 2022;29(5):1071-1073.
[doi:10.1038/s41418-022-00976-3](https://doi.org/10.1038/s41418-022-00976-3) · [PMID:35352025](https://pubmed.ncbi.nlm.nih.gov/35352025/)

Database release: The TP53 Database (R21, Jan 2025): <https://tp53.cancer.gov>

**Cancer Hotspots**
Bandlamudi C, et al. Cancer type-specific variation in patterns of driver alterations across 50,000 tumors. <https://www.cancerhotspots.org>

Chang MT, Asthana S, Gao SP, et al. Identifying recurrent mutations in cancer reveals widespread lineage diversity and mutational specificity. Nat Biotechnol. 2016;34(2):155-163.
[doi:10.1038/nbt.3391](https://doi.org/10.1038/nbt.3391) · [PMID:26619011](https://pubmed.ncbi.nlm.nih.gov/26619011/)

Chang MT, et al. Accelerating discovery of functional mutant alleles in cancer. Cancer Discov. 2018;8(2):174-183.
[doi:10.1158/2159-8290.CD-17-0321](https://doi.org/10.1158/2159-8290.CD-17-0321) · [PMID:29247016](https://pubmed.ncbi.nlm.nih.gov/29247016/)

**REVEL**
Ioannidis NM, Rothstein JH, Pejaver V, et al. REVEL: An Ensemble Method for Predicting the Pathogenicity of Rare Missense Variants. Am J Hum Genet. 2016;99(4):877-885.
[doi:10.1016/j.ajhg.2016.08.016](https://doi.org/10.1016/j.ajhg.2016.08.016) · [PMID:27666373](https://pubmed.ncbi.nlm.nih.gov/27666373/)

**PrimateAI-3D**
Gao H, Hamp T, Ede J, et al. The landscape of tolerated genetic variation in humans and other primates. Science. 2023;380(6648):eabn8197.
[doi:10.1126/science.abn8197](https://doi.org/10.1126/science.abn8197) · [PMID:37262156](https://pubmed.ncbi.nlm.nih.gov/37262156/)
