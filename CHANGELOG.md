# Changelog

All notable changes to OncoSieve are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.1.0] — 2026-07-28

Correctness release. A systematic per-source annotation audit surfaced 13
HIGH-severity and 11 MEDIUM-severity defects across all 7 parsers and the
display code. Every finding on the audit is addressed here. Numeric outputs
will change against v2.0 as a result — this release is intended to be run
against a fresh pipeline, not compared row-for-row against `v2.0` results.

### Fixed — clinical safety
- **ClinVar germline contamination.** `parse_clinvar.py` was looking for an
  `INFO/CLNORIGIN` field that does not exist in modern ClinVar VCFs; the
  actual field is `INFO/ORIGIN`. The regex never matched, so the somatic-
  origin filter was silently dead and hereditary germline pathogenic
  variants (BRCA1, MLH1, hereditary TP53 syndromes) entered the somatic
  whitelist whenever `CLNDN` mentioned a cancer term. Fixed the regex,
  clarified the docstring, and added an "unknown origin" (value 0) path
  that falls back to the CLNDN cancer-term gate.
- **TCGA `cancer_type` resolved from the wrong barcode field.** mc3
  barcodes are `TCGA-{TSS}-{Participant}-{Sample}-…` where `{TSS}` is a
  two-character Tissue Source Site code, not a study abbreviation.
  `_TCGA_PROJECT_MAP` never matched, so every TCGA row was labelled
  `TCGA_<TSS>` and each TSS was counted as its own cancer type. Embedded
  GDC's canonical tissueSourceSite table (830 entries) as `_TCGA_TSS_MAP`,
  added the six historical study abbreviations to `_TCGA_PROJECT_MAP`
  (`CNTL`, `COADREAD`, `FPPP`, `LCML`, `MISC`, `STES`), and rewrote
  `_cancer_type_from_barcode` to chain TSS → study → label. Unresolved
  barcodes now return `'unspecified'` (the aggregator sentinel) instead
  of a junk `TCGA_XX` string.
- **CancerHotspots inflated `n_cancer_types`.** The parser hard-coded
  `cancer_type = 'pan_cancer'`, so a variant seen in two real tumour
  types plus CancerHotspots counted as three distinct cancer types and
  cleared `tier1_min_cancer_types = 3` on artificial diversity.
  Emits `'unspecified'` instead; the pan-cancer semantics are already
  carried by `passes_hotspot`.
- **CancerHotspots + ClinVar consequence degraded to `'unknown'`.**
  Both parsers now route consequences through `map_consequence`, which
  handles VEP's `&`-joined multi-term consequences. ClinVar iterates
  every comma-separated `MC` entry and picks the most severe via an
  explicit severity ranking rather than always taking the first.

### Fixed — annotation correctness
- **`is_valid_allele` no longer accepts `-`.** MAF-style `-` alleles are
  not valid VCF and had been leaking into the output as `REF=-`,
  breaking `tabix` indexing and Mutect2 rescue key lookups. The
  aggregator now rejects them; parsers that emit `-` are responsible
  for VCF-normalising the row (anchor base) or dropping it.
- **COSMIC classification lookup is now required.** Missing the
  classification TSV previously left every phenotype ID as its own
  cancer type, promoting all COSMIC variants past the tier-1 threshold.
  A missing file now raises rather than warns. Failed lookups emit
  `'unspecified'` instead of the raw `cosoNNNNN` ID. `clean_allele` is
  applied to COSMIC REF/ALT for parity with the other parsers.
- **GENIE `FILTER` column is now read.** Non-PASS calls (artefacts the
  contributing centre had already rejected) no longer contribute to
  n_samples or tier promotion.
- **GENIE + TCGA HGVSc now carry the `ENST` transcript prefix.** Both
  parsers read `Transcript_ID` from the MAF and prepend it to bare
  `c.` HGVSc so the downstream MANE Select lookup can resolve to a
  RefSeq accession. GENIE also prefers the VEP `Consequence` column
  (multi-term `&`-joined) over the coarser MAF `Variant_Classification`.
- **TP53 indels no longer fabricate placeholder bases.** Deletions had
  been emitted with `REF == ALT` (the parser dropped the last-base
  step); insertions had used a literal `N` as the anchor; duplications
  sat at the wrong position. `parse_tp53` now accepts a
  `reference_fasta` parameter, opens it with `pyfaidx`, and normalises
  each indel to a proper VCF anchor at `pos-1`. When FASTA is
  unavailable (no config entry, no pyfaidx, or file missing) the row is
  DROPPED with a log warning rather than fabricated. SNVs and delins
  entries remain unchanged. `deletion` / `insertion` / `duplication`
  effect labels are now frame-derived from `(len(ref)-len(alt)) % 3`
  rather than blindly mapped to `frameshift`.
- **OncoKB `_map_effect` no longer labels every non-neutral entry as
  `missense`.** Truncating, fusion, amplification, splice, and indel
  alterations are inspected via the Alteration string and land in the
  right consequence bucket. `_hgvsp_to_oncokb` rejects strings with
  a suffix (`fs`, `ext`, `del`, `ins`, `dup`) so `p.Ala100Argfs*7`
  no longer submits `A100R` to the byProteinChange endpoint.
- **OncoKB cancer_type sentinel changed** from `'pan_cancer'` to
  `'unspecified'` for parity with the other sources.
- **ClinVar `CLNDN` is now split, sentinels filtered, all survivors
  pipe-joined.** `not_provided`, `not_specified`, `see_cases`, and
  empty entries no longer contaminate `n_cancer_types`.
- **CONSEQUENCE_MAP additions.** COSMIC's `Substitution - Silent`,
  `Substitution - coding silent`, `Complex - compound substitution`,
  and Splice donor/acceptor terms map to their intended standardised
  terms instead of degrading to `'unknown'`.

### Fixed — VCF/TSV output correctness
- **VCF `CANCER_TYPES` and `SOURCES` are now VCFv4.2-compliant arrays.**
  The `Number=.` declaration required a comma-separated array, but the
  writer emitted a single pipe-joined string, so
  `bcftools view -i 'INFO/SOURCES[*]="ClinVar"'` and similar array
  queries silently returned empty. The aggregator still stores pipe-
  joined values for TSV/Excel readability; the VCF writer splits and
  rejoins on commas with per-element escaping.
- **`_vcf_escape` extended** to percent-encode comma, CR, tab and
  newline (previously silently deleted), so control characters that
  slip in from a source file become visible rather than invisibly
  truncating a field. `GENE` and `CSQ` are now routed through the
  escape function alongside `HGVSc` / `HGVSp`.
- **TSV row order now karyotypic** (1, 2, …, 22, X, Y, MT) matching the
  VCF, using the same `_chrom_sort` map. Lexicographic order produced
  `1, 10, 11, …, 2` which surprised users diffing the two files.
- **`assign_tiers` logs the active thresholds** and warns on unknown
  or missing keys in the settings.yaml `tiering` block, so a typo no
  longer silently substitutes the built-in defaults.

### Fixed — HTML report
- **Dynamic column indices.** `_TIER_COL_IDX` and `_NSAMP_COL_IDX`
  were module-level constants derived from `DISPLAY_COLS`. If any
  upstream column was missing (e.g. no PrimateAI columns, no
  RefSeq ID), the report would sort by and colour the wrong column
  in the DataTable. Both indices are now computed inside
  `build_datatable` from the actual filtered column set and threaded
  through to the CSS and JS via `build_report`.
- **Report version is no longer hard-coded.** Introduced
  `ONCOSIEVE_VERSION` at module scope; the header `<h1>` and page
  `<title>` interpolate it instead of the literal "2.0".
- **Header meta strip is derived from the data.** REVEL and
  PrimateAI-3D annotations are advertised only when the relevant
  columns are present with at least one non-null value.
- **Tier 3 KPI badge added** next to Tier 1 / Tier 2 on both summary
  cards; the Excel summary sheet already counted Tier 3 so the KPI
  card is now consistent.
- **`cancer_types` added to `DISPLAY_COLS`** immediately after
  `n_cancer_types`; the previous surface let users filter on the
  count without seeing which cancer types contributed. A
  human-friendly `DISPLAY_COL_LABELS` map replaces raw snake_case
  identifiers in the DataTable header.
- **"Top Genes ... top 60" title corrected to "top 100"** to match
  the chart's actual `top_n = 100`.
- **TP53 chart colour** no longer collides with TCGA.
- **`database_versions.txt` passed through from `post_pipeline.py`.**
  Previously the pipeline-produced report auto-discovered the file
  only when `generate_report.py` was invoked directly; the
  post-pipeline entry point never passed the path, so the Sources
  versions table was silently absent from every clinician-facing
  report.

### Added
- **PrimateAI-3D tier promotion** (planned for v2, implemented here).
  A Tier 3 variant with `primateai_prediction == 'pathogenic'` and
  `primateai_percentile >= 0.9` is promoted to Tier 2. Logged with
  a count of promoted variants.
- **`_TCGA_TSS_MAP`** embedded in `parse_tcga.py` (830 entries,
  GDC-sourced).

### Changed
- **Repository slug** renamed from `Trethewey/OncoSieve` to
  `Trethewey/oncosieve` (lowercase). GitHub redirects the old URL;
  every in-repo URL reference now points at the lowercase slug.
- **Licence** changed from MIT to **AGPL-3.0-or-later**. Downstream
  distribution or network-accessible deployment of a modified
  version must offer the corresponding modified source under the
  same licence.
- **Packaging headers** no longer stamped `ONCOSIEVE v1.0` —
  `packaging/Dockerfile`, `packaging/Makefile`,
  `packaging/environment.yml`, and `packaging/docker-compose.yml`
  now read `OncoSieve`. `packaging/setup.py` version bumped to
  `2.1.0`.

### Notes
- `TP53` germline processing still requires the germline TSV
  configured in `config.yaml` (`GermlineVariants_GRCh38.csv`).
- Server-side pipeline must be re-run to produce a v2.1 whitelist —
  counts and tier assignments will differ from v2.0.

## [2.0.0] — 2026-05-19

### Added
- **PrimateAI-3D annotation** (`tools/annotate_primateai.py`). Polars lazy-scan over the 1.7 GB hg38 dataset; adds `primateai_score`, `primateai_percentile`, and `primateai_prediction` columns to the post-pipeline output. Pipeline auto-detects `data/PRIMATE_AI/PrimateAI-3D.hg38.txt.gz`.
- **PrimateAI-3D citation** in the references table (Gao et al., *Science* 2023).
- **Sources table** in the HTML report, auto-generated from `output/database_versions.txt`. Shows database, version/date, file vs API, and origin.
- **CSV download** button on the high-confidence variant table (DataTables Buttons).
- **Favicon** in the HTML report, embedded as a base64 data URI for self-contained portability.
- **Packaging directory** with `Dockerfile`, `docker-compose.yml`, conda `environment.yml`, `setup.py` (pip), and `Makefile`.
- **MIT licence** (`LICENSE`); previously labelled "Research Use Only" in `setup.py`.
- **CHANGELOG.md** (this file).

### Changed
- **Brand:** `ONCOSIEVE` → `OncoSieve` across code surfaces, HTML report, Excel sheet titles, and `setup.py`. File paths and the GitHub repo URL retain their casing.
- **Repository renamed** on GitHub from `ONCOSIEVE` to `OncoSieve`. Old URL continues to redirect.
- **HTML report header redesigned:** mark SVG instead of cropped raster; transparent background; new wordmark in `-apple-system/Segoe UI` matching the SVG logo; "2.0" version badge baseline-aligned beside the title; restructured subtitle (`Pan-Cancer Variant Whitelist – Generated DATE`); source line lists databases + REVEL + PrimateAI-3D.
- **HTML report variant table** now uses DataTables `data:` JSON-init instead of pre-rendered `<tr>` rows. ~40% smaller HTML (6.8 → 4.1 MB) and renders effectively instantly.
- **High-Confidence KPI card** recoloured from teal-green to logo red `#C0392B`.
- **References section** in the report: PrimateAI-3D entry added; text colour grey → white; font size reduced by 2 pt.
- **TP53 database version** pinned in `config.yaml` (`tp53.version: R21`); `build_whitelist.py` reads this in preference to the previous mtime fallback.
- **`run_oncosieve.sh`** auto-detects `data/PRIMATE_AI/PrimateAI-3D.hg38.txt.gz` and passes `--primateai` to `post_pipeline.py` when present.
- **`.gitignore`** reorganised into labelled sections; adds `settings.yaml`, `NOTES.md`, `nohup.out`.
- **README** rewritten: new MIT badge under the logo, Annotation Layers sub-table, updated source versions (ClinVar 2026-03-09, OncoKB v7.1 Apr 2026), new PrimateAI-3D data-prep section, new Packaging section, Licence section, output columns 23–25 documented, tiering section accurately states PrimateAI-3D is annotation-only in v2.
- **`settings.yaml`** untracked. The template is `settings.yaml.example` (tracked); users copy it to `settings.yaml` and add their OncoKB token locally.
- **VCF source header** simplified to `pan_cancer_whitelist_pipeline` (no version suffix).

### Removed
- **`packaging/install_and_run.bat`** — the Windows .bat installer already depended on Git Bash / WSL (for `tee` and bcftools). Docker is now the supported Windows distribution path.
- **Tool-version strings** (`ONCOSIEVE v1.0`, `pipeline_v1.0`) from code surfaces, banner, Excel title, and VCF source header. The version lives in the Git tag and the in-report "2.0" badge.
- **Tracked `settings.yaml`** — see Changed.

### Not implemented in this release
- **PrimateAI-3D tier promotion.** Documented as a future enhancement in `tools/annotate_primateai.py:326`. The scores are annotated but do not feed into the tier assignment.

## [1.0.0] — 2026-03-16

### Added
- Initial release of OncoSieve.
- Pan-cancer somatic variant whitelist built from seven curated databases:
  COSMIC v103, AACR GENIE v19.0, TCGA mc3 v0.2.8, ClinVar, OncoKB,
  TP53 Database, CancerHotspots v2 — ~46.4 million raw variants pre-deduplication.
- Tiered output (Tier 1 / 2 / 3) with VAF floors configurable per tier.
- MANE Select enrichment for transcript IDs.
- Mutect2 rescue (`mutect2_rescue.py`): applies the whitelist to FilterMutectCalls
  VCFs with tier-aware VAF floors.
- Post-pipeline (`tools/post_pipeline.py`): REVEL annotation, high-confidence
  filtering (drops ClinVar-only zero-sample rows), Excel exports.
- Interactive HTML report with Plotly charts and a DataTables variant browser.
- Orchestrator `run_oncosieve.sh` with `--from-intermediates`, `--skip-sources`,
  `--rescue` flags.
- Pre-run audit (`pre_check.py`) and validation harness (`test_pipeline.py`).
