# =============================================================================
# OncoSieve -- Makefile
#
# Usage:
#   make install          Install conda environment
#   make check            Run pre-flight checks
#   make run              Run full pipeline
#   make run-fast         Re-run from intermediates (skip parsing)
#   make postprocess      Run post-pipeline only (REVEL, Excel, HTML)
#   make rescue VCF=x.vcf.gz TUMOR=SAMPLE  Mutect2 rescue
#   make test             Validate outputs
#   make report           Generate HTML report only
#   make clean            Remove output and intermediate files
#   make clean-all        Remove everything including logs
# =============================================================================

SHELL := /bin/bash
.DEFAULT_GOAL := help

# Configuration
PYTHON       := python3
CONFIG       := config.yaml
OUTPUT_DIR   := output
PREFIX       := pan_cancer_whitelist_GRCh38
DATA_DIR     ?=
SKIP_SOURCES ?=

# Derived paths
WL_VCF       := $(OUTPUT_DIR)/$(PREFIX).vcf.gz
WL_TSV       := $(OUTPUT_DIR)/$(PREFIX).tsv.gz
WL_ANNOTATED := $(OUTPUT_DIR)/$(PREFIX)_annotated.tsv.gz
WL_FULL      := $(OUTPUT_DIR)/$(PREFIX)_full.tsv.gz
WL_HIGHCONF  := $(OUTPUT_DIR)/$(PREFIX)_highconf.tsv.gz
REVEL_FILE   := data/REVEL/revel_with_transcript_ids

# Build args
BUILD_ARGS := --config $(CONFIG)
ifdef DATA_DIR
  BUILD_ARGS += --data-dir $(DATA_DIR)
endif
ifdef SKIP_SOURCES
  BUILD_ARGS += --skip-sources $(SKIP_SOURCES)
endif

# =============================================================================

.PHONY: help install check run run-fast postprocess rescue test report clean clean-all

help:  ## Show this help
	@echo ""
	@echo "  OncoSieve -- Available targets:"
	@echo ""
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "    %-18s %s\n", $$1, $$2}'
	@echo ""
	@echo "  Variables:"
	@echo "    DATA_DIR=/path      Override data directory"
	@echo "    SKIP_SOURCES=x,y    Skip specific sources"
	@echo "    VCF=sample.vcf.gz   Mutect2 VCF for rescue"
	@echo "    TUMOR=SAMPLE_NAME   Tumor sample name for rescue"
	@echo ""

install:  ## Create conda environment and install dependencies
	conda env create -f packaging/environment.yml || conda env update -f packaging/environment.yml
	@echo ""
	@echo "  Environment created. Activate with:"
	@echo "    conda activate oncosieve"

check:  ## Run pre-flight checks
	$(PYTHON) pre_check.py --config $(CONFIG) $(if $(DATA_DIR),--data-dir $(DATA_DIR),) $(if $(SKIP_SOURCES),--skip-sources $(SKIP_SOURCES),)

run: check  ## Run full pipeline (parse + merge + filter + tier + output)
	@mkdir -p $(OUTPUT_DIR) intermediate logs
	$(PYTHON) build_whitelist.py $(BUILD_ARGS)
	tabix -f -p vcf $(WL_VCF)
	$(PYTHON) post_process_whitelist.py --whitelist $(WL_TSV) --out $(WL_ANNOTATED)
	@echo ""
	@echo "  Whitelist built: $(WL_VCF)"
	@echo "  Annotated TSV:   $(WL_ANNOTATED)"
	@if [ -f "$(REVEL_FILE)" ]; then \
		$(MAKE) postprocess; \
	else \
		echo "  REVEL file not found -- skipping post-pipeline."; \
	fi

run-fast:  ## Re-run from saved intermediates (skip parsing)
	$(PYTHON) build_whitelist.py $(BUILD_ARGS) --from-intermediates
	tabix -f -p vcf $(WL_VCF)
	$(PYTHON) post_process_whitelist.py --whitelist $(WL_TSV) --out $(WL_ANNOTATED)
	@if [ -f "$(REVEL_FILE)" ]; then $(MAKE) postprocess; fi

postprocess:  ## Run post-pipeline (REVEL annotation, Excel, HTML report)
	$(PYTHON) tools/post_pipeline.py \
		--whitelist $(WL_ANNOTATED) \
		--vcf $(WL_VCF) \
		--revel $(REVEL_FILE) \
		--out-dir $(OUTPUT_DIR)/

rescue:  ## Run Mutect2 rescue (requires VCF= and optional TUMOR=)
ifndef VCF
	$(error VCF is required. Usage: make rescue VCF=sample.vcf.gz TUMOR=SAMPLE_NAME)
endif
	$(PYTHON) mutect2_rescue.py \
		--input $(VCF) \
		--whitelist $(WL_VCF) \
		--output $(VCF:.vcf.gz=.rescued.vcf.gz) \
		--config $(CONFIG) \
		$(if $(TUMOR),--tumor-sample $(TUMOR),)
	tabix -f -p vcf $(VCF:.vcf.gz=.rescued.vcf.gz)

test:  ## Validate pipeline outputs
	$(PYTHON) test_pipeline.py

report:  ## Generate HTML report only
	$(PYTHON) tools/generate_report.py \
		--full $(WL_FULL) \
		--highconf $(WL_HIGHCONF) \
		--out-dir $(OUTPUT_DIR)/

clean:  ## Remove output and intermediate files
	rm -rf $(OUTPUT_DIR)/*.tsv.gz $(OUTPUT_DIR)/*.vcf.gz $(OUTPUT_DIR)/*.tbi
	rm -rf $(OUTPUT_DIR)/*.xlsx $(OUTPUT_DIR)/*.html
	rm -rf intermediate/*.tsv.gz

clean-all: clean  ## Remove everything including logs
	rm -rf logs/*
