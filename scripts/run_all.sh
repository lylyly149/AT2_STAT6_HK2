#!/bin/bash
# ==============================================================================
# AT2_STAT6_HK2: One-click reproducibility script
# ==============================================================================
#
# This script runs all analysis pipelines in sequence.
# Usage: ./scripts/run_all.sh
#
# Requirements:
#   - Conda with r_4_3 environment installed
#   - Raw data files configured in config/paths.yaml
#
# ==============================================================================

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Config file (prefer local override)
CONFIG_FILE="$PROJECT_DIR/config/paths.local.yaml"
if [ ! -f "$CONFIG_FILE" ]; then
    CONFIG_FILE="$PROJECT_DIR/config/paths.yaml"
fi
if [ ! -f "$CONFIG_FILE" ]; then
    echo -e "${RED}No config/paths.yaml or config/paths.local.yaml found.${NC}"
    exit 1
fi

# Lightweight YAML value extractor (single-line key: value)
yaml_val() {
    local key="$1"
    local file="$2"
    awk -F': *' -v k="$key" '
        $1 ~ "^[[:space:]]*" k "$" {
            val=$2
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", val)
            gsub(/^"|"$/, "", val)
            print val
            exit
        }
    ' "$file"
}

R_ENV_NAME="$(yaml_val 'r_env' "$CONFIG_FILE")"
R_ENV_NAME="${R_ENV_NAME:-r_4_3}"
PY_ENV_NAME="$(yaml_val 'py_env' "$CONFIG_FILE")"
PY_ENV_NAME="${PY_ENV_NAME:-py310}"
SCFEA_DIR="$(yaml_val 'scfea_dir' "$CONFIG_FILE")"

echo "=============================================="
echo "AT2_STAT6_HK2 Analysis Pipeline"
echo "=============================================="
echo ""
echo "Project directory: $PROJECT_DIR"
echo "Using config:      $CONFIG_FILE"
echo "R environment:     $R_ENV_NAME"
echo "Python env (scFEA): $PY_ENV_NAME"
echo ""

# Change to project directory
cd "$PROJECT_DIR"

# Create output directories
echo -e "${YELLOW}Creating output directories...${NC}"
mkdir -p data/processed
mkdir -p figures/figure1
mkdir -p figures/figures1
mkdir -p figures/figures2
mkdir -p figures/figures3
mkdir -p figures/figures6
mkdir -p figures/figures7
mkdir -p figures/figgures8
mkdir -p figures/figures9
mkdir -p figures/figures10
echo -e "${GREEN}Done.${NC}"
echo ""

# Check if conda environment exists (R)
echo -e "${YELLOW}Checking conda environment ${R_ENV_NAME}...${NC}"
if conda env list | awk '{print $1}' | grep -Fxq "${R_ENV_NAME}"; then
    echo -e "${GREEN}${R_ENV_NAME} environment found.${NC}"
else
    echo -e "${RED}${R_ENV_NAME} environment not found. Please create it first:${NC}"
    echo "  conda env create -f environment.yml"
    exit 1
fi
echo ""

# Activate conda environment
echo -e "${YELLOW}Activating ${R_ENV_NAME} environment...${NC}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${R_ENV_NAME}"
echo -e "${GREEN}Done.${NC}"
echo ""

# ==============================================================================
# Use existing scFEA outputs (skip rerun to avoid env crash)
# Copy counts/metadata/flux from original scFEA results into processed/
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Using existing scFEA outputs...${NC}"
echo "=============================================="
if [ -z "$SCFEA_DIR" ]; then
  echo -e "${RED}scfea_dir is not set in ${CONFIG_FILE}. Skipping scFEA copy step.${NC}"
else
  echo "scFEA source directory: $SCFEA_DIR"
fi

mkdir -p data/processed/scfea_input data/processed/scfea_input/ipf data/processed/scfea_flux_at2 data/processed/scfea_flux_ipf_at2

if [ -n "$SCFEA_DIR" ] && [ -d "$SCFEA_DIR" ]; then
  # Metadata and counts (if needed)
  for f in at2_counts.csv at2_metadata.csv; do
    if [ -f "$SCFEA_DIR/data/interim/scfea_input/$f" ]; then
      cp -f "$SCFEA_DIR/data/interim/scfea_input/$f" "data/processed/scfea_input/$f"
    fi
  done
  for f in ipf_at2_counts.csv ipf_at2_metadata.csv; do
    if [ -f "$SCFEA_DIR/data/interim/scfea_input/$f" ]; then
      cp -f "$SCFEA_DIR/data/interim/scfea_input/$f" "data/processed/scfea_input/$f"
    fi
  done

  # Flux outputs
  for f in adj_flux.csv adj_balance.csv; do
    if [ -f "$SCFEA_DIR/results/scfea_flux_at2/$f" ]; then
      cp -f "$SCFEA_DIR/results/scfea_flux_at2/$f" "data/processed/scfea_flux_at2/$f"
    fi
    if [ -f "$SCFEA_DIR/results/scfea_flux_ipf_at2/$f" ]; then
      cp -f "$SCFEA_DIR/results/scfea_flux_ipf_at2/$f" "data/processed/scfea_flux_ipf_at2/$f"
    fi
  done
  echo -e "${GREEN}Copied existing scFEA results into data/processed/.${NC}"
else
  echo -e "${YELLOW}scFEA directory not found; skipping copy (set raw_data.scfea_dir in config).${NC}"
fi
echo ""

# ==============================================================================
# Run Figure 1 analysis (scRNA-seq)
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Running Figure 1: scRNA-seq analysis${NC}"
echo "=============================================="
if [ -f "src/figure1/figure1.R" ]; then
    Rscript src/figure1/figure1.R
    echo -e "${GREEN}Figure 1 completed.${NC}"
else
    echo -e "${RED}Warning: src/figure1/figure1.R not found, skipping.${NC}"
fi
echo ""

# ==============================================================================
# Run Figure 2 analysis
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Running Figure 2: scRNA-seq epithelial analyses${NC}"
echo "=============================================="
if [ -f "src/figures2/figures2a.R" ]; then
    Rscript src/figures2/figures2a.R
    echo -e "${GREEN}Figure 2A completed.${NC}"
else
    echo -e "${RED}Warning: src/figures2/figures2a.R not found, skipping.${NC}"
fi

if [ -f "src/figures2/figures2bc.R" ]; then
    Rscript src/figures2/figures2bc.R
    echo -e "${GREEN}Figure 2B/C completed.${NC}"
else
    echo -e "${RED}Warning: src/figures2/figures2bc.R not found, skipping.${NC}"
fi

if [ -f "src/figures2/figures2d.R" ]; then
    Rscript src/figures2/figures2d.R
    echo -e "${GREEN}Figure 2D (DE + GO) completed.${NC}"
else
    echo -e "${RED}Warning: src/figures2/figures2d.R not found, skipping.${NC}"
fi
echo ""

# ==============================================================================
# Run Figure 3 analysis (macrophage & fibroblasts STAT6)
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Running Figure 3: macrophage/fibroblast STAT6 analyses${NC}"
echo "=============================================="
if [ -f "src/figures3/figures3a_b_e_f_g.R" ]; then
    Rscript src/figures3/figures3a_b_e_f_g.R
    echo -e "${GREEN}Figure 3A/B/E/F/G completed.${NC}"
else
    echo -e "${RED}Warning: src/figures3/figures3a_b_e_f_g.R not found, skipping.${NC}"
fi

if [ -f "src/figures3/figures3c.R" ]; then
    Rscript src/figures3/figures3c.R
    echo -e "${GREEN}Figure 3C completed.${NC}"
else
    echo -e "${RED}Warning: src/figures3/figures3c.R not found, skipping.${NC}"
fi
echo ""

# ==============================================================================
# Run Response Figure 1,3 analysis (Venn + bulk limma from raw counts / series)
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Running Response Fig 1,3: Venn diagrams${NC}"
echo "=============================================="
if [ -f "src/figures1/figures1b_run_limma_gse134692_from_rawcounts.R" ]; then
    Rscript src/figures1/figures1b_run_limma_gse134692_from_rawcounts.R
    echo -e "${GREEN}GSE134692 limma from raw counts completed.${NC}"
else
    echo -e "${RED}Warning: src/figures1/figures1b_run_limma_gse134692_from_rawcounts.R not found, skipping.${NC}"
fi
if [ -f "src/figures1/figures1d_run_limma_gse16846_from_series.R" ]; then
    Rscript src/figures1/figures1d_run_limma_gse16846_from_series.R
    echo -e "${GREEN}GSE16846 limma from series completed.${NC}"
else
    echo -e "${RED}Warning: src/figures1/figures1d_run_limma_gse16846_from_series.R not found, skipping.${NC}"
fi
if [ -f "src/figures1/figures1c_run_limma_gse32147_16w_from_series.R" ]; then
    Rscript src/figures1/figures1c_run_limma_gse32147_16w_from_series.R
    echo -e "${GREEN}GSE32147 (16w) limma from series completed.${NC}"
else
    echo -e "${RED}Warning: src/figures1/figures1c_run_limma_gse32147_16w_from_series.R not found, skipping.${NC}"
fi

if [ -f "src/figures1/figures1a-1.R" ]; then
    Rscript src/figures1/figures1a-1.R
    echo -e "${GREEN}Response Fig 1,3 completed.${NC}"
else
    echo -e "${RED}Warning: src/figures1/figures1a-1.R not found, skipping.${NC}"
fi

if [ -f "src/figures1/figures1j_k.R" ]; then
    Rscript src/figures1/figures1j_k.R
    echo -e "${GREEN}Response Fig 1J/K network analyses completed.${NC}"
else
    echo -e "${RED}Warning: src/figures1/figures1j_k.R not found, skipping.${NC}"
fi
echo ""

# ==============================================================================
# Run Figure 6 analysis (pseudotime + STAT6 alluvial)
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Running Figure 6: pseudotime heatmap and STAT6 alluvial${NC}"
echo "=============================================="
if [ -f "src/figures6/figures6.R" ]; then
    Rscript src/figures6/figures6.R
    echo -e "${GREEN}Figure 6 completed.${NC}"
else
    echo -e "${RED}Warning: src/figures6/figures6.R not found, skipping.${NC}"
fi
echo ""

# ==============================================================================
# Run Figure 7 analysis (spatial STAT6/EMT/glycolysis)
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Running Figure 7: spatial STAT6/EMT/glycolysis${NC}"
echo "=============================================="
if [ -f "src/figures7/figures7.R" ]; then
    Rscript src/figures7/figures7.R
    echo -e "${GREEN}Figure 7 completed.${NC}"
else
    echo -e "${RED}Warning: src/figures7/figures7.R not found, skipping.${NC}"
fi
echo ""

# ==============================================================================
# Run Comment 7 Figure 1 analysis (AT2 glycolysis)
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Running Comment7 Fig 1: AT2 glycolysis${NC}"
echo "=============================================="
if [ -f "src/figgures8/figgures8.R" ]; then
    Rscript src/figgures8/figgures8.R
    echo -e "${GREEN}Comment7 Fig 1 completed.${NC}"
else
    echo -e "${RED}Warning: src/figgures8/figgures8.R not found, skipping.${NC}"
fi
echo ""

# ==============================================================================
# Run Comment 7 Figure 2 analysis (IPF AT2 STAT6 glycolysis)
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Running Comment7 Fig 2: IPF AT2 STAT6 glycolysis${NC}"
echo "=============================================="
if [ -f "src/figures9/figures9.R" ]; then
    Rscript src/figures9/figures9.R
    echo -e "${GREEN}Comment7 Fig 2 completed.${NC}"
else
    echo -e "${RED}Warning: src/figures9/figures9.R not found, skipping.${NC}"
fi
echo ""

# ==============================================================================
# Run Comment 7 Figure 3 analysis (Metabolomics)
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Running Comment7 Fig 3: Metabolomics${NC}"
echo "=============================================="
if [ -f "src/figures10/figures10.R" ]; then
    Rscript src/figures10/figures10.R
    echo -e "${GREEN}Comment7 Fig 3 completed.${NC}"
else
    echo -e "${RED}Warning: src/figures10/figures10.R not found, skipping.${NC}"
fi
echo ""

# ==============================================================================
# Summary
# ==============================================================================
echo "=============================================="
echo -e "${GREEN}All analyses completed!${NC}"
echo "=============================================="
echo ""
echo "Output locations:"
echo "  - Figures:    $PROJECT_DIR/figures/"
echo "  - Response:   $PROJECT_DIR/figures/figures1/"
echo "  - Figure 2:   $PROJECT_DIR/figures/figures2/"
echo "  - Figure 3:   $PROJECT_DIR/figures/figures3/"
echo "  - Figure 6:   $PROJECT_DIR/figures/figures6/"
echo "  - Figure 7:   $PROJECT_DIR/figures/figures7/"
echo "  - Processed:  $PROJECT_DIR/data/processed/"
echo ""
echo "To view generated figures:"
echo "  ls -la figures/*/"
echo ""
