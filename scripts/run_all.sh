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

echo "=============================================="
echo "AT2_STAT6_HK2 Analysis Pipeline"
echo "=============================================="
echo ""
echo "Project directory: $PROJECT_DIR"
echo ""

# Change to project directory
cd "$PROJECT_DIR"

# Create output directories
echo -e "${YELLOW}Creating output directories...${NC}"
mkdir -p data/processed
mkdir -p figures/figure1
mkdir -p figures/response_fig1_3
mkdir -p figures/comment7_fig1
mkdir -p figures/comment7_fig2
mkdir -p figures/comment7_fig3
echo -e "${GREEN}Done.${NC}"
echo ""

# Check if conda environment exists (R)
echo -e "${YELLOW}Checking conda environment r_4_3...${NC}"
if conda env list | grep -q "r_4_3"; then
    echo -e "${GREEN}r_4_3 environment found.${NC}"
else
    echo -e "${RED}r_4_3 environment not found. Please create it first:${NC}"
    echo "  conda env create -f environment.yml"
    exit 1
fi
echo ""

# Activate conda environment
echo -e "${YELLOW}Activating r_4_3 environment...${NC}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate r_4_3
echo -e "${GREEN}Done.${NC}"
echo ""

# ==============================================================================
# Use existing scFEA outputs (skip rerun to avoid env crash)
# Copy counts/metadata/flux from original scFEA results into processed/
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Using existing scFEA outputs...${NC}"
echo "=============================================="
SCFEA_DIR="/mnt/d/researchproject/STAT6 and EMT/25_11_08_STAT6KO"
mkdir -p data/processed/scfea_input data/processed/scfea_input/ipf data/processed/scfea_flux_at2 data/processed/scfea_flux_ipf_at2

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
# Run Response Figure 1,3 analysis (Venn diagrams)
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Running Response Fig 1,3: Venn diagrams${NC}"
echo "=============================================="
if [ -f "src/response_fig1_3/response_fig1_3.R" ]; then
    Rscript src/response_fig1_3/response_fig1_3.R
    echo -e "${GREEN}Response Fig 1,3 completed.${NC}"
else
    echo -e "${RED}Warning: src/response_fig1_3/response_fig1_3.R not found, skipping.${NC}"
fi
echo ""

# ==============================================================================
# Run Comment 7 Figure 1 analysis (AT2 glycolysis)
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Running Comment7 Fig 1: AT2 glycolysis${NC}"
echo "=============================================="
if [ -f "src/comment7_fig1/comment7_fig1.R" ]; then
    Rscript src/comment7_fig1/comment7_fig1.R
    echo -e "${GREEN}Comment7 Fig 1 completed.${NC}"
else
    echo -e "${RED}Warning: src/comment7_fig1/comment7_fig1.R not found, skipping.${NC}"
fi
echo ""

# ==============================================================================
# Run Comment 7 Figure 2 analysis (IPF AT2 STAT6 glycolysis)
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Running Comment7 Fig 2: IPF AT2 STAT6 glycolysis${NC}"
echo "=============================================="
if [ -f "src/comment7_fig2/comment7_fig2.R" ]; then
    Rscript src/comment7_fig2/comment7_fig2.R
    echo -e "${GREEN}Comment7 Fig 2 completed.${NC}"
else
    echo -e "${RED}Warning: src/comment7_fig2/comment7_fig2.R not found, skipping.${NC}"
fi
echo ""

# ==============================================================================
# Run Comment 7 Figure 3 analysis (Metabolomics)
# ==============================================================================
echo "=============================================="
echo -e "${YELLOW}Running Comment7 Fig 3: Metabolomics${NC}"
echo "=============================================="
if [ -f "src/comment7_fig3/comment7_fig3.R" ]; then
    Rscript src/comment7_fig3/comment7_fig3.R
    echo -e "${GREEN}Comment7 Fig 3 completed.${NC}"
else
    echo -e "${RED}Warning: src/comment7_fig3/comment7_fig3.R not found, skipping.${NC}"
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
echo "  - Processed:  $PROJECT_DIR/data/processed/"
echo ""
echo "To view generated figures:"
echo "  ls -la figures/*/"
echo ""
