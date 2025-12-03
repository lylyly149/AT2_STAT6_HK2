# AT2_STAT6_HK2: STAT6 Pathway Analysis in IPF AT2 Cells

This repository contains the analysis code for studying STAT6 pathway activation and glycolysis metabolism in pulmonary fibrosis (PF) alveolar type II (AT2) cells.

## Project Description

Pulmonary fibrosis (PF) is a progressive fibrotic lung disease with limited therapeutic options. This study investigates the role of STAT6 signaling pathway in AT2 cell dysfunction and metabolic reprogramming during IPF progression. We integrated single-cell RNA sequencing (scRNA-seq), bulk RNA-seq, metabolic flux analysis (scFEA), and targeted metabolomics to comprehensively characterize STAT6-mediated glycolysis alterations in IPF.

Key findings include:
- Elevated STAT6 pathway activity in IPF AT2 cells
- Association between STAT6 activation and glycolysis metabolic flux
- Validation in STAT6 knockout mouse models

## Folder Structure

```
AT2_STAT6_HK2/
|
+-- README.md                    # Project documentation
+-- LICENSE                      # MIT License
+-- environment.yml              # Conda environment for R 4.3
+-- requirements.txt             # Python dependencies for scFEA (py310/py311)
|
+-- config/
|   +-- paths.yaml               # Parameterized data paths (commit-safe template)
|   +-- paths.local.yaml         # Local overrides (not tracked)
|
+-- src/
|   +-- figure1/                 # Main scRNA-seq figure
|   |   +-- figure1.R
|   +-- figures1/                # Response figures + PF/network analyses
|   |   +-- figures1a-1.R
|   |   +-- figures1b_run_limma_gse134692_from_rawcounts.R
|   |   +-- figures1c_run_limma_gse32147_16w_from_series.R
|   |   +-- figures1d_run_limma_gse16846_from_series.R
|   |   +-- figures1j_k.R        # PF-target heatmaps + network/TF analyses
|   +-- figures2/                # Epithelial re-analysis (dotplots/tSNE/DEG+GO)
|   |   +-- figures2a.R
|   |   +-- figures2bc.R
|   |   +-- figures2d.R
|   +-- figures3/                # Macrophage & fibroblast STAT6 violins
|   |   +-- figures3a_b_e_f_g.R
|   |   +-- figures3c.R
|   +-- figures6/                # Monocle pseudotime heatmap + STAT6 alluvial
|   |   +-- figures6.R
|   +-- figures7/                # Spatial STAT6/EMT/glycolysis maps
|   |   +-- figures7.R
|   +-- figgures8/               # Comment 7 figure 1 (scFEA glycolysis AT2)
|   |   +-- figgures8.R
|   +-- figures9/                # Comment 7 figure 2 (STAT6-stratified flux)
|   |   +-- figures9.R
|   +-- figures10/               # Comment 7 figure 3 (metabolomics STAT6 KO)
|   |   +-- figures10.R
|   +-- scfea/                   # scFEA exports/runs
|   |   +-- export_scfea_inputs.R
|   |   +-- run_scfea.py
|   +-- utils/Function_for_violin_plot.R
|
+-- scripts/
|   +-- run_all.sh               # One-click reproducibility script
|
+-- data/                        # User data workspace (git-ignored)
+-- figures/                     # Generated figures (git-ignored)
+-- logs/                        # Optional logs (git-ignored)
```

## Environment Setup

### R Environment (for scRNA-seq and statistical analysis)

Create the conda environment using:

```bash
conda env create -f environment.yml
conda activate r_4_3
```

Or manually install R 4.3+ with the following packages:
- Seurat (>= 5.0)
- GSVA
- AUCell
- monocle
- limma
- edgeR
- GEOquery
- ggpubr
- ggrastr
- ggsci
- ggalluvial
- patchwork
- msigdbr
- dynamicTreeCut
- WGCNA
- TFTF
- GOSemSim
- spatstat.geom / spatstat.core
- pheatmap
- VennDiagram
- biomaRt
- tidyverse
- ggplot2
- ggrepel
- pals
- viridis
- yaml

### Python Environment (for scFEA)

```bash
conda create -n py311 python=3.11  # or match config/paths.local.yaml: conda_envs.py_env
conda activate py311
pip install -r requirements.txt
```

## Reproducibility Instructions

### Step 1: Configure data paths

Fill `config/paths.local.yaml` with real paths (keep `paths.yaml` in sync for scripts that read it directly, e.g. `figures9.R`). Key entries:
- `raw_data.seurat_all`, `seurat_at2`, `seurat_epi`, `seurat_ipf`, `seurat_ipfat2`
- `raw_data.stat6_geneset`
- `raw_data.violin_plot_functions`
- `raw_data.scfea_dir` (used by `run_all.sh` to copy existing scFEA outputs)
- `raw_data.scfea_module_info`
- `raw_data.network.edge_file`
- `raw_data.spatial_sio2_56` (spatial dataset for Figure 7; falls back to `../spatial trans/dataset/SIO2_56`)
- Optional: `figures.figures7` if you want to place spatial outputs elsewhere
- PF target / GEO raw paths used in `figures1j_k.R` are currently absolute; adjust there if your layout differs.

### Step 2: Prepare scFEA outputs (Python)

- Install scFEA: `pip install -r requirements.txt && pip install git+https://github.com/changwn/scFEA.git`
- Export inputs (R): `Rscript src/scfea/export_scfea_inputs.R --config config/scfea_export_at2.yaml` and `..._ipf.yaml`
- Run scFEA (Python): `python src/scfea/run_scfea.py --config config/scfea_run_at2.yaml` and `..._ipf.yaml`
- If you already have flux outputs, place them in `raw_data.scfea_dir`; `scripts/run_all.sh` copies `data/interim/scfea_input/*` and `results/scfea_flux_*/*` into `data/processed/`.

### Step 3: Run all analyses (R)

```bash
chmod +x scripts/run_all.sh
./scripts/run_all.sh
```

The script reads `config/paths.local.yaml` (fallback `paths.yaml`) to pick the conda env name, scFEA source directory, and required inputs, then runs:
- Figure 1 (scRNA-seq overview)
- Figure 2 (epithelial re-analysis)
- Figure 3 (macrophage/fibroblast STAT6 violins)
- Response Figures 1/3 (limma reruns, Venns, PF-target heatmaps/network `figures1j_k.R`)
- Figure 6 (Monocle pseudotime heatmap + STAT6 alluvial)
- Figure 7 (spatial STAT6/EMT/glycolysis maps)
- Comment 7 Figures 1–3 (scFEA + metabolomics)

### Step 4: Run scripts individually (optional)

```bash
# Figure 1
Rscript src/figure1/figure1.R

# Figure 2
Rscript src/figures2/figures2a.R
Rscript src/figures2/figures2bc.R
Rscript src/figures2/figures2d.R

# Figure 3
Rscript src/figures3/figures3a_b_e_f_g.R
Rscript src/figures3/figures3c.R

# Response figures
Rscript src/figures1/figures1b_run_limma_gse134692_from_rawcounts.R
Rscript src/figures1/figures1d_run_limma_gse16846_from_series.R
Rscript src/figures1/figures1c_run_limma_gse32147_16w_from_series.R
Rscript src/figures1/figures1a-1.R
Rscript src/figures1/figures1j_k.R   # PF-target heatmaps + network/TF analyses

# Figure 6 / Figure 7
Rscript src/figures6/figures6.R
Rscript src/figures7/figures7.R

# Comment 7 figures
Rscript src/figgures8/figgures8.R
Rscript src/figures9/figures9.R
Rscript src/figures10/figures10.R
```

### Outputs

- `figures/figure1/` – main scRNA-seq figure panels.
- `figures/figures1/` – response figures, limma volcanoes, Venns, PF-target heatmaps, and `figures1jk_*` network outputs (TF scores, edge filters, module info, STAT6-annotated network).
- `figures/figures2/` – cell-type dotplots/tSNE plus `figures2d_deg_ipf_vs_nom_at2.csv` and GO results.
- `figures/figures3/` – macrophage UMAP/dotplot and STAT6 violins for macrophage subtypes and fibroblasts.
- `figures/figures6/` – `figure6_pseudotime_heatmap.pdf`, `state_DEGs_pseudotime.csv`, `figure6_alluvial_STAT6_positive_clusters.pdf`.
- `figures/figures7/` – spatial feature maps, STAT6/EMT/glycolysis scatter, and crossK/PCF plots.
- `figures/figgures8/`, `figures/figures9/`, `figures/figures10/` – comment 7 figure sets (glycolysis flux, STAT6-stratified flux, metabolomics).
- `data/processed/` – copied scFEA inputs/fluxes, scFEA differential flux tables (`scfea_flux_diff_ipf_at2_stat6p_all.csv`, `..._sig.csv`), and other intermediates.

## Data Availability

### Public Datasets

- **GSE122960**: scRNA-seq of IPF and healthy lung (GEO)
- **GSE134692**: Bulk RNA-seq of IPF lung tissue
- **GSE16846**: Bulk RNA-seq of bleomycin-induced pulmonary fibrosis (mouse)
- **GSE32147**: Bulk RNA-seq of crystalline silica–induced pulmonary fibrosis (rat)

### Study-specific Data

- Targeted metabolomics data from STAT6 KO mice: Available upon request to the corresponding contact
- Large scFEA inputs/flux matrices and Seurat objects are not tracked in git to keep the repository lightweight.
- This repository is for sharing analysis code for the scientific manuscript; fully processed Seurat objects (QC, normalization, dimensionality reduction, clustering, annotations) are available from the corresponding author upon reasonable request.

## Citation

If you use this code, please cite (manuscript under review):

> STAT6-mediated glycolysis reprogramming in AT2 cells contributes to idiopathic pulmonary fibrosis. [Journal]. [Year].

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions about the analysis code, please contact the corresponding author.
