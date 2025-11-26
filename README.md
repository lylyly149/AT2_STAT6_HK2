# AT2_STAT6_HK2: STAT6 Pathway Analysis in IPF AT2 Cells

This repository contains the analysis code for studying STAT6 pathway activation and glycolysis metabolism in idiopathic pulmonary fibrosis (IPF) alveolar type II (AT2) cells.

## Project Description

Idiopathic pulmonary fibrosis (IPF) is a progressive fibrotic lung disease with limited therapeutic options. This study investigates the role of STAT6 signaling pathway in AT2 cell dysfunction and metabolic reprogramming during IPF progression. We integrated single-cell RNA sequencing (scRNA-seq), bulk RNA-seq, metabolic flux analysis (scFEA), and targeted metabolomics to comprehensively characterize STAT6-mediated glycolysis alterations in IPF.

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
+-- requirements.txt             # Python dependencies for scFEA (py310)
|
+-- config/
|   +-- paths.yaml               # Parameterized data paths
|
+-- src/
|   +-- figure1/                 # Main figure 1 analysis
|   |   +-- figure1.R            # scRNA-seq analysis of STAT6 in IPF
|   |
|   +-- response_fig1_3/         # Reviewer response figures
|   |   +-- response_fig1_3.R    # Venn diagrams of DEG intersections
|   |
|   +-- comment7_fig1/           # Comment 7 figure 1
|   |   +-- comment7_fig1.R      # scFEA glycolysis in AT2 (IPF vs NOM)
|   |
|   +-- comment7_fig2/           # Comment 7 figure 2
|   |   +-- comment7_fig2.R      # scFEA glycolysis by STAT6 status
|   |
|   +-- comment7_fig3/           # Comment 7 figure 3
|       +-- comment7_fig3.R      # Metabolomics analysis (STAT6 KO)
|
+-- scripts/
|   +-- run_all.sh               # One-click reproducibility script
|
+-- data/                        # Data directory
|   +-- processed/               # Intermediate results (generated; large scFEA inputs/flux and .rds files are git-ignored)
|
+-- figures/                     # Output figures (generated locally; git-ignored)
    +-- figure1/                 # Created by scripts
    +-- response_fig1_3/         # Created by scripts
    +-- comment7_fig1/           # Created by scripts
    +-- comment7_fig2/           # Created by scripts
    +-- comment7_fig3/           # Created by scripts
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
- pheatmap
- VennDiagram
- biomaRt
- tidyverse
- ggplot2
- ggrepel
- patchwork
- pals
- viridis
- yaml

### Python Environment (for scFEA)

```bash
conda create -n py310 python=3.10
conda activate py310
pip install -r requirements.txt
```

## Reproducibility Instructions

### Step 1: Configure data paths

Edit `config/paths.yaml` to set the correct paths to your raw data files:
- scRNA-seq H5 files (GSE122960)
- Bulk RNA-seq differential expression results
- scFEA metabolic flux outputs
- Metabolomics data

### Step 2: Run scFEA (py310) — optional if you already have flux outputs

1) Install scFEA into the py310 env (once):

```bash
conda activate py310
pip install -r requirements.txt
pip install git+https://github.com/changwn/scFEA.git
```

2) Export scFEA inputs from Seurat objects (R env):

```bash
conda activate r_4_3
Rscript src/scfea/export_scfea_inputs.R --config config/scfea_export_at2.yaml
Rscript src/scfea/export_scfea_inputs.R --config config/scfea_export_ipf.yaml
```

3) Run scFEA (py310):

```bash
conda activate py310
python src/scfea/run_scfea.py --config config/scfea_run_at2.yaml
python src/scfea/run_scfea.py --config config/scfea_run_ipf.yaml
```

Flux and balance outputs are written under `data/processed/scfea_flux_at2/` and `data/processed/scfea_flux_ipf_at2/` (git-ignored). Adjust paths in the YAML configs as needed.

### Step 3: Run all analyses

```bash
chmod +x scripts/run_all.sh
./scripts/run_all.sh
```

Or run individual analysis scripts:

```bash
# Figure 1: scRNA-seq analysis
conda activate r_4_3
Rscript src/figure1/figure1.R

# Response figures: Venn diagrams
Rscript src/response_fig1_3/response_fig1_3.R

# Comment 7 figures: Glycolysis analysis
Rscript src/comment7_fig1/comment7_fig1.R
Rscript src/comment7_fig2/comment7_fig2.R
Rscript src/comment7_fig3/comment7_fig3.R
```

### Step 4: Outputs

Generated figures will be saved locally to the `figures/` directory (git-ignored, regenerate with the scripts):

| Figure | Description | Output Location |
|--------|-------------|-----------------|
| Figure 1A | UMAP of cell types | `figures/figure1/figure1a_*.pdf` |
| Figure 1B | Epithelial tSNE | `figures/figure1/figure1b_*.pdf` |
| Figure 1C | STAT6 GSVA violin | `figures/figure1/figure1c_*.pdf` |
| Figure 1D | STAT6 group ratio | `figures/figure1/figure1d_*.pdf` |
| Figure 1E | IPF epithelial clusters | `figures/figure1/figure1e_*.pdf` |
| Figure 1F | Marker gene dotplot | `figures/figure1/figure1f_*.pdf` |
| Figure 1G | STAT6 by cluster | `figures/figure1/figure1g_*.pdf` |
| Figure 1H | AUC signature | `figures/figure1/figure1h_*.pdf` |
| Figure 1I | Monocle pseudotime | `figures/figure1/figure1i_*.pdf` |
| Response Fig 1,3 | Venn diagrams | `figures/response_fig1_3/*.pdf` |
| Comment7 Fig 1A-C | AT2 glycolysis | `figures/comment7_fig1/*.pdf` |
| Comment7 Fig 2A-E | IPF AT2 STAT6 glycolysis | `figures/comment7_fig2/*.pdf` |
| Comment7 Fig 3A-C | Metabolomics | `figures/comment7_fig3/*.pdf` |

## Data Availability

### Public Datasets

- **GSE122960**: scRNA-seq of IPF and healthy lung (GEO)
- **GSE134692**: Bulk RNA-seq of IPF lung tissue
- **GSE16846**: Bulk RNA-seq of bleomycin-induced pulmonary fibrosis (mouse)
- **GSE32147**: Bulk RNA-seq of cigarette smoke-induced pulmonary fibrosis (rat)

### Study-specific Data

- Targeted metabolomics data from STAT6 KO mice: Available upon request to the corresponding contact
- Large scFEA inputs/flux matrices and Seurat objects are not tracked in git to keep the repository lightweight.
- Fully processed Seurat objects (QC, normalization, dimensionality reduction, clustering, annotations) are available from the corresponding author upon reasonable request.

## Citation

If you use this code, please cite:

> STAT6-mediated glycolysis reprogramming in AT2 cells contributes to idiopathic pulmonary fibrosis. [Journal]. [Year].

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions about the analysis code, please open an issue on GitHub or contact the corresponding contact.
