# Single-cell & spatial profiling of human sacrococcygeal teratomas
**Code repository for:** “Integrated single-nuclei and spatial transcriptomic profiling of human sacrococcygeal teratomas reveals heterogeneity in cellular composition and X-chromosome inactivation”

**Authors:** Ernesto J. Rojas, Krinio Giannikou, Benjamin J. Huang, Soo-Jin Cho, Marco A. Cordero, Deion Pena, Lan Vu, Aditya Bagrodia, S. Christopher Derderian, Tippi C. MacKenzie*, Diana J. Laird*  
*Co-corresponding authors: T.C. MacKenzie & D.J. Laird*  
[bioRxiv link](https://doi.org/10.1101/2025.07.21.665156)

---

## Data availability (read first)
- **Raw FASTQs** must email corresponding authors **diana.laird@ucsf.edu, tippi.mackenzie@ucsf.edu, ernesto.rojas@ucsf.edu,**.  We must have an MTA with UCSF and your institute.
- **Processed objects (Seurat/snRNA-seq & Visium), metadata, and helper tables**: [Zenodo DOI](https://doi.org/10.5281/zenodo.17469320).  
- **External references** (ovarian & ESC-derived teratomas): GEO: **GSE229343**, **GSE156170**, respectively.

> **Important:** This GitHub repository does **not** include patient-level or processed `.rds` objects. Download the Zenodo bundle and place files locally as described below.

---

## Folder layout (matches this repo)
```
2025_Rojas_SCTatlas/
├── analysis/
│   ├── 00_preprocessing/        # SoupX, QC, doublets
│   ├── 10_annotation/           # clustering + cell-type labels
│   └── 20_figures_main/         # figure notebooks (1–5)
├── data/
│   ├── raw/                     # Our emails are here for MTA request
│   ├── processed/               # Zenodo Link here too: https://doi.org/10.5281/zenodo.17469320
│   └── metadata_forProcessing/  # gene sets, XCI tables (text/xlsx/gmt)
├── scripts_xchr/                # variant calling helpers (cellsnp-lite, QC)
└── README.md
```

---

## Quick start (use processed data from Zenodo)
1) **Clone this repo**  
   ```bash
   git clone https://github.com/ejscience/2025_Rojas_SCTatlas
   cd 2025_Rojas_SCTatlas
   ```

2) **Install R dependencies** (R 4.3.0). If using `renv`, then later run `renv::restore()` once the lockfile is committed. For now please follow the list at the bottom of this `README.md` in "Software versions"

3) **Download processed objects from Zenodo** you expect:
   
   **Expected files (from Zenodo):**
   - `nucSamples_Harmony_LowResAdded.rds`
   - `SCT04edge.rds`, `SCT04middle.rds`, `SCT05.rds`
   - `MetadataFor_HiResNamedOrdered_SCT01renamed_chromaffinPosChanged_removedRedund_20250331.rds`
   - `metadataForXchr_fromGOGH_20241222.rds`

   **Already included in repo (metadata):**
   - `data/metadata/gencode_v27.-hg38-2.txt`
   - `data/metadata/msigdb.v2023.2.Hs.symbols.gmt`
   - `data/metadata/XLinkedGenes_StatusForXCI_fromTukiainen_2017.xlsx`

4) **Recreate all figure panels (Figures 1–5)**  
   - (run notebooks directly): see commands below.

---

## How to run the analysis (step-by-step)

### A) Preprocessing & QC (snRNA-seq)
Located in `analysis/00_preprocessing/`  
- `20240516_SoupX_Corrections.Rmd` — Ambient RNA removal with SoupX.  
- `20240516-2_Filtering_LowQualityCellsAndDoublets.Rmd` — QC thresholds & doublet removal.

**Output:** filtered Seurat objects (written to `data/processed/`, not tracked in Git).

### B) Cluster annotation & integration
Located in `analysis/10_annotation/`  
- `20240517_SCTsamplesSplit_DGE_Clustering.Rmd` — per-sample clustering & DE.  
- `20240518-1_AzimuthNaming.Rmd` — mapping to fetal atlas (Azimuth).  
- `20240518-2_GPTCelltype_SplitSCTs_NamingClusters.Rmd` — gene-set aided naming.  
- `20240518-3_NamingSplit_AddModuleScore_OtherTeratomas.Rmd` — reference scoring with ovarian/ESC-derived teratomas.  
- `20240522_SplitSCT_CheckingCellTypes.Rmd` — reconciliation across split objects.  
- `20240524_nucSamples_CorrectResolutionAndNamingCellTypes.Rmd` — final resolution & labels.

**Output:** integrated, annotated objects used downstream (these will be the same as `data/processed/` which are on Zenodo).

### C) Figure notebooks (reproduce main figures)
Located in `analysis/20_figures_main/`

**Figure 1 (atlas overview):**  
- `20250331_Figure1_NewObjectHere.Rmd`

**Figure 2 (pseudobulk PCA & stratification):**  
- `20250401_Figure2_PCA_andHeterogeneity.Rmd`

**Figure 3 (spatial deconvolution & mapping):**  
- `20250402_Figure3_SpatialAnalyses_KrinioCode.R`  
- `20250528_SpatialSCT_Pseudobulk_CorrelationAndPCAplot.Rmd`

**Figure 4 (PGC/pluripotency programs):**  
- `20250403_Figure4_PGCandPluripotentGeneExpression.Rmd`

**Figure 5 (X-chromosome inactivation & XaXa-like cells):**  
- `20250404_Figure5_Xchromosome.Rmd`  
- X-chr preprocessing helpers used by Figure 5:  
  - `analysis/20_figures_main/20250110_XChromosome_InferCNV_WITHxaRatios.Rmd` (inferCNV & X/A ratios)  
  - `analysis/20_figures_main/20250330_CallingSNPsFrom_scLinaX_preProcessing.Rmd` (cellsnp-lite preprocessing)


### D) X-chromosome variant calling helpers (these need to be run before Figure 5 when starting from scratch)
Located in `scripts/xchr/`  
- `RunningScriptsForVariantCalling.txt` — overview of the pipeline.  
- `run_cellsnp-lite.sh`, `move_barcodes_for_cellsnp_lite.sh` — CLI helpers.  
- `RCODE_process_cellsnp.r`, `RCODE_QC.r` — post-processing & QC.  

> Requires aligned BAMs and barcodes; these are not distributed here. See manuscript Methods for details.

---

## Software versions
Pinned in the paper; typical working versions below.

| Package / Tool          | Version |
|-------------------------|---------|
| R                       | 4.3.0   |
| SoupX                   | 1.6.2   |
| Seurat                  | 5.3.0   |
| scCustomize             | 3.0.1   |
| DoubletFinder           | 2.0.4   |
| scDblFinder             | 1.16.0  |
| SingleCellExperiment    | 1.24.0  |
| ACT                     | N/A     |
| Harmony                 | 1.2.0   |
| clustree                | 0.5.1   |
| infercnv                | 1.19.1  |
| scProportionTest        | 0.0.0.9 |
| perturbLM               | 1.0.0   |
| entropy                 | 1.3.1   |
| Semla                   | 1.1.6   |
| CARD                    | 1.1     |
| biomaRt                 | 2.58.2  |
| qvalue                  | 2.34.0  |
| msigdbr                 | 7.5.1   |
| clusterProfiler         | 4.10.1  |
| org.Hs.eg.db            | 3.18.0  |
| DESeq2                  | 1.42.1  |
| ggplot2                 | 3.5.1   |
| ggprism                 | 1.0.5   |
| circlize                | 0.4.16  |
| ComplexHeatmap          | 2.18.0  |
| ggpattern               | 1.1.1   |
| viridis                 | 0.6.5   |
| scales                  | 1.3.1   |
| ggrepel                 | 0.9.5   |
| CellRanger              | 7.0.1   |
| spaceranger             | 3.0.1   |
| LoupeBrowser            | 7.0.1   |
| cellsnp-lite            | 1.2.3   |
| Annovar                 | 2025Mar21 |

---

## Citation
Please cite:
Rojas EJ, Giannikou K, Huang BJ, Cho SJ, Cordero MA, Pena D, Vu L, Bagrodia A, Derderian SC, MacKenzie TC*, Laird DJ*. “Integrated single-nuclei and spatial transcriptomic profiling of human sacrococcygeal teratomas reveals heterogeneity in cellular composition and X-chromosome inactivation.” (2025).

## Contact
- ernesto.rojas@ucsf.edu  
- diana.laird@ucsf.edu  
- tippi.mackenzie@ucsf.edu

_Last updated: 2025-10-28_
