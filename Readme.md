# Spatial and Single-Cell Alternative Polyadenylation Analysis in Thyroid Cancer

![Status](https://img.shields.io/badge/status-in%20progress-yellow)

Pilot analysis pipeline studying alternative polyadenylation (APA) events in papillary thyroid carcinoma (PTC) using paired scRNA-seq and Visium spatial
transcriptomics data.

## Dataset

- **Source:** [GSE250521](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250521) (10x Genomics, Single Cell 3' v3.1 chemistry — scRNA-seq; 10x Visium — spatial)
- **Pilot design:** 3 Normal thyroid samples vs. 3 PTC samples (unpaired, independent patients: biological replicates, not patient-matched)
- **Scope note:** LPTC and ATC samples exist in the full dataset but are excluded from this pilot to keep the initial comparison clean (Normal vs. PTC only)

## Current status

This is an active work in progress. The scRNA-seq QC and integration pipeline is functional through initial clustering; malignant cell calling, APA quantification, and the spatial track are not yet implemented.

### Completed
- [x] Per-sample data loading (`Read10X`, GSM-prefixed barcodes preserved)
- [x] Per-sample QC via [scQCenrich](https://github.com/lemonlyy755/scQCenrich) (GMM-based adaptive thresholding, cluster-coherence-aware)
- [x] Per-sample doublet detection and removal (`scDblFinder`)
- [x] Merge into single Seurat object (doublet and low-quality-cell-free)
- [x] Normalization, variable feature selection (n=3000), scaling, PCA (30 PCs)
- [x] Batch/sample integration (Harmony)
- [x] Clustering (resolution 0.3 → 19 clusters) and `FindAllMarkers`

### In progress / next
- [ ] Cluster annotation (canonical thyroid + immune/stromal markers)
- [ ] Malignant cell calling on epithelial cluster(s) via CNV inference (inferCNV/CopyKAT)
- [ ] Cell-type composition comparison (Normal vs. PTC)
- [ ] Pseudobulk differential expression (patient-level aggregation, n=3 per group)
- [ ] FASTQ → BAM alignment via Galaxy/STARsolo (CB/UB-tagged BAM for APA calling)
- [ ] Single-cell APA quantification (Sierra / SCAPE)
- [ ] Visium spatial pipeline (Seurat spatial workflow, deconvolution via scRNA reference)
- [ ] Spatial APA quantification and Normal vs. PTC comparison


## Tools used

`Seurat` (v5) · `scQCenrich` · `scDblFinder` · `SoupX` · `harmony` · `dplyr` ·
(planned: `STARsolo`, `Sierra`/`SCAPE`, `inferCNV`/`CopyKAT`, spatial deconvolution TBD)

---
*This README reflects pipeline state as of the pilot analysis phase and will be updated as the malignant calling, APA, and spatial tracks are implemented.*