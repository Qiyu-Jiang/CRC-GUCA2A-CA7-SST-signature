# SUBMISSION_SUPP - supplementary package README

CRC three-gene signature (GUCA2A/CA7/SST) study - supplementary materials and code for journal upload.
Generated 2026-09-10; last updated 2026-10-08 (added Tables S36-S37 and the post-review sensitivity scripts 14/14b; all numbers verified against the project archives).

## Supplementary_Figures (17 files)
- Figure S1: PCA before/after ComBat batch correction (training cohort, top-2000 MAD genes)
- Figure S2: ANN in-sample ROC (reproduction of the reported AUC 0.987)
- Figure S3: ANN stratified 10-fold cross-validation ROC (AUC 0.907, 95% CI 0.830-0.984)
- Figure S4: ANN external validation ROC on GSE110224 (AUC 0.820) with logistic comparator (0.841)
- Figure S5: Feature-selection frequency over 100 repeats (80% stratified subsamples; LASSO + RF)
- Figure S6: ssGSEA boxplot, Senbabaoglu 20-cell signature (Control vs Tumour)
- Figure S7: Spearman heatmap, feature genes vs Senbabaoglu 20 cell types
- Figure S8: ssGSEA boxplot, Charoentong 28-cell signature (Control vs Tumour)
- Figure S9: Spearman heatmap, feature genes vs Charoentong 28 cell types
- Figure S10: Immune-checkpoint correlations, training cohort (15 checkpoint genes)
- Figure S11: Immune-checkpoint correlations, TCGA tumour cohort
- Figure S12: TCGA overall survival, KM curves for GUCA2A/CA7/SST (median split, n=242)
- Figure S13: TCGA non-synonymous mutations and GISTIC CNA summary barplot
- Figure S14: TCGA tumour expression boxplot (RSEM, three signature genes)
- Figure S15: CellChat rankNet comparison, EpiSigHigh vs EpiSigLow
- Figure S16: CellChat bubble plot, epithelial subsets as signal senders
- Figure S17: CellChat pathway-level communication difference (top 12)

## Supplementary_Tables (40 files; S1-S38 complete as CSVs)
- Table S38: HPA protein-level evidence for the three signature genes (antibody IDs, reliability, normal-colon cell-type staining, CRC IHC counts) - new, added 2026-10-19
- Table S23: Fully nested re-selection CV, held-out metrics (10 folds)
- Table S24: Features selected within each fold of the nested re-selection CV
- Table S28: Fixed-gene repeated CV summary (ANN 0.912 ± 0.011 / logistic 0.940 ± 0.003)
- Table S29: Fixed-gene repeated CV, per-partition metrics (5 seeds)
- Table S27: Full nested-CV held-out predictions (150 rows)
- Table S36: Checkpoint correlations with BH adjustment (both cohorts, 87 tests)
- Table S37: Frozen-z-score external sensitivity (ANN 0.706 / logistic 0.848)
- Note: Table_S33 = Table_S9 and Table_S33b = Table_S11 (identical files kept for numbering continuity)
- Table S25: Ten-seed consensus, full training set (43-candidate pool)
- Table S26: Ten-seed consensus under 80% downsampling
- Table S1: Dataset overview (GEO accessions, platforms, sample counts)
- Table S2: Leave-one-dataset-out sensitivity (3 genes, logFC/adjP per exclusion)
- Table S3: Downsampling stability - feature selection (3 seeds, fixed-RF version)
- Table S4: Downsampling stability - DEG of the three genes (3 seeds)
- Table S5: Seeded feature lists on the full training set (3 seeds)
- Table S6: ANN performance summary (repro/CV/external, AUC + DeLong CI)
- Table S7: ANN external predictions at sample level (GSE110224)
- Table S8: Feature-selection frequency (LASSO/RF, 100 repeats)
- Table S9/S10: Senbabaoglu ssGSEA scores / correlations
- Table S11/S12: Charoentong ssGSEA scores / correlations
- Table S13/S14: Immune-checkpoint correlations (training / TCGA)
- Table S15: TCGA non-synonymous mutations and GISTIC CNA (all zero for the three genes)
- Table S16: TCGA KM summary (log-rank p, HR)
- Table S17: TCGA GISTIC details (molecular-data endpoint output)
- Table S18/S19: UALCAN promoter-methylation boxstats / medians
- Table S33b: per-sample ssGSEA scores, Charoentong 28-cell (per-sample variant of S11)
- Table S20/S21: CellChat pathway diff / pathway probability summary

## Supplementary_Code (analysis scripts)
R 4.6.1 scripts (numbered in pipeline order) + Python helpers:
01-03: ANN validation, sensitivity analysis (LODO/downsampling/seeds), feature stability
04: selection frequency (100 repeats)
05/05b: ssGSEA with two immune signatures (GSVA)
06: checkpoint correlations
07: TCGA KM/mutation/CNA (cBioPortal API data)
08/08b: CellChat pipeline and merge comparison
09: UALCAN methylation boxplots from embedded Highcharts stats
10-12: Python helpers (cBioPortal API pull, Charoentong supplement fetch/parse)
main_pipeline_fixed__: the corrected original pipeline scripts (paths, seeds, input names)

## Environment
R 4.6.1 (limma, sva, WGCNA, glmnet, e1071, randomForest, neuralnet, pROC, GSVA, Seurat 5.5.1, CellChat 2.1.2, survival, pheatmap);
Python 3.13 (requests, openpyxl, pypdf). AutoDock Vina 1.2.7 standalone (download: github.com/ccsb-scripps/AutoDock-Vina/releases).

## Notes
- GDC (api.gdc.cancer.gov) returned persistent 500 errors during the retry window; the raw-count matrix was replaced by
  GEPIA2 screenshots (Figure panels) + cBioPortal data. Re-run 10_cbioportal_api_pull.py + a GDC manifest query when the server recovers.
- cellchat_diff.pdf is a placeholder (netVisual_diff is unavailable in CellChat 2.1.2); rankNet and bubble cover the comparison.
- Browser-verification screenshots (GEPIA2/UALCAN source pages) are kept under project archive\06.* and are not part of this package.

## Citation

If you use these materials, please cite the associated manuscript: "Loss of the enterocyte-differentiation genes GUCA2A, CA7 and SST as a diagnostic signature for colorectal cancer" (under review).

Zenodo DOI: to be minted via the Zenodo-GitHub integration after the first release (v1.0.0 tag is prepared for this purpose).
