# CPTAC Pan-Cancer PTMs

Repo for manuscript:
Geffen, Y., Anand, S., Akiyama, Y., Yaron, T., Song, Y., et al. "Pan-Cancer analysis of post-translational modifications reveals shared patterns of protein regulation." Cell (2023).

[[Link](https://www.cell.com/cell/fulltext/S0092-8674(23)00781-X)]

__Contact__: Shankara Anand - sanand@broadinstitute.org, Yo Akiyama - yakiyama@broadinstitute.org

---

**Repo Structure**:

`Fig1`:
* `01*_process_*`: processing files that take raw proteomic data from proteomic data commons and perform filtering, imputation for clustering, and PTM correction for protein abundance.
* `02_process_rna.ipynb`: notebook for processing transcriptomic data for all CPTAC cohorts. Transcriptomic data were processed using the [GTEx pipeline](https://github.com/broadinstitute/gtex-pipeline).
* `03_integrate_omics.ipynb`: notebook for integrating RNA, protein, and phosphoprotein for use in ARD-NMF clustering. These are then fed into a script to run SignatureAnalyzer, a bayesian variant of non-negative matrix factorization found [here](https://github.com/getzlab/CPTAC_PanCan_PTM_2023/blob/master/scripts/ardnmf.sh) in the **scripts** folder.
* `04_processing_clustering.ipynb`: processing results from SignatureAnalyzer.
* `05_pathways.ipynb`: pathway enrichment for each data modality in each signature.
* `Fig1ACD`: notebook for re-creating Figure 1 A,C,D.
  
`Fig2`:
* `01_cluster_samples.ipynb`: notebook for clustering CPTAC dataset based on signature weights.
* `Fig2A`: notebook for creating overview, signature figure.
* `Fig2F`: violin plots for comparison.

`Fig3` (visit directory for details):
* `01_mutational_signatures_extraction.ipynb`: Extracting mutational signatures
* `02_signature_thresholding_and_classification.ipynb`: Defining mutational signature thresholds for classifying DNA repair deficiencies
* `03_recluster_ddr_deficient.ipynb`: Clustering DNA repair deficient cancers with homologous recombination deficiency (HRD) and mismatch repair deficiency (MMRD)
* `04_HRD_Differential_Expression.ipynb`:  HRD differential expression analysis (global HRD vs HRP and acute vs chronic hypoxia HRD)
* `05_MMRD_Differential_Expression.ipynb`: MMRD vs MMRP differential expression analysis
* `06_MRN_complex_analysis_ccle.ipynb`: MRN complex proteomic and transcriptomic analysis in CCLE cell lines

`Fig4`: immune analyses

`Fig5`: histone analyses

`data`: see `./data/README.md`.

`downsampling`: folder with downsampling analyses to demonstrate robustness of clustering approach.

`funcs`: folder with python/R functions used throughout notebooks.

`scripts`: folder with scripts used throughout analyses.
