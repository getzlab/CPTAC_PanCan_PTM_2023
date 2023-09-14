# CPTAC Pan-Cancer PTMs

Repo for manuscript:
Geffen, Y., Anand, S., Akiyama, Y., Yaron, T., Song, Y., et al. "Pan-Cancer analysis of post-translational modifications reveals shared patterns of protein regulation." Cell (2023).

[[Link](https://www.cell.com/cell/fulltext/S0092-8674(23)00781-X)]

__Contact__: Shankara Anand - sanand@broadinstitute.org, Yo Akiyama - yakiyama@broadinstitute.org

---

**Repo Structure**:

`data` - see `./data/README.md`

`funcs` - folder with python/R functions and scripts

`Fig3` (visit directory for details)
* `01_mutational_signatures_extraction.ipynb`: Extracting mutational signatures
* `02_signature_thresholding_and_classification.ipynb`: Defining mutational signature thresholds for classifying DNA repair deficiencies
* `03_recluster_ddr_deficient.ipynb`: Clustering DNA repair deficient cancers with homologous recombination deficiency (HRD) and mismatch repair deficiency (MMRD)
* `04_HRD_Differential_Expression.ipynb`:  HRD differential expression analysis (global HRD vs HRP and acute vs chronic hypoxia HRD)
* `05_MMRD_Differential_Expression.ipynb`: MMRD vs MMRP differential expression analysis
* `06_MRN_complex_analysis_ccle.ipynb`: MRN complex proteomic and transcriptomic analysis in CCLE cell lines
