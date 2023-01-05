## Data for CPTAC Pan-Cancer PTM

Authors: Yo Akiyama (yakiyama@broadinstitute.org), Shankara Anand (sanad@broadinstitute.org)

**Geneset Dir** (`./genesets`)
* GMT files from [MSigDB](http://www.gsea-msigdb.org/gsea/msigdb/) and [Reactome] (https://reactome.org/)
* ReactomePathways_ddr.gmt: Reactome DNA repair pathways
* c2.cp.kegg.v7.4.symbols.gmt: KEGG genesets
* h.all.v7.4.symbols.gmt: Hallmark genesets
* full_geneset.json: Dictionary storing the CPTAC DNA Damage Response Geneset (from companion CPTAC DDR paper), kinases, phosphatases, acetyltransferases, and deacetylases

**Raw data** (`./raw`)
* tumor_rna_counts_pc_X.tsv: RNA read counts for protein coding genes (GENCODE v34 annotation)
* * tumor_rna_tpm_X.parquet: RNA TPM counts
* proteome_gene_centric.tsv.gz: Gene-level relative protein abundance
* phosphoproteome_X.tsv.gz: Relative phosphopeptide abundance
* acetylome_X.tsv.gz: Relative acetylpeptide abundance

**PTM corrected for protein abundance ("residuals")** (`./ptm_residuals`)
* phosphoproteome_raw_res_X.tsv
* acetylome_raw_res_X.tsv

**Other files** (`./`)
* union_somatic_maf_v1.1_hotspot_cptac2fix.maf: Pan-Cancer harmonized WXS mutation annotation formatted file

---

For details on CPTAC WXS harmonization, visit the [accompanying repo](https://github.com/getzlab/cptac_wxs_harmonize/tree/master)
* Authors: Qing Zhang (qzhang9@illumina.com), Yo Akiyama (yakiyama@broadinstitute.org)
