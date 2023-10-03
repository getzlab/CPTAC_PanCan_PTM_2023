## Data for CPTAC Pan-Cancer PTM

Authors: Yo Akiyama (yakiyama@broadinstitute.org), Shankara Anand (sanad@broadinstitute.org)

**Reference files** (`./reference`)
* `Readme_RefSeq.20180629.txt`: Readme file with information on how the reference fasta for this CPTAC project was created.
* `RefSeq.20180629_Human_ucsc_hg38_cpdbnr_mito_264contams_553smORFs.fasta`: Reference fasta used for the CPTAC project.
* `var_map_full_v4.tsv`: Processed features collected across all CPTAC data modalities for this project (RNA, proteome, PTM).

**Geneset Dir** (`./genesets`)
* GMT files from [MSigDB](http://www.gsea-msigdb.org/gsea/msigdb/) and [Reactome] (https://reactome.org/)
* ReactomePathways_ddr.gmt: Reactome DNA repair pathways
* c2.cp.kegg.v7.4.symbols.gmt: KEGG genesets
* h.all.v7.4.symbols.gmt: Hallmark genesets
* full_geneset.json: Dictionary storing the CPTAC DNA Damage Response Geneset (from companion CPTAC DDR paper), kinases, phosphatases, acetyltransferases, and deacetylases

**Other files** (`./`)
* union_somatic_maf_v1.1_hotspot_cptac2fix.maf: Pan-Cancer harmonized WXS mutation annotation formatted file

---

For details on CPTAC WXS harmonization, visit the [accompanying repo](https://github.com/getzlab/cptac_wxs_harmonize/tree/master)
* Authors: Qing Zhang (qzhang9@illumina.com), Yo Akiyama (yakiyama@broadinstitute.org)
