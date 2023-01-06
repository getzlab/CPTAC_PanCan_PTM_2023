#!/bin/bash

# -----------------------------
# Runs for Dendrogroups
# -----------------------------
for NUM in 1 2 3 4; do
  clumpsptm -i ../../../analysis/diffexp/061721_raw_dendro/full_de_cohort_cov.tsv \
            -w gsea_rank --maps ../ref_uniprot_072522/mapped_sites_to_pdbs.tsv --pdbstore ../../../../getzlab-CLUMPS2/clumps/db/ref/pdbs \
            --grouping $NUM --features acetylome phosphoproteome --threads 96 -v --subset positive --output_dir dendro_${NUM}_positive

  clumpsptm -i ../../../analysis/diffexp/061721_raw_dendro/full_de_cohort_cov.tsv \
            -w gsea_rank --maps ../ref_uniprot_072522/mapped_sites_to_pdbs.tsv --pdbstore ../../../../getzlab-CLUMPS2/clumps/db/ref/pdbs \
            --grouping $NUM --features acetylome phosphoproteome --threads 96 -v --subset negative --output_dir dendro_${NUM}_negative
done

# -----------------------------
# Runs for Immune Groups
# -----------------------------
for GROUP in ImmuneHot ImmuneWarm ImmuneCool ImmuneCold; do
  clumpsptm -i ../../../analysis/Fig_immuno_metabolism/diffexp/061721_raw_res_immune_ciber/full_de_cohort_cov.parquet \
            -w gsea_rank --maps ../ref_uniprot_072522/mapped_sites_to_pdbs.tsv --pdbstore ../../../../getzlab-CLUMPS2/clumps/db/ref/pdbs \
            --grouping $GROUP --features acetylome phosphoproteome --threads 96 -v --subset positive --output_dir ${GROUP}_ciber_positive

  clumpsptm -i ../../../analysis/Fig_immuno_metabolism/diffexp/061721_raw_res_immune_ciber/full_de_cohort_cov.parquet\
            -w gsea_rank --maps ../ref_uniprot_072522/mapped_sites_to_pdbs.tsv --pdbstore ../../../../getzlab-CLUMPS2/clumps/db/ref/pdbs \
            --grouping $GROUP --features acetylome phosphoproteome --threads 96 -v --subset negative --output_dir ${GROUP}_ciber_negative
done
