#!bin/bash

# https://github.com/getzlab/SignatureAnalyzer/tree/master
# pip3 install signatureanalyzer
# Using SignatureAnalyzer v0.0.8
BASEDIR=../data/processed/061721/imputed_res/

# Full wtih RNA
RUNNAME=061721_imputed_res_reg
signatureanalyzer -o ../analysis/signatures/${RUNNAME} \
                  -t matrix \
                  -n 100 \
                  --K0 50 \
                  --max_iter 20000 \
                  --object gaussian \
                  --prior_on_H L2 \
                  --prior_on_W L2 \
                  --random_seed 42 \
                  ${BASEDIR}/pan_reg_X.parquet
