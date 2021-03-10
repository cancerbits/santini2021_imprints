#!/usr/bin/env Rscript
#
# Master script executing the majority of the analysis for the paper in sequence.
# Preceding steps that are not covered in this script:
# 1. Primary processing of WGBS reads: For this we used the pypiper/looper framework (http://looper.databio.org/en/latest/, http://pypiper.databio.org/en/latest/) with the configuration file in metadata/config.yaml
# 2. Primary processing of asRNA-seq reads: For this we used the Allelome.Pro pipeline (https://academic.oup.com/nar/article/43/21/e146/2468102) and the respective code is in src/rnaseq/*
#
# run("imprinting")

run("imprinting", "init")
run("imprinting", "gold_standard")
run("imprinting", "load") 
run("imprinting", "define_regions")
run("imprinting", "region_enrichments")
run("imprinting", "motifs")
run("imprinting", "annotate")
run("imprinting", "imprint_control_regs")
run("imprinting", "expression")
run("imprinting", "export")
run("imprinting", "external")

