#!/usr/bin/env Rscript
#
# Master script executing the entire analysis for the paper in sequence.
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

