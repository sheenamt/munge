Changes for Munge
================
v3.4
====
 * specify polyhunter outfile with -o
 * specify CADD file from pipeline repo rather than parsing from annovar files
 * test for all variant types in varscan parser
 * update tests for monoseq/polyt
 * change polyt top-level to use samplesheet order for output order
 * change headers in polyt output
 * remove exception from xlsxmaker
 * update tests to ensure sorting of loadlist by well works
 * update code and tests to use manifest for sort order, get rid of exception in xlsxmaker
 * update for parsing correct number of controls in samplesheet creation 

v3.3
====
 * test for all variant types in varscan parser
 * parse multiple calls at same position in genotyper, use qX depth

v2.0.1
====
 * munge quality_metrics updated:
  * One quality parsing script works for hs metrics and/or quality metrics file
  * Print git version of pipeline in quality_analysis file
 * munge annovar_summary updated:
  * Add splice site effect prediction by AdaBoost and Random Forest from dbscSNV version 1.1
  * UW frequency now prints only for applicable assay-machine combo of sample 
  * Variant_type is now a list of variant types
  * Use “--separate” flag to print all transcripts in annovar
  * Read count printed from gatk, varscan and varscanSNP, when available
 * munge amplicon_metrics added:
  * Print amplicon coverage metrics
