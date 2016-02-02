Changes for Munge
================

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
