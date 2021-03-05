# Data and code for Miller-ter Kuile et al. 2021

This repository includes the data and code to reproduce the analyses in a manuscript by Miller-ter Kuile et al. examining the roles of predator body size and hunting strategy on prey size selection based on DNA metabarcoding data. 

# Explanation of folders and contents

## Drafts

Drafts of manuscript, supplementary materials, and notes from meetings with collaborators on this project

## art

A folder of watercolor artwork produced for this project by Miller-ter Kuile. 

## code

The code scripts for reproducing the bioinformatics, data tidying, and data analysis steps in this paper. 

### 1_denoising

A subfolder including three sub-folders of the steps taken to perform the bioinformatics steps for this project

#### 1a_trim_sequences

*'cutadapt_primer_trimming.txt'*: a text file of command line script to cut primer adapters from all DNA sequencing data

#### 1b_dada2

*'1b1_dada2_separate_cross_run_error.R'*: R script for performing the error generating step in the dada2 algorithm to compare cross-sequencing run error rates. These error rates should be similar in order to combine multiple sequencing runs to perform the full dada2 denoising algorithm

*'1b2_cross_run_error.R'*: R script for producing the figures comparing cross-sequencing run errors

*'1b3_dada2_combined_runs.R'*: R script for performing the dada2 algorithm to produce the final denoised ASV list and ASV-by-sample tables used in subsequent analyses.

#### 1c_unoise

*'unoise_code.txt'*: command line script used to perform the unoise3 algorithm on the sequencing data. The input had already been run through the primer-cutting step.

#### 1d_blast_code

*'knot.blast.bash.script.txt'*: bash script submitted to UCSB computing cluster for assigning DNA sequences to taxonomic identifiers based on the nucleotide GenBank database.

### 2_qc

*2a_prelim_sequencing_depth_check.R*: R script for rarefaction curve for our initial MiSeq NANO run with 4 samples to determine how many samples to run per sequencing run to achieve sufficient sequencing depth. 

*2b_seq_depth_dadavsunoise.R*: R script comparing sequencing depth per sample for the dada2 and unoise3 algorithms

*2c_cross-run_sample_comparisons.R*: R script of analyses of the 19 samples re-run on all four sequencing runs to report run-to-run variability in sequencing depth and resulting prey ASV assignments.

*2d_multiple_individuals.R*: R script of analyses comparing samples that consisted of only one versus multiple individuals to verify that samples with multiple individuals were not over-represented in interaction data.

*2e_control_samples.R*: R script examining ASV assignment to positive and negative controls as an estimate of ASV assignment specificity and sequencing error due to sequence jumping. 

### 3_data_prep

*3a_low_sampling_depth_cutoff.R*: R script used to determine the quantile cutoff for removing samples with insufficient sequencing depth to be compared to the rest of the samples.

*3b_master_taxonomy_list.R*: R script combining GenBank and BOLD taxonomic assignments and merging these into a master list. 

*3c_rarefy_samples.R*: R script used to rarefy the samples for subsequent analyses.

*3d_subset_diet_taxonomies.R*: R script for removing all ASVs either not assigned a taxonomic identification or which did not match to prey taxonomies (e.g. predator DNA, environmental DNA such as fungi, etc.).

*3e_master_body_size_list.R*: R script combining all sources of body size data for predators and prey

*3f_predator_mass_length_models.R*: R script used to create mixed effects models of predator mass-length to give mass values to predators for which these data were not collected.

*3g_interactions_w_body_size.R*: R script that combines the predator-prey interactions to the body size dataset and produces the final dataset used in statstical analyses. 

*3h_palmyra_missing_bs_data_feb21.R*: R script used to help determine future lab work priority for body size measurements.

### 4_analyses

*4a_pred_prey_size.R*: R script used to run the predator prey mass linear mixed effects models and model selection process. Includes code for creating figures. 

*4b_feeding_mode_ratios.R*: R Script of the feeding mode predator-prey size ratio statistical models and model selection process. Includes code for creating figures.

### 5_supp_figures

*data_summary_tables.R*: R script to produce summary tables of sample sizes by predator species

*prey_phylogeny.R*: R script to produce the supplementary prey phylogeny figure.










