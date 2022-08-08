
**code repo for the manuscript:** 


_Gene transcriptional profiles in gonads of Bacillus taxa (Phasmida) with different cytological mechanisms of automictic parthenogenesis._


---


![alt text](https://upload.wikimedia.org/wikipedia/commons/0/0b/Bacillus_rossius_Livorno.jpg)

This repository contains the code used to study the molecular groundplan of automixis.

Experiment reads have been deposited under the NCBI BioProject [**PRJNA578804**](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA578804).

The manuscript will be made available [here]() when published.


---


The repo is organized in four folders:

- scripts_snakefiles

	- snakefile_orthology_inference
	- snakefile_transcriptome_annotation
	- snakefile_transcriptome_filter_aa
	- snakefile_transcriptome_goterms

- scripts_shell
	
	- assign_GOterms_to_orthogroups.sh
	- differential_expression.sh
	- gblocks.sh
	- msa_align.sh
	- omega_pos_sel_formatting.sh
	- orthogroups-differential-expression.sh
	- orthologs_assignment_exp_matrixes.sh
	- pep_2_cds_orthologs.sh
	- phylostratigraphy_BAT.sh
	- phylostratigraphy_BGM.sh
	- phylostratigraphy_BRO.sh
	- reformat_orthogroups.sh
	- taxonomy_assignement.sh

- [scripts_R](https://github.com/for-giobbe/gene-transcriptional-profiles-in-automictic-taxa/tree/main/scripts_R)

	- gene_enrichment.R
	- normalized_counts_PCA.R
	- phylostratigraphy_and_sankyeys.R
	- phylostratigraphy_z-test.R

- intermediate_files

	- a zipped folder containing all the intermediate files necessary to generate the figures using the scripts in the scripts_R folder.
