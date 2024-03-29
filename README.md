
**code repo for the manuscript:** 


_Gene transcriptional profiles in gonads of Bacillus taxa (Phasmida) with different cytological mechanisms of automictic parthenogenesis._


![alt text](https://upload.wikimedia.org/wikipedia/commons/0/0b/Bacillus_rossius_Livorno.jpg)


This repository contains the code used to study the molecular groundplan of automixis in Bacillus stick insects. 


Experiment reads have been deposited under the NCBI BioProject [**PRJNA578804**](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA578804).


The manuscript is available [here](https://zoologicalletters.biomedcentral.com/articles/10.1186/s40851-022-00197-z).


The scripts names are largely self-explanatory and the repo is organized in four folders:


- [scripts_snakefiles](https://github.com/for-giobbe/gene-transcriptional-profiles-in-automictic-taxa/tree/main/scripts_snakefiles)

	- snakefile_orthology_inference 	
	- snakefile_transcriptome_annotation	
	- snakefile_transcriptome_filter_aa
	- snakefile_transcriptome_goterms


- [scripts_shell](https://github.com/for-giobbe/gene-transcriptional-profiles-in-automictic-taxa/tree/main/scripts_shell)
	
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


- [intermediate_files](https://github.com/for-giobbe/gene-transcriptional-profiles-in-automictic-taxa/tree/main/intermediate_files.zip)

	- folder containing intermediate files necessary for  paper figures, generated using the scripts in [scripts_R](https://github.com/for-giobbe/gene-transcriptional-profiles-in-automictic-taxa/tree/main/scripts_R).


*NB*: several abbreviations are used, including:

- BAT for _Bacillus atticus_
- BRO for _Bacillus rossius_
- BGM for _Bacillus grandii_
- BP for Biological Processes
- MF for Molecular Functions
