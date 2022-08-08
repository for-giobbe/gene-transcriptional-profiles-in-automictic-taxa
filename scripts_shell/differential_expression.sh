####################################################################################################################################### map reads

for i in $(ll ../raw_reads/*.gz | awk -F "/" '{print $3}' | awk -F "_" '{print $1"_"$2"_"$3"_"$4}' | sort -u);

do species=$(echo $i | awk -F "_" '{print $1}')

../../Miniconda3_x86_64/opt/trinity-2.4.0/util/align_and_estimate_abundance.pl --seqType fq --left "../raw_reads/"$i"_1.fastq.gz" --right "../raw_reads/"$i"_2.fastq.gz" \
--transcripts $species".def.nuc.fa" --est_method RSEM  --aln_method bowtie2 --trinity_mode --output_dir $i --prep_reference;

done

####################################################################################################################################### DE BAT

../../Miniconda3_x86_64/opt/trinity-2.4.0/util/abundance_estimates_to_matrix.pl  --est_method RSEM  --out_prefix BAT_RSEM --name_sample_by_basedir BAT_F_G_RNA03/RSEM.genes.results BAT_F_L_RNA03/RSEM.genes.results \
BAT_F_G_RNA05/RSEM.genes.results BAT_F_L_RNA05/RSEM.genes.results BAT_F_G_RNA06/RSEM.genes.results BAT_F_L_RNA06/RSEM.genes.results BAT_F_G_RNA07/RSEM.genes.results BAT_F_L_RNA07/RSEM.genes.results \
BAT_F_G_RNA11/RSEM.genes.results BAT_F_L_RNA11/RSEM.genes.results BAT_F_G_RNA12/RSEM.genes.results BAT_F_L_RNA12/RSEM.genes.results

../../Miniconda3_x86_64/opt/trinity-2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix BAT_RSEM.counts.matrix --method DESeq2 --samples_file BAT_samples_file --output DE_DESeq2_BAT

####################################################################################################################################### DE BRO

../../Miniconda3_x86_64/opt/trinity-2.4.0/util/abundance_estimates_to_matrix.pl  --est_method RSEM  --out_prefix BRO_RSEM --name_sample_by_basedir BRO_F_G_RNA01/RSEM.genes.results BRO_F_L_RNA01/RSEM.genes.results \
BRO_F_G_RNA03/RSEM.genes.results BRO_F_L_RNA03/RSEM.genes.results BRO_F_G_RNA06/RSEM.genes.results BRO_F_L_RNA06/RSEM.genes.results BRO_F_G_RNA07/RSEM.genes.results BRO_F_L_RNA07/RSEM.genes.results \
BRO_F_G_RNA08/RSEM.genes.results BRO_F_L_RNA08/RSEM.genes.results BRO_F_G_RNA10/RSEM.genes.results BRO_F_L_RNA10/RSEM.genes.results

../../Miniconda3_x86_64/opt/trinity-2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix BRO_RSEM.counts.matrix --method DESeq2 --samples_file BRO_samples_file --output DE_DESeq2_BRO


####################################################################################################################################### DE BGM f


../../Miniconda3_x86_64/opt/trinity-2.4.0/util/abundance_estimates_to_matrix.pl  --est_method RSEM  --out_prefix BGM_RSEM_f_only --name_sample_by_basedir BGM_F_G_RNA10/RSEM.genes.results BGM_F_L_RNA10/RSEM.genes.results \
BGM_F_G_RNA12/RSEM.genes.results BGM_F_L_RNA12/RSEM.genes.results BGM_F_G_RNA14/RSEM.genes.results BGM_F_L_RNA14/RSEM.genes.results BGM_F_G_RNA15/RSEM.genes.results BGM_F_L_RNA15/RSEM.genes.results \
BGM_F_G_RNA17/RSEM.genes.results BGM_F_L_RNA17/RSEM.genes.results BGM_F_G_RNA19/RSEM.genes.results BGM_F_L_RNA19/RSEM.genes.results

../../Miniconda3_x86_64/opt/trinity-2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix BGM_RSEM_f_only.counts.matrix --method DESeq2 --samples_file BGM_f_only_samples_file --output DE_DESeq2_BGM_f_only


####################################################################################################################################### DE BGM m


../../Miniconda3_x86_64/opt/trinity-2.4.0/util/abundance_estimates_to_matrix.pl  --est_method RSEM  --out_prefix BGM_RSEM_m_only --name_sample_by_basedir BGM_M_G_RNA04/RSEM.genes.results BGM_M_L_RNA04/RSEM.genes.results \
BGM_M_G_RNA06/RSEM.genes.results BGM_M_L_RNA06/RSEM.genes.results BGM_M_G_RNA08/RSEM.genes.results BGM_M_L_RNA08/RSEM.genes.results BGM_M_G_RNA20/RSEM.genes.results BGM_M_L_RNA20/RSEM.genes.results \
BGM_M_G_RNA21/RSEM.genes.results BGM_M_L_RNA21/RSEM.genes.results BGM_M_G_RNA22/RSEM.genes.results BGM_M_L_RNA22/RSEM.genes.results

../../Miniconda3_x86_64/opt/trinity-2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix BGM_RSEM_m_only.counts.matrix --method DESeq2 --samples_file BGM_m_only_samples_file --output DE_DESeq2_BGM_m_only



