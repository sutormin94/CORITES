#!bin/bash

PWD=/home/niagara/Storage/D_Sutormin/Metagenomes/N_Rusanova/Spongy/2022_nanopore/I_palmata/I_palmata_2016_2022_core_assebmlies/
Reads_path=${PWD}Bridging_reads_Shared_contigs_four_binners_bact_1_mapq_30_br_strength_500.fastq

for i in `ls -a ${PWD}trycycler_cluster/`
do
 echo 'Working with' ${i}
 trycycler reconcile --reads ${Reads_path} --cluster_dir ${PWD}trycycler_cluster/${i}
done