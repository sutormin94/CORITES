###############################################
##Dmitry Sutormin, 2023##
### Collect reads uniquely aligned to a set of contigs for re-assembly.
############

import os
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
import math
import pandas as pd
import copy
import numpy as np
import time
import subprocess


############
# Read core contigs list.
############

def read_core_contigs_list(Core_contigs_path):
    
    Core_contigs_list=[]
    
    filein=open(Core_contigs_path, 'r')
    
    for line in filein:
        line=line.rstrip()
        Core_contigs_list.append(line)
        
    filein.close()
    
    print(f'Number of core contigs retrieved: {len(Core_contigs_list)}')
    
    return Core_contigs_list


############
# Read sam file. Find uniquely mapped reads.
############

def read_sam_get_unique_reads(sam_file_path, mapq_threshold, multimap_thr, Core_contigs_list):
    
    sam_filein=open(sam_file_path, 'r')
    
    Long_reads_dict={}
    
    count_aligs=0
    
    for line in sam_filein:
        if line[0]!="@":
            line=line.rstrip().split('\t')
            read_id=line[0]
            sam_flag=int(line[1])
            contig_name=line[2]
            mapq=int(line[3])
            if (sam_flag!=4) and (mapq>mapq_threshold):
                if read_id not in Long_reads_dict:
                    Long_reads_dict[read_id]=[contig_name]
                else:
                    if contig_name not in Long_reads_dict[read_id]:
                        Long_reads_dict[read_id].append(contig_name)
        count_aligs+=1
        
        if count_aligs%1000000==0:
            print(f'Processed {count_aligs} alignments')
            
    sam_filein.close()
    
    Read_IDs_list=[]
    
    for read_id, contig_ar in Long_reads_dict.items():
        if len(contig_ar)<=multimap_thr:
            wrong_target=0
            for contig_name in contig_ar:
                if contig_name not in Core_contigs_list:
                    wrong_target=1
            if wrong_target==0:
                Read_IDs_list.append(read_id)                 
    
    print(f'{len(Read_IDs_list)} uniquely mapped reads found.')
        
    return Read_IDs_list


############
# Write list of uniquely mapped reads IDs.
############

def write_read_ids(dataset_name, read_type, Read_IDs_list, bridge_strength_threshold, mapq_threshold, output_path):
    
    reads_ids_out_path=f'{output_path}Uniquely_mapped_reads_{dataset_name}_mapq_{mapq_threshold}_br_strength_{bridge_strength_threshold}.list'
    
    fileout=open(reads_ids_out_path, 'w')
    
    for read_id in Read_IDs_list:
        if read_type=='nanopore':
            fileout.write(f'{read_id}\n')
        elif read_type=='illumina_paired':
            fileout.write(f'{read_id}/1\n{read_id}/2\n')
    
    fileout.close()
    
    return


############
# Write uniquely mapped reads sequences.
############

def get_un_map_reads_seqs(dataset_name, reads_path, Read_IDs_list, 
                          bridge_strength_threshold, mapq_threshold, mode, seqtk_path, output_path):
    
    um_reads_out_path=f'{output_path}Uniquely_mapped_reads_{dataset_name}_mapq_{mapq_threshold}_br_strength_{bridge_strength_threshold}_{mode}.fastq'
    
    if mode=="biopython":
        
        UM_reads_ar=[]
        
        read_num=0
        
        for record in SeqIO.parse(reads_path, "fastq"):
            record_id=record.id
            if record_id in Read_IDs_list:
                UM_reads_ar.append(record)
                
            read_num+=1
            
            if read_num%10000==0:
                print(f'{read_num} reads processed..')
    
        SeqIO.write(UM_reads_ar, um_reads_out_path, 'fastq')
        
    elif mode=="manual":
        
        filein_reads=open(reads_path, 'r')
        fileout_um_reads=open(um_reads_out_path, 'w')
        
        read_num=0
        
        for line in filein_reads:
            if line[0]=='@':
                
                read_num+=1
                
                if read_num%10000==0:
                    print(f'{read_num} reads processed..')                
                
                read_id=line[1:].split(' ')[0]
                
                if read_id in Read_IDs_list:
                    fileout_um_reads.write(f'{line.rstrip()}\n')
                    write_the_read=1
                else:
                    write_the_read=0
            else:
                if write_the_read==1:
                    fileout_um_reads.write(f'{line.rstrip()}\n')
                
                elif write_the_read==0:
                    continue                      
        
        filein_reads.close()
        fileout_um_reads.close() 
        
    elif mode=="seqtk":
        
        reads_ids_file_path=f'{output_path}Uniquely_mapped_reads_{dataset_name}_mapq_{mapq_threshold}_br_strength_{bridge_strength_threshold}.list'
        seqtk_command=f'{seqtk_path} subseq {reads_path} {reads_ids_file_path} > {um_reads_out_path}'
        
        subprocess.run(seqtk_command, shell=True)
    
    return


############
# Wrapper function.
############

def Collect_uniq_reads_wrapper(dataset_name, mode, read_type, seqtk_path, core_contigs_path, bridge_strength_threshold, 
                               mapq_threshold, multimap_thr, reads_path, sam_file_path, output_path):
    
    # Read core contigs list.
    Core_contigs_list=read_core_contigs_list(core_contigs_path)  
    
    # Read sam file identify uniquely mapped reads.
    Read_IDs_list=read_sam_get_unique_reads(sam_file_path, mapq_threshold, multimap_thr, Core_contigs_list)
    
    # Write list of uniquely mapped reads IDs.
    write_read_ids(dataset_name, read_type, Read_IDs_list, bridge_strength_threshold, mapq_threshold, output_path)
    
    # Get uniquely mapped reads sequences.
    get_un_map_reads_seqs(dataset_name, reads_path, Read_IDs_list, 
                          bridge_strength_threshold, mapq_threshold, mode, seqtk_path, output_path)
    
    return


starttime=time.time()

# Dataset name.
Dataset_name="Refined_contigs_core_bin_24_3_and_003_illumina"

# Mode of fastq file reading and filtering: biopython, manual or seqtk.
Mode="seqtk"

# Path to seqtk executable. Required only in the seqtk mode.
Seqtk_path="/home/lam16/Programs/seqtk/seqtk"

# Read types: 'nanopore', 'illumina_paired'
Read_type="illumina_paired"

# Mapq threshold.
Mapq_threshold=30

# Bridge strength threshold.
Bridge_strength_threshold=100

# Multimapping threshold - an allowed number of sites for a read to be aligned to.
Multimap_thr=4

# Path to core contigs file.
Core_contigs_path=f'/home/lam16/Sponges_reads/I_palmata/refine_for_core_assebmly/Contig_to_contig_bridges/I_palmata_2016_bact_1/Refined_contigs_core_bin_24_3_and_003_mapq_{Mapq_threshold}_br_strength_{Bridge_strength_threshold}__iteration_5.txt'

# Path to fastq.
Reads_path="/home/lam16/Sponges_reads/I_palmata/refine_for_core_assebmly/I_palmata_2016_bact1_R12_illumina.fastq"

# Path to sam file.
Sam_file_path="/home/lam16/Sponges_reads/I_palmata/refine_for_core_assebmly/I_palmata_2016_bact1_illumina.sam"

# Output path.
Output_path="/home/lam16/Sponges_reads/I_palmata/refine_for_core_assebmly/Contig_to_contig_bridges/I_palmata_2016_bact_1/"

Collect_uniq_reads_wrapper(Dataset_name, Mode, Read_type, Seqtk_path, Core_contigs_path, Bridge_strength_threshold, 
                           Mapq_threshold, Multimap_thr, Reads_path, Sam_file_path, Output_path)

print(f'Total working time, min: {(time.time()-starttime)/60}')
