###############################################
##Dmitry Sutormin, 2023##
### Collect core contigs and bridge reads for re-assembly.
############

import argparse
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


parser=argparse.ArgumentParser(description="Get sequences of core contigs and sequences of bridging reads.")
parser.add_argument("-dn", "--dataset_name", required=True, help="Dataset name.")
parser.add_argument("-fcl", "--fin_contigs_list", required=True, help="Path to file with names of contigs comprising the final expanded core set.")
parser.add_argument("-fbf", "--fin_bridge_file", required=True, help="Path to file containing information about contig-to-contig bridges.")
parser.add_argument("-mqt", "--mapq_thr", required=True, help="Mapq threshold to filter alignments.")
parser.add_argument("-bst", "--bridge_str_thr", required=True, help="Contig-to-contig bridge strength threshold to filter contig links.")
parser.add_argument("-af", "--assembly_file", required=True, help="Assembly file (fasta) to take contig sequences from.")
parser.add_argument("-rf", "--long_read_file", required=True, help="Long read file (fastq) to take read sequences from.")
parser.add_argument("-rsm", "--read_search_mode", required=True, default='seqtk', choices=['seqtk', 'manual'], help="Algorithm to use for retrival of reads sequences.")
parser.add_argument("-sp", "--seqtk_path", required=True, default='seqtk', help="Path to seqtk executable.")
parser.add_argument("-o", "--output", required=True, help="Path to output directory.")
args, unknown=parser.parse_known_args()

Dataset_name=args.dataset_name
Core_contigs_path=args.fin_contigs_list
Bridges_path=args.fin_bridge_file
Bridge_strength_threshold=int(args.bridge_str_thr)
Mapq_threshold=int(args.mapq_thr)
Assembly_path=args.assembly_file
Long_reads_path=args.long_read_file
Mode=args.read_search_mode
Seqtk_path=args.seqtk_path
Output_path=args.output


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
# Reads bridge file, get bridging read IDs.
############

def read_bridge_file_get_reads(Dataset_name, Bridges_path, mapq_threshold, bridge_strength_threshold, output_path):
    
    Read_IDs_list=[]
    
    filein=open(Bridges_path, 'r')
    
    for line in filein:
        if line[0]!='#':
            line=line.rstrip().split('\t') 
            bridge_strength=int(line[2])
            if bridge_strength>=bridge_strength_threshold:
                reads_id_ar=line[4].split('!|!')
                Read_IDs_list+=reads_id_ar
                
    filein.close()
    
    Read_IDs_list_u=list(set(Read_IDs_list))
    
    # Write read IDs.
    bridging_reads_ids_path=os.path.join(output_path, f'Bridges_reads_IDs_{Dataset_name}_mapq_{mapq_threshold}_br_strength_{bridge_strength_threshold}.list')
    with open(bridging_reads_ids_path, "w") as good_out:
        for Read_ID in Read_IDs_list_u:
            good_out.write("{}\n".format(Read_ID))   
    good_out.close()      
    
    print(f'Number of unique bridging reads retrieved: {len(Read_IDs_list_u)}')
    
    return Read_IDs_list_u


############
# Write core contigs sequences.
############

def get_core_contig_seqs(Dataset_name, Assembly_path, Core_contigs_list, 
                         bridge_strength_threshold, mapq_threshold, output_path):
    
    core_cont_out_path=os.path.join(output_path, f'Expanded_core_contigs_{Dataset_name}_mapq_{mapq_threshold}_br_strength_{bridge_strength_threshold}.fasta')
    
    fileout=open(core_cont_out_path, 'w')
    
    for record in SeqIO.parse(Assembly_path, "fasta"):
        record_id=record.id
        if record_id in Core_contigs_list:
            record_seq=str(record.seq)
            fileout.write(f'>{record_id}\n{record_seq}\n')
            
    fileout.close()
    
    return


############
# Write bridging reads sequences.
############

def get_br_reads_seqs(Dataset_name, Long_reads_path, Read_IDs_list, 
                      bridge_strength_threshold, mapq_threshold, output_path):
    
    br_reads_out_path=os.path.join(output_path, f'Bridging_reads_{Dataset_name}_mapq_{mapq_threshold}_br_strength_{bridge_strength_threshold}.fastq')
    
    Bridging_reads_ar=[]
    
    read_num=0
    
    for record in SeqIO.parse(Long_reads_path, "fastq"):
        record_id=record.id
        if record_id in Read_IDs_list:
            Bridging_reads_ar.append(record)
            
        read_num+=1
        
        if read_num%10000==0:
            print(f'{read_num} reads processed..')

    SeqIO.write(Bridging_reads_ar, br_reads_out_path, 'fastq')
    
    return


############
# Wrapper function.
############

def Collect_contigs_reads_wrapper(Dataset_name, Core_contigs_path, Bridges_path, bridge_strength_threshold, 
                                  mapq_threshold, Assembly_path, reads_path, mode, seqtk_path, output_path):
    
    # Read core contigs list.
    Core_contigs_list=read_core_contigs_list(Core_contigs_path)
    
    # Read bridges file, get read IDs.
    Read_IDs_list=read_bridge_file_get_reads(Dataset_name, Bridges_path, mapq_threshold, bridge_strength_threshold, output_path)
    
    # Get core contigs sequences.
    get_core_contig_seqs(Dataset_name, Assembly_path, Core_contigs_list, bridge_strength_threshold, mapq_threshold, output_path)
    
    # Get bridging reads sequences.
    if mode=="manual":
        get_br_reads_seqs(Dataset_name, reads_path, Read_IDs_list, bridge_strength_threshold, mapq_threshold, output_path)
        
    elif mode=="seqtk":
        reads_ids_file_path=os.path.join(output_path, f'Bridges_reads_IDs_{Dataset_name}_mapq_{mapq_threshold}_br_strength_{bridge_strength_threshold}.list')
        br_reads_out_path=os.path.join(output_path, f'Bridging_reads_{Dataset_name}_mapq_{mapq_threshold}_br_strength_{bridge_strength_threshold}.fastq')
        seqtk_command=f'{seqtk_path} subseq {reads_path} {reads_ids_file_path} > {br_reads_out_path}'
        
        subprocess.run(seqtk_command, shell=True)        
    
    return



starttime=time.time()


Collect_contigs_reads_wrapper(Dataset_name, Core_contigs_path, Bridges_path, Bridge_strength_threshold, 
                              Mapq_threshold, Assembly_path, Long_reads_path, Mode, Seqtk_path, Output_path)

print(f'Total working time, min: {(time.time()-starttime)/60}')
