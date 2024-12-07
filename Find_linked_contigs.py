###############################################
##Dmitry Sutormin, 2023##
### Find contigs linked to a core set of contigs using the information from long reads.
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

parser=argparse.ArgumentParser(description="Find contigs linked with long reads. Iteratively expand core contigs set.")
parser.add_argument("-dn", "--dataset_name", required=True, help="Dataset name.")
parser.add_argument("-sam", "--sam_file", required=True, help="Path to sam file with long reads vs contigs alignment.")
parser.add_argument("-mqt", "--mapq_thr", required=True, help="Mapq threshold to filter alignments.")
parser.add_argument("-bst", "--bridge_str_thr", required=True, help="Contig-to-contig bridge strength threshold to filter contig links.")
parser.add_argument("-mi", "--max_iterations", required=True, help="Maximal number of iterations to expand core contigs set.")
parser.add_argument("-cl", "--contigs_list", required=True, help="Path to file with names of contigs comprising the core set.")
parser.add_argument("-samtools", "--samtools_path", required=True, default='samtools', help="Path to samtools executable.")
parser.add_argument("-o", "--output", required=True, help="Path to output directory.")
args, unknown=parser.parse_known_args()

Dataset_name=args.dataset_name
Sam_file_path=args.sam_file
Mapq_threshold=int(args.mapq_thr)
Bridge_strength_thr=int(args.bridge_str_thr)
Max_num_iterations=int(args.max_iterations)
Path_core_cont_list=args.contigs_list
Samtools_path=args.samtools_path
Output_path=args.output


starttime=time.time()
Logfile_path=os.path.join(Output_path, 'Find_linked_contigs.log')
logfile=open(Logfile_path, 'w', buffering=1)

c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Expansion of core contigs started.', file=sys.stdout)
logfile.write(f'[{current_time}] Expansion of core contigs started.\n')


############
# Read core contigs list.
############

def read_core_contigs(path_core_cont_list):
    
    core_contigs_list=[]
    
    filein=open(path_core_cont_list, 'r')
    
    for line in filein:
        line=line.rstrip()
        core_contigs_list.append(line)
    
    filein.close()
    
    print(f'{len(core_contigs_list)} core contigs to be used as anchors.')
    
    return core_contigs_list


############
# Read sam file. Build contig to contig graph.
############

def read_sam_file(sam_file_path, mapq_threshold):
    
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
    
    contig_to_contig_dict={}
    
    for read_id, contig_list in Long_reads_dict.items():
        if len(contig_list)>0:
            for i in range(len(contig_list)):
                for j in range(len(contig_list)):
                    if j>i:
                        contig_pair=f'{contig_list[i]}___{contig_list[j]}'
                        contig_pair_rev=f'{contig_list[j]}___{contig_list[i]}'
                        if contig_pair in contig_to_contig_dict:
                            contig_to_contig_dict[contig_pair][0]+=1
                            contig_to_contig_dict[contig_pair][1].append(read_id)
                        elif contig_pair_rev in contig_to_contig_dict:
                            contig_to_contig_dict[contig_pair_rev][0]+=1
                            contig_to_contig_dict[contig_pair_rev][1].append(read_id)  
                        else:
                            contig_to_contig_dict[contig_pair]=[1, [read_id]]                    
    
    print(f'{len(contig_to_contig_dict)} contig to contig bridges found.')
        
    return contig_to_contig_dict


############
# Calculate distribution of core contigs coverage depth.
############

def core_contigs_cov_depth(dataset_name, coverage_file_path, core_contigs_list, output_path):
    
    cov_filein=open(coverage_file_path, 'r')
    Cov_dict={}
    for line in cov_filein:
        if line[0]!="#":
            line=line.rstrip().split('\t')
            Contig_name=line[0]
            Mean_cov_depth=float(line[6])
            Cov_dict[Contig_name]=Mean_cov_depth
    
    cov_filein.close()
    
    core_contigs_cov_depth_list=[]
    for contig_name in core_contigs_list:
        if contig_name in Cov_dict:
            cov_depth=Cov_dict[contig_name]
            core_contigs_cov_depth_list.append(cov_depth)
            
    mean_cov_depth=np.mean(core_contigs_cov_depth_list)
    c_time=time.localtime()
    current_time=time.strftime("%H:%M:%S", c_time) 
    print(f'[{current_time}] Mean coverage depth of a core set contigs: {mean_cov_depth}.')
    logfile.write(f'[{current_time}] Mean coverage depth of a core set contigs: {mean_cov_depth}.\n')    
            
    fig=plt.figure(figsize=(3, 3), dpi=100)
    plt.hist(core_contigs_cov_depth_list, bins=20, facecolor='blue', edgecolor='black', alpha=0.5)
    plt.axvline(x=mean_cov_depth, color='r', linestyle='dashed', linewidth=2)
    plt.xlabel('Coverage depth')
    plt.ylabel('Number of contigs')
    plt.tight_layout()
    hist_filename=f'{dataset_name}_core_contigs_coverage_depth_distribution.png'
    plot_outpath=os.path.join(output_path, hist_filename)
    plt.savefig(plot_outpath, dpi=300, figsize=(4, 4))
    plt.close()
    
    return


############
# Calculate core contigs set statistics.
############

def core_contigs_stat(core_contigs_list, alg_iteration):
    
    Number_of_core_contigs=len(core_contigs_list)
    print(f'###Iteration {alg_iteration} ###\nTotal number of core contigs for algorithm iteration {alg_iteration}: {Number_of_core_contigs}')
    
    Total_contig_len=[]
    for contig_name in core_contigs_list:
        contig_len=int(contig_name.split('_')[3])
        Total_contig_len.append(contig_len)
        
    Tot_length_of_contigs=np.sum(Total_contig_len)
    print(f'Total length of core contigs for algorithm iteration {alg_iteration}: {Tot_length_of_contigs} bp')
    Med_length_of_contigs=np.median(Total_contig_len)
    print(f'Median length of core contigs for algorithm iteration {alg_iteration}: {Med_length_of_contigs} bp')
    
    return Number_of_core_contigs, Tot_length_of_contigs, Med_length_of_contigs


############
# Pull contigs by anchor contigs.
############

def pull_contigs_by_anchors(dataset_name, contig_to_contig_dict, core_contigs_list, bridge_strength_thr, mapq_threshold, alg_iteration, output_path):
    
    fileout=open(os.path.join(output_path, f'Bridges_info_mapq_{mapq_threshold}_br_strength_{bridge_strength_thr}__iteration_{alg_iteration}.tab'), 'w')
    
    # Core contig flag=1 if core contig name is in the first column.
    # Core contig flag=12 if core contigs names are both in the first and second columns - indication of an inter-core-contig bridge.
    fileout.write('#Contig 1\tContig 2\tBridge strength\tCore contig flag\tSupporting reads\n') 
    
    bridge_strength_list=[]
    
    for contig_pair, bridge_info in contig_to_contig_dict.items():
        contig_pair_ar=contig_pair.split('___')
        contig_1=contig_pair_ar[0]
        contig_2=contig_pair_ar[1]
        
        if (contig_1 in core_contigs_list) and (contig_2 in core_contigs_list):
            #print(f'Inter-core-contig bridge detected between {contig_1} and {contig_2}: {bridge_info}')
            fileout.write(f'{contig_1}\t{contig_2}\t{bridge_info[0]}\t12\t')
            for read_name in bridge_info[1]:
                fileout.write(f'{read_name}!|!')
            fileout.write('\n')
            
            bridge_strength_list.append(float(bridge_info[0]))
            
        elif (contig_1 in core_contigs_list) and (contig_2 not in core_contigs_list):
            #print(f'Bridge found between core contig {contig_1} and {contig_2}: {bridge_info}')
            fileout.write(f'{contig_1}\t{contig_2}\t{bridge_info[0]}\t1\t')
            for read_name in bridge_info[1]:
                fileout.write(f'{read_name}!|!')
            fileout.write('\n')
            
            bridge_strength_list.append(float(bridge_info[0]))
            
            if bridge_info[0]>=bridge_strength_thr:
                core_contigs_list.append(contig_2)
            
        elif (contig_1 not in core_contigs_list) and (contig_2 in core_contigs_list):
            #print(f'Bridge found between core contig {contig_2} and {contig_1}: {bridge_info}')
            fileout.write(f'{contig_2}\t{contig_1}\t{bridge_info[0]}\t1\t')
            for read_name in bridge_info[1]:
                fileout.write(f'{read_name}!|!')
            fileout.write('\n')
            
            bridge_strength_list.append(float(bridge_info[0]))
            
            if bridge_info[0]>=bridge_strength_thr:
                core_contigs_list.append(contig_1)            
            
        else:
            continue
        
    fileout.close()
    
    if alg_iteration==1:
        
        
        bridge_strength_list_sorted=sorted(bridge_strength_list)
        bs_percentile_index=int(len(bridge_strength_list)*0.9)
        bs_percentile=bridge_strength_list_sorted[bs_percentile_index]
        
        c_time=time.localtime()
        current_time=time.strftime("%H:%M:%S", c_time) 
        print(f'[{current_time}] 90% percentile for core contigs bridge strength: {bs_percentile}.')
        logfile.write(f'[{current_time}] 90% percentile for core contigs bridge strength: {bs_percentile}.\n')           
        
        fig=plt.figure(figsize=(3, 3), dpi=100)
        plt.hist(bridge_strength_list, bins=20, facecolor='red', edgecolor='black', alpha=0.5)
        plt.axvline(x=bs_percentile, color='blue', linestyle='dashed', linewidth=2)
        plt.xlabel('Bridge strength')
        plt.ylabel('Number of contig to contig links')
        plt.tight_layout()
        hist_filename=f'{dataset_name}_core_contigs_bs_distribution.png'
        plot_outpath=os.path.join(output_path, hist_filename)
        plt.savefig(plot_outpath, dpi=300, figsize=(4, 4))
        plt.close()        
            
    return core_contigs_list


############
# Write core contigs list.
############

def write_contig_list(dataset_name, core_contigs_list, bridge_strength_thr, mapq_threshold, alg_iteration, output_path):
    
    fileout=open(os.path.join(output_path, f'{dataset_name}_mapq_{mapq_threshold}_br_strength_{bridge_strength_thr}__iteration_{alg_iteration}.txt'), 'w')
    
    for contig_name in core_contigs_list:
        fileout.write(f'{contig_name}\n')
        
    fileout.close()
    
    return


############
# Wrapper function.
############

def linked_contig_wrapper(dataset_name, sam_file_path, path_core_cont_list, mapq_threshold, 
                          bridge_strength_thr, max_num_iterations, samtools_bin_path, output_path):
    
    # Prepare alignment to estimate coverage depth.
    coverage_file_path=sam_file_path[:-4]+".coverage.tab"
    if not os.path.isfile(coverage_file_path):
        
        c_time=time.localtime()
        current_time=time.strftime("%H:%M:%S", c_time) 
        print(f'[{current_time}] Preparation of coverage depth file started.')
        logfile.write(f'[{current_time}] Preparation of coverage depth file started.\n')
        
        bamfile_path=sam_file_path[:-4]+".bam"
        if not os.path.isfile(bamfile_path):
            
            c_time=time.localtime()
            current_time=time.strftime("%H:%M:%S", c_time) 
            print(f'[{current_time}] Bam file preparation started.')
            logfile.write(f'[{current_time}] Bam file preparation started.\n')
            
            samtools_view_command=f'{samtools_bin_path} view -@ 16 -b {sam_file_path} > {bamfile_path}'
            subprocess.run(samtools_view_command, shell=True)
            
        if not os.path.isfile(f'{bamfile_path}.sorted'):
            
            c_time=time.localtime()
            current_time=time.strftime("%H:%M:%S", c_time) 
            print(f'[{current_time}] Bam file sorting started.')
            logfile.write(f'[{current_time}] Bam file sorting started.\n')
            
            samtools_bam_sort=f'{samtools_bin_path} sort -@ 16 {bamfile_path} > {bamfile_path}.sorted'
            subprocess.run(samtools_bam_sort, shell=True)
            
        if not os.path.isfile(f'{bamfile_path}.sorted.bai'):
            
            current_time=time.strftime("%H:%M:%S", c_time) 
            print(f'[{current_time}] Bam file indexing started.')
            logfile.write(f'[{current_time}] Bam file indexing started.\n')
            
            samtools_index_command=f'{samtools_bin_path} index -@ 16 {bamfile_path}.sorted'
            subprocess.run(samtools_index_command, shell=True)
               
        c_time=time.localtime()
        current_time=time.strftime("%H:%M:%S", c_time) 
        print(f'[{current_time}] Coverage depth calculation started.')
        logfile.write(f'[{current_time}] Coverage depth calculation started.\n')    
        
        samtools_coverage_command=f'{samtools_bin_path} coverage {bamfile_path}.sorted > {coverage_file_path}'
        subprocess.run(samtools_coverage_command, shell=True)
        
        c_time=time.localtime()
        current_time=time.strftime("%H:%M:%S", c_time) 
        print(f'[{current_time}] Removing temporary files.')
        logfile.write(f'[{current_time}] Removing temporary files.\n')  
        
        rm_command=f'rm {bamfile_path} {bamfile_path}.sorted {bamfile_path}.sorted.bai'
        subprocess.run(rm_command, shell=True)
        
        c_time=time.localtime()
        current_time=time.strftime("%H:%M:%S", c_time) 
        print(f'[{current_time}] Preparation of coverage depth file finished.')
        logfile.write(f'[{current_time}] Preparation of coverage depth file finished.\n') 
        
    else:
        
        c_time=time.localtime()
        current_time=time.strftime("%H:%M:%S", c_time) 
        print(f'[{current_time}] Coverage depth iformation found: {coverage_file_path}')
        logfile.write(f'[{current_time}] Coverage depth iformation found: {coverage_file_path}\n')         
    
    # Read sam file. Build contig to contig graph.
    contig_to_contig_dict=read_sam_file(sam_file_path, mapq_threshold)
    
    # Read the initial core contigs list.
    core_contigs_list=read_core_contigs(path_core_cont_list)
    
    # Get coverage depth for core contigs.
    core_contigs_cov_depth(dataset_name, coverage_file_path, core_contigs_list, output_path)
    
    # Algorithm iterations indicator.
    alg_iteration=0
    
    # Keep core contig set statistics.
    Num_cont_ar=[]
    Tot_len_cont_ar=[]
    Med_len_cont_ar=[]
    
    # Calculate initial core contigs set statistics.
    Number_of_core_contigs, Tot_length_of_contigs, Med_length_of_contigs=core_contigs_stat(core_contigs_list, alg_iteration)
    Num_cont_ar.append(Number_of_core_contigs)
    Tot_len_cont_ar.append(Tot_length_of_contigs)
    Med_len_cont_ar.append(Med_length_of_contigs)
    
    while alg_iteration<=max_num_iterations:
        
        alg_iteration+=1
    
        # Pull contigs by anchor contigs.
        core_contigs_list=pull_contigs_by_anchors(dataset_name, contig_to_contig_dict, core_contigs_list, bridge_strength_thr, mapq_threshold, alg_iteration, output_path)
        
        # Calculate new core contigs set statistics.
        Number_of_core_contigs, Tot_length_of_contigs, Med_length_of_contigs=core_contigs_stat(core_contigs_list, alg_iteration)
        
        if Number_of_core_contigs>Num_cont_ar[-1]:
            Num_cont_ar.append(Number_of_core_contigs)
            Tot_len_cont_ar.append(Tot_length_of_contigs)
            Med_len_cont_ar.append(Med_length_of_contigs)
            
            # Write new core contigs list.
            write_contig_list(dataset_name, core_contigs_list, bridge_strength_thr, mapq_threshold, alg_iteration, output_path)
        
        else:
            Num_cont_ar.append(Number_of_core_contigs)
            Tot_len_cont_ar.append(Tot_length_of_contigs)
            Med_len_cont_ar.append(Med_length_of_contigs)
            
            # Write final core contigs list.
            write_contig_list(dataset_name, core_contigs_list, bridge_strength_thr, mapq_threshold, alg_iteration, output_path)            
            
            print(f'Algorithm converged at iteration {alg_iteration}.')
            
            break
    
    return


linked_contig_wrapper(Dataset_name, Sam_file_path, Path_core_cont_list, Mapq_threshold, 
                      Bridge_strength_thr, Max_num_iterations, Samtools_path, Output_path)


c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time) 
print(f'[{current_time}] Core contigs set expansion finished.')
logfile.write(f'[{current_time}] Core contigs set expansion finished.\n')
print(f'[{current_time}] Total working time for core contigs set expansion, min: {(time.time()-starttime)/60}\n')
logfile.write(f'[{current_time}] Total working time for core contigs set expansion, min: {(time.time()-starttime)/60}\n')
logfile.close()
