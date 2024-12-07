###############################################
##Anastasiia Rusanova, Dmitry Sutormin, 2023##
### Find contigs shared between several metagenomic bins.
############

import sys
import argparse
import os
import time
import numpy as np
from Bio import SeqIO

parser=argparse.ArgumentParser(description="Find contigs shared between metagenomic bins.")
parser.add_argument("-dn", "--dataset_name", required=True, help="Dataset name.")
parser.add_argument("-bf", "--files", required=True, nargs='+', help="List of files with contigs.")
parser.add_argument("-o", "--output", help="Path to output directory.")
args, unknown=parser.parse_known_args()

Dataset_name=args.dataset_name
cont_files_list=args.files
output_path=args.output

Logfile_path=os.path.join(output_path, 'Core_contigs_set_construction.log')
logfile=open(Logfile_path, 'w', buffering=1)

c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Collecting shared contigs, construction of core set of contigs.', file=sys.stdout)
logfile.write(f'[{current_time}] Construction of core set of contigs started.\n')

################
# Collect initial core contigs.
################

def list_of_shared_contigs(Dataset_name, cont_files_list, output_path):
    
    contigs_dict={}
    
    for infile_path in cont_files_list:
        
        for record in SeqIO.parse(infile_path, "fasta"):
            record_id=record.id        
            if record_id not in contigs_dict:
                contigs_dict[record_id]=[record]
            else:
                contigs_dict[record_id].append(record)
    
    core_contigs_list_outpath=os.path.join(output_path, f'{Dataset_name}_initial_core_contigs.list')
    core_contigs_fasta_outpath=os.path.join(output_path, f'{Dataset_name}_initial_core_contigs.fasta')
    list_outfile=open(core_contigs_list_outpath, "w")
    fasta_outfile=open(core_contigs_fasta_outpath, "w")
    core_contigs_list_stat=[]
    for contig_name, records_ar in contigs_dict.items():
        if len(records_ar)==len(cont_files_list):
            list_outfile.write("{}\n".format(contig_name))
            fasta_outfile.write(f'>{records_ar[0].id}\n{str(records_ar[0].seq)}\n')
            core_contigs_list_stat.append(len(str(records_ar[0].seq)))
    
    list_outfile.close()
    fasta_outfile.close()
    
    Tot_len=np.sum(core_contigs_list_stat)
    Med_len=np.median(core_contigs_list_stat)
    
    c_time=time.localtime()
    current_time=time.strftime("%H:%M:%S", c_time)    
    print(f'[{current_time}] Core contigs set total length: {Tot_len} bp')
    logfile.write(f'[{current_time}] Core contigs set total length: {Tot_len} bp\n')
    c_time=time.localtime()
    current_time=time.strftime("%H:%M:%S", c_time)    
    print(f'[{current_time}] Core contigs set median length: {Med_len} bp')
    logfile.write(f'[{current_time}] Core contigs set median length: {Med_len} bp\n')    
    
    return

starttime=time.time()

list_of_shared_contigs(Dataset_name, cont_files_list, output_path)

c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time) 
print(f'[{current_time}] Core contigs set construction finished.')
logfile.write(f'[{current_time}] Core contigs set construction finished.\n')
print(f'[{current_time}] Total working time for core contigs set construction, min: {(time.time()-starttime)/60}\n')
logfile.write(f'[{current_time}] Total working time for core contigs set construction, min: {(time.time()-starttime)/60}\n')
logfile.close()