###############################################
##Dmitry Sutormin, 2023##
### Split contigs into groups based on initial source assembly.
############

import argparse
import os
import sys
from Bio import SeqIO
import time

parser=argparse.ArgumentParser(description="Split expanded core contigs by original assemblies.")
parser.add_argument("-dn", "--dataset_name", required=True, help="Dataset name.")
parser.add_argument("-ecf", "--expanded_contigs", required=True, help="Fasta file with expanded core set of contigs.")
parser.add_argument("-mqt", "--mapq_thr", required=True, help="Mapq threshold to filter alignments.")
parser.add_argument("-bst", "--bridge_str_thr", required=True, help="Contig-to-contig bridge strength threshold to filter contig links.")
parser.add_argument("-af", "--assembly_files", nargs='+', help="List of assemblies from which bins were derived.")
parser.add_argument("-o", "--output", help="Path to output directory.")
args, unknown=parser.parse_known_args()

Dataset_name=args.dataset_name
Contig_file_path=args.expanded_contigs
Mapq_threshold=int(args.mapq_thr)
Bridge_strength_thr=int(args.bridge_str_thr)
Assembly_path_list=args.assembly_files
Output_path=args.output



def split_contigs_wrapper(dataset_name, contig_file_path, mapq_threshold, 
                          bridge_strength_thr, assembly_path_list, output_path):    
    
    # Create output directory.
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
        
    # Dictionary of initias source assemblies.
    assembly_path_dict={}
    
    for assembly_file_path in assembly_path_list:
        assembly_file_name=assembly_file_path.split('/')[-1][:-6]
        assembly_path_dict[assembly_file_name]=assembly_file_path
    
    # Get composition of initial assemblies.
    Assembly_compos_dict={}
    
    for assembly_name, assembly_path in assembly_path_dict.items():
        for record in SeqIO.parse(assembly_path, "fasta"):
            record_id=record.id
            Assembly_compos_dict[record_id]=assembly_name
            
    # Split contigs into groups.
    Cont_groups_dir={}
    
    for record in SeqIO.parse(contig_file_path, "fasta"):
        record_id=record.id 
        record_seq=str(record.seq)
        Source_assembly=Assembly_compos_dict[record_id]
        if Source_assembly in Cont_groups_dir:
            Cont_groups_dir[Source_assembly].append(record)
        else:
            Cont_groups_dir[Source_assembly]=[record]
     
    # Write contigs groups into separate fasta files.
    for source_assembly, list_of_contigs in Cont_groups_dir.items():
        SeqIO.write(list_of_contigs, os.path.join(output_path, f'{dataset_name}_{source_assembly}_core_contigs_br_thr_{bridge_strength_thr}_mapq_thr_{mapq_threshold}.fasta'), 'fasta')
          
    return


starttime=time.time()

split_contigs_wrapper(Dataset_name, Contig_file_path, Mapq_threshold, 
                      Bridge_strength_thr, Assembly_path_list, Output_path)

print(f'Total working time, min: {(time.time()-starttime)/60}')

