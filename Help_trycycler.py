##Dmitry Sutormin, 2023##
### Help trycycler.
############

import os
import subprocess
import re
from Bio import SeqIO
import numpy as np

# Trycycler directory.
Trycycler_path="/home/niagara/Storage/D_Sutormin/Metagenomes/N_Rusanova/Spongy/2022_nanopore/I_palmata/I_palmata_2016_2022_core_assebmlies/trycycler_cluster/"


def mark_not_reconciled(trycycler_path):

    Clusters_ar=os.listdir(trycycler_path)
    
    count_rec=0
    count_not_rec=0
    
    for cluster_dir in Clusters_ar:
        reconc_file=os.path.join(trycycler_path, cluster_dir, '2_all_seqs.fasta')
        if not os.path.isfile(reconc_file):
            mv_command=f'mv {os.path.join(trycycler_path, cluster_dir)} {os.path.join(trycycler_path, "Not_reconcile_"+cluster_dir)}'
            subprocess.run(mv_command, shell=True)
            count_not_rec+=1
        else:
            count_rec+=1
            
    print(f'Number of clusters succesfully reconciled: {count_rec}')
    print(f'Number of clusters NOT reconciled: {count_not_rec}')
    
    return

#mark_not_reconciled(Trycycler_path)

def mark_single(trycycler_path):

    Clusters_ar=os.listdir(trycycler_path)
    
    count_not_rec_sing=0
    count_not_rec_mult=0
    
    for cluster_dir in Clusters_ar:
        if 'Not_reconcile_' in cluster_dir:
            Init_contigs_clust=os.listdir(os.path.join(trycycler_path, cluster_dir, '1_contigs'))
            if len(Init_contigs_clust)==1:
                mv_command=f'mv {os.path.join(trycycler_path, cluster_dir)} {os.path.join(trycycler_path, "Single_"+cluster_dir)}'
                subprocess.run(mv_command, shell=True)
                count_not_rec_sing+=1
            else:
                count_not_rec_mult+=1
            
    print(f'Number of clusters NOT clustered NOT reconciled: {count_not_rec_sing}')
    print(f'Number of clusters clustered NOT reconciled: {count_not_rec_mult}')
    
    return

#mark_single(Trycycler_path)


def collect_Trycycler_output(trycycler_path):
    
    fasta_list=[]
    fasta_len_ar=[]
    
    Clusters_ar=os.listdir(trycycler_path)
    
    for cluster_dir in Clusters_ar:
        if re.match('^cluster_', cluster_dir):    
            res_fasta=os.path.join(trycycler_path, cluster_dir, '7_final_consensus.fasta')
            for record in SeqIO.parse(res_fasta, "fasta"):
                fasta_list.append(record)
                fasta_len_ar.append(len(str(record.seq)))
        elif re.match('^Single_Not_', cluster_dir): 
            for contig_file in os.listdir(os.path.join(trycycler_path, cluster_dir, '1_contigs')):
                if 'fasta' in contig_file:
                    res_fasta=os.path.join(trycycler_path, cluster_dir, '1_contigs', contig_file)
                    print(res_fasta)
                    for record in SeqIO.parse(res_fasta, "fasta"):
                        fasta_list.append(record)
                        fasta_len_ar.append(len(str(record.seq)))
        elif re.match('^Not_reconcile_', cluster_dir): 
            longest_record=0
            longest_record_length=0
            for contig_file in os.listdir(os.path.join(trycycler_path, cluster_dir, '1_contigs')):
                if 'fasta' in contig_file:
                    res_fasta=os.path.join(trycycler_path, cluster_dir, '1_contigs', contig_file)
                    for record in SeqIO.parse(res_fasta, "fasta"): 
                        if len(str(record.seq))>longest_record_length:
                            longest_record_length=len(str(record.seq))
                            longest_record=record
            fasta_list.append(longest_record)
            fasta_len_ar.append(longest_record_length)
    
    print(f'Number of contigs collapsed with Trycycler: {len(fasta_len_ar)}')  
    print(f'Total length of contigs collapsed with Trycycler: {np.sum(fasta_len_ar)} bp')
    print(f'Median length of contigs collapsed with Trycycler: {np.median(fasta_len_ar)} bp')
    
    SeqIO.write(fasta_list, os.path.join(Trycycler_path, 'Trycycler_final_contigs.fasta'), 'fasta')
    
    return

collect_Trycycler_output(Trycycler_path)

