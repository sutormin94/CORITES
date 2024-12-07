###############################################
##Anastasiia Rusanova, Dmitry Sutormin, 2024##
### CORITES pipeline wrapper script: identifies core contigs and expand the core using information from long reads.
############

import sys
import argparse
import os
import subprocess
import time

parser=argparse.ArgumentParser(description="Main script for bin improvement pipeline.")
parser.add_argument("-dn", "--dataset_name", required=True, help="Dataset name.")
parser.add_argument("-bf", "--bin_files", required=True, nargs='+', help="List of files with contigs.")
parser.add_argument("-af", "--assembly_files", nargs='+', help="List of assemblies from which bins were derived.")
parser.add_argument("-rf", "--long_read_files", nargs='+', help="List of fastq files with long reads to be used for bridging and scaffolding.")
parser.add_argument("-mm2", "--minimap2", default="minimap2", help="Path to minimap2 read aligner.")
parser.add_argument("-rn", "--run_alignment", action='store_true', help="Run alignment of all long reads vs all contigs or not.")
parser.add_argument("-mqt", "--mapq_thr", required=True, help="Mapq threshold to filter alignments.")
parser.add_argument("-bst", "--bridge_str_thr", required=True, help="Contig-to-contig bridge strength threshold to filter contig links.")
parser.add_argument("-mi", "--max_iterations", required=True, help="Maximal number of iterations to expand core contigs set.")
parser.add_argument("-rsm", "--read_search_mode", required=True, default='seqtk', choices=['seqtk', 'manual'], help="Algorithm to use for retrival of reads sequences.")
parser.add_argument("-sp", "--seqtk_path", required=True, default='seqtk', help="Path to seqtk executable.")
parser.add_argument("-samtools", "--samtools_path", required=True, default='samtools', help="Path to samtools executable.")
parser.add_argument("-ChM2", "--CheckM2_env", required=True, default='checkm2', help="Conda environment with CheckM2 installed and avaliable in path.")
parser.add_argument("-qst", "--quast_path", required=True, default='quast.py', help="Path to quast executable.")
parser.add_argument("-o", "--output", help="Path to output directory.")
args, unknown=parser.parse_known_args()

Dataset_name=args.dataset_name
cont_files_list=args.bin_files
if args.assembly_files:
    assembly_files_list=args.assembly_files
if args.long_read_files:
    long_read_files_list=args.long_read_files
Minimap2_path=args.minimap2
if args.run_alignment:
    run_universal_alignment=True
Mapq_threshold=int(args.mapq_thr)
Bridge_strength_thr=int(args.bridge_str_thr)
Max_num_iterations=int(args.max_iterations)
Mode=args.read_search_mode
Seqtk_path=args.seqtk_path
Samtools_path=args.samtools_path
checkM_env=args.CheckM2_env
Quast_path=args.quast_path
output_path=args.output

starttime=time.time()
Logfile_path=os.path.join(output_path, 'Core_contigs_set_expansion_pipeline.log')
logfile=open(Logfile_path, 'w', buffering=1)

c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Core contigs set expansion pipeline started.', file=sys.stdout)
logfile.write(f'[{current_time}] Core contigs set expansion pipeline started.\n')


current_script_directory=os.path.dirname(os.path.abspath(sys.argv[0]))


# Gets core contigs.
c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Construction of core contigs set started.', file=sys.stdout)
logfile.write(f'[{current_time}] Construction of core contigs set started.\n')

core_contigs_fasta_outpath=os.path.join(output_path, f'{Dataset_name}_initial_core_contigs.fasta')

if not os.path.isfile(core_contigs_fasta_outpath):
    core_prep_script=os.path.join(current_script_directory, 'Find_shared_contigs_from_different_binners.py')
    Collect_core_command=f'python3 {core_prep_script} --dataset_name {Dataset_name} --output {output_path} --files'
    for file_path in cont_files_list:
        Collect_core_command+=f' {file_path}'
        
    subprocess.run(Collect_core_command, shell=True)


Fin_bins_path=os.path.join(output_path, 'Final_expanded_bins')
if not os.path.isdir(Fin_bins_path):
    os.mkdir(Fin_bins_path)
    
core_contigs_fasta_fin_outpath=os.path.join(output_path, 'Final_expanded_bins', f'{Dataset_name}_initial_core_contigs.fasta')

if not os.path.isfile(core_contigs_fasta_fin_outpath):
    subprocess.run(f'cp {core_contigs_fasta_outpath} {Fin_bins_path}', shell=True)

c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Construction of core contigs set finished.', file=sys.stdout)
logfile.write(f'[{current_time}] Construction of core contigs set finished.\n')


# Prepares universal alignment. Concatenates fasta files together and fastq files together, then aligns with minimap2.
c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Merging of metagenomic assemblies and long read dataset started.', file=sys.stdout)
logfile.write(f'[{current_time}] Merging of metagenomic assemblies and long read dataset started.\n')

if args.assembly_files:
    
    seed_assembly=assembly_files_list[0]
    
    c_time=time.localtime()
    current_time=time.strftime("%H:%M:%S", c_time)
    print(f'[{current_time}] Processing {seed_assembly}.', file=sys.stdout)
    logfile.write(f'[{current_time}] Processing {seed_assembly}.\n')
    
    universal_assembly_path=os.path.join(output_path, f'{Dataset_name}_all_contigs.fasta')
    
    if not os.path.isfile(universal_assembly_path):
        subprocess.run(f'cp {seed_assembly} {universal_assembly_path}', shell=True)
        
        if len(assembly_files_list)>1:
            
            for i in range(len(assembly_files_list)-1):
                
                current_assembly=assembly_files_list[i+1]
                
                c_time=time.localtime()
                current_time=time.strftime("%H:%M:%S", c_time)
                print(f'[{current_time}] Processing {current_assembly}.', file=sys.stdout)
                logfile.write(f'[{current_time}] Processing {current_assembly}.\n')
                
                subprocess.run(f'cat {current_assembly} >> {universal_assembly_path}', shell=True)
else:
    
    universal_assembly_path=os.path.join(output_path, f'{Dataset_name}_all_contigs.fasta')
    
if args.long_read_files:
    
    seed_fastq=long_read_files_list[0]
    
    c_time=time.localtime()
    current_time=time.strftime("%H:%M:%S", c_time)
    print(f'[{current_time}] Processing {seed_fastq}.', file=sys.stdout)
    logfile.write(f'[{current_time}] Processing {seed_fastq}.\n')
    
    universal_fastq_path=os.path.join(output_path, f'{Dataset_name}_all_long_reads.fastq')
    
    if not os.path.isfile(universal_fastq_path):
    
        if seed_fastq[-3:]=='.gz':
            subprocess.run(f'gunzip {seed_fastq} -c > {os.path.join(output_path, seed_fastq.split("/")[-1][:-3])}', shell=True)
            seed_fastq=os.path.join(output_path, seed_fastq.split("/")[-1][:-3])
        
        subprocess.run(f'cp {seed_fastq} {universal_fastq_path}', shell=True)
        
        if len(long_read_files_list)>1:
            
            for i in range(len(long_read_files_list)-1):
                current_fastq=long_read_files_list[i+1]
                
                c_time=time.localtime()
                current_time=time.strftime("%H:%M:%S", c_time)
                print(f'[{current_time}] Processing {current_fastq}.', file=sys.stdout)
                logfile.write(f'[{current_time}] Processing {current_fastq}.\n')
                
                if current_fastq[-3:]=='.gz':
                    subprocess.run(f'gunzip {current_fastq} -c > {os.path.join(output_path, current_fastq.split("/")[-1][:-3])}', shell=True)
                    current_fastq=os.path.join(output_path, current_fastq.split("/")[-1][:-3])            
                subprocess.run(f'cat {current_fastq} >> {universal_fastq_path}', shell=True)    
else:
    
    universal_fastq_path=os.path.join(output_path, f'{Dataset_name}_all_long_reads.fastq')  
    
c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Merging of metagenomic assemblies and long read dataset finished.', file=sys.stdout)
logfile.write(f'[{current_time}] Merging of metagenomic assemblies and long read dataset finished.\n')

if 'run_universal_alignment' in locals():
    
    c_time=time.localtime()
    current_time=time.strftime("%H:%M:%S", c_time)
    print(f'[{current_time}] Alignment of merged long reads to merged assemblies started.', file=sys.stdout)
    logfile.write(f'[{current_time}] Alignment of merged long reads to merged assemblies started.\n')
    
    universal_sam_path=os.path.join(output_path, f'{Dataset_name}_all_long_reads_vs_all_contigs.sam')
    minimap2_command=f'{Minimap2_path} -ax map-ont {universal_assembly_path} {universal_fastq_path} > {universal_sam_path}'
    subprocess.run(minimap2_command, shell=True)
    
    c_time=time.localtime()
    current_time=time.strftime("%H:%M:%S", c_time)
    print(f'[{current_time}] Alignment of merged long reads to merged assemblies finished.', file=sys.stdout)
    logfile.write(f'[{current_time}] Alignment of merged long reads to merged assemblies finished.\n')    

    
# Iteratively expands core contigs set using long read alignment.
c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Iterative expansion of a core set of contigs started.', file=sys.stdout)
logfile.write(f'[{current_time}] Iterative expansion of a core set of contigs started.\n')

universal_sam_path=os.path.join(output_path, f'{Dataset_name}_all_long_reads_vs_all_contigs.sam')
core_contigs_path=os.path.join(output_path, f'{Dataset_name}_initial_core_contigs.list')
expand_core_script=os.path.join(current_script_directory, 'Find_linked_contigs.py')

Expand_core_path=os.path.join(output_path, 'Expand_contigs_core')
if not os.path.isdir(Expand_core_path):
    os.mkdir(Expand_core_path)

expand_core_command=f'python3 {expand_core_script} --dataset_name {Dataset_name} --sam_file {universal_sam_path} --mapq_thr {Mapq_threshold} --bridge_str_thr {Bridge_strength_thr} --max_iterations {Max_num_iterations} --contigs_list {core_contigs_path} --samtools_path {Samtools_path} --output {Expand_core_path}'
subprocess.run(expand_core_command, shell=True)

c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Iterative expansion of a core set of contigs finished.', file=sys.stdout)
logfile.write(f'[{current_time}] Iterative expansion of a core set of contigs finished.\n')


# Gets sequences for expanded core contigs & sequences of bridging reads.
c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Retrival of contigs and reads sequences for the expanded core contigs set started.', file=sys.stdout)
logfile.write(f'[{current_time}] Retrival of contigs and reads sequences for the expanded core contigs set started.\n')

final_alg_iteration=0
list_of_exp_files=os.listdir(Expand_core_path)
for exp_file in list_of_exp_files:
    if f'Bridges_info_mapq_{Mapq_threshold}_br_strength_{Bridge_strength_thr}__iteration' in exp_file:
        iter_number=int(exp_file.split('__iteration_')[1][:-4])
        if iter_number>final_alg_iteration:
            final_alg_iteration=iter_number

fin_br_file=os.path.join(Expand_core_path, f'Bridges_info_mapq_{Mapq_threshold}_br_strength_{Bridge_strength_thr}__iteration_{final_alg_iteration}.tab')
fin_cont_list=os.path.join(Expand_core_path, f'{Dataset_name}_mapq_{Mapq_threshold}_br_strength_{Bridge_strength_thr}__iteration_{final_alg_iteration}.txt')
get_sequences_script=os.path.join(current_script_directory, 'Collect_core_contigs_and_bridge_reads.py')

get_sequences_command=f'python3 {get_sequences_script} --dataset_name {Dataset_name} --fin_contigs_list {fin_cont_list} --fin_bridge_file {fin_br_file} --mapq_thr {Mapq_threshold} --bridge_str_thr {Bridge_strength_thr} --assembly_file {universal_assembly_path} --long_read_file {universal_fastq_path} --read_search_mode {Mode} --seqtk_path {Seqtk_path} --output {Expand_core_path}'
subprocess.run(get_sequences_command, shell=True)

c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Retrival of contigs and reads sequences for the expanded core contigs set finished.', file=sys.stdout)
logfile.write(f'[{current_time}] Retrival of contigs and reads sequences for the expanded core contigs set finished.\n')


# Splits expanded contigs sets by original assemblies.
c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Splitting of expanded core contigs set by original assemblies started.', file=sys.stdout)
logfile.write(f'[{current_time}] Splitting of expanded core contigs set by original assemblies started.\n')

expanded_contigs_file=os.path.join(Expand_core_path, f'Expanded_core_contigs_{Dataset_name}_mapq_{Mapq_threshold}_br_strength_{Bridge_strength_thr}.fasta')
split_contigs_script=os.path.join(current_script_directory, 'Split_contigs_into_groups.py')

split_contigs_command=f'python3 {split_contigs_script} --dataset_name {Dataset_name} --expanded_contigs {expanded_contigs_file} --mapq_thr {Mapq_threshold} --bridge_str_thr {Bridge_strength_thr} --assembly_files'
if args.assembly_files:
    for file_path in assembly_files_list:
        split_contigs_command+=f' {file_path}'
split_contigs_command+=f' --output {Fin_bins_path}'

subprocess.run(split_contigs_command, shell=True)

c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Splitting of expanded core contigs set by original assemblies finished.', file=sys.stdout)
logfile.write(f'[{current_time}] Splitting of expanded core contigs set by original assemblies finished.\n')


# Runs checkM for assessment of contig sets - original and expanded.
c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Assessment of core contigs sets quality by CheckM started.', file=sys.stdout)
logfile.write(f'[{current_time}] Assessment of core contigs sets quality by CheckM started.\n')

CheckM_outdir_path=os.path.join(Fin_bins_path, 'CheckM')

checkM_command=f'conda run -n {checkM_env} checkm2 predict --input {Fin_bins_path}/*.fasta --output-directory {CheckM_outdir_path} --force --threads 20'
subprocess.run(checkM_command, shell=True)

c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Assessment of core contigs sets quality by CheckM finished.', file=sys.stdout)
logfile.write(f'[{current_time}] Assessment of core contigs sets quality by CheckM finished.\n')


# Runs quast for assessment of contig sets - original and expanded.
c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Assessment of core contigs sets quality by Quast started.', file=sys.stdout)
logfile.write(f'[{current_time}] Assessment of core contigs sets quality by Quast started.\n')

Quast_outdir_path=os.path.join(Fin_bins_path, 'Quast')
if not os.path.isdir(Quast_outdir_path):
    os.mkdir(Quast_outdir_path)
    
quast_command=f'python2 {Quast_path} --output-dir {Quast_outdir_path} --threads 20 {Fin_bins_path}/*.fasta'
subprocess.run(quast_command, shell=True)

c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time)
print(f'[{current_time}] Assessment of core contigs sets quality by Quast finished.', file=sys.stdout)
logfile.write(f'[{current_time}] Assessment of core contigs sets quality by Quast finished.\n')


c_time=time.localtime()
current_time=time.strftime("%H:%M:%S", c_time) 
print(f'[{current_time}] Core contigs set expansion pipeline finished.')
logfile.write(f'[{current_time}] Core contigs set expansion pipeline finished.\n')
print(f'[{current_time}] Total working time for core contigs set expansion, min: {(time.time()-starttime)/60}\n')
logfile.write(f'[{current_time}] Total working time for core contigs set expansion, min: {(time.time()-starttime)/60}\n')
logfile.close()