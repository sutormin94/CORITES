# CORITES
 Core contigs iterative expansion and scaffolding pipeline
 
## Description 
The CORITES algorithm is developed for refinment of metagenomic bins. 

**Step 0:** Metagenome binning

CORITES takes bins prepared by different binners from the same metagenome which represent the same genome.
E.g., these bins are classified identically (e.g., by GTDB) or sharing significant similarity.

**Step 1:** Identification of a core set of contigs

A core set of contigs is built from contigs shared between the selected bins. 

**Step 2:** Iterative expansion of the core set of contigs

The core set is iteratively expanded by other contigs derived from the binned assembly, using information from long-reads aligned to the assembly with minimap2. 
A new contig is included in the core set if it was bridged with any of the core set contigs by the number of long reads exceeding the threshold value specified with the bridge_strength parameter. 
Usually, the iterative expansion module converged after 3-8 iterations without dramatic inflation of the initial core set. Inflation is observed only with low bridge_strength parameter values. 

**Step 3:** Scaffolding with longstitch

The expanded core set is scaffolded with longstitch using long reads. 

**Step 4:** Scaffolding with assembly graph

Obtained final bins are additionally polished through the alignment with the assembly graph generated by SPAdes. 
A pair of contigs from the bin was combined if their order can be unambiguously determined using the topology of the assembly graph. 
If such ordering exists, the two contigs were merged, and the gap between them was filled with the highest covered path in the assembly graph.


## Usage

### Iterative expansion module

```
>python Main_pipeline.py --help
usage: Main_pipeline.py [-h] -dn DATASET_NAME -bf BIN_FILES [BIN_FILES ...] [-af ASSEMBLY_FILES [ASSEMBLY_FILES ...]] [-rf LONG_READ_FILES [LONG_READ_FILES ...]] [-mm2 MINIMAP2] [-rn]
                        -mqt MAPQ_THR -bst BRIDGE_STR_THR -mi MAX_ITERATIONS -rsm {seqtk,manual} -sp SEQTK_PATH -samtools SAMTOOLS_PATH -ChM2 CHECKM2_ENV -qst QUAST_PATH [-o OUTPUT]

CORITES pipeline main script.

optional arguments:
  -h, --help            show this help message and exit
  -dn DATASET_NAME, --dataset_name DATASET_NAME
                        Dataset name.
  -bf BIN_FILES [BIN_FILES ...], --bin_files BIN_FILES [BIN_FILES ...]
                        List of files with contigs (fasta).
  -af ASSEMBLY_FILES [ASSEMBLY_FILES ...], --assembly_files ASSEMBLY_FILES [ASSEMBLY_FILES ...]
                        List of assemblies from which bins were derived (fasta).
  -rf LONG_READ_FILES [LONG_READ_FILES ...], --long_read_files LONG_READ_FILES [LONG_READ_FILES ...]
                        List of fastq files with long reads to be used for bridging and scaffolding (fastq).
  -mm2 MINIMAP2, --minimap2 MINIMAP2
                        Path to minimap2 read aligner executable.
  -rn, --run_alignment  Run alignment of all long reads vs all contigs or not.
  -mqt MAPQ_THR, --mapq_thr MAPQ_THR
                        Mapq threshold to filter alignments.
  -bst BRIDGE_STR_THR, --bridge_str_thr BRIDGE_STR_THR
                        Contig-to-contig bridge strength threshold to filter contig links.
  -mi MAX_ITERATIONS, --max_iterations MAX_ITERATIONS
                        Maximal number of iterations to expand core contigs set.
  -rsm {seqtk,manual}, --read_search_mode {seqtk,manual}
                        Algorithm to use for retrival of reads sequences.
  -sp SEQTK_PATH, --seqtk_path SEQTK_PATH
                        Path to seqtk executable.
  -samtools SAMTOOLS_PATH, --samtools_path SAMTOOLS_PATH
                        Path to samtools executable.
  -ChM2 CHECKM2_ENV, --CheckM2_env CHECKM2_ENV
                        Conda environment with CheckM2 installed and avaliable in path.
  -qst QUAST_PATH, --quast_path QUAST_PATH
                        Path to quast executable.
  -o OUTPUT, --output OUTPUT
                        Path to output directory.
```


**Requirements:**

For the initial run, metagenomic bins (in fasta), full metagenome assembly(s) (in fasta), and long reads (in fastq) should be provided as well as the run_alignment option to generate the reads alignment.
For further runs (e.g., when adjusting the mapq_thr, bridge_str_thr, and max_iterations paremeters), the run_alignment can be skipped.


