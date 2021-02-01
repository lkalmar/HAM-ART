# HAM-ART
Hi-C metagenomics pipeline for tracking antimicrobial resistance genes in complex microbial communities
# Introduction
HAM-ART is an easy to use complete bioinformatics pipeline for generating metagenome assembled genomes (MAGs) from the combination of traditional metagenomics and Hi-C metagenomics sequencing. The pipeline uses its own project management system to make the workflow user friendly (while it runs in a Linux command line environment, it doesn't require advanced command line experience). The pipeline can handle pooled Hi-C sequencing libraries, and use them in combination with multiple metagenomics libraries.
# Installation and configuration
## Installing dependecies
HAM-ART requires a generic Linux command line environment with PERL programming language interpreter installed (installed together with the operating system in most of the Linux distributions). Before using the HAM-ART pipeline you need to install a few software the pipeline uses. Some of these software doesn't provide a conventional installer (e.g. conda), but only a download option. Downloaded software are not accessible anywhere in your filesystem (except if you put them in your PATH), but the HAM-ART dependency configuration file can deal with the problem (see below). The list of required software:
* clumpify.sh from the [BBMap package](https://sourceforge.net/projects/bbmap/)
* [Bowtie2 aligner](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [FLASh](https://ccb.jhu.edu/software/FLASH/)
* metaSPAdes from the [SPAdes software package](https://cab.spbu.ru/software/spades/)
* [ABRicate](https://github.com/tseemann/abricate)
* [Samtools](http://www.htslib.org)
* [louvain](https://sourceforge.net/projects/louvain/) - Louvain method for community detection in large graphs
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [MASH](https://mash.readthedocs.io/en/latest/index.html)
* [GTDB-tk](https://github.com/Ecogenomics/GTDBTk)
* [Julia programming language interpeter](https://julialang.org)
## Setting up the environment
1. Create a directory in your Linux filesystem that will serve as the root directory for the pipeline (_pipeline root folder_). Please make sure you have plenty of space allocated to the mountpoint the folder in (at least hundreds of GB, but for multiple projects several TB will be needed). All data related to your Hi-C metagenomics projects (from raw data to results) will be stored in the subfolders of this _pipeline root folder_.
1. Create the following sub-folders in your _pipeline root folder_: ```raw_seq_data```, ```pre_assembly```, ```assembly```, ```post_assembly```, ```results```, ```project_defs```, ```host_ref_seqs```
1. Download the HAM-ART repository and copy the following files to your _pipeline root folder_: ```ham_art.pl```, ```deps.txt```, ```ham_art_prj_template.txt```
1. Make sure the ```ham_art.pl``` script is executable, if not you can change it by using the command line command ```chmod u+x ham_art.pl```
1. Edit the ```deps.txt``` file, to set up the dependencies. Lines starting with ```#``` are comments only, read them to understand the syntax. For each software dependencies, you need to set up the path to the executable/folder.
   1. __This file has to be edited only once, if you are not familiar with Linux filesystem and command line environment, you may ask someone with more experience for help.__
   1. Don't change the first field in the line, those are the internal identifiers for the software dependencies
   1. The second fields are defining the amount of CPU core and/or memory to be used when running the executable. Please read the comment lines to get the information. The amount set in this file are the optimal CPU core / memory amounts, usually only a little performace increase can be achieved by increasing them (except the metaspades.py, if you have more than 32, increase it). If you have less CPU cores in your Linux server, decrease the amount to that value where needed
   1. If the executable path is in your PATH variable (you can execute them by only calling the application name even from the pipeline root folder), you don't need to edit the 3rd field in the line
   1. If the executable path is not in your PATH variable, or you want to use a specific version for the pipeline (e.g. if you download a new version for SPAdes), you can include the full application path in the 3rd field of the line. Simply replace the application name with its full path (e.g. ```metaspades.py``` -> ```/home/user/SPAdes-3.13.0-Linux/bin/metaspades.py```)
   1. Please be aware that in a few cases the line defines a path to a folder, not to the executable, please read the comment lines carefully!
1. If you move the ```ham_art.pl``` main script file to a different location, always make sure, that you move it together with the ```deps.txt``` file, these two files has to be in the same directory.
1. Download and index host genome(s). The HAM-ART pipeline cleans up the metagenomics sequencing data by removing those reads that are aligning to the host genome (e.g. if you are working with a human faecal sample, the potentially sequenced gut epithelial cell DNA). Make sure that you have the host organism whole genome sequence downloaded, and put the indexed (bowtie2 indexed) files in the ```host_ref_seqs/``` folder. For easier use, use simple prefix for the indexed files, you will need to use that prefix later in the project definition file (for an example workflow, see the example project below). 
# Running the HAM-ART pipeline
The HAM-ART pipeline has its own project management system to make the pipeline steps and analysing multiple samples easy. For the analysis, you have to define a project identifier for each sample, define the input files, the project category, and the host organism.
## Setting up a project category
The project category in the HAM-ART pipeline is usually a well identifiable name for a research study that involves multiple samples. During the pipeline run, all raw data, temporary files and final results will be stored in project category specific folders within the pipeline sub-folders. It is also very important to note, that the clade-refinement step ot the pipeline is always performed within one project category. To set up a project category, you simply create a sub-folder in the ```raw_seq_data/``` folder, and copy your raw sequencing data files there. Please do not use whitespace (space, tabulator) in the folder and file names, you can use underscore instead (for an example of project category setup, see the example workfolw below). You don't need to create the project category specific sub-folders in other folders (e.g. in ```results/```), those will be created automatically.
## Setting up a project
In the HAM-ART pipeline individual samples are analysed as projects, each project with its own unique identifier. The unique project identifier is used during the analysis steps, and in the combination of multiple samples during the clade-refinement. To set up a new project, you need the following steps (for an example setup, see the example workflow below):
1. Copy the ```ham_art_prj_template.txt``` to a new filename within the ```project_defs/``` folder (ideally by naming it with the unique identifier)
1. Edit the new file and define the project identifier, the raw data file names, the project folder and the host organism prefix in the appropriate lines
1. After saving the new file, initiate the project by being in the ```project_defs/``` folder and executing the command ```../ham_art.pl initialise <project_definition_filename>``` (where you replace the spaceholder with your filename).
## Running the pipeline steps
The HAM-ART pipeline 3 main steps and one optional step, that you have to execute from the _pipeline root folder_ are the following (in each case ```<project_ID(s)``` can be a single ID or multiple IDs separated by comma (see the example workflow below)):
* __Pre-assembly.__ De-duplicate paired-end raw reads, remove host organism matching raw reads, quality filtering and merge overlapping paired-end reads. Run this step by:
```
./ham_art.pl pre_assembly <project_ID(s)>
```
* __Assembly.__ Re-fragment Hi-C reads based on restriction site sequence pattern (Please note, that the HAM-ART pipeline can only handle single enzyme digested Hi-C libraries by the _HpyCH4IV_ restriction enzyme) and do the _de novo_ assembly. Run this step by:
```
./ham_art.pl assembly <project_ID(s)>
```
* __Post-assembly.__ Re-align Hi-C reads, extract inter-contig contacts, build up the contact network, find communities in the network, split potentially mixed communities of contigs (based on contig coverage), perform basic MAG extension and create final MAGs. Run this step by:
```
./ham_art.pl post_assembly <project_ID(s)>
```
* __Clade-refinement.__ If you have multiple samples / projects within the same project category, you can significantly increase the number and quality of you final MAGs by using the HAM-ART pipeline specific clade-refinement step. For the detailed description of this step, please look at the original research paper (cited below). For this step, you need more than 1 project to run:
```
./ham_art.pl clade_refinement <project_IDs>
```
## Output files and folders
During the HAM-ART pipeline run most of the temporary files and folders are deleted, but some of them (together with the results) are preserved for potential future analysis:
* Final MAGs for a single project before clade refinement are in the folder ```results/<project_category>/<project_ID>_clusters/```
* Final MAGs for the involved projects after clade refinement are in their clade specific folder within ```results/<project_category>/combined_results_.../refined_CC_files/allbig_by_clade_extended/```
* Information on final MAGs and their AMR gene associations are in the file ```results/<project_category>/combined_results_.../refined_CC_files/analysis_results/full_amr_table_filtered.txt```
* Clade specific taxonomy, MAG quality, distance trees, iTOL annotations and AMR association information are in the folder ```results/<project_category>/combined_results_.../refined_CC_files/analysis_results/clade_summary_files```
* Project specific _de novo_ assembly results (metaSPAdes contigs.fasta outputs), binary inter-contig information, initial consensus cluster definitions and AMR gene containing contigs are in the folder ```post_assembly/<project_category>/```
* Filtered, deduplicated, dehosted and merged raw read libraries (```LIB_PEL1_...``` and ```LIB_PEL2_...``` for paired-end reads, ```LIB_SEL_...``` files for merged reads and re-fragmentes Hi-C reads) are in the folder ```pre_assembly/<project_category>/```
# Example workflow
Coming soon...
# Citation
Kalmar L, Gupta S, Kean IRL, Ba X, Hadjirin N, Lay E, Vries S, Bateman M, Bartlett H, Hernandez-Garcia J, Tucker AW, Restiff O, Stevens MP, Wood J, Maskell DJ, Grant AJ, Holmes M.
HAM-ART: An optimised culture-free Hi-C metagenomics pipeline for tracking antimicrobial resistance genes in complex microbial communities.
(Under submission)
