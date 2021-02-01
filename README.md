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
1. Edit the ```deps.txt``` file, to set up the dependencies. Lines starting with ```#``` are comments only, read them to understand the syntax. For each software dependencies, you need to set up the path to the executable/folder.
   1. __This file has to be edited only once, if you are not familiar with Linux filesystem and command line environment, you may ask someone with more experience for help.__
   1. Don't change the first field in the line, those are the internal identifiers for the software dependencies
   1. The second fields are defining the amount of CPU core and/or memory to be used when running the executable. Please read the comment lines to get the information. The amount set in this file are the optimal CPU core / memory amounts, usually only a little performace increase can be achieved by increasing them (except the metaspades.py, if you have more than 32, increase it). If you have less CPU cores in your Linux server, decrease the amount to that value where needed
   1. If the executable path is in your PATH variable (you can execute them by only calling the application name even from the pipeline root folder), you don't need to edit the 3rd field in the line
   1. If the executable path is not in your PATH variable, or you want to use a specific version for the pipeline (e.g. if you download a new version for SPAdes), you can include the full application path in the 3rd field of the line. Simply replace the application name with its full path (e.g. ```metaspades.py``` -> ```/home/user/SPAdes-3.13.0-Linux/bin/metaspades.py```)
   1. Please be aware that in a few cases the line defines a path to a folder, not to the executable, please read the comment lines carefully!
1. If you move the ```ham_art.pl``` main script file to a different location, always make sure, that you move it together with the ```deps.txt``` file, these two files has to be in the same directory.
1. Download and index host genome(s). The HAM-ART pipeline cleans up the metagenomics sequencing data by removing those reads that are aligning to the host genome (e.g. if you are working with a human faecal sample, the potentially sequenced gut epithelial cell DNA). Make sure that you have the host organism whole genome sequence downloaded, and put the indexed (bowtie2 indexed) files in the ```host_ref_seqs``` folder. For easier use, use simple prefix for the indexed files, you will need to use that prefix later in the project definition file (for an example workflow, see the example project below). 
# Running the HAM-ART pipeline
The HAM-ART pipeline has its own project management system to make the pipeline steps and analysing multiple samples easy. For the analysis, you have to define a project identifier for each sample, define the input files, the project category, and the host organism.
## Setting up a project category
The project category in the HAM-ART pipeline is usually a well identifiable name for a research study that involves multiple samples. During the pipeline run, all raw data, temporary files and final results will be stored in project category specific folders within the pipeline sub-folders. It is also very important to note, that the clade-refinement step ot the pipeline is always performed within one project category. To set up a project category, you simply create a sub-folder in the ```raw_seq_data``` folder, and copy your raw sequencing data files there. Please do not use whitespace (space, tabulator) in the folder and file names, you can use underscore instead (for an example of project category setup, see the example workfolw below).


