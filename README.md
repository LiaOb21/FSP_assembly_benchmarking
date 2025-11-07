# Snakemake workflow: `FSP_assembly_benchmarking`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/snakemake-workflows/snakemake-workflow-template/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/snakemake-workflows/snakemake-workflow-template/actions/workflows/main.yml)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/LiaOb21/FSP_assembly_benchmarking)
[![DOI](https://zenodo.org/badge/998768660.svg)](https://doi.org/10.5281/zenodo.17454052)

A Snakemake workflow for *de novo* genome assembly using Illumina reads. 

This workflow was developed for the [Fungarium Sequencing Project (FSP)](https://www.kew.org/science/our-science/projects/sequencing-kews-fungarium) at Royal Botanic Gardens, Kew. 

<img width="1920" height="1080" alt="FSP_assembly_benchmarking_DAG" src="https://github.com/user-attachments/assets/b908d9c6-d2f0-44ab-8f05-fb689d156374" />


The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/LiaOb21/FSP_assembly_benchmarking).

Detailed information about workflow configuration can also be found in the [`config/README.md`](config/README.md).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

- [Snakemake workflow: `FSP_assembly_benchmarking`](#snakemake-workflow-fsp_assembly_benchmarking)
  - [Overview](#overview)
  - [Usage](#usage)
    - [Clone the repo](#clone-the-repo)
    - [Prepare the required inputs](#prepare-the-required-inputs)
      - [1. Your input data directory structure](#1-your-input-data-directory-structure)
      - [2. Download BUSCO databases](#2-download-busco-databases)
      - [3. Prepare the taxonomy file](#3-prepare-the-taxonomy-file)
    - [Running the workflow locally](#running-the-workflow-locally)
      - [1. Set up Snakemake environment](#1-set-up-snakemake-environment)
      - [2. Run the workflow](#2-run-the-workflow)
    - [Running the workflow on the cluster](#running-the-workflow-on-the-cluster)
      - [1. Set up Snakemake environment](#1-set-up-snakemake-environment-1)
      - [2. Create submission script](#2-create-submission-script)
      - [3. Submit and monitor](#3-submit-and-monitor)
  - [Deployment options (needs update)](#deployment-options-needs-update)
  - [Authors](#authors)
  - [References (needs update)](#references-needs-update)
  - [TODO](#todo)

## Overview


The workflow was developed to benchmark the performance of different short reads assemblers using hDNA (historical DNA).

It may be useful to other projects that deal with difficult samples as the FSP, and need to find the best short reads assembler(s) for their own case.

The workflow was designed to process several samples in parallel. 

Each sample is assembled using the following assemblers:

- [SPAdes](https://github.com/ablab/spades) - multi k-mer assembler
- [MEGAHIT](https://github.com/voutcn/megahit) - multi k-mer assembler
- [AbySS](https://github.com/bcgsc/abyss) - single k-mer assembler
- [SparseAssembler](https://github.com/yechengxi/SparseAssembler) - single k-mer assembler
- [Minia](https://github.com/GATB/minia?tab=readme-ov-file) - single k-mer assembler
- [MaSuRCA](https://github.com/alekseyzimin/masurca/tree/master) - single k-mer assembler

Note: the workflow was designed to use previously pre-processed reads, although it works also with raw reads (but keep in mind that this doesn't pre-process reads automatically).

The user can choose between three different strategies to set the k-mer size for the genome assembly process:
1) `manual`: the k-mer size is set manually for each assembler in [`config/README.md`](config/README.md).
2) `kmergenie`: [KMerGenie](https://github.com/movingpictures83/KMerGenie) is used to estimate the best k-mer size for genome assembly for each library. Note that KMerGenie should be used with haploid or diploid genomes only.
3) `reads_length`: [SeqKit](https://bioinf.shenwei.me/seqkit/) is used to calculate the median reads length and the k-mer size is set to be 2/3rds of the median.

The assemblies produced by each assembler for each sample are then quality inspected with the following tools:

- [BUSCO](https://busco.ezlab.org/busco_userguide.html) - evaluates each produced assembly quality in terms of expected gene content. It is run twice for each sample, once using a general dataset, and once using a more closely related dataset. See [`config/README.md`](config/README.md) for more details.
- [QUAST](https://quast.sourceforge.net/docs/manual.html) - computes assembly statistics.
- [MerquryFK](https://github.com/thegenemyers/MERQURY.FK) - computes k-mer analysis for each assembly and compares its content with the k-mer computed for raw reads by [FastK](https://github.com/thegenemyers/FASTK).


The best assembly for each sample is then selected. The selection is based on the highest complete and single BUSCO content (absolute number, not percentage). If more than one assembly have the highest complete and single copy BUSCO score, the assembly with the highest aUN (calculated by QUAST) among those is selected.

The best assembly is then aligned to the reads using [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) and [samtools](https://www.htslib.org/doc/samtools.html) and handed over to [Pilon](https://github.com/broadinstitute/pilon), to improve the assembly by performing error correction and gap filling.

The improved best assembly is then quality assessed again with BUSCO, QUAST, and MerquryFK, and aligned back to the reads with bwa-mem2 and samtools, which is also used to analyse the genome coverage.

## Usage

### Clone the repo

```
git clone https://github.com/LiaOb21/FSP_assembly_benchmarking.git
```

### Prepare the required inputs

#### 1. Your input data directory structure

The directory that contains your input samples (e.g. `data`) must be structured in the following way:
```
data/
├── 048ds
│   ├── 048ds_merge.fq.gz
│   ├── 048ds_trimmed.R1.fq.gz
│   ├── 048ds_trimmed.R2.fq.gz
│   ├── 048ds_unmerged.R1.fq.gz
│   └── 048ds_unmerged.R2.fq.gz
└── 048ss
    ├── 048ss_merge.fq.gz
    ├── 048ss_trimmed.R1.fq.gz
    ├── 048ss_trimmed.R2.fq.gz
    ├── 048ss_unmerged.R1.fq.gz
    └── 048ss_unmerged.R2.fq.gz
```

In this example `048ds` and `048ss` are the only two samples present in the input directory. Note how each subdirectory is named after the sample, and the files inside each sample subdirectory have a standardised name. This is crucial for Snakemake to work properly.

#### 2. Download BUSCO databases

This workflow runs BUSCO on each produced assembly twice using two different databases. One database more "general" (e.g. fungi_odb), and one as closely related as possible to the organisms analysed (e.g. agaricales_odb).

The recommendation here is to store all the BUSCO databases in a single directory, to allow the workflow to automatically use the closest database to the organism for the second busco run.

First of all, obtain the full list of busco lineages and use the awk command to extract the lineages of the group of interest (e.g. fungi). The `_odb12` suffix refers to the version of the busco databases. The following commands can be used also with other groups of organisms or database versions, modifying the relevant parts of the code.

```
cd to/where/you/want/to/store/busco/databases

conda activate busco 

busco --list > busco_lineages.txt

awk '/fungi_odb12/{flag=1; indent=length($0)-length(ltrim($0)); print "fungi_odb12"; next} 
     flag && /- [a-z_]*_odb12/ {
         current_indent=length($0)-length(ltrim($0))
         if(current_indent <= indent) flag=0
         else print gensub(/.*- ([a-z_]*_odb12).*/, "\\1", "g")
     }
     function ltrim(s) { sub(/^[ \t\r\n]+/, "", s); return s }' busco_lineages.txt > fungi_busco_lineages.txt
```

`fungi_busco_lineages.txt` will be the input for this field in the config file:
```
busco:
  database_list: "to/where/you/want/to/store/busco/databases/fungi_busco_lineages.txt"
```

You will need to download all the databases you need for your analysis before hands. You can easily do this in this way:

```
for i in $(cat fungi_busco_lineages.txt); do
  echo "downloading $i database"
  busco --download_path . --download $i
done
```

As you will download dabases before starting the workflow, you should set `"--offline": True` in the `config.yaml` for the BUSCO parameters. For more details, please refer to [`config/README.md`](config/README.md).

#### 3. Prepare the taxonomy file

We need a tab separated taxonomy file containing at least the following columns: 'Sample', 'Family', 'Order', 'Class', 'Phylum'. Name the columns exactly as shown here (with capital letters).

Example:

<img width="901" height="267" alt="image" src="https://github.com/user-attachments/assets/fc163f59-b988-4e30-a0b0-251b8df171ca" />

Here the column 'Genus' will be ignored (as BUSCO databases are available up to family level), but it's okay if it is in the file.

### Running the workflow locally

:mushroom::mushroom::mushroom: **Note: before running the workflow you should set up your `config.yml`. You can find it in [`config/config.yml`](config/config.yml). The [`config/README.md`](config/README.md) explains how to set up the `config.yml`.** :mushroom::mushroom::mushroom:

#### 1. Set up Snakemake environment

First install snakemake using conda (be sure to install version 9+):
```
conda create -n snakemake snakemake
conda activate snakemake
```

If you want to monitor the workflow through a console you can install the [`snkmt` plugin](https://github.com/cademirch/snkmt) within the snakemake environment (this is not mandatory):
```
pip install snakemake-logger-plugin-snkmt # used for monitoring resources
```
#### 2. Run the workflow

:mushroom::mushroom::mushroom: **IMPORTANT:** Snakemake workflows must be always run from the cloned directory. :mushroom::mushroom::mushroom:

To run the workflow, you can use the following command for using the regular inputs reads forward and reverse (R1 and R2):
```
cd ~/FSP_assembly_benchmarking
snakemake --configfile config/config.yml --software-deployment-method conda --snakefile workflow/Snakefile --cores 8
```

If you want to use `snkmt` add this flag to the previous command: `--logger snkmt`.

When you use `workflow/Snakefile`, the only files you need in [your input data directory](README.md#1-your-input-data-directory-structure) are `*_trimmed.R1.fq.gz` and `*_trimmed.R2.fq.gz`.

If you want to use merged reads (that you have to merge previously) you can run the following command:

```
cd ~/FSP_assembly_benchmarking
snakemake --configfile config/config.yml --software-deployment-method conda --snakefile workflow/Snakefile_merged --cores 8
```

Note that for the workflow to run properly when using `workflow/Snakefile_merged` the input directory must contain all the files showed in [Your input data directory structure](README.md#1-your-input-data-directory-structure), including `*_unmerged.R1.fq.gz`, `*_unmerged.R2.fq.gz`, `*_trimmed.R1.fq.gz`, and `*_trimmed.R2.fq.gz`, as some of the assemblers needs these file in order to run, and the reads alignment is always performed using `*_trimmed.R1.fq.gz` and `*_trimmed.R2.fq.gz`.

Even in this case you can add the `--logger snkmt` flag.

Be sure to always name the files as shown in `1. Your data directory structure` section.


### Running the workflow on the cluster

#### 1. Set up Snakemake environment

```
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
pip install snakemake-executor-plugin-cluster-generic # to handle SLURM job submission
pip install snakemake-logger-plugin-snkmt # used for monitoring resources (optional)
```

#### 2. Create submission script

:mushroom::mushroom::mushroom: **IMPORTANT:** Snakemake workflows must be always run from the cloned directory (i.e. `FSP_assembly_benchmarking`). Therefore, the following script must be placed in that directory. :mushroom::mushroom::mushroom:

Create `run_snakemake.sh`:

```
#!/bin/bash
#SBATCH --job-name=snakemake-fsp
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=snakemake-%j.out
#SBATCH --error=snakemake-%j.err

# Create logs directory
mkdir -p logs

# Load conda environment
source ~/.bashrc
conda activate snakemake

# Run Snakemake (individual jobs will be distributed to appropriate partitions)
echo "Starting Snakemake workflow at $(date)"
snakemake --profile profile/ --logger snkmt
echo "Workflow completed at $(date)"
```

With this sbatch script example, Snakemake will look for the default snakefil, which is `workflow/Snakefile`, and will run the workflow using that file. The only files you need in [your input data directory](README.md#1-your-input-data-directory-structure) in this case are `*_trimmed.R1.fq.gz` and `*_trimmed.R2.fq.gz`.


If you want to use (previously) merged reads, you should modify the above script to use `workflow/Snakefile_merged` instead:

```
#!/bin/bash
#SBATCH --job-name=snakemake-fsp
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=snakemake-%j.out
#SBATCH --error=snakemake-%j.err

# Create logs directory
mkdir -p logs

# Load conda environment
source ~/.bashrc
conda activate snakemake

# Run Snakemake (individual jobs will be distributed to appropriate partitions)
echo "Starting Snakemake workflow at $(date)"
snakemake --snakefile workflow/Snakefile_merged --profile profile/ --logger snkmt
echo "Workflow completed at $(date)"
```

Note that for the workflow to run properly when using `workflow/Snakefile_merged` the input directory must contain all the files showed in [Your input data directory structure](README.md#1-your-input-data-directory-structure), including `*_unmerged.R1.fq.gz`, `*_unmerged.R2.fq.gz`, `*_trimmed.R1.fq.gz`, and `*_trimmed.R2.fq.gz`, as some of the assemblers needs these file in order to run, and the reads alignment is always performed using `*_trimmed.R1.fq.gz` and `*_trimmed.R2.fq.gz`.

:mushroom::mushroom::mushroom: **IMPORTANT**: These are examples of sbatch script. Make sure to set up your sbatch script according to your system settings. :mushroom::mushroom::mushroom:

In alternative, you can launch Snakemake in a `screen` shell using `srun` for an interactive run.

#### 3. Submit and monitor

Submit the workflow throgh the sbatch script:

```
sbatch run_snakemake.sh
```

Monitor the workflow:

```
# Monitor using snkmt console
conda activate snakemake
snkmt console

# Monitor with SLURM
squeue -u $USER
squeue -o "%.18i %.100j %.8u %.2t %.10M %.6D %R" -u $USER      # to see full name of jobs including sample name

# Monitor through the logs
less snakemake-JOBID.err      # to have a look to the log
grep -A 20 "Error in rule" snakemake-JOBID.err       # to screen the log for possible errors
```

## Deployment options (needs update)

To run the workflow from command line, change the working directory.

```bash
cd path/to/snakemake-workflow-name
```

Adjust options in the default config file [`config/config.yml`](config/config.yml).
Before running the complete workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

To run the workflow with test files using **conda**:

```bash
snakemake --cores 2 --sdm conda --directory .test
```

To run the workflow with **apptainer** / **singularity**, add a link to a container registry in the `Snakefile`, for example `container: "oras://ghcr.io/<user>/<repository>:<version>"` for Github's container registry.
Run the workflow with:

```bash
snakemake --cores 2 --sdm conda apptainer --directory .test
```

## Authors

- Lia Obinu
  - Royal Botanic Gardens, Kew
  - [ORCID profile](https://orcid.org/0000-0002-3208-323X)

The other members of the FSP bioinformatics team, [Wu Huang](https://github.com/Hazelhuangup) and [Niall Garvey](https://github.com/NiallG1), and [George Mears](https://github.com/George-Mears) also contributed to the development and testing of this part of the workflow (which is only a part of the full pipeline developed for the FSP).

## References (needs update)

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

## TODO




- Update the [deployment](#deployment-options) and [references](#references) sections.
- Update the `README.md` badges. Add or remove badges for `conda`/`singularity`/`apptainer` usage depending on the workflow's [deployment](#deployment-options) options.

