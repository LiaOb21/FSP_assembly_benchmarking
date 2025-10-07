# Snakemake workflow: `FSP_assembly_benchmarking`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/snakemake-workflows/snakemake-workflow-template/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/snakemake-workflows/snakemake-workflow-template/actions/workflows/main.yml)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/LiaOb21/FSP_assembly_benchmarking)

A Snakemake workflow for *de novo* genome assembly using Illumina reads. 

This workflow was developed for the Fungarium Sequencing Project (FSP) at Royal Botanic Gardens, Kew. 



- [Snakemake workflow: `FSP_assembly_benchmarking`](#snakemake-workflow-fsp_assembly_benchmarking)
  - [Overview](#overview)
  - [Usage](#usage)
    - [Download the repo](#download-the-repo)
    - [Prepare the required inputs](#prepare-the-required-inputs)
    - [Running the workflow locally:](#running-the-workflow-locally)
    - [Running the workflow on the cluster](#running-the-workflow-on-the-cluster)
  - [Note for kmergenie usage](#note-for-kmergenie-usage)
  - [Note for abyss](#note-fro-abyss)
  - [obtain the list of busco databases](#obtain-the-list-of-busco-databases)
  - [Deployment options](#deployment-options)
  - [Authors](#authors)
  - [References](#references)
  - [TODO](#todo)

## Overview


The workflow was developed to benchmark the performance of different short reads assemblers using aDNA (ancient DNA). 

It may be useful to other projects that deal with difficult samples as the FSP, and need to find the best short reads assembler(s) for their own case.

The workflow was designed to process several samples in parallel. The flow chart above illustrates how the workflow looks like when processing a single sample (example: 048ds).

Each sample is assembled using the following assemblers:

- [SPAdes](https://github.com/ablab/spades)
- [MEGAHIT](https://github.com/voutcn/megahit)
- [AbySS](https://github.com/bcgsc/abyss)
- [SparseAssembler](https://github.com/yechengxi/SparseAssembler)
- [Minia](https://github.com/GATB/minia?tab=readme-ov-file)
- [MaSuRCA](https://github.com/alekseyzimin/masurca/tree/master)


The assemblies produced by each assembler for each sample are then quality inspected with the following tools

- [BUSCO](https://busco.ezlab.org/busco_userguide.html) - evaluates each produced assembly quality in terms of expected gene content. It is run twice for each sample, once using a general dataset, and once using a more closely related dataset. See [`config/README.md`](config/README.md) for more details.
- [QUAST](https://quast.sourceforge.net/docs/manual.html) - computes assembly statistics.
- [MerquryFK](https://github.com/thegenemyers/MERQURY.FK) - computes k-mer analysis for each assembly and compares its content with the k-mer computed for raw reads by [FastK](https://github.com/thegenemyers/FASTK).

## Usage

### Download the repo

```
git clone https://github.com/LiaOb21/FSP_assembly_benchmarking.git
```

### Prepare the required inputs

#### 1. Your data directory structure

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

The reccommendation here is to store all the BUSCO databases in a single directory, to allow the workflow to automatically use the closest database to the organism for the second busco run.

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
  database_list: path/to/fungi_busco_lineages.txt
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

We need a tab separated taxonomy file containing the following columns: 'Sample', 'Family', 'Order', 'Class', 'Phylum'. Name the columns exactly as shown here (with capital letters).

Example:



### Running the workflow locally:

**Note: before running the workflow you should set up your `config.yml`. You can find it in [`config/config.yml`](config/config.yml). The [`config/README.md`](config/README.md) explains how to set up the `config.yml`.**

#### 1. Set up Snakemake environment

First install snakemake using conda (be sure to install version 9+):
```
conda create -n snakemake snakemake
conda activate snakemake
```

If you want to monitor the workflow through a console you can install the `snkmt` plugin within the snakemake environment (this is not mandatory):
```
pip install snakemake-logger-plugin-snkmt # used for monitoring resources
```
#### 2. Run the workflow

To run the workflow, you can use the following command for using the regular inputs reads forward and reverse (R1 and R2):
```
snakemake --configfile config/config.yml --software-deployment-method conda --snakefile workflow/Snakefile --cores 8
```

If you want to use `snkmt` add this flag to the previous command: `--logger snkmt`.

If you want to use merged reads (that you have to merge previously) you can run the following command:

```
snakemake --configfile config/config.yml --software-deployment-method conda --snakefile workflow/Snakefile_merged --cores 8
```

Even in this case you can add the `--logger snkmt` flag.

Be sure to always name the files as shown in `1. Your data directory structure` section.


### Running the workflow on the cluster

#### 1. Set up Snakemake environment

```
conda create -c conda-forge -c bioconda -n snakemake snakemake  #Use the full version of the commands, as otherwise an older version may be downloaded.
conda activate snakemake
pip install snakemake-executor-plugin-cluster-generic # to handle SLURM job submission
pip install snakemake-logger-plugin-snkmt # used for monitoring resources (optional)
```

#### 2. Run the workflow


```
screen -S snakemake_test
conda activate snakemake
snakemake --profile profile/
```


To set up the config file you must at least indicate the following:
```
nano config/config.yml
---
input_dir: "/home/lobinu/test_data/clean_data/00_test" #use absolute paths and do not add `/` at the end
output_dir: "/home/lobinu/scratch/FSP_assembly_benchmarking/results"  #use absolute paths and do not add `/` at the end
busco:
  lineage1: "resources/fungi_odb12"
  lineage2: "resources/basidiomycota_odb12"
```



## Note for kmergenie usage

The first time that the workflow is run using kmergenie, the kmer size for the assemblers will be considered a "to be determined (TBD)" parameter by Snakemake. The second time that exactly the same workflow is executed using kmergenie, the kmer size is changed from TBD to the value estimated by kmergenie. This will cause the trigger of all of the downstream rules, even if only a parameter for let's say QUAST was changed in the second run. In this case, to avoid re-running the whole workflow, please add `--rerun-triggers mtime` to your Snakemake command. This is recommended only if you are sure that no other changes were made to the workflow.

## Note fro abyss

Pay attention to B parameter, that can trigger OOM errors.


## obtain the list of busco databases

This is an attempt to process samples with different taxonomy together, rather than processing them by groups. 

First of all, obtain the full list of busco lineages and use the awk command to extract the lineages of the group of interest (e.g. fungi). The `_odb12` suffix refers to the version of the busco database.

```
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

You will need to download all the database you need for your analysis. You can easily do this in this way:

```
for i in $(cat fungi_busco_lineages.txt); do
  echo "downloading $i database"
  busco --download_path . --download $i
done
```

`fungi_busco_lineages.txt` will be the input for this in the config file:
```
busco:
  database_list: path/to/fungi_busco_lineages.txt
```

To achieve this, we also need a tab separated taxonomy file containing the following columns: 'Sample', 'Family', 'Order', 'Class', 'Phylum'. Name them exactly as shown here (with capital letters).

-----

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/LiaOb21/FSP_assembly_benchmarking).

Detailed information about input data and workflow configuration can also be found in the [`config/README.md`](config/README.md).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

## Deployment options

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

- Firstname Lastname
  - Affiliation
  - ORCID profile
  - home page

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

## TODO

- Replace `<owner>` and `<repo>` everywhere in the template with the correct user name/organization, and the repository name. The workflow will be automatically added to the [snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/index.html) once it is publicly available on Github.
- Replace `<name>` with the workflow name (can be the same as `<repo>`).
- Replace `<description>` with a description of what the workflow does.
- Update the [deployment](#deployment-options), [authors](#authors) and [references](#references) sections.
- Update the `README.md` badges. Add or remove badges for `conda`/`singularity`/`apptainer` usage depending on the workflow's [deployment](#deployment-options) options.
- Do not forget to also adjust the configuration-specific `config/README.md` file.
