:mushroom::mushroom::mushroom: ATTENTION PLEASE :mushroom::mushroom::mushroom:

Read this document carefully to set up the `config.yml` for your workflow!

Plese notice that the indentation in `config.yml` is important.

## Resources settings

```
# Set memory and threads for high demanding rules
high:
  mem_mb: 45000 # memory in MB
  t: 32 # number of threads
  partition: "himem"

# Set memory and threads for medium demanding rules
medium_high:
  mem_mb: 16000 # memory in MB
  t: 16 # number of threads
  partition: "long" # partition to use for medium memory jobs

# Set memory and threads for medium demanding rules
medium:
  mem_mb: 10000 # memory in MB
  t: 8 # number of threads
  partition: "medium" # partition to use for medium memory jobs

# Set memory and threads for low demanding rules
low:
  mem_mb: 4000 # memory in MB
  t: 4 # number of threads
  partition: "short" # partition to use for low memory jobs

# Set memory and threads for very low demanding rules
very_low:
  mem_mb: 500 # memory in MB
  t: 4 # number of threads
  partition: "short" # partition to use for low memory jobs
```

This section of the `config.yml` allows to set memory, threads, and partition for the tools included in the snakemake workflow. The rules are divided in high, medium, and low based on empirical observations on memory requirements. Depending on the samples characteristics (e.g. number of reads, genome size) and HPC requirements, you may need to adjust these values. Note that this part of the config is not relevant if you are running Snakemake locally, in which case you should always use the `--cores` flag in the initial command, which allocates a number of cores to the whole workflow.

| Rule                 | High (memory / CPUs / partition )                           | Medium-High (memory / CPUs / partition )                    | Medium (memory / CPUs / partition )                         | Low (memory / CPUs / partition )                            | Very Low (memory / CPUs / partition )                       |
| -------------------- | ----------------------------------------------------------- | ----------------------------------------------------------- | ----------------------------------------------------------- | ----------------------------------------------------------- | ----------------------------------------------------------- |
| kmergenie            |                                                             |                                                             | :heavy_check_mark: / :heavy_check_mark: /:heavy_check_mark: |                                                             |                                                             |
| seqkit               |                                                             |                                                             |                                                             |                                                             | :heavy_check_mark: / :heavy_check_mark: /:heavy_check_mark: |
| decompress           |                                                             |                                                             |                                                             |                                                             | :heavy_check_mark: / 1 CPU /:heavy_check_mark:              |
| masurca_config       |                                                             |                                                             |                                                             |                                                             | :heavy_check_mark: / 1 CPU /:heavy_check_mark:              |
| fastk                |                                                             |                                                             |                                                             | :heavy_check_mark: / :heavy_check_mark: /:heavy_check_mark: |                                                             |
| spades               | :heavy_check_mark: / :heavy_check_mark: /:heavy_check_mark: |                                                             |                                                             |                                                             |                                                             |
| megahit              |                                                             | :heavy_check_mark: / :heavy_check_mark: /:heavy_check_mark: |                                                             |                                                             |                                                             |
| abyss                | :heavy_check_mark: / :x: / :heavy_check_mark:               |                                                             | :x: / :heavy_check_mark: /:x:                               |                                                             |                                                             |
| masurca              | :heavy_check_mark: / :heavy_check_mark: /:heavy_check_mark: |                                                             |                                                             |                                                             |                                                             |
| minia                |                                                             |                                                             | :heavy_check_mark: / :heavy_check_mark: /:heavy_check_mark: |                                                             |                                                             |
| sparseassembler      |                                                             | :x:/ 1CPU / :heavy_check_mark:                              | :heavy_check_mark: / 1 CPU /:x:                             |                                                             |                                                             |
| get_busco_db         |                                                             |                                                             |                                                             |                                                             | :heavy_check_mark: / 1 CPU /:heavy_check_mark:              |
| busco                |                                                             |                                                             | :heavy_check_mark: / :heavy_check_mark: /:heavy_check_mark: |                                                             |                                                             |
| quast                |                                                             |                                                             |                                                             | :heavy_check_mark: / :heavy_check_mark: /:heavy_check_mark: |                                                             |
| merquryfk            |                                                             |                                                             |                                                             | :heavy_check_mark: / :heavy_check_mark: /:heavy_check_mark: |                                                             |
| select_best_assembly |                                                             |                                                             |                                                             |                                                             | :heavy_check_mark: / 1 CPU /:heavy_check_mark:              |
| bwa                  |                                                             | :heavy_check_mark: / :heavy_check_mark: /:heavy_check_mark: |                                                             |                                                             |                                                             |
| pilon                | :heavy_check_mark: / :x: / :heavy_check_mark:               |                                                             | :x: / :heavy_check_mark: /:x:                               |                                                             |                                                             |
| coverage_viz         |                                                             |                                                             |                                                             |                                                             | :heavy_check_mark: / 1 CPU /:heavy_check_mark:              |

Note that the number of threads allocated in this section of the config file is only used when the tool allows multithreading. If not, the number of threads is set to 1 internally (shown in the table as "1 CPU").

Also, these are the starting values. If a rule does not complete successfully due to insufficient memory allocation, the workflow will automatically attempt again increasing the allocated memory and threads for a maximum of three attempts:

- first attempt: memory set in config; threads set in config
- second attempt: memory set in config x 2; threads set in config + 4
- third attempt: memory set in config x 3; threads set in config + 8

Memory caps internally at 500GB, threads at 128.

The values shown here in the `config/README.md` are suitable values for fungal samples. If dealing with organisms that typically have larger genomes (e.g. plants), or smaller genomes (e.g. bacteria), these values should be changed for a more efficient use of resources on HPC.

## Input and output directories

```
# Set the path to the input directory (this directory must contain
# one subdirectory for each sample, and the subdirectories must be named after the sample)
input_dir: "/home/lia/git_repos/FSP_assembly_benchmarking/test_data/" # use absolute path and add trailing slash

# Customise the results directory
output_dir: "/home/lia/git_repos/FSP_assembly_benchmarking/results/" # use absolute path and add trailing slash
```

This is simple as it is. You have to set up these directories in your `config.yml`. You only need to be careful about two things:

1. make sure the input directory follows the structure described in [Your input data directory structure](../README.md#1-your-input-data-directory-structure) in the main README.md
2. make sure to add the trailing slash to your paths

## K-mer strategy

```
# K-mer size determination strategy
# Note that if you set this to "manual", you must also set the k-mer size for spades, megahit, abyss, sparseassembler, minia, and masurca.
# If you set this to "kmergenie", kmergenie will be used to determine the optimal k-mer size and the values defined below for each rule will be ignored.
# If you set this to "reads_length", seqkit will be used to calculate the median length of the reads, and the kmer will be set to 2/3rds of this value. Note that the k-mer values defined below for each rule will be ignored.
# For spades when kmer_strategy is set to "kmergenie" or "reads_length" the values found by these two methods will be added to this list: 21, 33, 55,77.
# For megahit when kmer_strategy is set to "kmergenie" or "reads_length" the values found by these two methods will be added to this list:21, 29, 39, 59, 79, 99, 119, 141
# For the other assemblers only the predicted k-mer will be used when kmer_strategy is set to "kmergenie" or "reads_length".

kmer_strategy:
  mode: "manual" # "kmergenie" for KmerGenie, "reads_length" to calculate the k-mer size based on the median reads length, "manual" for user-defined k-mer sizes.

# Customise kmergenie parameters
kmergenie:
  "k": 121 # largest k-mer size to consider (default: 121)
  "l": 15 # smallest k-mer size to consider (default: 15)
  optional_params: {} # if you do not want to use any optional parameters, add `{}`
#    "--diploid": True
```

As explained in the comments of the `config.yml` itself, this section allows the user to select between three different strategies to determine the k-mer size used for the genome assembly process:

1. `mode: "manual"`, allows the user to set the k-mer size manually for each assembler. This means that when this strategy is chosen, the user must set the k-mer size in each assembler's section of the `config.yml`.
2. `mode: "kmergenie"`, when this is set the workflow uses kmergenie to estimate an ideal k-mer size based on library characteristics. Note that this tool must be used only when dealing with haploid or diploid samples.
3. `mode: reads_length`, uses Seqkit to calculate the median read length for each library and the k-mer size is set to be 2/3rds of this value.

### Notes for kmergenie usage

As explained above, this tool should be used only in case of haploid or diploid samples. When using dipoids, the flag `--diploid` must be enabled.

### Other notes for k-mer strategy

The first time that the workflow is run using `mode:kmergenie` or `mode:reads_length`, the kmer size for the assemblers will be considered a "to be determined (TBD)" parameter by Snakemake. The second time that exactly the same workflow is executed using the same k-mer strategy, the k-mer size is changed from TBD to the estimated value. This will cause the trigger of all of the downstream rules, even if only a parameter for let's say QUAST was changed in the second run. In this case, to avoid re-running the whole workflow, please add `--rerun-triggers mtime` to your Snakemake command. This is recommended only if you are sure that no other changes were made to the workflow. In any case, it's recommended to run a snakemake command using the flag `--dry-run` to check which rules will be triggered again before launching the real command.

## Assemblers selection

```
# Choose which assemblers to run
assemblers:
  spades: True
  megahit: True
  abyss: True
  sparseassembler: True
  minia: True
  masurca: True
```

This section of the workflow allows the user which assemblers to use in their workflow, simply setting True or False for each of them. This feature is particularly useful when a workflow does not complete due to the failure of a single assembler. In such cases, the user can start again the workflow changing settings in this section, and re-start the workflow from where it stopped.

## Spades

```
# Customise spades parameters
# please do not use the `-o` option here, it is set into the workflow
spades:
  "k": "auto" # list of k-mer sizes (must be odd and less than 128) [default: 'auto']
  optional_params: # if you do not want to use any optional parameters, add `{}` here
    "--only-assembler": True # use this to skip the read error correction step (reads are already cleaned)
```

This is the section of the `config.yml` where the user can set parameters for spades assembler. The `-o` flag is set internally in the workflow and the user must not use this in `optional_params`. Note that when using `kmergenie` or `reads_length` k-mer strategies the `k` parameter here is ignored. For further information about using spades parameters, refer to the [official documentation](https://ablab.github.io/spades/index.html).

## Megahit

```
# Customise megahit parameters
megahit:
  "k": "31,51,71,91,99" #comma-separated list of kmer size all must be odd, in the range 15-255, increment <= 28 [21,29,39,59,79,99,119,141]
  optional_params: # if you do not want to use any optional parameters, add `{}` here
    "--no-mercy": True
    "--min-count": 3 # suggested value for single genome assembly
```

This is the section of the `config.yml` where the user can set parameters for megahit assembler. The `-o` flag is set internally in the workflow and the user must not use this in `optional_params`. Note that when using `kmergenie` or `reads_length` k-mer strategies the `k` parameter here is ignored. The optional parameters `--no-mercy` and `--min-count 3` are suggested for single genome assembly with coverage above 30x, and should not be used for heavily contaminated samples. For further information about using megahit parameters, refer to the [official documentation](https://github.com/voutcn/megahit/wiki).

## Abyss

```
# Customise abyss parameters
# please do not use the `-j` and `-C` options here, they are set into the workflow
abyss:
  "k": 21
  "B": "8G" # note that this must be lower than "mem_mb", otherwise abyss fails. See docs at https://github.com/bcgsc/abyss?tab=readme-ov-file#bloom-filter-mode
  optional_params: # if you do not want to use any optional parameters, add `{}` here
    "kc": 2
```

This is the section of the `config.yml` where the user can set parameters for abyss assembler. The `-C` and `name` flags are set internally in the workflow and the user must not use these in `optional_params`. Note that when using `kmergenie` or `reads_length` k-mer strategies the `k` parameter here is ignored. For further information about using abyss parameters, refer to the [official documentation](https://github.com/bcgsc/abyss). Please, pay attention to the `B` parameter and ensure it is always below the memory allocated to high demanding rules in [this part](#resources-settings) of the `config.yml`. Setting inadequate values for this flag can trigger OOM errors. The official documentation provides guidelines about how to set this parameter according to the genome size of the organisms to analyse.

## SparseAssembler

```
# Customise sparseassembler parameters
sparseassembler:
  "k": 21
  "GS": 50000000 # genome size estimation in bp (used for memory pre-allocation), suggest a large value if possible.(e.g. ~ 3x genome size)
  "Scaffold": 1 # for scaffolding with paired-end reads, set to 1
  "ExpCov": 10 # expected coverage (for scaffolding)
  optional_params: # if you do not want to use any optional parameters, add `{}` here
    "g": 10
    "LD": 0
    "NodeCovTh": 1
    "EdgeCovTh": 0
```

This is the section of the `config.yml` where the user can set parameters for sparseassembler assembler. Note that when using `kmergenie` or `reads_length` k-mer strategies the `k` parameter here is ignored. For further information about using sparseassembler parameters, refer to the [official documentation](https://github.com/yechengxi/SparseAssembler). We believe this tools is no longer maintained.

## Minia

```
# Customise minia parameters
minia:
  "k": 21
  optional_params: {} # if you do not want to use any optional parameters, add `{}` here
#    "-histo2D": 1 # enable 2D histogram output
```

This is the section of the `config.yml` where the user can set parameters for minia assembler. The `-out` flag is set internally in the workflow and the user must not use this in `optional_params`. Note that when using `kmergenie` or `reads_length` k-mer strategies the `k` parameter here is ignored. For further information about using minia parameters, refer to the [official documentation](https://github.com/GATB/minia).

## Masurca

```
# Customise masurca parameters
masurca:
  fragment_mean: 500 # mean insert size for your PE library (adjust to your data)
  fragment_stdev: 50 # stdev for your PE library (adjust to your data)
  k: auto # values between 25 and 127 are supported, auto will compute the optimal size based on the read data and GC content
  jf_size: 100000000 # this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*20
  ca_parameters: cgwErrorRate=0.15 # set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
# Note that the following parameters are always set as follows:
# USE_GRID=0
# CLOSE_GAPS=1
# SOAP_ASSEMBLY=0
# MEGA_READS_ONE_PASS=0
# other parameters for long-reads, jump reads, mate-pairs, and slurm related options (customisable from masurca config) are ignored (for slurm, snakemake will take care of it).
```

This is the section of the `config.yml` where the user can set parameters for masurca assembler. The comments in the `config.yml` are already explaining which parameters are set internally in the workflow. Since these values are going to be used for setting up the config file for masurca, here it's not possible to set any optional parameter. Note that when using `kmergenie` or `reads_length` k-mer strategies the `k` parameter here is ignored. For further information about using masurca parameters, refer to the [official documentation](https://github.com/alekseyzimin/masurca).

## BUSCO

```
# Customise busco parameters
busco:
  taxonomy_file: "/home/lia/git_repos/FSP_assembly_benchmarking/resources/taxonomy.csv" # path to the file containing the taxonomy information for the samples to process (see Readme.md for instructions on how to format it)
  database_list: "/home/lia/git_repos/FSP_assembly_benchmarking/resources/fungi_busco_lineages.txt" # path to the file containing the list of downloaded busco databases (see Readme.md for instructions on how to create it)
  path_to_busco_bds: "/home/lia/git_repos/FSP_assembly_benchmarking/resources/lineages/" # path to the directory containing the busco databases. Please add trailing slash (see Readme.md for instructions on how to download them)
  lineage_general: "fungi_odb12" # chose a "general" dataset, in our case "fungi_odb12"
  extension: "_odb12" # extension of the busco datasets, in our case "_odb12"
  optional_params:
    "--mode": "genome" # always keep this flag
    #    "--metaeuk": True
    "--offline": True # always keep this flag
```

This section of the `config.yml` allows the user to set the parameters for BUSCO. The file contaning the information about the taxonomy of the samples must be created by the user as explained [here](../README.md#3-prepare-the-taxonomy-file) and its path provided to the workflow as input for `taxonomy_file`. Before running the workflow, the user must download all the lineages that could be useful for analysing their samples as explained [here](../README.md#2-download-busco-databases). The input for `database_list` is the path to a text file containing all the downloaded busco lineages, that can be obtained as explained [here](../README.md#2-download-busco-databases). `path_to_busco_dbs` refers to the directory containing the downloaded busco lineages. Please note that if you use [this method](../README.md#2-download-busco-databases) normally a subdirectory named `lineages` is automatically created within the directory where the download command is executed. Do not forget the trailing slash in `path_to_busco_dbs`.
This workflow runs BUSCO on each produced assembly twice using two different databases. One database more "general" (e.g. fungi_odb), and its must be set `config.yml` busco section under `lineage_general` parameter. The second database is automatically selected by the workflow based on the samples' taxonomy to be as closely related as possible to the organisms analysed (e.g. agaricales_odb). The `extension` parameter has been added to allow compatibility with other versions of busco databases, as they are periodically updated (examples: `_odb10`, `_odb12`, etc.).

The flags set internally that the user must not set from the `config.yml` are the following: `--out_path`, `-f`.

The methods used here are suitable for organisms different from fungi, including higher taxonomic ranks (e.g. eukarya, viridiplantae, bacteria). If you whish to try this workflow with other organisms, just change the relevant parts of the code [here](../README.md#2-download-busco-databases) to download the right databases to analyse your samples.

For further information about using BUSCO parameters, refer to the [official documentation](https://busco.ezlab.org/).

## QUAST

```
# Customise quast parameters
quast:
  optional_params: {} # if you do not want to use any optional parameters, add `{}` here
#    "--circos": True # produces a Circos plot
```

This section of the `config.yml` allows the user to set the parameters for QUAST. The `-o` flag is set internally in the workflow and the user must not use this in `optional_params`. We do not recommend using a reference genome for QUAST analysis, unless all the samples can be compared with the same reference genome, as the user is only able to set this parameter once from the `config.yml` adding the `-r` flag in the optional parameters.

Note that QUAST in this workflow is always run with the flag `--min-contig 250`, which is set internally. This means that the user should not use this flag, and that the statstics produced by QUAST do not take into account contigs shorter than 250 bp.

For further information about using QUAST parameters, refer to the [official documentation](https://quast.sourceforge.net/docs/manual.html#sec2.3).

## FastK

```
# Customise fastK parameters
fastk:
  k: 17 # k-mer size, adjust to your data
  t: 1 # include kmer frequency >= t in the output
```

This section of the `config.yml` allows the user to set the parameters for FastK. To calculate the optimal k-mer size for this step, we recommend using [`best_k.sh`](https://github.com/marbl/merqury/blob/master/best_k.sh) from the original [Merqury](https://github.com/marbl/merqury) suite. For fungal genomes, k-mers between 17 and 21 are usually fine (calculated on typical fungal genome size).

For further information about FastK, refer to the [official documentation](https://github.com/thegenemyers/FASTK).

## MerquryFK

```
# Customise merquryfk parameters
merquryfk:
  optional_params: # if you do not want to use any optional parameters, add `{}` here
    "-z": True # outputs assembly plots in addition to spectra copy number
#    "-pdf": True # outputs results in PDF format (default is png)
```

This section of the `config.yml` allows the user to set the parameters for MerquryFK. The `-lfs` flag is set internally in the workflow and the user must not use this in `optional_params`. The user can play with other visualisation settings (e.g. X and Y axes) through `optional_params`.

For further information about using MerquryFK parameters, refer to the [official documentation](https://github.com/thegenemyers/MERQURY.FK).

## Pilon

```
pilon:
  optional_params: # if you do not want to use any optional parameters, add `{}` here
    "--K": 47 # Kmer size used by internal assembler (default 47).
    "--fix": "all" # A comma-separated list of categories of issues to try to fix: "all": all of the above (default);
    "--changes": True # Output a file with the changes made to the assembly.
```

This section of the `config.yml` allows the user to set the parameters for Pilon. In our testing we noticed that using `--fix snps,indels,gaps` rather than `--fix all` significantly speeds up Pilon.

For further information about using Pilon parameters, refer to the [official documentation](https://github.com/broadinstitute/pilon/wiki).
