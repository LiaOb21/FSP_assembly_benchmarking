:mushroom::mushroom::mushroom: ATTENTION PLEASE :mushroom::mushroom::mushroom:

Read this document carefully to set up the `config.yml` for your workflow! 

## Resources settings

```
# Set memory and threads for high demanding rules
high:
  mem_mb: 45000 # memory in MB
  t: 16 # number of threads
  partition: "himem" # partition to use for high memory jobs

# Set memory and threads for medium demanding rules
medium:
  mem_mb: 16000 # memory in MB
  t: 16 # number of threads
  partition: "medium" # partition to use for medium memory jobs

# Set memory and threads for low demanding rules
low:
  mem_mb: 5000 # memory in MB
  t: 8 # number of threads
  partition: "short" # partition to use for low memory jobs
```

This section of the `config.yml` allows to set memory, threads, and partition for the tools included in the snakemake workflow. The rules are divided in high, medium, and low based on empirical observations on memory requirements. Depending on the samples characteristics (e.g. number of reads, genome size) and HPC requirements, you may need to adjust these values. Note that this part of the config is not relevant if you are running Snakemake locally, in which case you should always use the `--cores` flag in the initial command, which allocates a number of cores to the whole workflow.

| Rule | High | Medium | Low |
|---|---|---|---|
| kmergenie |  | :heavy_check_mark: | |
| seqkit  | | :heavy_check_mark: | |
| decompress | | | :heavy_check_mark:
| masurca_config | | | :heavy_check_mark:
| fastk  | | :heavy_check_mark: | |
| spades | :heavy_check_mark: | | |
| megahit | | :heavy_check_mark: | |
| abyss | :heavy_check_mark: | | |
| masurca | :heavy_check_mark: | | |
| minia |  | :heavy_check_mark: | |
| sparseassembler |  | :heavy_check_mark: | |
| busco |  | :heavy_check_mark: | |
| quast | | | :heavy_check_mark:
| merquryfk | | | :heavy_check_mark:
| select_best_assembly | | | :heavy_check_mark:
| bwa | :heavy_check_mark: | | |
| pilon | :heavy_check_mark: | | |
| coverage_viz | | | :heavy_check_mark: |

Note that the number of threads allocated in this section of the config file is only used when the tool allows multithreading. If not, the number of threads is set to 1 internally.

Also, these are the starting values. If a rule does not complete successfully due to insufficient memory allocation, the workflow will automatically attempt again increasing the allocated memory and threads for a maximum of three attempts:
- first attempt: memory set in config; threads set in config
- second attempt: memory set in config x 2; threads set in config + 4
- third attempt: memory set in config x 2; threads set in config + 8

Memory caps internally at 500GB, threads at 128.


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
```

As explained in the comments of the `config.yml` itself, this section allows the user to select between three different strategies to determine the k-mer size used for the genome assembly process:
1. `mode: "manual"`, allows the user to set the k-mer size manually for each assembler. This means that when this strategy is chosen, the user must set the k-mer size in each assembler's section of the `config.yml`.
2. `mode: "kmergenie"`, when this is set