:mushroom::mushroom::mushroom: ATTENTION PLEASE :mushroom::mushroom::mushroom:

Read this document carefully to set up the `config.yml` for your workflow! 

### Resources settings

```
# Set memory and threads for high demanding rules
high:
  mem_mb: 45000 # memory in MB
  t: 32 # number of threads
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

This section of the `config.yml` allows to set memory, threads, and partition for the tools included in the snakemake workflow. The rules are divided in high, medium, and low based on empirical observations. 