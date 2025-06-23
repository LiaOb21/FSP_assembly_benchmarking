# This rule runs minia assembler with the parameters specified in the config file.
import glob
import os


rule minia:
    input:
        merged_in= f"{input_dir}" + "/{sample}/{sample}_merge.fq.gz",
    output:
        scaffolds = f"{output_dir}" + "/{sample}/minia/{sample}.contigs.fa",
    params:
        k = config["minia"]["k"],
        result_prefix = f"{output_dir}" + "/{sample}/minia/{sample}",
        optional_params = " ".join(
            k for k, v in config["minia"]["optional_params"].items() if v is True
        ),
    threads: config["threads"],  # access threads from config
    log:
        "logs/{sample}/minia.log",
    resources:
        mem_mb=config["mem_mb"],  # access memory from config
    conda:
        "../envs/minia.yaml"
    shell:
        """
        minia -in {input.merged_in} -kmer-size {params.k} -out {params.result_prefix} -nb-cores {threads} {params.optional_params} >> {log} 2>&1
        """