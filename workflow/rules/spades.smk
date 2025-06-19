# This rule runs spades assembler with the parameters specified in the config file.
import glob
import os


rule spades:
    input:
        forward_in= f"{input_dir}" + "/{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in= f"{input_dir}" + "/{sample}/{sample}_trimmed.R2.fq.gz",
    output:
        result_dir = directory(f"{output_dir}" + "/{sample}"),
        scaffolds = f"{output_dir}" + "/{sample}/scaffolds.fasta",
    params:
        optional_params = " ".join(
            k for k, v in config["spades"]["optional_params"].items() if v is True
        ),
    threads: config["threads"],  # access threads from config
    log:
        "logs/{sample}/spades.log",
    resources:
        mem_mb=config["mem_mb"],  # access memory from config
    conda:
        "../envs/spades.yaml"
    shell:
        """
        spades.py -t {threads} -1 {input.forward_in} -2 {input.reverse_in} -o {output.result_dir} {params.optional_params} >> {log} 2>&1
        """