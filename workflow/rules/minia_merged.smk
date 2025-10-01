# This rule runs minia assembler with the parameters specified in the config file.
import glob
import os


rule minia:
    input:
        merged_in=f"{input_dir}" + "{sample}/{sample}_merge.fq.gz",
        kmergenie_result=get_kmergenie_dependency,
    output:
        scaffolds=f"{output_dir}" + "{sample}/minia/{sample}.contigs.fa",
        link_assembly=f"{output_dir}" + "{sample}/assemblies/{sample}_minia.fa",
    params:
        k=lambda wildcards: get_single_kmer(wildcards, "minia", "k"),
        result_prefix=lambda wildcards, output: os.path.splitext(output.scaffolds)[
            0
        ].replace(".contigs", ""),
        optional_params=" ".join(
            f"{k}" if v is True else f"{k} {v}"
            for k, v in config["minia"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: get_medium_threads
    resources:
        mem_mb=get_medium_mem,
        partition=config["medium"]["partition"],
    log:
        "logs/{sample}/minia.log",
    benchmark:
        "benchmark/{sample}/minia.txt"
    conda:
        "../envs/minia.yaml"
    shell:
        """
        echo "Running minia with the following command:" >> {log} 2>&1
        echo "minia -in {input.merged_in} -kmer-size {params.k} -out {params.result_prefix} -nb-cores {threads} {params.optional_params}" >> {log} 2>&1
        minia -in {input.merged_in} -kmer-size {params.k} -out {params.result_prefix} -nb-cores {threads} {params.optional_params} >> {log} 2>&1

        ln -srn {output.scaffolds} {output.link_assembly}
        """
