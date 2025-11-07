# This rule runs megahit assembler with the parameters specified in the config file.
import glob
import os


rule megahit_merged:
    input:
        merged_in=f"{input_dir}" + "{sample}/{sample}_merge.fq.gz",
        kmergenie_result=get_kmergenie_dependency,
    output:
        result_dir=directory(f"{output_dir}" + "{sample}/megahit"),
        link_assembly=f"{output_dir}" + "{sample}/assemblies/{sample}_megahit.fa",
    params:
        k=lambda wildcards: get_kmer_list(wildcards, "megahit", "k"),
        scaffolds=f"{output_dir}" + "{sample}/megahit/final.contigs.fa",
        optional_params=" ".join(
            f"{k}" if v is True else f"{k} {v}"
            for k, v in config["megahit"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: get_medium_high_threads
    resources:
        mem_mb=get_medium_high_mem,
        partition=config["medium_high"]["partition"],
    log:
        "logs/{sample}/megahit.log",
    benchmark:
        "benchmark/{sample}/megahit.txt"
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        echo "Running megahit with the following command:" >> {log} 2>&1
        echo "megahit -t {threads} -r {input.merged_in} -o {output.result_dir} --k-list {params.k} {params.optional_params}" >> {log} 2>&1
        megahit -t {threads} -r {input.merged_in} -o {output.result_dir} --k-list {params.k} {params.optional_params} >> {log} 2>&1

        ln -srn {params.scaffolds} {output.link_assembly}
        """
