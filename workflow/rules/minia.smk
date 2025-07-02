# This rule runs minia assembler with the parameters specified in the config file.
import glob
import os


rule minia:
    input:
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
    output:
        scaffolds=f"{output_dir}" + "{sample}/minia/{sample}.contigs.fa",
        link_assembly=f"{output_dir}" + "assemblies/{sample}/{sample}_minia.fa",
    params:
        k=config["minia"]["k"],
        result_prefix=lambda wildcards, output: os.path.splitext(output.scaffolds)[
            0
        ].replace(".contigs", ""),
        file_list=lambda wildcards, output: os.path.join(
            os.path.dirname(output.scaffolds), f"{wildcards.sample}_files.txt"
        ),
        optional_params=" ".join(
            k for k, v in config["minia"]["optional_params"].items() if v is True
        ),
    threads: config["threads"]  # access threads from config
    log:
        "logs/{sample}/minia.log",
    resources:
        mem_mb=config["mem_mb"],  # access memory from config
    conda:
        "../envs/minia.yaml"
    shell:
        """
        # Create file list for minia
        echo "{input.forward_in}" > {params.file_list}
        echo "{input.reverse_in}" >> {params.file_list}

        minia -in {params.file_list} -kmer-size {params.k} -out {params.result_prefix} -nb-cores {threads} {params.optional_params} >> {log} 2>&1

        ln -srn {output.scaffolds} {output.link_assembly}
        """
