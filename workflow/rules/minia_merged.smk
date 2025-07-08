# This rule runs minia assembler with the parameters specified in the config file.
import glob
import os


rule minia:
    input:
        merged_in=f"{input_dir}" + "{sample}/{sample}_merge.fq.gz",
    output:
        scaffolds=f"{output_dir}" + "{sample}/minia/{sample}.contigs.fa",
        link_assembly=f"{output_dir}" + "assemblies/{sample}/{sample}_minia.fa",
    params:
        k=config["minia"]["k"],
        result_prefix=lambda wildcards, output: os.path.splitext(output.scaffolds)[
            0
        ].replace(".contigs", ""),
        optional_params=" ".join(
            k for k, v in config["minia"]["optional_params"].items() if v is True
        ),
    threads: get_scaled_threads  # Use scaling function
    log:
        "logs/{sample}/minia.log",
    benchmark:
        "benchmark/{sample}/minia.txt"
    resources:
        mem_mb=get_scaled_mem,  # Use scaling function
    conda:
        "../envs/minia.yaml"
    shell:
        """
        echo "Running minia with the following command:" >> {log} 2>&1
        echo "minia -in {input.merged_in} -kmer-size {params.k} -out {params.result_prefix} -nb-cores {threads} {params.optional_params}" >> {log} 2>&1
        minia -in {input.merged_in} -kmer-size {params.k} -out {params.result_prefix} -nb-cores {threads} {params.optional_params} >> {log} 2>&1

        ln -srn {output.scaffolds} {output.link_assembly}
        """
