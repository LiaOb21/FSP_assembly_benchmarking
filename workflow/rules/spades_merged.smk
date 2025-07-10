# This rule runs spades assembler with the parameters specified in the config file.
import glob
import os


rule spades:
    input:
        merged_in=f"{input_dir}" + "{sample}/{sample}_merge.fq.gz",
        unmerged_r1=f"{input_dir}" + "{sample}/{sample}_unmerged.R1.fq.gz",
        unmerged_r2=f"{input_dir}" + "{sample}/{sample}_unmerged.R2.fq.gz",
    output:
        result_dir=directory(f"{output_dir}" + "{sample}/spades"),
        scaffolds=f"{output_dir}" + "{sample}/spades/scaffolds.fasta",
        link_assembly=f"{output_dir}" + "assemblies/{sample}/{sample}_spades.fa",
    params:
        optional_params=" ".join(
            f"{k} {v}" if v is not True else k
            for k, v in config["spades"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: get_scaled_threads  # Use scaling function
    log:
        "logs/{sample}/spades.log",
    benchmark:
        "benchmark/{sample}/spades.txt"
    resources:
        mem_mb=get_scaled_mem,  # Use scaling function
    conda:
        "../envs/spades.yaml"
    shell:
        """
        echo "Running spades with the following command:" >> {log} 2>&1
        echo "spades.py -t {threads} --merged {input.merged_in} -1 {input.unmerged_r1} -2 {input.unmerged_r2} -o {output.result_dir} {params.optional_params}" >> {log} 2>&1
        spades.py -t {threads} --merged {input.merged_in} -1 {input.unmerged_r1} -2 {input.unmerged_r2} -o {output.result_dir} {params.optional_params} >> {log} 2>&1
        ln -srn {output.scaffolds} {output.link_assembly}
        """
