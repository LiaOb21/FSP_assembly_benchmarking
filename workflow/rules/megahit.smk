# This rule runs megahit assembler with the parameters specified in the config file.
import glob
import os


rule megahit:
    input:
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
    output:
        result_dir=directory(f"{output_dir}" + "{sample}/megahit"),
        link_assembly=f"{output_dir}" + "assemblies/{sample}/{sample}_megahit.fa",
    params:
        scaffolds=f"{output_dir}" + "{sample}/megahit/final.contigs.fa",
        optional_params=" ".join(
            k for k, v in config["megahit"]["optional_params"].items() if v is True
        ),
    threads: get_scaled_threads  # Use scaling function
    log:
        "logs/{sample}/megahit.log",
    benchmark:
        "benchmark/{sample}/megahit.txt"
    resources:
        mem_mb=get_scaled_mem,  # Use scaling function
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        echo "Running megahit with the following command:" >> {log} 2>&1
        echo "megahit -t {threads} -1 {input.forward_in} -2 {input.reverse_in} -o {output.result_dir} {params.optional_params}" >> {log} 2>&1
        megahit -t {threads} -1 {input.forward_in} -2 {input.reverse_in} -o {output.result_dir} {params.optional_params} >> {log} 2>&1

        ln -srn {params.scaffolds} {output.link_assembly}
        """
