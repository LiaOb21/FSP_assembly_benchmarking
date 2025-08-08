# This rule runs pilon on the draft assemblies
import glob
import os


rule pilon:
    input:
        assembly=f"{output_dir}" + "assemblies/{sample}/{sample}_{assembler}.fa",
        sorted_bam=f"{output_dir}" + "{sample}/bwa_mem2_samtools/{assembler}/{sample}_{assembler}_sorted.bam",
    output:
        pilon_fasta=f"{output_dir}" + "{sample}/pilon/{assembler}/{sample}_{assembler}_pilon.fasta",
        link_pilon_assembly=f"{output_dir}" + "assemblies/{sample}/{sample}_{assembler}_pilon.fa",
    params:
        results_prefix=lambda wildcards, output: os.path.splitext(output.pilon_fasta)[
            0
        ].replace(".fasta", ""),
        optional_params=" ".join(
            f"{key}" if value is True else f"{key} {value}"
            for key, value in config["pilon"]["optional_params"].items()
            if value and value is not False and value != ""
        ),
    threads: get_scaled_threads  # Use scaling function
    log:
        "logs/{sample}/pilon_{assembler}.log",
    benchmark:
        "benchmark/{sample}/pilon_{assembler}.txt"
    resources:
        mem_mb=get_scaled_mem,  # Use scaling function
    conda:
        "../envs/pilon.yaml"
    shell:
        """
        echo "Running pilon with the following command:" >> {log} 2>&1
        echo "pilon --genome {input.assembly} --bam {input.sorted_bam} --output {params.results_prefix} --threads {threads} {params.optional_params} --verbose" >> {log} 2>&1
        pilon --genome {input.assembly} --bam {input.sorted_bam} --output {params.results_prefix} --threads {threads} {params.optional_params} --verbose >> {log} 2>&1
        ln -srn {output.pilon_fasta} {output.link_pilon_assembly}
        """
