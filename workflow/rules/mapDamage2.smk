# This rule runs mapDamage2 on the draft assemblies
import glob
import os


rule mapDamage2:
    input:
        assembly=f"{output_dir}" + "assemblies/{sample}/{sample}_{assembler}.fa",
        sorted_bam=f"{output_dir}" + "{sample}/bwa_mem2_samtools/{assembler}/{sample}_{assembler}_sorted.bam",
    output:
        results_dir=directory(f"{output_dir}" + "{sample}/mapDamage2/{assembler}"),
    params:
        optional_params=" ".join(
            f"{key}" if value is True else f"{key} {value}"
            for key, value in config["mapDamage2"]["optional_params"].items()
            if value and value is not False and value != ""
        ),
    threads: 1
    log:
        "logs/{sample}/mapDamage2_{assembler}.log",
    benchmark:
        "benchmark/{sample}/mapDamage2_{assembler}.txt"
    resources:
        mem_mb=get_scaled_mem,  # Use scaling function
    conda:
        "../envs/mapDamage2.yaml"
    shell:
        """
        echo "Running mapDamage with the following command:" >> {log} 2>&1
        echo "mapDamage -i {input.sorted_bam} -r {input.assembly} -d {output.results_dir} {params.optional_params}" >> {log} 2>&1
        mapDamage -i {input.sorted_bam} -r {input.assembly} -d {output.results_dir} {params.optional_params} >> {log} 2>&1
        """
