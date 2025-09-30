# This rule runs abyss assembler with the parameters specified in the config file.
import glob
import os


rule abyss:
    input:
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
        kmergenie_result=get_kmergenie_dependency,
    output:
        result_dir=directory(f"{output_dir}" + "{sample}/abyss"),
        scaffolds=f"{output_dir}" + "{sample}/abyss/abyss-scaffolds.fa",
        link_assembly=f"{output_dir}" + "{sample}/assemblies/{sample}_abyss.fa",
    params:
        k=lambda wildcards: get_single_kmer(wildcards, "abyss", "k"),
        B=config["abyss"]["B"],
        optional_params=" ".join(
            (
                key
                if value is True
                else f"{key} {value}" if key.startswith("-") else f"{key}={value}"
            )
            for key, value in config["abyss"]["optional_params"].items()
            if value and value is not False and value != ""
        ),
    threads: get_high_threads
    resources:
        mem_mb=get_high_mem,
    log:
        "logs/{sample}/abyss.log",
    benchmark:
        "benchmark/{sample}/abyss.txt"
    conda:
        "../envs/abyss.yaml"
    shell:
        """
        echo "Running abyss with the following command:" >> {log} 2>&1
        echo "abyss-pe -C {output.result_dir} name=abyss k={params.k} B={params.B} in='{input.forward_in} {input.reverse_in}' -j {threads} {params.optional_params}" >> {log} 2>&1
        abyss-pe -C {output.result_dir} name=abyss k={params.k} B={params.B} in='{input.forward_in} {input.reverse_in}' -j {threads} {params.optional_params} >> {log} 2>&1

        ln -srn {output.scaffolds} {output.link_assembly}
        """
