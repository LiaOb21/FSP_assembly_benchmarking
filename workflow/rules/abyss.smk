# This rule runs abyss assembler with the parameters specified in the config file.


rule abyss:
    wildcard_constraints:
        reads_type="R1R2"
    input:
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
        kmergenie_result=get_kmergenie_dependency,
    output:
        result_dir=directory(f"{output_dir}" + "{reads_type}/{strategy}/{sample}/abyss"),
        scaffolds=f"{output_dir}" + "{reads_type}/{strategy}/{sample}/abyss/abyss-scaffolds.fa",
        link_assembly=f"{output_dir}" + "assemblies/{sample}/{sample}_{reads_type}_{strategy}_abyss.fa",
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
    threads: get_medium_threads
    resources:
        mem_mb=get_high_mem,
        partition=config["high"]["partition"],
    log:
        "logs/{sample}/abyss_{reads_type}_{strategy}.log",
    benchmark:
        "benchmark/{sample}/abyss_{reads_type}_{strategy}.txt"
    conda:
        "../envs/abyss.yaml"
    container:
        "docker://quay.io/biocontainers/abyss:2.3.10--hf316886_2"
    shell:
        """
        echo "Running abyss with the following command:" >> {log} 2>&1
        echo "abyss-pe -C {output.result_dir} name=abyss k={params.k} B={params.B} in='{input.forward_in} {input.reverse_in}' -j {threads} {params.optional_params}" >> {log} 2>&1
        abyss-pe -C {output.result_dir} name=abyss k={params.k} B={params.B} in='{input.forward_in} {input.reverse_in}' -j {threads} {params.optional_params} >> {log} 2>&1

        ln -srn {output.scaffolds} {output.link_assembly}
        """
