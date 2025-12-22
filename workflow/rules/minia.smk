# This rule runs minia assembler with the parameters specified in the config file.


rule minia:
    wildcard_constraints:
        reads_type="R1R2"
    input:
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
        kmergenie_result=get_kmergenie_dependency,
    output:
        scaffolds=f"{output_dir}" + "{reads_type}/{strategy}/{sample}/minia/{sample}.contigs.fa",
        link_assembly=f"{output_dir}" + "assemblies/{sample}/{sample}_{reads_type}_{strategy}_minia.fa",
    params:
        k=lambda wildcards: get_single_kmer(wildcards, "minia", "k"),
        result_prefix=lambda wildcards, output: os.path.splitext(output.scaffolds)[
            0
        ].replace(".contigs", ""),
        file_list=lambda wildcards, output: os.path.join(
            os.path.dirname(output.scaffolds), f"{wildcards.sample}_files.txt"
        ),
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
        "logs/{sample}/minia_{reads_type}_{strategy}.log",
    benchmark:
        "benchmark/{sample}/minia_{reads_type}_{strategy}.txt"
    conda:
        "../envs/minia.yaml"
    container:
        "docker://quay.io/biocontainers/minia:3.2.6--h22625ea_5"
    shell:
        """
        # Create file list for minia
        echo "{input.forward_in}" > {params.file_list}
        echo "{input.reverse_in}" >> {params.file_list}

        echo "Running minia with the following command:" >> {log} 2>&1
        echo "minia -in {params.file_list} -kmer-size {params.k} -out {params.result_prefix} -nb-cores {threads} {params.optional_params}" >> {log} 2>&1
        minia -in {params.file_list} -kmer-size {params.k} -out {params.result_prefix} -nb-cores {threads} {params.optional_params} >> {log} 2>&1

        ln -srn {output.scaffolds} {output.link_assembly}
        """
