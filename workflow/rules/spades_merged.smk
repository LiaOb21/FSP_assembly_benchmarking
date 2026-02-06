# This rule runs spades assembler with the parameters specified in the config file.


rule spades_merged:
    wildcard_constraints:
        reads_type="merged",
    input:
        merged_in=f"{input_dir}" + "{sample}/{sample}_merge.fq.gz",
        unmerged_r1=f"{input_dir}" + "{sample}/{sample}_unmerged.R1.fq.gz",
        unmerged_r2=f"{input_dir}" + "{sample}/{sample}_unmerged.R2.fq.gz",
        kmergenie_result=get_kmergenie_dependency,
    output:
        result_dir=directory(
            f"{output_dir}" + "{reads_type}/{strategy}/{sample}/spades"
        ),
        scaffolds=f"{output_dir}"
        + "{reads_type}/{strategy}/{sample}/spades/scaffolds.fasta",
        link_assembly=f"{output_dir}"
        + "assemblies/{sample}/{sample}_{reads_type}_{strategy}_spades.fa",
    params:
        k=lambda wildcards: get_kmer_list(wildcards, "spades", "k"),
        optional_params=" ".join(
            f"{k} {v}" if v is not True else k
            for k, v in config["spades"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: get_high_threads
    resources:
        mem_mb=get_high_mem,
        partition=config["high"]["partition"],
    log:
        "logs/{sample}/spades_{reads_type}_{strategy}.log",
    benchmark:
        "benchmark/{sample}/spades_{reads_type}_{strategy}.txt"
    conda:
        "../envs/spades.yaml"
    container:
        "docker://quay.io/biocontainers/spades:4.2.0--h8d6e82b_2"
    shell:
        """
        echo "Running spades with the following command:" >> {log} 2>&1
        echo "spades.py -t {threads} --merged {input.merged_in} -1 {input.unmerged_r1} -2 {input.unmerged_r2} -o {output.result_dir} -k {params.k} {params.optional_params}" >> {log} 2>&1
        spades.py -t {threads} --merged {input.merged_in} -1 {input.unmerged_r1} -2 {input.unmerged_r2} -o {output.result_dir} -k {params.k} {params.optional_params} >> {log} 2>&1
        ln -srn {output.scaffolds} {output.link_assembly}
        """
