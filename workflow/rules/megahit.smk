# This rule runs megahit assembler with the parameters specified in the config file.


rule megahit:
    input:
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
        kmergenie_result=get_kmergenie_dependency,
    output:
        result_dir=directory(f"{output_dir}" + "{sample}/megahit"),
        link_assembly=f"{output_dir}" + "{sample}/assemblies/{sample}_megahit.fa",
    params:
        k=lambda wildcards: get_kmer_list(wildcards, "megahit", "k"),
        scaffolds=f"{output_dir}" + "{sample}/megahit/final.contigs.fa",
        optional_params=" ".join(
            f"{k}" if v is True else f"{k} {v}"
            for k, v in config["megahit"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: get_medium_high_threads
    resources:
        mem_mb=get_medium_high_mem,
        partition=config["medium_high"]["partition"],
    log:
        "logs/{sample}/megahit.log",
    benchmark:
        "benchmark/{sample}/megahit.txt"
    conda:
        "../envs/megahit.yaml"
    container:
        "docker://quay.io/biocontainers/megahit:1.2.9--haf24da9_8"
    shell:
        """
        echo "Running megahit with the following command:" >> {log} 2>&1
        echo "megahit -t {threads} -1 {input.forward_in} -2 {input.reverse_in} -o {output.result_dir} --k-list {params.k} {params.optional_params}" >> {log} 2>&1
        megahit -t {threads} -1 {input.forward_in} -2 {input.reverse_in} -o {output.result_dir} --k-list {params.k} {params.optional_params} >> {log} 2>&1

        ln -srn {params.scaffolds} {output.link_assembly}
        """
