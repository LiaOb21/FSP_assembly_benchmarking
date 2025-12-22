# This rule runs SparseAssembler assembler with the parameters specified in the config file.


rule sparseassembler_merged:
    wildcard_constraints:
        reads_type="merged"
    input:
        merged_in=f"{output_dir}" + "{reads_type}/fqreads/{sample}/{sample}_merged.fq",
        kmergenie_result=get_kmergenie_dependency,
    output:
        scaffolds=f"{output_dir}" + "{reads_type}/{strategy}/{sample}/sparseassembler/SuperContigs.txt",
        link_assembly=f"{output_dir}"
        + "assemblies/{sample}/{sample}_{reads_type}_{strategy}_sparseassembler.fa",
    params:
        result_dir=lambda wildcards, output: os.path.dirname(output.scaffolds),
        k=lambda wildcards: get_single_kmer(wildcards, "sparseassembler", "k"),
        GS=config["sparseassembler"]["GS"],
        Scaffold=config["sparseassembler"]["Scaffold"],
        ExpCov=config["sparseassembler"]["ExpCov"],
        optional_params=" ".join(
            f"{k}" if v is True else f"{k} {v}"
            for k, v in config["sparseassembler"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: 1
    resources:
        mem_mb=get_medium_mem,
        partition=config["medium_high"]["partition"],
    log:
        "logs/{sample}/sparseassembler_{reads_type}_{strategy}.log",
    benchmark:
        "benchmark/{sample}/sparseassembler_{reads_type}_{strategy}.txt"
    conda:
        "../envs/sparseassembler.yaml"
    container:
        "docker://quay.io/biocontainers/sparseassembler:20160205--h9948957_11"
    shell:
        """
        cd {params.result_dir}
        echo "Running SparseAssembler with the following command:" >> sparseassembler.log 2>&1
        echo "SparseAssembler k {params.k} GS {params.GS} f {input.merged_in} Scaffold {params.Scaffold} ExpCov {params.ExpCov} {params.optional_params}" >> sparseassembler.log 2>&1
        SparseAssembler k {params.k} GS {params.GS} f {input.merged_in} Scaffold {params.Scaffold} ExpCov {params.ExpCov} {params.optional_params} >> sparseassembler.log 2>&1

        cd -

        cp {params.result_dir}/sparseassembler.log {log}
        
        ln -srn {output.scaffolds} {output.link_assembly}
        """
