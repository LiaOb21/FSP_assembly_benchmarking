# This rule runs SparseAssembler assembler with the parameters specified in the config file.
import glob
import os


rule sparseassembler:
    input:
        forward_in=f"{output_dir}" + "fqreads/{sample}/{sample}_trimmed.R1.fq",
        reverse_in=f"{output_dir}" + "fqreads/{sample}/{sample}_trimmed.R2.fq",
        kmergenie_result=get_kmergenie_dependency,
    output:
        scaffolds=f"{output_dir}" + "{sample}/sparseassembler/SuperContigs.txt",
        link_assembly=f"{output_dir}"
        + "{sample}/assemblies/{sample}_sparseassembler.fa",
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
    threads: get_scaled_threads  # Use scaling function
    log:
        "logs/{sample}/sparseassembler.log",
    benchmark:
        "benchmark/{sample}/sparseassembler.txt"
    resources:
        mem_mb=get_scaled_mem,  # Use scaling function
    conda:
        "../envs/sparseassembler.yaml"
    shell:
        """
        cd {params.result_dir}

        echo "Running SparseAssembler with the following command:" >> sparseassembler.log 2>&1
        echo "SparseAssembler k {params.k} GS {params.GS} p1 {input.forward_in} p2 {input.reverse_in} Scaffold {params.Scaffold} ExpCov {params.ExpCov} {params.optional_params}" >> sparseassembler.log 2>&1
        SparseAssembler k {params.k} GS {params.GS} p1 {input.forward_in} p2 {input.reverse_in} Scaffold {params.Scaffold} ExpCov {params.ExpCov} {params.optional_params} >> sparseassembler.log 2>&1

        cd -

        cp {params.result_dir}/sparseassembler.log {log}
        
        ln -srn {output.scaffolds} {output.link_assembly}
        """
