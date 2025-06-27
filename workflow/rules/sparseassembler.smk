# This rule runs SparseAssembler assembler with the parameters specified in the config file.
import glob
import os


rule sparseassembler:
    input:
        forward_in = f"{output_dir}" + "/fqreads/{sample}/{sample}_trimmed.R1.fq",
        reverse_in = f"{output_dir}" + "/fqreads/{sample}/{sample}_trimmed.R2.fq",
    output:
        scaffolds = f"{output_dir}" + "/{sample}/sparseassembler/SuperContigs.txt",
        link_assembly = f"{output_dir}" + "/assemblies/{sample}/{sample}_sparseassembler.fa"
    params:
        result_dir = directory(f"{output_dir}" + "/{sample}/sparseassembler"),
        k = config["sparseassembler"]["k"],
        GS = config["sparseassembler"]["GS"],
        Scaffold = config["sparseassembler"]["Scaffold"],
        ExpCov = config["sparseassembler"]["ExpCov"],
        optional_params = " ".join(
            k for k, v in config["sparseassembler"]["optional_params"].items() if v is True
        ),
    threads: config["threads"],  # access threads from config
    log:
        "logs/{sample}/sparseassembler.log",
    resources:
        mem_mb=config["mem_mb"],  # access memory from config
    conda:
        "../envs/sparseassembler.yaml"
    shell:
        """
        cd {params.result_dir}
        SparseAssembler k {params.k} GS {params.GS} p1 {input.forward_in} p2 {input.reverse_in} Scaffold {params.Scaffold} ExpCov {params.ExpCov} {params.optional_params} >> sparseassembler.log 2>&1

        cd -

        cp {params.result_dir}/sparseassembler.log {log}
        
        ln -srn {output.scaffolds} {output.link_assembly}
        """