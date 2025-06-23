# This rule runs SparseAssembler assembler with the parameters specified in the config file.
import glob
import os


rule sparseassembler:
    input:
        forward_in = f"{output_dir}" + "/fqreads/{sample}/{sample}_trimmed.R1.fq",
        reverse_in = f"{output_dir}" + "/fqreads/{sample}/{sample}_trimmed.R2.fq",
    output:
        result_dir = directory(f"{output_dir}" + "/{sample}/sparseassembler"),
        scaffolds = f"{output_dir}" + "/{sample}/sparseassembler/SuperContigs.txt",
    params:
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
        SparseAssembler k {params.k} GS {params.GS} p1 {input.forward_in} p2 {input.reverse_in} Scaffold {params.Scaffold} ExpCov {params.ExpCov} {params.optional_params} >> {log} 2>&1
        mkdir -p $(dirname {output.scaffolds})
        mv *txt {output.result_dir}
        mv *HT* {output.result_dir}
        """