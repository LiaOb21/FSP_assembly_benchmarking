# This rule runs abyss assembler with the parameters specified in the config file.
import glob
import os


rule abyss:
    input:
        forward_in= f"{input_dir}" + "/{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in= f"{input_dir}" + "/{sample}/{sample}_trimmed.R2.fq.gz",
    output:
        result_dir = directory(f"{output_dir}" + "/{sample}/abyss"),
        scaffolds = f"{output_dir}" + "/{sample}/abyss/abyss-scaffolds.fa",
        link_assembly = f"{output_dir}" + "/assemblies/{sample}/{sample}_abyss.fa"
    params:
        k = config["abyss"]["k"],
        B = config["abyss"]["B"],
        optional_params = " ".join(
            k for k, v in config["abyss"]["optional_params"].items() if v is True
        ),
    threads: config["threads"],  # access threads from config
    log:
        "logs/{sample}/abyss.log",
    resources:
        mem_mb=config["mem_mb"],  # access memory from config
    conda:
        "../envs/abyss.yaml"
    shell:
        """
        abyss-pe -C {output.result_dir} name=abyss k={params.k} B={params.B} in='{input.forward_in} {input.reverse_in}' -j {threads} {params.optional_params} >> {log} 2>&1

        ln -srn {output.scaffolds} {output.link_assembly}
        """