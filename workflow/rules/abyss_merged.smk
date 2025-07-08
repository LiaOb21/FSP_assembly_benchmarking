# This rule runs abyss assembler with the parameters specified in the config file.
import glob
import os


rule abyss:
    input:
        merged_in=f"{input_dir}" + "{sample}/{sample}_merge.fq.gz",
    output:
        result_dir=directory(f"{output_dir}" + "{sample}/abyss"),
        scaffolds=f"{output_dir}" + "{sample}/abyss/abyss-scaffolds.fa",
        link_assembly=f"{output_dir}" + "assemblies/{sample}/{sample}_abyss.fa",
    params:
        k=config["abyss"]["k"],
        B=config["abyss"]["B"],
        optional_params=" ".join(
            f"{k}={v}" if v is not True else k
            for k, v in config["abyss"]["optional_params"].items()
            if v
        ),
    threads: config["threads"]  # access threads from config
    log:
        "logs/{sample}/abyss.log",
    resources:
        mem_mb=config["mem_mb"],  # access memory from config
    conda:
        "../envs/abyss.yaml"
    shell:
        """
        echo "Running abyss with the following command:" >> {log} 2>&1
        echo "abyss-pe -C {output.result_dir} name=abyss k={params.k} B={params.B} se='{input.merged_in}' -j {threads} {params.optional_params}" >> {log} 2>&1
        abyss-pe -C {output.result_dir} name=abyss k={params.k} B={params.B} se='{input.merged_in}' -j {threads} {params.optional_params} >> {log} 2>&1

        ln -srn {output.scaffolds} {output.link_assembly}
        """
