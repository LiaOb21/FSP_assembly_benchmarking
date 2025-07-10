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
            (
                key
                if value is True
                else f"{key} {value}" if key.startswith("-") else f"{key}={value}"
            )
            for key, value in config["abyss"]["optional_params"].items()
            if value and value is not False and value != ""
        ),
    threads: get_scaled_threads  # Use scaling function
    log:
        "logs/{sample}/abyss.log",
    benchmark:
        "benchmark/{sample}/abyss.txt"
    resources:
        mem_mb=get_scaled_mem,  # Use scaling function
    conda:
        "../envs/abyss.yaml"
    shell:
        """
        echo "Running abyss with the following command:" >> {log} 2>&1
        echo "abyss-pe -C {output.result_dir} name=abyss k={params.k} B={params.B} se='{input.merged_in}' -j {threads} {params.optional_params}" >> {log} 2>&1
        abyss-pe -C {output.result_dir} name=abyss k={params.k} B={params.B} se='{input.merged_in}' -j {threads} {params.optional_params} >> {log} 2>&1

        ln -srn {output.scaffolds} {output.link_assembly}
        """
