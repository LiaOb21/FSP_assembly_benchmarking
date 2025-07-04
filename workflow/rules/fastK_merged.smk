# This rule generates a FastK table needed to run merquryFK
import glob
import os


rule fastk:
    input:
        merged_in=f"{input_dir}" + "/{sample}/{sample}_merge.fq.gz",
    output:
        ktab=f"{output_dir}" + "{sample}/fastk/fastk_table.ktab",
        hist=f"{output_dir}" + "{sample}/fastk/fastk_table.hist",
    params:
        k=config["fastk"]["k"],
        result_prefix=lambda wildcards, output: os.path.splitext(output.ktab)[0],
        t=config["fastk"]["t"],
        memory_gb=lambda wildcards, resources: resources.mem_mb // 1024,
        temp_dir=lambda wildcards, output: os.path.join(
            os.path.dirname(output.ktab), "temp"
        ),
        optional_params=" ".join(
            k for k, v in config["fastk"]["optional_params"].items() if v is True
        ),
    threads: config["threads"]  # access threads from config
    log:
        "logs/{sample}/fastk.log",
    resources:
        mem_mb=config["mem_mb"],  # access memory from config
    conda:
        "../envs/merquryFK.yaml"
    shell:
        """
        mkdir -p {params.temp_dir}
        FastK -v -T{threads} -k{params.k} -M{params.memory_gb} {input.merged_in} -N{params.result_prefix} -t{params.t} -P{params.temp_dir} {params.optional_params} >> {log} 2>&1
        rm -rf {params.temp_dir}  # clean up temporary directory
        """
