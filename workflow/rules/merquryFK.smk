# This rule generates a runs merquryFK to evaluate assembly quality based on kmer
import glob
import os


rule merquryfk:
    input:
        ktab=f"{output_dir}" + "{sample}/fastk/fastk_table.ktab",
        assembly=f"{output_dir}" + "{sample}/assemblies/{sample}_{assembler}.fa",
    output:
        stats=f"{output_dir}"
        + "{sample}/merquryfk/{assembler}/merquryfk.completeness.stats",
        qv=f"{output_dir}" + "{sample}/merquryfk/{assembler}/merquryfk.qv",
    params:
        result_prefix=lambda wildcards, output: os.path.splitext(output.qv)[0],
        temp_dir=lambda wildcards, output: os.path.join(
            os.path.dirname(output.qv), "temp"
        ),
        optional_params=" ".join(
            f"{k}" if v is True else f"{k} {v}"
            for k, v in config["merquryfk"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: get_low_threads
    resources:
        mem_mb=get_low_mem,
    log:
        "logs/{sample}/merquryfk_{sample}_{assembler}.log",
    benchmark:
        "benchmark/{sample}/merquryFK_{sample}_{assembler}.txt"
    conda:
        "../envs/merquryFK.yaml"
    shell:
        """
        mkdir -p {params.temp_dir}

        echo "Running merquryFK with the following command:" >> {log} 2>&1
        echo "MerquryFK -lfs -v -T{threads} -P{params.temp_dir} {input.ktab} {input.assembly} {params.result_prefix} {params.optional_params}" >> {log} 2>&1
        MerquryFK -lfs -v -T{threads} -P{params.temp_dir} {input.ktab} {input.assembly} {params.result_prefix} {params.optional_params} >> {log} 2>&1
        rm -rf {params.temp_dir}
        """
