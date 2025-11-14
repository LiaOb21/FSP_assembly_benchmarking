# This rule generates a FastK table needed to run merquryFK


rule fastk:
    input:
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
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
    threads: get_low_threads
    resources:
        mem_mb=get_low_mem,
        partition=config["low"]["partition"],
    log:
        "logs/{sample}/fastk.log",
    benchmark:
        "benchmark/{sample}/fastk.txt"
    conda:
        "../envs/merquryFK.yaml"
    container:
        "docker://quay.io/biocontainers/merquryfk:1.1.3--h71df26d_0"
    shell:
        """
        mkdir -p {params.temp_dir}

        echo "Running FastK with the following command:" >> {log} 2>&1
        echo "FastK -v -T{threads} -k{params.k} -M{params.memory_gb} {input.forward_in} {input.reverse_in} -N{params.result_prefix} -t{params.t} -P{params.temp_dir}" >> {log} 2>&1
        FastK -v -T{threads} -k{params.k} -M{params.memory_gb} {input.forward_in} {input.reverse_in} -N{params.result_prefix} -t{params.t} -P{params.temp_dir}>> {log} 2>&1
        rm -rf {params.temp_dir}  # clean up temporary directory
        """
