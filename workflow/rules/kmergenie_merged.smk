# This rule runs kmergenie to find the best kmer value for subsequent assembly if the config file specifies kmer_strategy: mode: auto.


rule kmergenie_merged:
    wildcard_constraints:
        reads_type="merged",
    input:
        merged_in=f"{input_dir}" + "{sample}/{sample}_merge.fq.gz",
    output:
        kmergenie_report=f"{output_dir}"
        + "{reads_type}/{strategy}/{sample}/kmergenie/{sample}_report.html",
        best_kmer=f"{output_dir}"
        + "{reads_type}/{strategy}/{sample}/kmergenie/{sample}_best_kmer.txt",
    params:
        k=config["kmergenie"]["k"],
        l=config["kmergenie"]["l"],
        result_prefix=lambda wildcards, output: os.path.splitext(
            output.kmergenie_report
        )[0].replace("_report", ""),
        optional_params=" ".join(
            f"{k}" if v is True else f"{k} {v}"
            for k, v in config["kmergenie"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: get_medium_threads
    resources:
        mem_mb=get_medium_mem,
        partition=config["medium"]["partition"],
    log:
        "logs/{sample}/kmergenie_{reads_type}_{strategy}.log",
    benchmark:
        "benchmark/{sample}/kmergenie_{reads_type}_{strategy}.txt"
    conda:
        "../envs/kmergenie.yaml"
    container:
        "docker://quay.io/biocontainers/kmergenie:1.7051--py27r40h077b44d_11"
    shell:
        """
        echo "Running kmergenie with the following command:" >> {log} 2>&1
        echo "kmergenie -o {params.result_prefix} -l {params.l} -k {params.k} -t {threads} {params.optional_params} {input.merged_in}" >> {log} 2>&1
        kmergenie -o {params.result_prefix} -l {params.l} -k {params.k} -t {threads} {params.optional_params} {input.merged_in} >> {log} 2>&1

        # Extract best k-mer value from HTML report using awk
        awk -f workflow/scripts/extract_best_k.awk < {output.kmergenie_report} > {output.best_kmer}
        """
