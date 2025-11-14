# This rule runs seqkit to calculate the read median length and the kmer value as 2/3rds of the median if the config file specifies kmer_strategy: mode: reads_length.


rule seqkit:
    input:
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
    output:
        seqkit_results=f"{output_dir}" + "{sample}/seqkit/{sample}_seqkit.txt",
        kmer=f"{output_dir}" + "{sample}/seqkit/{sample}_kmer_value.txt",
    threads: get_very_low_threads
    resources:
        mem_mb=get_very_low_mem,
        partition=config["very_low"]["partition"],
    log:
        "logs/{sample}/seqkit.log",
    benchmark:
        "benchmark/{sample}/seqkit.txt"
    conda:
        "../envs/seqkit.yaml"
    container:
        "docker://quay.io/biocontainers/seqkit:2.10.1--he881be0_0"
    shell:
        """
        echo "Running seqkit with the following command:" >> {log} 2>&1
        echo "seqkit stats -a {input.forward_in} -j {threads} > {output.seqkit_results}" >> {log} 2>&1
        seqkit stats -a {input.forward_in} -j {threads} > {output.seqkit_results} >> {log} 2>&1
        # calculate kmer as 2/3rds of median read length
        awk '{{print $10}}' {output.seqkit_results} | tail -n1 | awk '{{result = int($1 * 2 / 3); print (result % 2 == 0) ? result - 1 : result}}' > {output.kmer}
        echo "Calculated kmer value (2/3rds of median read length):" >> {log} 2>&1
        cat {output.kmer} >> {log} 2>&1
        """
