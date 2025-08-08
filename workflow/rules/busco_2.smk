rule busco_2:
    input:
        assembly=f"{output_dir}" + "assemblies/{sample}/{sample}_{assembler}_pilon.fa",
    output:
        general_dir=directory(f"{output_dir}" + "{sample}/busco_general_pilon/{assembler}"),
        specific_dir=directory(f"{output_dir}" + "{sample}/busco_specific_pilon/{assembler}"),
    params:
        lineage_general=config["busco"]["lineage_general"],
        lineage_specific=config["busco"]["lineage_specific"],
        optional_params=" ".join(
            f"{k}" if v is True else f"{k} {v}"
            for k, v in config["busco"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: get_scaled_threads  # Use scaling function
    log:
        "logs/{sample}/busco_{sample}_{assembler}_pilon.log",
    benchmark:
        "benchmark/{sample}/busco_{sample}_{assembler}_pilon.txt"
    resources:
        mem_mb=get_scaled_mem,  # Use scaling function
    conda:
        "../envs/busco.yaml"
    shell:
        """
        echo "Running busco with the following command:" >> {log} 2>&1
        echo "busco -i {input.assembly} --out_path {output.general_dir} -l {params.lineage_general} -f -c {threads} {params.optional_params}" >> {log} 2>&1
        busco -i {input.assembly} --out_path {output.general_dir} -l {params.lineage_general} -f -c {threads} {params.optional_params} >> {log} 2>&1
        
        echo "Running busco with the following command:" >> {log} 2>&1
        echo "busco -i {input.assembly} --out_path {output.specific_dir} -l {params.lineage_specific} -f -c {threads} {params.optional_params}" >> {log} 2>&1
        busco -i {input.assembly} --out_path {output.specific_dir} -l {params.lineage_specific} -f -c {threads} {params.optional_params} >> {log} 2>&1
        """