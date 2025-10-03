
rule get_busco_db:
    input:
        taxonomy=config["busco"]["taxonomy_file"],
        busco_lineages=config["busco"]["database_list"],
    output:
        specific_db=f"{output_dir}" + "{sample}/busco_specific/busco_db.txt",
    params:
        fallback=config["busco"]["lineage_general"],
        extension=config["busco"]["extension"],
    threads: 1
    resources:
        mem_mb=get_low_mem,
        partition=config["low"]["partition"],
    log:
        "logs/{sample}/get_busco_db.log",
    benchmark:
        "benchmark/{sample}/get_busco_db.txt"
    conda:
        "../envs/data_viz.yaml"
    shell:
        """
        python workflow/scripts/map_busco_db_to_sample.py -t {input.taxonomy} \
            -b {input.busco_lineages} \
            -o {output.specific_db} \
            -s {wildcards.sample} \
            -f {params.fallback} \
            -e {params.extension} \
            >> {log} 2>&1
        """
