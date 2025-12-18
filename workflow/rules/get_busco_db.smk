
rule get_busco_db:
    input:
        taxonomy=config["busco"]["taxonomy_file"],
        busco_lineages=config["busco"]["database_list"],
    output:
        specific_db=f"{output_dir}" + "busco_db/{sample}/busco_db.txt",
    params:
        fallback=config["busco"]["lineage_general"],
        extension=config["busco"]["extension"],
    threads: 1
    resources:
        mem_mb=get_very_low_mem,
        partition=config["very_low"]["partition"],
    log:
        "logs/busco_db/{sample}/get_busco_db.log",
    benchmark:
        "benchmark/busco_db/{sample}/get_busco_db.txt"
    conda:
        "../envs/data_viz.yaml"
    container:
        "docker://lobinu21/data_viz_py:latest"
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
