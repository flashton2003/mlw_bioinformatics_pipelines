configfile:'config.yaml'

## expand statement goes at the end (bottom) of each path in the dag
rule all:
    input:
        expand('{root_dir}/{sample}.kraken_report.txt', sample = config["samples"], root_dir = config["root"])


# Check sed & proper sample name dir

rule kraken2:
   input:
        read = '{root_dir}/samples/{sample}',
        db = '/home/ubuntu/data/belson/kraken/' #replce with from config
   output:
       '{root_dir}/{sample}.kraken_report.txt'
   threads: 8
   shell:
        'conda activate /home/ubuntu/data/belson/kraken'
        'kraken2 --gzip-compressed --use-names --output {output}  --db {input.db} --report {output} --threads {threads} --confidence 0.9 --memory-mapping {input.read}'