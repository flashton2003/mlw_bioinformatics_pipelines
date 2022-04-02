configfile:'config.yaml'

## expand statement goes at the end (bottom) of each path in the dag
rule all:
    input:
        expand('{root_dir}/{sample}.kraken_report.txt', sample = config["samples"], root_dir = config["root"])


# Check sed & proper sample name dir

rule kraken2:
   input:
        read = '{root_dir}/samples/{sample}',
        db = 'kraken' #replce with from config
   output:
       out = '{root_dir}/{sample}.kraken_report.txt',
       report = '{root_dir}/{sample}.kraken'
   threads: 8
   shell:
        '''
        conda activate /home/ubuntu/data/belson/kraken
        kraken2 --use-names --threads 4 --db kraken --report {output.out} --gzip-compressed {input.read}> {output.report}
        '''