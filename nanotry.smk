configfile:'config_nano.yaml'

rule all:
    input:
        expand('{root_dir}/{sample}/assembly_stat/{sample}_reads.assembly_stats.tsv', sample = config["samples"], root_dir = config["root"])

rule assembly_stats:
    input:
        read = '{root_dir}/{sample}/{sample}.fastq.gz'
    output:
        stats = '{root_dir}/{sample}/assembly_stat/{sample}_reads.assembly_stats.tsv'
    shell:
        '''
        mkdir -p $( dirname {input.read} )
        conda activate assembly-stat
        assembly-stats <(gzip -cd {input.read}) | tee {output.stats}
        '''