configfile: 'config_illumina.yaml'

## expand statement goes at the end (bottom) of each path in the dag
rule all:
    input:
        expand('{root_dir}/{sample}/{sample}_1_fastqc.zip',sample = config['sample'], root_dir = config['root_dir'])
        expand('{root_dir}/{sample}/{sample}_2_.fastqc.zip', sample = config['sample'], root_dir = config['root_dir'])

rule pre_trimming:
    input:
        r1 = lambda wildcards : config[wildcards.sample]["r1"],
		r2 = lambda wildcards : config[wildcards.sample]["r2"]
    output:
        ['{root_dir}/{sample}/{sample}_1_fastqc.zip', '{root_dir}/{sample}/{sample}_2_fastqc.zip']
    # conda:
    #     '../../envs/fastqc.yaml'
    shell:
        '''
        conda activate fastqc
        fastqc {wildcards.root_dir}/{wildcards.sample}/{input.r1}
        fastqc {wildcards.root_dir}/{wildcards.sample}/{input.r2}
        '''