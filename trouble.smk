configfile:'config_nano.yaml'

rule all:
    input:
        expand('{root_dir}/{sample}/assembly_stat/{sample}_reads.assembly_stats.tsv',sample = config["samples"], root_dir = config["root"]),
        expand('{root_dir}/{sample}/kraken2/{sample}.kraken_report.txt',sample = config["samples"], root_dir = config["root"]),
        expand('{root_dir}/{sample}/Refseeker/{sample}.Refseeker.txt',sample = config["samples"], root_dir = config["root"]),
        expand('{root_dir}/{sample}/read_depth/{sample}_coverage.txt',sample = config["samples"], root_dir = config["root"]),
        expand('{root_dir}/{sample}/Flye',sample = config["samples"], root_dir = config["root"]),
        expand('{root_dir}/{sample}/Medaka',sample = config["samples"], root_dir = config["root"]),
        expand('{root_dir}/{sample}/Bakta',sample = config["samples"], root_dir = config["root"])

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

rule kraken2:
   input:
        read = '{root_dir}/{sample}/{sample}.fastq.gz',
        db = config['kraken_db']
   output:
       out1 = '{root_dir}/{sample}/kraken2/{sample}.kraken_report.txt',
       out2 = '{root_dir}/{sample}/kraken2/{sample}.kraken'
   threads: 8
   shell:
        '''
        mkdir -p $( dirname {input.read} )
        conda activate kraken2
        kraken2 --use-names --threads {threads} --db {input.db} --report {output.out1} --gzip-compressed {input.read}> {output.out2}
        '''

rule minimap:
	input:
		read = rules.assembly_stats.input.read,
		ref = config['ref_genome'] 
	output:
		temp('{root_dir}/{sample}/{sample}.sam')
	shell:
		'minimap2 -ax map-ont {input.ref} {input.read} > {output}'

rule sam2bam:
	input:
		rules.minimap.output
	output:
		temp('{root_dir}/{sample}/{sample}.bam')
	shell:
        'samtools sort -@ 8 -o {output} {input}'

rule depth_calc:
	input:
		rules.sam2bam.output
	output:
		'{root_dir}/{sample}/read_depth/{sample}_coverage.txt'
	shell:
        '''
        mkdir -p {wildcards.root_dir}/{wildcards.sample}/read_depth
		samtools depth -aa {input} | awk '{{sum+=$3}} END {{print \"Average = \",sum/NR}}' > {output}
        '''

rule flye:
	input:
		rules.assembly_stats.input.read
	output:
		directory('{root_dir}/{sample}/Flye')
	conda:
		'/home/ubuntu/data/belson/Guppy5_guppy3_comparison/napa/scripts/envs/flye.yml'
	shell:
		'flye --nano-hq {input} -g 5m -o {output} -t 8 --plasmids'

rule racon:
	input:
		genome = rules.flye.output,
		nano = rules.assembly_stats.input.read
	output:
		racon = temp('{root_dir}/{sample}/{sample}racon.fasta'),
		paf = temp('{root_dir}/{sample}/{sample}.racon.paf')
	shell:
        '''
		minimap2 -x map-ont {input.genome}/assembly.fasta {input.nano} > {output.paf}
        racon -t 4 {input.nano} {output.paf} {input.genome}/assembly.fasta > {output.racon}
        '''

rule medaka:
	input:
		genome = rules.racon.output.racon,
		nano = rules.assembly_stats.input.read,
        model = config['medaka_model']
	output:
		directory('{root_dir}/{sample}/Medaka')
	shell:
        '''
        conda activate medaka
		medaka_consensus -i {input.nano} -d {input.genome} -t 8  -m {input.model} -o {output}
        '''

rule bakta:
    input:
        rules.medaka.output,
        config['bakta_db']
    output:
        directory('{root_dir}/{sample}/Bakta')
    shell:
        'bakta --db {input[1]} {input[0]}/consensus.fasta -o {output}'