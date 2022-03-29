import os

def read_todo_list(todo_list):
    with open(todo_list) as fi:
        lines = fi.readlines()
        lines = [x.strip() for x in lines]
    return lines

todo_list = read_todo_list(config['todo_list'])
root_dir = config['root_dir']

kraken_threads = 8

assert os.path.exists(root_dir)
if not os.path.exists(qc_results_dir):  
    os.makedirs(qc_results_dir)

## expand statement goes at the end (bottom) of each path in the dag
rule all:
    input:
        expand('{root_dir}/{sample}/{sample}_contigs.assembly_stats.tsv', sample = todo_list, root_dir = root_dir),
        expand('{root_dir}/{sample}/kraken2/{sample}.kraken_report.txt', sample = todo_list, root_dir = root_dir),
        expand('{root_dir}/{sample}/Refseeker/{sample}.Refseeker.txt', sample = todo_list, root_dir = root_dir),
        expand('{root_dir}/{sample}/read_depth/{sample}.read_depth.txt', sample = todo_list, root_dir = root_dir)

rule assembly_stats:
    input:
        read = expand('{root_dir}/{sample}/{sample}.fastq.gz', sample = todo_list, root_dir = root_dir)
    output:
        stats = '{root_dir}/{sample}/{sample}_contigs.assembly_stats.tsv'
    conda:
        '../../envs/assembly_stats.yaml'
    shell:
        'assembly-stats -t {input.read} > {output.stats}'

rule kraken2:
   input:
        rules.assembly_stats.input.read
   output:
       kraken_report = '{root_dir}/{sample}/kraken2/{sample}.kraken_report.txt'
   threads: kraken_threads
   conda:
       '../../envs/kraken2.yaml'
   shell:
       'kraken2 --gzip-compressed --use-names --output {output}  --db {DB_PATH} --report {output.kraken_report} --threads {threads} --confidence 0.9 --memory-mapping {input}'
 
rule ref_seeker:
    input:
        specie = rules.assembly_stats.input.read,
        db = '...'
    output:
        refSeeker_results = '{root_dir}/{sample}/Refseeker/{sample}.Refseeker.txt'
    conda:
        '../../envs/refseeker.yaml'
    shell:
        'referenceseeker {input.db} {input.specie} | tee {output}'

##Reference seeker? If no reference genome given.
##Conditional logic within snakemake
##Mapping against a reference and getting depth
rule minimap:
	input:
		specie = rules.assembly_stats.input.read,
		ref = expand('/home/ubuntu/data/belson/reference/2021.04.01/{ref}_contigs.fa',ref=refs)
	output:
		temp('{root_dir}/{sample}/minimap/{sample}.sam')
	shell:
		'minimap2 -ax map-ont {input.ref} {input.specie} > {output}'

rule sam2bam:
	input:
		rules.minimap.output
	output:
		temp('{root_dir}/{sample}/minimap/{sample}.bam')
	shell:
        'samtools sort -@ 8 -o {output} {input}'

rule depth_calc:
	input:
		rules.sort_bam.output
	output:
		'{root_dir}/{sample}/minimap/{sample}_coverage.txt'
	shell:
		"samtools depth -aa {input} | awk '{{sum+=$3}} END {{print \"Average = \",sum/NR}}' > {output}"

