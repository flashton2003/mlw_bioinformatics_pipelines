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
        expand('{root_dir}/{sample}/{sample}.fastq.gz', sample = todo_list, root_dir = root_dir)
    output:
        stats = '{root_dir}/{sample}/{sample}_contigs.assembly_stats.tsv'
    conda:
        '../../envs/assembly_stats.yaml'
    shell:
        'assembly-stats -t {input} > {output.stats}'

rule kraken2:
   input:
        expand('{root_dir}/{sample}/{sample}.fastq.gz', sample = todo_list, root_dir = root_dir)
   output:
       kraken_report = '{root_dir}/{sample}/kraken2/{sample}.kraken_report.txt'
   threads: kraken_threads
   conda:
       '../../envs/kraken2.yaml'
   shell:
       'kraken2 --gzip-compressed --use-names --output {output}  --db {DB_PATH} --report {output.kraken_report} --threads {threads} --confidence 0.9 --memory-mapping {input}'
 
rule ref_seeker:
    input:
        specie = expand('{root_dir}/{sample}/{sample}.fastq.gz', sample = todo_list, root_dir = root_dir),
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
		sam = rules.reads_stat.input,
		ref = expand('/home/ubuntu/data/belson/reference/2021.04.01/{ref}_contigs.fa',ref=refs)
	output:
		temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.08/{sample}.sam')
	shell:
		'minimap2 -ax map-ont {input.ref} {input.sam} > {output}'

rule sam2bam:
	input:
		rules.minimap.output
	output:
		temp('/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.08/{sample}.bam')
	shell:
		'samtools view -b {input} -o {output}'

rule sort_bam:
	input:
		rules.sam2bam.output
	output:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.08/{sample}.sorted.bam'
	shell:
		'samtools sort {input} -o {output}'

rule depth_calc:
	input:
		rules.sort_bam.output
	output:
		'/home/ubuntu/data/belson/isangi_nanopore/qc/results/2021.01.08/{sample}_coverage.txt'
	shell:
		"samtools depth -aa {input} | awk '{{sum+=$3}} END {{print \"Average = \",sum/NR}}' > {output}"

