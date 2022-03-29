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

    Nanopore "basic"
        Read length distribution
        Read count
        Classification using e.g. Kraken or something.
            Uses loads of RAM
        Maybe a different script?
            Reference seeker? If no reference genome given.
                Conditional logic within snakemake?
            Mapping against a reference and getting depth
## expand statement goes at the end (bottom) of each path in the dag
rule all:
    input:
        expand('{root_dir}/{sample}/{sample}_contigs.assembly_stats.tsv', sample = todo_list, root_dir = root_dir),
        expand('{root_dir}/{sample}/kraken2/{sample}.kraken_report.txt', sample = todo_list, root_dir = root_dir),
        expand('{root_dir}/{sample}/Refseeker/{sample}.Refseeker.txt', sample = todo_list, root_dir = root_dir),
        expand('{root_dir}/{sample}/read_depth/{sample}.read_depth.txt', sample = todo_list, root_dir = root_dir)

rule assembly_stats:
    input:
        assembly = rules.move_shovill_output.output.final
    output:
        stats = '{root_dir}/{sample}/shovill_bbduk/{sample}_contigs.assembly_stats.tsv',
    conda:
        '../../envs/assembly_stats.yaml'
    shell:
        'assembly-stats -t {input.assembly} > {output.stats}'

rule mlst:
    input:
        assembly = rules.move_shovill_output.output.final
    output:
        mlst_results = '{root_dir}/{sample}/mlst/{sample}.mlst.tsv'
    conda:
        '../../envs/mlst.yaml'
    shell:
        'mlst --scheme senterica --nopath {input.assembly} > {output.mlst_results}'

rule sistr:
    input:
        assembly = rules.move_shovill_output.output.final
    output:
        sistr_results = '{root_dir}/{sample}/sistr/{sample}.sistr.tab'
    conda:
        '../../envs/sistr.yaml'
    shell:
        'sistr --qc -f tab -t 4 -o {output.sistr_results} {input.assembly}'

rule star_amr:
    input:
        rules.move_shovill_output.output.final
    output:
        star_amr_output = = '{root_dir}/{sample}/star_amr/{sample}.star_amr.tsv'
    conda:
        '../../envs/staramr.yaml'
    shell:
        'amrfinder -u'

rule amr_finder_plus:
    input:
        assembly = rules.move_shovill_output.output.final
    output:
        amr_finder_plus_results = '{root_dir}/{sample}/amr_finder_plus/{sample}.amr_finder_plus.tsv'
    conda:
        '../../envs/amrfinderplus.yaml'
    shell:
        'amrfinder -n {input.assembly} -O Salmonella --output {output.amr_finder_plus_results} --threads 4 --name {wildcards.sample} -d {amrfinder_db}'

#rule snippy:
#    input:
#        r1 = rules.bbduk.output.r1,
#        r2 = rules.bbduk.output.r2
#    output:
#        '{root_dir}/{sample}/snippy_bbduk/{sample}.consensus.subs.fa'
#    conda:
#        '../../envs/snippy.yaml'
#    shell:
#        'snippy --outdir {root_dir}/{wildcards.sample}/snippy_bbduk --reference {ref_genome} --R1 {input.r1} --R2 {input.r2} --cpus 8 --force --prefix {wildcards.sample}'

#rule kraken2:
#    input:
#        r1 = rules.bbduk.output.r1,
#        r2 = rules.bbduk.output.r2
#    output:
#        kraken_report = '{root_dir}/{sample}/kraken2/{sample}.kraken_report.txt'
#    threads: kraken_threads
#    conda:
#        '../../envs/kraken2.yaml'
#    shell:
#        'kraken2 --gzip-compressed --use-names --output - --db /home/ubuntu/external_tb/kraken2/database/2020.09.28/gtdb_r89_54k_kraken2_full --report {output.kraken_report} --threads {threads} --confidence 0.9 --memory-mapping --paired {input.r1} {input.r2}'
 
