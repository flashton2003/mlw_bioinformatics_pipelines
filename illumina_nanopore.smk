Illumina basic
    Trimming
    Fastqc post trimming
    Combine with MultiQC

Illumina mapping and SNP calling

    Mapping and SNP calling with PHEnix
    Output basic stats about sequencing depth, number of SNPs, number of indels, number of Ns, proportion of reads mapped, etc.
    Output consensus genome

Illumina assembly

    SKESA (?)
    Assembly stats

import os

def read_todo_list(todo_list):
    with open(todo_list) as fi:
        lines = fi.readlines()
        lines = [x.strip() for x in lines]
    return lines

#configfile: "/home/ubuntu/scripts/snakemake/configs/salmonella.yaml"
#configfile: "/home/ubuntu/.config/snakemake/salmonella_slurm/config.yaml"

todo_list = read_todo_list(config['todo_list'])
root_dir = config['root_dir']
amrfinder_db = '/home/ubuntu/hao_shigella/salmonella/reference_genomes/2021.02.10/2020-12-17.1'
qc_results_dir = config['qc_results_dir']

kraken_threads = 8

assert os.path.exists(root_dir)
if not os.path.exists(qc_results_dir):  
    os.makedirs(qc_results_dir)

## expand statement goes at the end (bottom) of each path in the dag
rule all:
    input:
        expand(['{root_dir}/{sample}/{sample}_bbduk_1.fastq.gz', '{root_dir}/{sample}/{sample}_bbduk_2.fastq.gz'], sample = todo_list, root_dir = root_dir),
        expand('{root_dir}/{sample}/mlst/{sample}.mlst.tsv', sample = todo_list, root_dir = root_dir),
        expand('{root_dir}/{sample}/sistr/{sample}.sistr.tab', sample = todo_list, root_dir = root_dir),
        expand('{root_dir}/{sample}/shovill_bbduk/{sample}_contigs.assembly_stats.tsv', sample = todo_list, root_dir = root_dir),
        expand('{root_dir}/{sample}/amr_finder_plus/{sample}.amr_finder_plus.tsv', sample = todo_list, root_dir = root_dir)

rule pre_trimming:
    input:
        ['{root_dir}/{sample}/{sample}_1.fastq.gz', '{root_dir}/{sample}/{sample}_2.fastq.gz']
    output:
        ['{root_dir}/{sample}/{sample}_1_fastqc.zip', '{root_dir}/{sample}/{sample}_2_fastqc.zip', '{root_dir}/{sample}/{sample}_1_fastqc.html', '{root_dir}/{sample}/{sample}_2_fastqc.html']
    conda:
        '../../envs/fastqc.yaml'
    shell:
        'fastqc {input}'

rule move_fastqc_output:
    input:
        rules.fastqc.output
    output:
        ['{root_dir}/{sample}/fastqc/{sample}_1_fastqc.zip', '{root_dir}/{sample}/fastqc/{sample}_2_fastqc.zip', '{root_dir}/{sample}/fastqc/{sample}_1_fastqc.html', '{root_dir}/{sample}/fastqc/{sample}_2_fastqc.html']
    run:
        cmds = zip(input, output)
        dirname = os.path.dirname(input[0])
        if not os.path.exists(dirname):
            shell('mkdir -p {dirname}')
        for c in cmds:
            print(c[0], c[1])
            shell('mv {c[0]} {c[1]}')

rule bbduk:
    input:
        r1 = '{root_dir}/{sample}/{sample}_1.fastq.gz',
        r2 = '{root_dir}/{sample}/{sample}_2.fastq.gz'

    output:
        r1 = '{root_dir}/{sample}/{sample}_bbduk_1.fastq.gz', 
        r2 = '{root_dir}/{sample}/{sample}_bbduk_2.fastq.gz'
    conda:
        '../../envs/bbmap.yaml'
    shell:
        'bbduk.sh threads=8 ref=/home/ubuntu/external_tb/references/2019.04.22/adapters.fa in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} ktrim=r k=23 mink=11 hdist=1 tbo tpe qtrim=r trimq=20 minlength=50'

rule post_trimming:
    input:
        ['{root_dir}/{sample}/{sample}_1.fastq.gz', '{root_dir}/{sample}/{sample}_2.fastq.gz']
    output:
        ['{root_dir}/{sample}/{sample}_1_fastqc.zip', '{root_dir}/{sample}/{sample}_2_fastqc.zip', '{root_dir}/{sample}/{sample}_1_fastqc.html', '{root_dir}/{sample}/{sample}_2_fastqc.html']
    conda:
        '../../envs/fastqc.yaml'
    shell:
        'fastqc {input}'
## do the expand bit in the multiqc as this is the last section which requires all these, and snakemake works by 'pulling'
rule multiqc:
    input:
        expand(['{root_dir}/{sample}/fastqc/{sample}_1_fastqc.zip', '{root_dir}/{sample}/fastqc/{sample}_2_fastqc.zip', '{root_dir}/{sample}/fastqc/{sample}_1_fastqc.html', '{root_dir}/{sample}/fastqc/{sample}_2_fastqc.html'], sample = todo_list, root_dir = root_dir)
    output:
        '{qc_results_dir}/multiqc_report.html'
    #conda:
    #    '../../envs/multiqc.yaml'
    run:
        shell('conda activate multiqc')
        shell('multiqc -o {qc_results_dir} {input}')


rule shovill:
    params:
        threads = 8,
        ram = 32
    input:
        r1 = rules.bbduk.output.r1,
        r2 = rules.bbduk.output.r2
    output:
        final = '{root_dir}/{sample}/shovill_bbduk/contigs.fa',
        graph = '{root_dir}/{sample}/shovill_bbduk/contigs.gfa',
        spades = '{root_dir}/{sample}/shovill_bbduk/spades.fasta'
    conda:
        '../../envs/shovill.yaml'
    shell:
        'shovill --outdir {root_dir}/{wildcards.sample}/shovill_bbduk -R1 {input.r1} -R2 {input.r2} --cpus {params.threads} --ram {params.ram} --force'

rule move_shovill_output:
    input:
        final = rules.shovill.output.final,
        graph = rules.shovill.output.graph,
        spades = rules.shovill.output.spades
    output:
        final = '{root_dir}/{sample}/shovill_bbduk/{sample}_contigs.fa',
        graph = '{root_dir}/{sample}/shovill_bbduk/{sample}_contigs.gfa',
        spades = '{root_dir}/{sample}/shovill_bbduk/{sample}_spades.fasta'

    shell:
        '''mv {input.final} {output.final}
        mv {input.graph} {output.graph}
        mv {input.spades} {output.spades}
        '''

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
 

