configfile: 'config_illumina.yaml'

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
    # conda:
    #     '../../envs/fastqc.yaml'
    shell:
        '''
        conda activate fastqc
        fastqc {input}
        '''

rule move_fastqc_output:
    input:
        rules.pre_trimming.output
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
    # conda:
    #     '../../envs/bbmap.yaml'
    shell:
        '''
        conda activate bbduk
        bbduk.sh threads=8 ref=/home/ubuntu/external_tb/references/2019.04.22/adapters.fa in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} ktrim=r k=23 mink=11 hdist=1 tbo tpe qtrim=r trimq=20 minlength=50
        '''

rule post_trimming:
    input:
        ['{root_dir}/{sample}/{sample}_1.fastq.gz', '{root_dir}/{sample}/{sample}_2.fastq.gz']
    output:
        ['{root_dir}/{sample}/{sample}_1_fastqc.zip', '{root_dir}/{sample}/{sample}_2_fastqc.zip', '{root_dir}/{sample}/{sample}_1_fastqc.html', '{root_dir}/{sample}/{sample}_2_fastqc.html']
    conda:
        '../../envs/fastqc.yaml'
    shell:
        '''
        conda activate fastqc
        fastqc {input}
        '''
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


rule phenix:
    params:
        threads = 8
    input:
        r1 = rules.bbduk.output.r1,
        r2 = rules.bbduk.output.r2,
        ref = '...'
    output:
        '{root_dir}/{sample}/phenix'
    # conda:
    #     '../../envs/phenix.yaml'
    shell:
        '''
        conda activate phenix
        phenix.py run_snp_pipeline.py -r1 {input.r1} -r2 {input.r2} -r {input.ref} --sample-name {wildcards.sample} --mapper bwa --variant gatk --filters min_depth:5,mq_score:30
        '''

rule skesa:
    input:
        r1 = rules.bbduk.output.r1,
        r2 = rules.bbduk.output.r2
    output:
        '{root_dir}/{sample}/skesa/{sample}_skesa.fa'
    # conda:
    #     '../../envs/skesa.yaml'
    shell:
        '''
        conda activate skesa
        skesa --reads {input.r1},{input.r2} --cores 4 --memory 48 > {output}
        '''

rule assembly_stats:
    input:
        rules.skesa.output
    output:
        '{root_dir}/{sample}/skesa/{sample}_contigs.assembly_stats.tsv'
    conda:
        '../../envs/assembly_stats.yaml'
    shell:
        '''
        conda activate assembly-stat
        assembly-stats -t {input} > {output}
        '''



