'''
Author: Allen Zhu
Affiliation: Meren Lab at the University of Chicago
Aim: A simple Snakemake workflow to automate the processing of database files into anvi-profiles.
Date: Wed Sep 27
Run: snakemake (-s Snakefile) 
'''

BASE_DIR = "/Users/allenzhu/Google Drive/[Rotation_Meren-and-Pan]/tRNA-profiles/"

SAMPLES = ['tRNA-DM-HF-M01']
''', 'tRNA-DM-HF-M02', 'tRNA-DM-HF-M03', 'tRNA-DM-HF-M01', 'tRNA-DM-HF-M02', 'tRNA-DM-HF-M03',\
           'tRNA-HF-M01', 'tRNA-HF-M02', 'tRNA-HF-M03', 'tRNA-LF-M01', 'tRNA-LF-M02', 'tRNA-LF-M03']
'''

# tRNA-seq-tools needs to be installed in order to use this snakemake file.
# prior to running this Snakefile, we had installed a virtual environment of snakemake within the virtual environment of tRNA-seq-tools.
TRNA_SEQ_TOOLS_ACTIVATE_PATH = '~/virtual-envs/tRNA-seq-tools/bin/activate'

# Print functions to test a working snakemake file.
print('First, activating tRNA-seq-tools.')
print(expand("hello_{sample}", sample=SAMPLES))

# Rules

rule all:
    input:
        expand('../workflow_snakemake/anvi_bam_{sample}.bam-ANVIO_PROFILE/PROFILE.db', sample=SAMPLES)

rule db_to_fasta:
    input:
        '../tRNA-db-profiles/{sample}.db'
    output:
        full_length_fasta = '../workflow_snakemake/full_length_reads_{sample}.fa',    
        all_read_fasta = '../workflow_snakemake/all_reads_{sample}.fa'
    shell:
        'trna-get-sequences -p {input} -o {output.full_length_fasta} --full-length-only; trna-get-sequences -p {input} -o {output.all_read_fasta}'

rule reformat_fasta_for_anvi:
    input:
        '../workflow_snakemake/full_length_reads_{sample}.fa'
    output:
        '../workflow_snakemake/fixed_full_length_reads_{sample}.fa'
    shell:
        'anvi-script-reformat-fasta {input} -o {output} -l 0 --simplify-names'

rule bowtie_build_index:
    input:
        full_length_fasta = '../workflow_snakemake/fixed_full_length_reads_{sample}.fa',
        all_read_fasta = '../workflow_snakemake/all_reads_{sample}.fa'
    output:
        all_reads_sam = '../workflow_snakemake/sam_file_{sample}.sam'
    params:
        seed_prefix = '../workflow_snakemake/seed_index_{sample}'
    threads: 4
    shell:
        'bowtie2-build {input[0]} {params.seed_prefix}; bowtie2 --threads {threads} -x {params.seed_prefix} -f {input[1]} -S {output}'

rule sam_to_bam:
    input:
        '../workflow_snakemake/sam_file_{sample}.sam'
    output:
        '../workflow_snakemake/raw_bam_{sample}.bam'
    shell:
        'samtools view -F 4 -bS {input} > {output}'

rule anvi_init_bam:
    input:
        '../workflow_snakemake/raw_bam_{sample}.bam'
    output:
        '../workflow_snakemake/anvi_bam_{sample}.bam'
    shell:
        'anvi-init-bam {input} -o {output}'


rule anvi_gen_contigs_db:
    input: 
        '../workflow_snakemake/fixed_full_length_reads_{sample}.fa'
    output:
        '../workflow_snakemake/full_length_reads_{sample}.db'
    shell:
        'anvi-gen-contigs-database -f {input} -o {output}'

rule anvi_profile:
    input:
        '../workflow_snakemake/anvi_bam_{sample}.bam',
        '../workflow_snakemake/full_length_reads_{sample}.db'
    output:
        '../workflow_snakemake/anvi_bam_{sample}.bam-ANVIO_PROFILE/PROFILE.db'
    params:
        min_contig_length = 0,
        output_dir = '../anvi-profiles',
        sample_name = 'anvi_profile_dobby'
    shell:
        'anvi-profile -i {input[0]} -c {input[1]} -M {params.min_contig_length} -W' #--output-dir {params.output_dir} --sample-name {params.sample_name}'

'''
Below, this is not yet implemented, mainly because they are optional for the snakemake workflow and can be dealt with later.
rule remove_unncessary_files_maybeoptional:
    input: 
        NAMEOFINDEXHERE
        nameofrawbamfile
        nameofsamfile
    shell:
        'rm samfile rawbamfile indexfiles'

rule anvi_run_hmms_highly_recommended:
    input:
        'tRNA-db-profiles/full_length_reads_{sample}.db'
    shell:
        'anvi-run-hmms -c {input}'

rule anvi_merge:
I need to ask Meren something about anvi_merge.
'''