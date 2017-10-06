'''
    Author: Allen Zhu
    Affiliation: Meren Lab at the University of Chicago
    Aim: A simple Snakemake pipeline to automate the processing of reading tRNA databases into anvi'o.
    Date of Creation: Wed Sep 27, 2017
    Last Updated: Thursday Oct 5, 2017

    This is a snakemake file for a tRNA metagenomics workflow in the Meren Lab.
    It uses tRNA-seq-tools, bowtie, and anvi'o.

    Here is an outline of the following steps:
        Convert sample database (.db) files to FASTA (.fa) files with all the reads, and a full-length reference contig FASTA file.
        Reformatting the FASTA files for anvi'o.
        Mapping the sequences to the assemblies using bowtie2.
        Generating an anvi'o contigs database. (and run hmm profile)
        Generating anvi'o profile database.

    The following must be in the working directory:
        config.json - contains the configuration information necessary for the pipeline to run specified input.

        samples.txt - a TAB-delimited file to dsecribe where the samples are.
            The header line should be "sample", "file". 
            The first column should have sample name. The second column should have file path. 

    Run Example: 
        snakemake (-s Snakefile) 
'''

import os
import pandas as pd


# The configuration file with the essential configurations for the workflow
configfile: "config.json"
# Setting the names of the directories
dir_list = ["LOGS_DIR", "FASTA_DIR", "BOWTIE_DIR", "BAM_DIR", "ANVI_BAM_DIR", "CONTIGS_DB_DIR", "PROFILE_DIR", "MERGE_DIR"]
dir_names = ["00_LOGS", "01_QC", "02_BOWTIE", "03_BAM", "04_ANVIBAM", "05_CONTIGSDB", "06_ANVIO_PROFILE", "07_MERGED"]
dirs_dict = dict(zip(dir_list, dir_names))

os.makedirs(dirs_dict["LOGS_DIR"], exist_ok=True)

def A(_list, d, default_value = ""):
    '''
        A helper function to make sense of config details.
        string_list is a list of strings (or a single string)
        d is a dictionary

        this function checks if the strings in x are nested values in y.
        For example if x = ['a','b','c'] then this function checkes if the
        value y['a']['b']['c'] exists, if it does then it is returned

        You give the function a _list. If it's not a list, it converts it into a list.
        While the list is not empty, it pops out the first element. If the first element is in the dictionary d provided in the input,
        the dictionry is updated to that first element. Then it pops out the next element to see if that next element is embedded in the first element.
        Works well with dictionaries from the .json file or embedded lists. 
        If the final element is found, it will return the valuue corresponding to the key in that diciontary.
        -Allen Zhu
    '''
    if type(_list) is not list:
        # converting to list for the cases of only one item
        _list = [_list]
    while _list:
        a = _list.pop(0)
        if a in d:
            d = d[a]
        else:
            return default_value
    return d

def T(rule_name, N=1): return A([rule_name,'threads'], config, default_value=N)
'''
    A helper function to get the user-defined number of threads for a rule.
    This function takes in a rule name, which gets put in a list along with 'threads'. It then uses the A helper function preceding this function.
    The A helper function checks if the rule_name is in the dictionary, and then looks for a 'thread' key inside this rule_name dictionary.
    -Allen Zhu
'''

# for key in A("output_dirs", config):
#     if key not in dir_list:
#         raise ConfigError("You define the name for directory '%s' in your "\
#                         "config file, but the only available folders are: "\
#                          "%s" % (key, dir_list))

    # this updates the 
#    dirs_dict[key] = A(key, config["output_dirs"])


# Load sample file, with default file "samples.txt"
samples_txt_file = A("samples_txt", config, default_value="samples.txt")
samples_info = pd.read_csv(samples_txt_file, sep='\t', index_col=False)
SAMPLES = list(samples_info['sample'])
FILES = list(samples_info['file'])
MIN_CONTIG_LENGTH = A(['anvi_profile', 'min_contig_length'], config, default_value="0")
#group_names = sample_names


BASE_DIR = "/Users/allenzhu/Google Drive/[Rotation_Meren-and-Pan]/tRNA-profiles/"

# Print functions to test a working snakemake file.
print('First, activating tRNA-seq-tools.')
print(expand("hello_{sample}", sample=SAMPLES))

# Rules

rule all:
    input:
        expand(dirs_dict["PROFILE_DIR"] + "/{sample}/PROFILE.db", sample=SAMPLES)
        #expand('../workflow_snakemake/anvi_bam_{sample}.bam-ANVIO_PROFILE/PROFILE.db', sample=SAMPLES)

rule db_to_fasta:
    input:
        expand("{file}", file=FILES)
    output:
        full_length_fasta = dirs_dict["FASTA_DIR"] + '/full_length_reads_{sample}.fa',    
        all_read_fasta = dirs_dict["FASTA_DIR"] + '/all_reads_{sample}.fa'
    shell:
        'trna-get-sequences -p {input} -o {output.full_length_fasta} --full-length-only; trna-get-sequences -p {input} -o {output.all_read_fasta}'

rule reformat_fasta_for_anvi:
    input:
        rules.db_to_fasta.output.full_length_fasta
        #'../workflow_snakemake/full_length_reads_{sample}.fa'
    output:
        dirs_dict["FASTA_DIR"] + '/fixed_full_length_reads_{sample}.fa'
        #'../workflow_snakemake/fixed_full_length_reads_{sample}.fa'
    shell:
        'anvi-script-reformat-fasta {input} -o {output} -l 0 --simplify-names'

rule bowtie_build_index:
    input:
        rules.reformat_fasta_for_anvi.output, 
        rules.db_to_fasta.output.all_read_fasta
        #full_length_fasta = '../workflow_snakemake/fixed_full_length_reads_{sample}.fa',
        #all_read_fasta = '../workflow_snakemake/all_reads_{sample}.fa'
    output:
        #all_reads_sam = 
        dirs_dict['BOWTIE_DIR'] + '/sam_file_{sample}.sam'
    params:
        seed_prefix = dirs_dict['BOWTIE_DIR'] + '/seed_index_{sample}'
    threads: 4
    shell:
        'bowtie2-build {input[0]} {params.seed_prefix}; bowtie2 --threads {threads} -x {params.seed_prefix} -f {input[1]} -S {output}'

rule sam_to_bam:
    input:
        rules.bowtie_build_index.output
        #'../workflow_snakemake/sam_file_{sample}.sam'
    output:
        dirs_dict['BAM_DIR'] + '/raw_bam_{sample}.bam'
        #'../workflow_snakemake/raw_bam_{sample}.bam'
    shell:
        'samtools view -F 4 -bS {input} > {output}'

rule anvi_init_bam:
    input:
        rules.sam_to_bam.output
        #'../workflow_snakemake/raw_bam_{sample}.bam'
    output:
        dirs_dict['ANVI_BAM_DIR'] + '/anvi_bam_{sample}.bam'
        #'../workflow_snakemake/anvi_bam_{sample}.bam'
    shell:
        'anvi-init-bam {input} -o {output}'


rule anvi_gen_contigs_db:
    input: 
        rules.reformat_fasta_for_anvi.output
        #'../workflow_snakemake/fixed_full_length_reads_{sample}.fa'
    output:
        dirs_dict['CONTIGS_DB_DIR'] + '/full_length_reads_{sample}.db'
        #'../workflow_snakemake/full_length_reads_{sample}.db'
    shell:
        'anvi-gen-contigs-database -f {input} -o {output}'

rule anvi_profile:
    input:
        rules.anvi_init_bam.output, 
        rules.anvi_gen_contigs_db.output
        #'../workflow_snakemake/anvi_bam_{sample}.bam',
        #'../workflow_snakemake/full_length_reads_{sample}.db'
    output:
        dirs_dict['PROFILE_DIR'] + '/{sample}/PROFILE.db'
        #'../workflow_snakemake/anvi_bam_{sample}.bam-ANVIO_PROFILE/PROFILE.db'
    params:
        min_contig_length = MIN_CONTIG_LENGTH,
        output_dir = dirs_dict['PROFILE_DIR'] + '/{sample}',
        sample_name = 'anvi_profile_dobby'
    shell:
        'anvi-profile -i {input[0]} -c {input[1]} -M {params.min_contig_length} -W --output-dir {params.output_dir}' #--sample-name {params.sample_name}'

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