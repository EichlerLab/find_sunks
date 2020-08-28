import glob
import re
import os
import sys
import pandas as pd

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.prefix('source %s/env.cfg; ')

manifest_df = pd.read_csv('sunk_manifest.tab', sep='\t', header=None, names=['sample', 'asm'], index_col=0)

def findRef(wildcards):
	return manifest_df.at[wildcards.sample, 'asm']

rule all:
	input:
		expand("sunks/{sample}.unique_kmers.fasta", sample=manifest_df.index)

rule run_jellyfish:
	input: 
		ref = findRef
	output: 
		counts = temp('jellyfish/{sample}.mer_counts.jf')
	params: 
		sge_opts="-pe serial 8 -l mfree=4G"
	resources:
		mem=4
	threads: 8
	shell:
		'''
		jellyfish count -m 30 -s 400000000 -t {threads} -L 1 -U 1 -o {output.counts} {input.ref}
		'''


rule run_jellyfish_dump:
	input: 
		count = rules.run_jellyfish.output.counts
	output: 
		tab = temp("jellyfish/{sample}.unique_kmers.tab")
	params: 
		sge_opts="-l mfree=20G"
	resources:
		mem=32
	threads: 1
	shell:
		'''
		jellyfish dump -L 1 -U 1 -c -t -o {output.tab} {input.count}
		'''

rule sunk_convert_kmers_into_fasta:
	input: 
		tab = rules.run_jellyfish_dump.output.tab
	output: 
		sunks = "sunks/{sample}.unique_kmers.fasta"
	resources:
		mem = 8
	threads: 1
	shell:
		'''
		awk '{{ print \">\"$1\"\\n\"$1 }}' {input.tab} > {output.sunks}
		'''

