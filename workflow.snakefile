import yaml
import os
import pandas as pd
import re

####	Set up config

configfile: "config.yaml"

input_samples=config["sample_fastqs"]


####	End of config


### Process

rule all:
	input:
		expand("results/{sample}/fastq_process/{sample}_r1_processed.barcoded.fastq.gz", sample=input_samples.keys()),
		expand("results/{sample}/fastq_process/{sample}_r2_processed.barcoded.fastq.gz", sample=input_samples.keys()),
		expand("results/{sample}/{sample}_frags.sorted.bed.gz", sample=input_samples.keys()),
		expand("results/{sample}/{sample}_frags.sorted.bed.gz.tbi", sample=input_samples.keys()),
		expand("results/{sample}/align/{sample}.sorted.dedup.bam", sample=input_samples.keys()),
		expand("results/{sample}/align/{sample}.sorted.dedup.bam.bai", sample=input_samples.keys()),
		expand("results/{sample}/macs2_callpeak/{sample}_peaks.broadPeak", sample=input_samples.keys()),
		expand("results/{sample}/raw_peak_bc_matrix/barcodes.tsv", sample=input_samples.keys()),
		expand("results/{sample}/raw_peak_bc_matrix/peaks.bed", sample=input_samples.keys()),
		expand("results/{sample}/raw_peak_bc_matrix/matrix.mtx", sample=input_samples.keys()),
		expand("results/{sample}/filtered_peak_bc_matrix/barcodes.tsv", sample=input_samples.keys()),
		expand("results/{sample}/filtered_peak_bc_matrix/peaks.bed", sample=input_samples.keys()),
		expand("results/{sample}/filtered_peak_bc_matrix/matrix.mtx", sample=input_samples.keys()),
#		expand("results/{sample}/{sample}.sorted.dedup.bam", sample=input_samples.keys()),
#		expand("results/{sample}/{sample}_edit_distance.tsv", sample=input_samples.keys()),
#		expand("results/{sample}/{sample}_per_umi_per_position.tsv", sample=input_samples.keys()),
#		expand("results/{sample}/{sample}_per_umi.tsv", sample=input_samples.keys()),
		expand("results/{sample}/fastq_process/fastqc/{sample}_r1_trimmed_fastqc.html", sample=input_samples.keys()),
		expand("results/{sample}/fastq_process/fastqc/{sample}_r2_trimmed_fastqc.html", sample=input_samples.keys()),
		expand("results/{sample}/fastq_process/cutadapt/{sample}_r1_discarded.fastq.gz", sample=input_samples.keys()),
		expand("results/{sample}/fastq_process/cutadapt/{sample}_r2_discarded.fastq.gz", sample=input_samples.keys()),


"""
Trim fastqs to remove degenerate sequences (e.g. homopolymers)
Filter out any read pairs where R1 is trimmed, as these reads are missing UMI and/or spatial barcode
Note: This may not actually be necessary if we don't care about UMI
"""
rule fastq_trim:
	input:
		read1 = lambda wildcards: config['sample_fastqs'][wildcards.sample]['R1'],
		read2 = lambda wildcards: config['sample_fastqs'][wildcards.sample]['R2']
	output:
		read1_trimmed = "results/{sample}/fastq_process/{sample}_r1_trimmed.fastq.gz",
		read2_trimmed = "results/{sample}/fastq_process/{sample}_r2_trimmed.fastq.gz",
		read1_discarded = "results/{sample}/fastq_process/cutadapt/{sample}_r1_discarded.fastq.gz",
		read2_discarded = "results/{sample}/fastq_process/cutadapt/{sample}_r2_discarded.fastq.gz",
		discarded_read_names = "results/{sample}/fastq_process/cutadapt/discarded_reads.lst"
	message:
		"Trimming R1 and R2"
	log:
		"logs/{sample}/fastq_trim.log"
	threads: 8
	run:
#		shell("cutadapt --minimum-length 28 --cores={threads} -b 'AAAAAA;o=6' -b 'TTTTTT;o=6' -b 'CCCCCC;o=6' -b 'GGGGGG;o=6' -B 'AAAAAA;o=6' -B 'TTTTTT;o=6' -B 'CCCCCC;o=6' -B 'GGGGGG;o=6' --times 15 -o {output.read1_trimmed} -p {output.read2_trimmed} {input.read1} {input.read2} &> {log}")
		shell("cutadapt --minimum-length 28 --cores={threads} -b 'AAAAAAAAAAAAAA;o=12' -b 'CCCCCCCCCCCCCC;o=12' -b 'GGGGGGGGGGGGGG;o=12' -B 'AAAAAAAAAAAAAA;o=12' -B 'CCCCCCCCCCCCCC;o=12' -B 'GGGGGGGGGGGGGG;o=12' -q 15,10 --times 15 -o {output.read1_trimmed} -p {output.read2_trimmed} {input.read1} {input.read2} &> {log}")
		shell("comm -23 <( zcat {input.read1} | paste - - - - | grep -Po '^@\K\S+' | sort ) <( zcat {output.read1_trimmed} | paste - - - - | grep -Po '^@\K\S+' | sort ) > {output.discarded_read_names}")
		shell("seqtk subseq {input.read1} {output.discarded_read_names} | gzip > {output.read1_discarded}")
		shell("seqtk subseq {input.read2} {output.discarded_read_names} | gzip > {output.read2_discarded}")

"""
Generate fastqc report for trimmed reads
"""
rule fastqc:
	input:
		read1_trimmed = "results/{sample}/fastq_process/{sample}_r1_trimmed.fastq.gz",
		read2_trimmed = "results/{sample}/fastq_process/{sample}_r2_trimmed.fastq.gz"
	output:
		read1_report = "results/{sample}/fastq_process/fastqc/{sample}_r1_trimmed_fastqc.html",
		read2_report = "results/{sample}/fastq_process/fastqc/{sample}_r2_trimmed_fastqc.html"
	message:
		"Generating FastQC read quality report"
	log:
		"logs/{sample}/fastqc.log"
	run:
		shell("fastqc -o results/{wildcards.sample}/fastq_process/fastqc -f fastq {input.read1_trimmed} {input.read2_trimmed} &> {log}")

"""
Filter spatial barcodes and perform error correction by previously calculated whitelist
"""
rule process_umi_bc:
	input:
		read1 = "results/{sample}/fastq_process/{sample}_r1_trimmed.fastq.gz",
		read2 = "results/{sample}/fastq_process/{sample}_r2_trimmed.fastq.gz"
	output:
		r1_temp = temporary("results/{sample}/fastq_process/{sample}_r1_processed.barcoded.fastq"),
		r1_fastq="results/{sample}/fastq_process/{sample}_r1_processed.barcoded.fastq.gz",
		r2_fastq="results/{sample}/fastq_process/{sample}_r2_processed.barcoded.fastq.gz"
	message:
		"Filtering barcodes and performing UMI deduplication"
	params:
	# Uses lambda function to defer evaluation of wildcards; for some reason snakemake thinks there's a wildcard in the params
		whitelist=lambda wildcards: config["whitelist"],
		r1_format=lambda wildcards: config["r1_format"]
	log: 
		umitools_log="logs/{sample}/process_umi_bc_umitools.log",
	run:
		shell("zcat {input.read1} | umi_tools extract --extract-method=regex --whitelist={params.whitelist} --bc-pattern='{params.r1_format}' --read2-in={input.read2} --read2-out={output.r2_fastq} --log {log.umitools_log} >> {output.r1_temp} 2>/dev/null")
		shell("gzip -c {output.r1_temp} > {output.r1_fastq}")



"""
Align with bowtie2
Move UMI and spatial barcode from read name to tags in the bam file
We expect the read name to end in _{16bp_spatial_barcode}_{12bp_umi}.
"""
rule align:
	input:
		r1_fastq = "results/{sample}/fastq_process/{sample}_r1_processed.barcoded.fastq.gz",
		r2_fastq = "results/{sample}/fastq_process/{sample}_r2_processed.barcoded.fastq.gz"
	output:
		nobc_bam = temporary("results/{sample}/align/{sample}_nobc.sorted.bam"),
		nobc_index = temporary("results/{sample}/align/{sample}_nobc.sorted.bam.bai"),
		unfiltered_bam = temporary("results/{sample}/align/{sample}_unfiltered.sorted.bam"),
		unfiltered_index= temporary("results/{sample}/align/{sample}_unfiltered.sorted.bam.bai"),
		bam = "results/{sample}/align/{sample}.sorted.bam",
		index = "results/{sample}/align/{sample}.sorted.bam.bai"
	message:
		"Aligning to specified reference genome with bowtie2"
	params:
		ref=config["ref"],
		mode=config['mode'],
		remove_artifacts=config['remove_artifacts'],
		remove_artifacts_method=config['remove_artifacts_method'],
		blacklist_file=config['blacklist_file']
	log: 
		align = "logs/{sample}/align_bowtie2.log",
		nametotag = "logs/{sample}/align_nametotag.log"
	threads: 16
	run:
#		shell("bowtie2 --very-sensitive -q --phred33 --end-to-end -t -x {params.ref} -U {input.fastq} -p {threads} | samtools view -bS - | samtools sort - -o {output.nobc_bam} -@ {threads} &>{log.align}")
		if params.mode=='se':
			shell("bowtie2 --very-sensitive-local -q --phred33 --local -t -x {params.ref} -U {input.r2_fastq} -p {threads} | samtools view -bS - | samtools sort - -o {output.nobc_bam} -@ {threads} &>{log.align}")
		elif params.mode=='pe':
			shell("bowtie2 --very-sensitive-local -q --phred33 --local -t -x {params.ref} -1 {input.r1_fastq} -2 {input.r2_fastq} -p {threads} | samtools view -bS - | samtools sort - -o {output.nobc_bam} -@ {threads} &>{log.align}")
		shell("samtools index {output.nobc_bam} -o {output.nobc_index} &> {log.align}")
		if params.remove_artifacts and pararms.remove_artifacts_method == 'blacklist':
			shell("bedtools intersect -abam {output.nobc_bam} -b {params.blacklist_file} -v > {output.unfiltered_bam}")
			shell("samtools index {output.unfiltered_bam} -o {output.unfiltered_index}")
			shell("sinto nametotag -b {output.unfiltered_bam} --barcode_regex '(?<=\_)[ACTG]*(?=\_)' --tag CB | sinto nametotag -b - --barcode_regex '[ACTG]*$' --tag UM | samtools view -b - > {output.bam} 2> {log.nametotag}")
		else:	
			shell("samtools view -h {output.nobc_bam} | sinto nametotag -b - --barcode_regex '(?<=\_)[ACTG]*(?=\_)' --tag CB | sinto nametotag -b - --barcode_regex '[ACTG]*$' --tag UM | samtools view -b - > {output.bam} 2> {log.nametotag}")
		shell("samtools index {output.bam} -o {output.index}")

"""
Also perform read deduplication
Note: We prefer to use genomic position rather than UMI for deduplication as some reads
were likely amplified by IVT before UMIs were attached
"""
rule dedup:
	input:
		bam = "results/{sample}/align/{sample}.sorted.bam",
		index = "results/{sample}/align/{sample}.sorted.bam.bai"
	output:
#		statsfile1 = "results/{sample}/{sample}_edit_distance.tsv",
#		statsfile2 = "results/{sample}/{sample}_per_umi_per_position.tsv",
#		statsfile3 = "results/{sample}/{sample}_per_umi.tsv",
		dedup_bam = "results/{sample}/align/{sample}.sorted.dedup.bam",
		dedup_index = "results/{sample}/align/{sample}.sorted.dedup.bam.bai"
	message:
		"Performing UMI deduplication with umi_tools and bam file barcode assignment"
	log: 
		umitools_log = "logs/{sample}/dedup_umitools.log"
	params:
		use_umi = config['umi_dedup']
	run:
		if params.use_umi:
			shell("umi_tools dedup --stdin={input.bam} --log={log.umitools_log} --extract-umi-method=tag --umi-tag=UM --cell-tag=CB --per-cell --output-stats=results/{wildcards.sample}/{wildcards.sample} > {output.dedup_bam}")
		else:
			shell("umi_tools dedup --stdin={input.bam} --log={log.umitools_log} --per-cell --ignore-umi --extract-umi-method=tag --umi-tag=UM --cell-tag=CB > {output.dedup_bam}")
		shell("samtools index {output.dedup_bam} -o {output.dedup_index}")


"""
Note: our data is single end so we fake a fragments file by just converting our bam file to a correctly formatted tsv
according to https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments

We should also consider using e.g. MACS2 estimated fragment length to generate a more correct fragments file
"""
rule build_fragments_file:
	input:
		bam = "results/{sample}/align/{sample}.sorted.dedup.bam"
	output:
		frag_file_unsorted = temporary("results/{sample}/{sample}_frags.bed"),
		frag_file_sorted = temporary("results/{sample}/{sample}_frags.sorted.bed"),
		frag_file_sorted_gzip = "results/{sample}/{sample}_frags.sorted.bed.gz",
		frag_file_sorted_index = "results/{sample}/{sample}_frags.sorted.bed.gz.tbi"
	log: "logs/{sample}/build_fragments_file.log"
	message: "Creating fragments file"
	run:
		shell("python ./scripts/create_fragments.py {input.bam} {output.frag_file_unsorted} 2>{log}")
		shell("sort -k1,1 -k2,2n {output.frag_file_unsorted} 1>{output.frag_file_sorted} 2>{log}")
		shell("bgzip {output.frag_file_sorted} -c 1>{output.frag_file_sorted_gzip} 2>{log}")
		shell("tabix -p bed {output.frag_file_sorted_gzip} 2>{log}")


"""
Peak calling with MACS2
Note: macs2 estimates fragment size(?) - we might want to think about using that
"""
rule call_peaks:
	input:
		bam = "results/{sample}/align/{sample}.sorted.dedup.bam"
	output:
		peaks = "results/{sample}/macs2_callpeak/{sample}_peaks.broadPeak"
	message: "Calling peaks with MACS2"
	log: "logs/{sample}/call_peaks.log"
	params:
		outdir = "results/{sample}/macs2_callpeak",
		fdr = config['macs2_fdr'],
		genome_size = config['macs2_genomesize'],
		name = '{sample}'
	shell:
		"""	
		macs2 callpeak --SPMR -B -q {params.fdr} --keep-dup all \
           -g {params.genome_size} -f BAM --slocal 5000 --llocal 50000 --broad --max-gap 1000 \
           -t {input.bam} --outdir {params.outdir} -n {params.name} > {log} 2>&1
        """

rule create_feature_mtx:
	input:
		frag_file = "results/{sample}/{sample}_frags.sorted.bed.gz",
		peaks = "results/{sample}/macs2_callpeak/{sample}_peaks.broadPeak"
	output:
		"results/{sample}/raw_peak_bc_matrix/barcodes.tsv",
		"results/{sample}/raw_peak_bc_matrix/peaks.bed",
		"results/{sample}/raw_peak_bc_matrix/matrix.mtx"
	message: "Creating feature by spot matrix"
	log: "logs/{sample}/create_feature_mtx.log"
	params:
		output_dir = "results/{sample}/raw_peak_bc_matrix",
		remove_artifact = config["remove_artifacts"],
		remove_artifact_method = config["remove_artifacts_method"],
		run_type = config["run_type"],
		genome_fasta = config["genome_fasta"],
		window = config["window_length"],
		p_len = config["poly_a_t_threshold"]
	run:
		if params.remove_artifact and params.remove_artifact_method == 'poly_a_t':
			print('Removing artifactual peaks...')
			if not os.path.exists(params.output_dir):
				os.mkdir(params.output_dir)
			shell("Rscript ./scripts/process_artifacts.R {params.genome_fasta} {input.peaks} {params.run_type} {params.window} {params.p_len} {params.output_dir} &> {log}")
			if (params.run_type == 'debug'):
				shell("Rscript ./scripts/create_feature_mtx.R {input.frag_file} {input.peaks} {params.output_dir} &> {log}")
			elif (params.run_type == 'processing'):
				shell("Rscript ./scripts/create_feature_mtx.R {input.frag_file} results/{wildcards.sample}/raw_peak_bc_matrix/filtered_peaks.broadPeak {params.output_dir} &> {log}")
		else:
			shell("Rscript ./scripts/create_feature_mtx.R {input.frag_file} {input.peaks} {params.output_dir} &> {log}")


rule filter_feature_mtx:
	input:
		barcodes = "results/{sample}/raw_peak_bc_matrix/barcodes.tsv",
		peaks = "results/{sample}/raw_peak_bc_matrix/peaks.bed",
		matrix = "results/{sample}/raw_peak_bc_matrix/matrix.mtx",
	output:
		barcodes = "results/{sample}/filtered_peak_bc_matrix/barcodes.tsv",
		peaks = "results/{sample}/filtered_peak_bc_matrix/peaks.bed",
		matrix = "results/{sample}/filtered_peak_bc_matrix/matrix.mtx"
	message: "Filtering spots by aligned tissue image"
	log: "logs/{sample}/filter_feature_mtx.log"
	params:
		json = lambda wildcards: config["sample_alignment_jsons"][wildcards.sample],
		output_dir="results/{sample}/filtered_peak_bc_matrix",
		spot_coords=config["spot_coords"]
	run:
		if (params.json != 'none'):
			shell('Rscript ./scripts/filter_feature_mtx.R {params.json} {input.barcodes} {input.peaks} {input.matrix} {params.spot_coords} {params.output_dir} &> {log}')
		else:
			shell('python -u ./scripts/align_image.py')
			shell('Rscript ./scripts/filter_feature_mtx.R results/{sample}/fiducial.json {input.barcodes} {input.peaks} {input.matrix} {params.spot_coords} {params.output_dir} &> {log}')