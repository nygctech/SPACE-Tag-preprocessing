# Dummy rule to generate all output files
rule all_spatial_output:
  input:
    expand("results/{sample}/raw_peak_bc_matrix/barcodes.tsv", sample=input_samples.keys()),
		expand("results/{sample}/raw_peak_bc_matrix/peaks.bed", sample=input_samples.keys()),
		expand("results/{sample}/raw_peak_bc_matrix/matrix.mtx", sample=input_samples.keys()),
		expand("results/{sample}/filtered_peak_bc_matrix/barcodes.tsv", sample=input_samples.keys()),
		expand("results/{sample}/filtered_peak_bc_matrix/peaks.bed", sample=input_samples.keys()),
		expand("results/{sample}/filtered_peak_bc_matrix/matrix.mtx", sample=input_samples.keys()),
    expand("results/{sample}/{sample}_frags.sorted.bed.gz", sample=input_samples.keys()),
		expand("results/{sample}/{sample}_frags.sorted.bed.gz.tbi", sample=input_samples.keys())

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

rule create_feature_mtx:
	input:
		frag_file = "results/{sample}/{sample}_frags.sorted.bed.gz",
		peaks = "results/{sample}/macs2_callpeak/{sample}_peaks.narrowPeak"
	output:
		"results/{sample}/raw_peak_bc_matrix/barcodes.tsv",
		"results/{sample}/raw_peak_bc_matrix/peaks.bed",
		"results/{sample}/raw_peak_bc_matrix/matrix.mtx"
	message: "Creating feature by spot matrix"
	log: "logs/{sample}/create_feature_mtx.log"
	params:
		output_dir = "results/{sample}/raw_peak_bc_matrix",
		remove_artifact = config["remove_artifacts"],
		genome_fasta = config["genome_fasta"],
		window = config["window_length"],
		p_len = config["poly_a_t_threshold"]
	run:
		if (params.remove_artifacts):
			shell("Rscript ./scripts/process_artifacts.R {params.genome_fasta} {input.peaks} not_debug {params.window} {params.p_len} {params.output_dir} &> {log}")
      shell("Rscript ./scripts/create_feature_mtx.R {input.frag_file} results/{wildcards.sample}/raw_peak_bc_matrix/filtered_peaks.narrowPeak {params.output_dir} &> {log}")
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

