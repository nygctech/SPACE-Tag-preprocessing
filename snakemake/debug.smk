rule all_debug_output:
  input:
    expand("results/{sample}/debug_plots/peak_nucleotide_freq.pdf", sample=input_samples.keys()),
    expand("results/{sample}/debug_plots/artifact_peak_frequency_plot.pdf", sample=input_samples.keys())


# Generate nucleotide frequency plot and artifact peak frequency plot
rule debug_artifact_peaks:
  input:
    expand("results/{sample}/macs2_callpeak/{sample}_peaks.narrowPeak")
  output:
    expand("results/{sample}/debug_plots/peak_nucleotide_freq.pdf", sample=input_samples.keys()),
    expand("results/{sample}/debug_plots/artifact_peak_frequency_plot.pdf", sample=input_samples.keys())
  message:
		"Generating debug plots for artifactual pA/pT peaks"
  params:
    genome_fasta = config["genome_fasta"],
    window = config["window_length"],
		p_len = config["poly_a_t_threshold"],
    output_dir = = "results/{sample}/debug_plots"
	log: 
		"logs/{sample}/debug_artifact_peaks.log"
	threads: 8
	run:
    shell("Rscript ./scripts/process_artifacts.R {params.genome_fasta} {input.peaks} not_debug {params.window} {params.p_len} {params.output_dir} &> {log}")

# Given a reference peak set or sample name, 
def get_tornado_ref(param):
  ref_peaks = param
  # If param is a sample name, return the macs2 called peaks
  # Else we assume it is the path to a peak file
  if param in sample_fastqs.keys():
    ref_peaks = f"results/{param}/macs2_callpeak/{param}_peaks.narrowPeak"
  return ref_peaks

rule generate_tornado_plot:
  input:
    expand()
  output:
    "results/tornado_plot.pdf"
  message:
    "Generating tornado plot"
  params:
    reference_peaks = get_tornado_ref(config["tornado_plot_reference"]),
    files = lambda wildcards 
  log:
    "logs/tornado_plot.log"
  threads: 4
  run:
    shell("./scripts/create_tornado_plot.sh {params.reference_peaks}")
  
  
  