configfile: 'analyze_seekdeep.yaml'

output_root=config['output_folder']+'/'

rule all:
	input:
		samp_heatmap_plot=output_root+'samp_read_heatmap.html',
		samp_hap_heatmap_plot=output_root+'samp_COI_heatmap.html',
		rep_heatmap_plot=output_root+'rep_read_heatmap.html',
		rep_hap_heatmap_plot=output_root+'rep_COI_heatmap.html',
		config_duplicate=output_root+'copied_config_files/analyze_seekdeep.yaml',
		yaml_database=output_root+'read_dict.yaml',
		final_file=output_root+config['final_output']

rule copy_files:
	input:
		original_config='analyze_seekdeep.yaml',
		original_snakemake='analyze_seekdeep.smk',
		original_scripts='input_files/scripts'
	output:
		config_duplicate=output_root+'copied_config_files/analyze_seekdeep.yaml',
		snakemake_duplicate=output_root+'copied_config_files/analyze_seekdeep.smk',
		scripts_duplicate=directory(output_root+'copied_config_files/scripts')
	shell:
		'''
		cp {input.original_config} {output.config_duplicate}
		cp {input.original_snakemake} {output.snakemake_duplicate}
		cp -r {input.original_scripts} {output.scripts_duplicate}
		'''

rule make_database:
	'''
	examines the final summary table to figure out how many reads exist for each
	sample, amplicon, and rep.
	'''
	input:
		seekdeep_results=config['skdp_folder']
	params:
		seekdeep_subdir=config['seekdeep_subdir']
	output:
		yaml_database=output_root+'read_dict.yaml'
	script:
		'input_files/scripts/make_database.py'

rule get_stats:
	'''
	gets failure stats for primers. Options include:
	primer failed for all samples, sample failed for all primers, sample failed
	for a subset of primers, replicate of a sample failed for all primers,
	replicate of a sample failed for a subset of primers
	'''
	input:
		yaml_database=output_root+'read_dict.yaml'
	output:
		empty_primers=output_root+'empty_items/empty_primers.txt',
		empty_reps=output_root+'empty_items/empty_reps.txt',
		empty_samples=output_root+'empty_items/empty_samples.txt'
	script:
		'input_files/scripts/get_stats.py'

rule plot_rep_heatmap:
	'''
	plots a heatmap and a table with read counts
	'''
	input:
		yaml_database=output_root+'read_dict.yaml',
	params:
		sorted_reps=expand('{reps}', reps=config['sorted_reps'])
	output:
		count_heatmap_plot=output_root+'rep_read_heatmap.html',
		hap_heatmap_plot=output_root+'rep_COI_heatmap.html',
		read_counts=output_root+'rep_read_counts.tsv',
		hap_counts=output_root+'rep_hap_counts.tsv',
#		pickle_sorted_parasitemia=output_root+config['output_folder']+'_sorted_parasitemia.pkl'
	script:
		'input_files/scripts/plot_rep_heatmap.py'

rule plot_samp_heatmap:
	'''
	plots a heatmap and a table with read counts
	'''
	input:
		yaml_database=output_root+'read_dict.yaml',
	output:
		count_heatmap_plot=output_root+'samp_read_heatmap.html',
		hap_heatmap_plot=output_root+'samp_COI_heatmap.html',
		read_counts=output_root+'samp_read_counts.tsv',
		hap_counts=output_root+'samp_hap_counts.tsv',
	script:
		'input_files/scripts/plot_samp_heatmap.py'

rule make_aa_database:
	'''
	examines the final summary table to figure out prevalences of individual
	amino acids for each sample, amplicon, and rep.
	'''
	input:
		seekdeep_results=config['skdp_folder'],
		sites_of_interest='input_files/sites_of_interest.tsv'
	params:
		seekdeep_subdir=config['seekdeep_subdir']
	output:
#		pickle_aa_db=output_root+'aa_dict.pkl',
		yaml_aa_db=output_root+'aa_dict.yaml'
	script:
		'input_files/scripts/make_aa_database.py'

rule plot_rep_aa_freqs:
	input:
#		pickle_aa_db=output_root+'aa_dict.pkl',
		yaml_aa_db=output_root+'aa_dict.yaml',
		yaml_main_db=output_root+'read_dict.yaml',
		#pickle_sorted_parasitemia=output_root+config['output_folder']+'_sorted_parasitemia.pkl'
	params:
		sorted_reps=expand('{reps}', reps=config['sorted_reps'])
	output:
		aa_heatmap=output_root+'rep_amino_acid_heatmap.html',
		aa_tsv=output_root+'rep_amino_acid_fracs.tsv',
		#count_pickle=output_root+config['output_folder']+'_amino_acid_counts.pkl',
		count_yaml=output_root+'rep_amino_acid_counts.yaml',
		count_table=output_root+'rep_count_table.tsv'
	script:
		'input_files/scripts/plot_rep_amino_acid_heatmap.py'

rule plot_samp_aa_freqs:
	input:
#		pickle_aa_db=output_root+'aa_dict.pkl',
		aa_heatmap=output_root+'rep_amino_acid_heatmap.html',
		yaml_aa_db=output_root+'aa_dict.yaml',
		yaml_main_db=output_root+'read_dict.yaml'
		#pickle_sorted_parasitemia=output_root+config['output_folder']+'_sorted_parasitemia.pkl'
	output:
		aa_heatmap=output_root+'samp_amino_acid_heatmap.html',
		aa_tsv=output_root+'samp_amino_acid_fracs.tsv',
	script:
		'input_files/scripts/plot_samp_amino_acid_heatmap.py'
