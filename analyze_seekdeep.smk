configfile: 'input_files/analyze_seekdeep.yaml'

output_root='output_files/'+config['output_folder']+'/'

rule all:
	input:
		heatmap_plot=output_root+config['output_folder']+'_read_heatmap.html',
		hap_heatmap_plot=output_root+config['output_folder']+'_COI_heatmap.html',
		config_duplicate=output_root+'config.yaml',
		yaml_database=output_root+'read_dict.yaml',
		aa_heatmap=output_root+config['output_folder']+'_amino_acid_heatmap.html'

rule copy_files:
	input:
		original_config='analyze_seekdeep.yaml',
		original_snakemake='analyze_seekdeep.smk',
		original_scripts='input_files/scripts'
	output:
		config_duplicate=output_root+'copied_config_files/analyze_seekdeep.yaml',
		snakemake_duplicate=output_root+'copied_config_files/analyze_seekdeep.smk',
		scripts_duplicate=output_root+'copied_config_files/scripts'
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
		zipped_status=config['zipped_status'],
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

rule plot_heatmap:
	'''
	plots a heatmap and a table with read counts
	'''
	input:
		yaml_database=output_root+'read_dict.yaml',
		parasitemia_path='input_files/sample_parasitemia.csv'
	output:
		count_heatmap_plot=output_root+config['output_folder']+'_read_heatmap.html',
		hap_heatmap_plot=output_root+config['output_folder']+'_COI_heatmap.html',
		read_counts=output_root+config['output_folder']+'_read_counts.tsv',
		hap_counts=output_root+config['output_folder']+'_hap_counts.tsv',
		yaml_sorted_parasitemia=output_root+config['output_folder']+'_sorted_parasitemia.yaml',
#		pickle_sorted_parasitemia=output_root+config['output_folder']+'_sorted_parasitemia.pkl'
	script:
		'input_files/scripts/plot_heatmap.py'

rule make_aa_database:
	'''
	examines the final summary table to figure out prevalences of individual
	amino acids for each sample, amplicon, and rep.
	'''
	input:
		seekdeep_results=config['skdp_folder'],
		sites_of_interest='input_files/sites_of_interest.tsv'
	params:
		zipped_status=config['zipped_status'],
		seekdeep_subdir=config['seekdeep_subdir']
	output:
#		pickle_aa_db=output_root+'aa_dict.pkl',
		yaml_aa_db=output_root+'aa_dict.yaml'
	script:
		'input_files/scripts/make_aa_database.py'

rule plot_aa_freqs:
	input:
#		pickle_aa_db=output_root+'aa_dict.pkl',
		yaml_aa_db=output_root+'aa_dict.yaml',
		yaml_main_db=output_root+'read_dict.yaml',
		#pickle_sorted_parasitemia=output_root+config['output_folder']+'_sorted_parasitemia.pkl'
		yaml_sorted_parasitemia=output_root+config['output_folder']+'_sorted_parasitemia.yaml'
	output:
		aa_heatmap=output_root+config['output_folder']+'_amino_acid_heatmap.html',
		aa_tsv=output_root+config['output_folder']+'_amino_acid_fracs.tsv',
		#count_pickle=output_root+config['output_folder']+'_amino_acid_counts.pkl',
		count_yaml=output_root+config['output_folder']+'_amino_acid_counts.yaml',
		count_table=output_root+config['output_folder']+'_count_table.tsv'
	script:
		'input_files/scripts/plot_amino_acid_heatmap.py'
