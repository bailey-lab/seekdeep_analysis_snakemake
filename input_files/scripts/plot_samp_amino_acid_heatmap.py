import plotly.express as px
#import pickle
import yaml
import pandas as pd
import math
import statistics

#main_db=pickle.load(open(snakemake.input['pickle_main_db'], 'rb'))
main_db=yaml.safe_load(open(snakemake.input['yaml_main_db']))
#pickle_sorted_parasitemia=pickle.load(open(snakemake.input['sorted_parasitemia'], 'rb'))
aa_db=yaml.safe_load(open(snakemake.input['yaml_aa_db'], 'rb'))
aa_heatmap=snakemake.output['aa_heatmap']
aa_tsv=open(snakemake.output['aa_tsv'], 'w')
#count_pickle=snakemake.output['count_pickle']

def get_fractions(main_db, aa_db):
	all_sites, count_dict,fraction_dict, amino_acids, all_samples={},{},{},[],{}
	for primer in main_db:
		for sample in main_db[primer]:
			for replicate in main_db[primer][sample]:
				all_samples.setdefault(sample, {})
				all_samples[sample].setdefault(replicate, {})
				hap_db=main_db[primer][sample][replicate]
				for haplotype in hap_db:
					count=main_db[primer][sample][replicate][haplotype]
					for amino_acid in aa_db[haplotype]:
						site=amino_acid[:-1]
						amino_acid=f'{primer}_{amino_acid}'
						amino_acids.append(amino_acid)
						all_sites.setdefault(f'{primer}_{site}', []).append(amino_acid)
						count_dict.setdefault(replicate, {})
						count_dict[replicate].setdefault(primer, {})
						count_dict[replicate][primer].setdefault(site, {})
						count_dict[replicate][primer][site].setdefault(amino_acid, 0)
						count_dict[replicate][primer][site][amino_acid]+=count
						count_dict[replicate][primer][site].setdefault('total', 0)
						count_dict[replicate][primer][site]['total']+=count
	#for replicate in count_dict:
	for sample in sorted(list(all_samples.keys())):
		for replicate in sorted(list(all_samples[sample].keys())):
			#for primer in count_dict[replicate]:
				#for site in count_dict[replicate][primer]:
			for site in all_sites:
					#"site" is e.g. ama1_162
					#total=count_dict[replicate][primer][site]['total']
				primer, primer_site=site.split('_')
				for amino_acid in all_sites[site]:
					#"amino_acid" is e.g. 162N
					#for amino_acid in count_dict[replicate][primer][site]:
					if replicate in count_dict and primer in count_dict[replicate] and primer_site in count_dict[replicate][primer]:
						if amino_acid in count_dict[replicate][primer][primer_site]:
							total=count_dict[replicate][primer][primer_site]['total']
							count=count_dict[replicate][primer][primer_site][amino_acid]
							all_samples[sample][replicate].setdefault(amino_acid, count/total)
						else: #no counts associated with this allele
							all_samples[sample][replicate].setdefault(amino_acid, 0)
					else: #either replicate has no haplotypes in this dataset or primer doesn't exist or locus doesn't exist
						all_samples[sample][replicate].setdefault(amino_acid, 'NA')
	fraction_dict={}
	for sample in all_samples:
		for replicate in all_samples[sample]:
			for amino_acid in all_samples[sample][replicate]:
				fraction_dict.setdefault(sample, {})
				fraction_dict[sample].setdefault(amino_acid, {})
				fraction_dict[sample][amino_acid].setdefault('reps', [])
				if all_samples[sample][replicate][amino_acid]!='NA':
					fraction_dict[sample][amino_acid]['reps'].append(all_samples[sample][replicate][amino_acid])
		for amino_acid in fraction_dict[sample]:
			if len(fraction_dict[sample][amino_acid]['reps'])>0:
				fraction_dict[sample][amino_acid]=statistics.mean(fraction_dict[sample][amino_acid]['reps'])
			else:
				fraction_dict[sample][amino_acid]='NA'
	return count_dict, fraction_dict, all_samples, sorted(list(set(amino_acids)))

def populate_graph_list(fraction_dict, amino_acids):
	aa_graphing_list=[]
	y_values=sorted(fraction_dict.keys())
	x_values=amino_acids
	aa_tsv.write('replicate'+'\t'+'\t'.join(x_values)+'\n')
	for sample in sorted(fraction_dict.keys()):
		sample_list=[]
		for amino_acid in amino_acids:
			sample_list.append(fraction_dict[sample][amino_acid])
		aa_graphing_list.append(sample_list)
		aa_tsv.write(sample+'\t'+'\t'.join(map(str, sample_list))+'\n')
	return aa_graphing_list, x_values, y_values

def plot_heatmap(graphing_list, x_values, y_values, x_title, y_title, count_title, output_path, width=2000, height=4000):
	#print(graphing_list)
	fig = px.imshow(graphing_list, aspect='auto', labels=dict(x=x_title, y=y_title,
	color=count_title), x=x_values, y=y_values)
	fig.update_xaxes(side="top")
	fig.update_layout(width=width, height=height, autosize=False)
	fig.write_html(output_path)

count_dict, fraction_dict, all_samples, amino_acids=get_fractions(main_db, aa_db)
aa_graphing_list, x_values, y_values=populate_graph_list(fraction_dict, amino_acids)
plot_heatmap(aa_graphing_list, x_values, y_values, 'amino acids', 'replicates', 'fraction of site', aa_heatmap)