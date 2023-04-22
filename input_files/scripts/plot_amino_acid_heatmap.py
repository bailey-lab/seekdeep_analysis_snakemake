import plotly.express as px
#import pickle
import yaml
import pandas as pd
import math

#main_db=pickle.load(open(snakemake.input['pickle_main_db'], 'rb'))
main_db=yaml.safe_load(open(snakemake.input['yaml_main_db']))
#pickle_sorted_parasitemia=pickle.load(open(snakemake.input['sorted_parasitemia'], 'rb'))
yaml_sorted_parasitemia=yaml.safe_load(open(snakemake.input['yaml_sorted_parasitemia']))
aa_db=yaml.safe_load(open(snakemake.input['yaml_aa_db'], 'rb'))
aa_heatmap=snakemake.output['aa_heatmap']
aa_tsv=open(snakemake.output['aa_tsv'], 'w')
#count_pickle=snakemake.output['count_pickle']
count_yaml=snakemake.output['count_yaml']
count_table=snakemake.output['count_table']

def get_fractions(main_db, aa_db):
	count_dict,fraction_dict={},{}
	for primer in main_db:
		for sample in main_db[primer]:
			for replicate in main_db[primer][sample]:
				hap_db=main_db[primer][sample][replicate]
				for haplotype in hap_db:
					count=main_db[primer][sample][replicate][haplotype]
					for amino_acid in aa_db[haplotype]:
						site=amino_acid[:-1]
						amino_acid=f'{primer}_{amino_acid}'
						count_dict.setdefault(replicate, {})
						count_dict[replicate].setdefault(primer, {})
						count_dict[replicate][primer].setdefault(site, {})
						count_dict[replicate][primer][site].setdefault(amino_acid, 0)
						count_dict[replicate][primer][site][amino_acid]+=count
						count_dict[replicate][primer][site].setdefault('total', 0)
						count_dict[replicate][primer][site]['total']+=count
	amino_acids=set([])
	for replicate in count_dict:
		for primer in count_dict[replicate]:
			for site in count_dict[replicate][primer]:
				total=count_dict[replicate][primer][site]['total']
				for amino_acid in count_dict[replicate][primer][site]:
					if amino_acid!='total':
						amino_acids.add(amino_acid)
						count=count_dict[replicate][primer][site][amino_acid]
						fraction_dict.setdefault(replicate, {})
						fraction_dict[replicate].setdefault(amino_acid, count/total)
	return count_dict, fraction_dict, sorted(list(amino_acids))

def populate_graph_list(fraction_dict, sorted_parasitemia, amino_acids):
	aa_graphing_list=[]
	y_values=sorted_parasitemia
	x_values=amino_acids
	aa_tsv.write('replicate'+'\t'+'\t'.join(x_values)+'\n')
	for replicate in sorted_parasitemia:
		replicate_list=[]
		sample_replicate=replicate.split('_')[0]
		for amino_acid in amino_acids:
			if sample_replicate in fraction_dict and amino_acid in fraction_dict[sample_replicate]:
				replicate_list.append(fraction_dict[sample_replicate][amino_acid])
			else:
				replicate_list.append(0)
		aa_graphing_list.append(replicate_list)
		aa_tsv.write(replicate+'\t'+'\t'.join(map(str, replicate_list))+'\n')
	return aa_graphing_list, x_values, y_values

def plot_heatmap(graphing_list, x_values, y_values, x_title, y_title, count_title, output_path, width=2000, height=4000):
	print(graphing_list)
	fig = px.imshow(graphing_list, aspect='auto', labels=dict(x=x_title, y=y_title,
	color=count_title), x=x_values, y=y_values)
	fig.update_xaxes(side="top")
	fig.update_layout(width=width, height=height, autosize=False)
	fig.write_html(output_path)

def make_count_table(count_dict, sorted_parasitemia, output_file):
	'''
	plots the counts associated with every amino acid of every gene in every
	replicate. Replicates are rows, with columns for gene name, amino acid
	position, mutation, mutation read count, and total read count
	'''
	output_file.write(f'rep\tgene\tposition\tmutation\tmutation_count\ttotal_count\n')
	for rep in sorted_parasitemia:
		rep=rep.split('_')[0]
		if rep in count_dict:
			for gene in count_dict[rep]:
				for pos in count_dict[rep][gene]:
					total=count_dict[rep][gene][pos]['total']
					for mut in count_dict[rep][gene][pos]:
						if mut!='total':
							count=count_dict[rep][gene][pos][mut]
							output_file.write(f'{rep}\t{gene}\t{pos}\t{mut}\t{count}\t{total}\n')
count_dict, fraction_dict, amino_acids=get_fractions(main_db, aa_db)
make_count_table(count_dict, yaml_sorted_parasitemia, open(count_table, 'w'))
aa_graphing_list, x_values, y_values=populate_graph_list(fraction_dict, yaml_sorted_parasitemia, amino_acids)
plot_heatmap(aa_graphing_list, x_values, y_values, 'amino acids', 'replicates', 'fraction of site', aa_heatmap)
#pickle.dump(count_dict, open(count_pickle, 'wb'))
yaml.dump(count_dict, open(count_yaml, 'w'))