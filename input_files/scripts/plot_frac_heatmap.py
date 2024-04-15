'''
based on plot_heatmap, but this version plots read counts and read fractions for
individual haplotypes (rather than the sum across all haplotypes).
'''
import plotly.express as px
import pickle
import yaml
import pandas as pd
import math
sorted_replicates=snakemake.params['sorted_reps']
main_db=yaml.safe_load(open(snakemake.input['yaml_database']))
count_tsv_file=open(snakemake.output.read_counts, 'w')
frac_tsv_file=open(snakemake.output.frac_counts, 'w')
count_heatmap_plot=snakemake.output.count_heatmap_plot
frac_heatmap_plot=snakemake.output.frac_heatmap_plot

def reorganize_db(main_db):
	reorganized_db={}
	for primer in main_db:
		for sample in main_db[primer]:
			for replicate in main_db[primer][sample]:
				total=0
				for hap in main_db[primer][sample][replicate]:
					total+=main_db[primer][sample][replicate][hap]
				for hap in main_db[primer][sample][replicate]:
					count=main_db[primer][sample][replicate][hap]
					new_replicate=replicate+'_'+primer+'_'+hap
					if new_replicate not in reorganized_db:
						reorganized_db[new_replicate]={}
					if primer not in reorganized_db[new_replicate]:
						reorganized_db[new_replicate][primer]=[count, round(count/total, 3)]
	return reorganized_db

def populate_graph_list(reorganized_db):
	count_graphing_list=[]
	frac_graphing_list=[]
	y_values=[]
	x_values=sorted(list(main_db.keys()))
	count_tsv_file.write('replicate_hap'+'\t'+'\t'.join(x_values)+'\n')
	frac_tsv_file.write('replicate_hap'+'\t'+'\t'.join(x_values)+'\n')
	for replicate in sorted(list(reorganized_db.keys())): #this could be from a pre-sorted snakemake yaml file
		count_rep_list=[]
		frac_rep_list=[]
		y_values.append(replicate) #this could be passed as additional list arguments for each replicate
		for primer in x_values:
			print('replicate is', replicate, 'primer is', primer)
			if replicate in reorganized_db and primer in reorganized_db[replicate]:
#				count_rep_list.append(math.log(reorganized_db[replicate][primer][0]+1, 2))
				count_rep_list.append(reorganized_db[replicate][primer][0])
				frac_rep_list.append(reorganized_db[replicate][primer][1])
			else:
				count_rep_list.append(0)
				frac_rep_list.append(0)
		count_graphing_list.append(count_rep_list)
		frac_graphing_list.append(frac_rep_list)
		count_tsv_file.write(replicate+'\t'+'\t'.join(map(str, count_rep_list))+'\n')
		frac_tsv_file.write(replicate+'\t'+'\t'.join(map(str, frac_rep_list))+'\n')
	return count_graphing_list, frac_graphing_list, x_values, y_values

#def get_totals(graphing_list, x_values, y_values, extraction_dict):

def plot_heatmap(graphing_list, x_values, y_values, x_title, y_title, count_title, output_path, width=2000, height=4000):
#	print(graphing_list)
	fig = px.imshow(graphing_list, aspect='auto', labels=dict(x=x_title, y=y_title,
	color=count_title), x=x_values, y=y_values)
	fig.update_xaxes(side="top")
	fig.update_layout(width=width, height=height, autosize=False)
	fig.write_html(output_path)

reorganized_db=reorganize_db(main_db)
count_graphing_list, frac_graphing_list, x_values, y_values=populate_graph_list(reorganized_db)
plot_heatmap(count_graphing_list, x_values, y_values, 'amplicons', 'replicate_plus_hap', 'log2(read cnts+1)', count_heatmap_plot)
plot_heatmap(frac_graphing_list, x_values, y_values, 'amplicons', 'replicate_plus_hap', 'read fraction', frac_heatmap_plot)
