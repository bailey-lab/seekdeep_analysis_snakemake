import plotly.express as px
import pickle
import yaml
import pandas as pd
import math
sorted_replicates=snakemake.params['sorted_reps']
main_db=yaml.safe_load(open(snakemake.input['yaml_database']))
count_tsv_file=open(snakemake.output['read_counts'], 'w')
hap_tsv_file=open(snakemake.output['hap_counts'], 'w')
count_heatmap_plot=snakemake.output['count_heatmap_plot']
hap_heatmap_plot=snakemake.output['hap_heatmap_plot']
reorganized_db={}

def summarize_haps(hap_dict):
	counts, haps=0,0
	for hap in hap_dict:
		counts+=hap_dict[hap]
		haps+=1
	return counts, haps

def reorganize_db(main_db):
	for primer in main_db:
		for sample in main_db[primer]:
			for replicate in main_db[primer][sample]:
				counts, haps=summarize_haps(main_db[primer][sample][replicate])
				if replicate not in reorganized_db:
					reorganized_db[replicate]={}
				if primer not in reorganized_db[replicate]:
					reorganized_db[replicate][primer]=counts, haps
	return reorganized_db

def populate_graph_list(reorganized_db):
	count_graphing_list=[]
	hap_graphing_list=[]
	y_values=[]
	x_values=sorted(list(main_db.keys()))
	count_tsv_file.write('replicate'+'\t'+'\t'.join(x_values)+'\n')
	hap_tsv_file.write('replicate'+'\t'+'\t'.join(x_values)+'\n')
	for replicate in sorted_replicates: #this could be from a pre-sorted snakemake yaml file
		replicate_list=replicate.split(',')
		replicate=replicate_list[0]
		count_rep_list=[]
		hap_rep_list=[]
		y_values.append('_'.join(replicate_list)) #this could be passed as additional list arguments for each replicate
		for primer in x_values:
			print('replicate is', replicate, 'primer is', primer)
			if replicate in reorganized_db and primer in reorganized_db[replicate]:
				count_rep_list.append(math.log(reorganized_db[replicate][primer][0]+1, 2))
				hap_rep_list.append(reorganized_db[replicate][primer][1])
			else:
				count_rep_list.append(math.log(1, 2))
				hap_rep_list.append(0)
		count_graphing_list.append(count_rep_list)
		hap_graphing_list.append(hap_rep_list)
		count_tsv_file.write(replicate+'\t'+'\t'.join(map(str, count_rep_list))+'\n')
		hap_tsv_file.write(replicate+'\t'+'\t'.join(map(str, hap_rep_list))+'\n')
	return count_graphing_list, hap_graphing_list, x_values, y_values

#def get_totals(graphing_list, x_values, y_values, extraction_dict):

def plot_heatmap(graphing_list, x_values, y_values, x_title, y_title, count_title, output_path, width=2000, height=4000):
#	print(graphing_list)
	fig = px.imshow(graphing_list, aspect='auto', labels=dict(x=x_title, y=y_title,
	color=count_title), x=x_values, y=y_values)
	fig.update_xaxes(side="top")
	fig.update_layout(width=width, height=height, autosize=False)
	fig.write_html(output_path)

reorganized_db=reorganize_db(main_db)
count_graphing_list, hap_graphing_list, x_values, y_values=populate_graph_list(reorganized_db)
plot_heatmap(count_graphing_list, x_values, y_values, 'amplicons', 'replicate_plus_metadata', 'log2(read cnts+1)', count_heatmap_plot)
plot_heatmap(hap_graphing_list, x_values, y_values, 'amplicons', 'replicate_plus_metadata', 'haplotype count', hap_heatmap_plot)
