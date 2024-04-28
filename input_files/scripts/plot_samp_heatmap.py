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
			reps, count_sum, hap_sum=0,0,0
			for replicate in main_db[primer][sample]:
				counts, haps=summarize_haps(main_db[primer][sample][replicate])
				reps+=1
				count_sum+=counts
				hap_sum+=haps
				if sample not in reorganized_db:
					reorganized_db[sample]={}
				if primer not in reorganized_db[sample]:
					reorganized_db[sample][primer]=count_sum, reps, hap_sum
	return reorganized_db

def populate_graph_list(reorganized_db):
	count_graphing_list=[]
	hap_graphing_list=[]
	y_values=[]
	x_values=sorted(list(main_db.keys()))
	count_tsv_file.write('sample'+'\t'+'\t'.join(x_values)+'\n')
	hap_tsv_file.write('sample'+'\t'+'\t'.join(x_values)+'\n')
	for sample in sorted(list(reorganized_db.keys())):
		count_samp_list=[]
		hap_samp_list=[]
		y_values.append(sample)
		for primer in x_values:
			counts, reps, haps=reorganized_db[sample][primer]
			print('sample is', sample, 'primer is', primer)
			count_samp_list.append(math.log(count_sum/reps+1, 2))
			hap_samp_list.append(hap_sum/reps)
		count_graphing_list.append(count_samp_list)
		hap_graphing_list.append(hap_samp_list)
		count_tsv_file.write(sample+'\t'+'\t'.join(map(str, count_samp_list))+'\n')
		hap_tsv_file.write(replicate+'\t'+'\t'.join(map(str, hap_samp_list))+'\n')
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
plot_heatmap(count_graphing_list, x_values, y_values, 'amplicons', 'samples', 'log2(read cnts+1)', count_heatmap_plot)
plot_heatmap(hap_graphing_list, x_values, y_values, 'amplicons', 'samples', 'haplotype count', hap_heatmap_plot)
