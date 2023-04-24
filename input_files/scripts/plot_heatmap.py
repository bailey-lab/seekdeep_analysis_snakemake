import plotly.express as px
import pickle
import yaml
import pandas as pd
import math
#main_db=pickle.load(open(snakemake.input['pickle_database'], 'rb'))
main_db=yaml.safe_load(open(snakemake.input['yaml_database']))
count_tsv_file=open(snakemake.output['read_counts'], 'w')
hap_tsv_file=open(snakemake.output['hap_counts'], 'w')
parasitemia_path=snakemake.input['parasitemia_path']
count_heatmap_plot=snakemake.output['count_heatmap_plot']
hap_heatmap_plot=snakemake.output['hap_heatmap_plot']
#pickle_sorted_parasitemia_path=snakemake.output['pickle_sorted_parasitemia']
yaml_sorted_parasitemia_path=snakemake.output['yaml_sorted_parasitemia']
reorganized_db={}

def calculate_replicate_parasitemia(parasitemia_path):
	'''
	This function doesn't calculate replicate parasitemia as much as parse it
	from an input file.
	'''
	para_dict={}
	for line_number, line in enumerate(open(parasitemia_path)):
		if line_number>0:
			line=line.strip().split(',')
			to_parse, sample=line[1].split('_')
			level=to_parse[:3]
			rep=to_parse[-1:]
			replicate=f'{level}-{sample}-Rep{rep}'
			para_dict[replicate]=float(line[-1])
	return para_dict

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
	para_dict=calculate_replicate_parasitemia(parasitemia_path)
	parasitemia_sorted=sorted([[para_dict[replicate], replicate] for replicate in reorganized_db])
	sorted_replicates=[replicate[1] for replicate in parasitemia_sorted]
	for replicate in sorted_replicates: #this could be from a pre-sorted snakemake yaml file
		count_rep_list=[]
		hap_rep_list=[]
		y_values.append(replicate+'_'+str(para_dict[replicate])) #this could be passed as additional list arguments for each replicate
		for primer in x_values:
			count_rep_list.append(math.log(reorganized_db[replicate][primer][0]+1, 2))
			hap_rep_list.append(reorganized_db[replicate][primer][1])
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
plot_heatmap(count_graphing_list, x_values, y_values, 'amplicons', 'replicate_parasitemia', 'log2(read cnts+1)', count_heatmap_plot)
plot_heatmap(hap_graphing_list, x_values, y_values, 'amplicons', 'replicate_parasitemia', 'haplotype count', hap_heatmap_plot)

#pickle.dump(y_values, open(pickle_sorted_parasitemia_path, 'wb'))
yaml.dump(y_values, open(yaml_sorted_parasitemia_path, 'w'))
