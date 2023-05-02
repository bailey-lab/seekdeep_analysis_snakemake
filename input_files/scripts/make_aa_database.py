'''
creates a database of amino acid types for every haplotype of every primer of a seekdeep run.
'''
import csv
import pickle
import yaml
import gzip
import os

snakemake_folder=snakemake.input['seekdeep_results']
#database_pickle=snakemake.output['pickle_aa_db']
database_yaml=snakemake.output['yaml_aa_db']
seekdeep_subdir=snakemake.params['seekdeep_subdir']

sample_file=f'{snakemake_folder}/{seekdeep_subdir}/info/sampNames.tab.txt'
sites_file=snakemake.input['sites_of_interest']


primer_list=[line.strip().split()[0] for line in open(sample_file)][1:]

def make_sites_dict(sites_file):
	'''
	makes a dictionary of amino acid sites of interest (so we can filter out
	site calls from seekdeep that we aren't interested in)
	'''
	sites_dict={}
	for line_number, line in enumerate(open(sites_file)):
		line=line.strip().split('\t')
		if line_number>0:
			sites_dict.setdefault(line[3], []).append(line[7])
	return sites_dict

def make_aa_db(primer_list, sites_dict):
	aa_dict={}
	for primer in primer_list:
		zip_test=f'{snakemake_folder}/{seekdeep_subdir}/popClustering/{primer}/analysis/selectedClustersInfo.tab.txt.gz'
		unzip_test=f'{snakemake_folder}/{seekdeep_subdir}/popClustering/{primer}/analysis/selectedClustersInfo.tab.txt'
		if os.path.isfile(zip_test):
			file_handle=gzip.open(zip_test, mode='rt')
		elif os.path.isfile(unzip_test):
			file_handle=open(unzip_test)
		else:
			continue
		summary_reader=csv.DictReader(file_handle, delimiter='\t')
		for row in summary_reader:
			haplotype=row['h_popUID']
			AA_unparsed=row['h_AATyped'].split('--')
			if len(AA_unparsed)>1:
				aa_sites=AA_unparsed[1].split(':')
			else:
				aa_sites=[]
			final_list=[]
			for aa in aa_sites:
				pos=aa[:-1]
#				if pos in sites_dict:
#					final_list.append(aa)
				final_list.append(aa)
			aa_dict[haplotype]=final_list
	return aa_dict

sites_dict=make_sites_dict(sites_file)
aa_dict=make_aa_db(primer_list, sites_dict)
#pickle.dump(aa_dict, open(database_pickle, 'wb'))
yaml.dump(aa_dict, open(database_yaml, 'w'))
