'''
creates a database of read counts for every haplotype of every replicate of
every sample of every primer of a seekdeep run.
'''
import csv
import pickle
import gzip
import yaml
import os

snakemake_folder=snakemake.input['seekdeep_results']
seekdeep_subdir=snakemake.params['seekdeep_subdir']
#database_pickle=snakemake.output['pickle_database']
database_yaml=snakemake.output['yaml_database']
sample_file=f'{snakemake_folder}/{seekdeep_subdir}/info/sampNames.tab.txt'

#def get_extraction_totals:
	

def make_empty_db(sample_file):
	empty_db={}
	for line_number, line in enumerate(open(sample_file)):
		line=line.strip().split()
		for entry_number, entry in enumerate(line):
			if entry.startswith('MID'):
				line[entry_number]=entry[3:]
		primer, sample, reps=line[0], line[1], line[2:]
		if primer not in empty_db:
			empty_db[primer]={}
		if sample not in empty_db[primer]:
			empty_db[primer][sample]={}
		for rep in reps:
			if rep not in empty_db[primer][sample]:
				empty_db[primer][sample][rep]={}
	return empty_db

def populate_db(main_db, snakemake_folder):
	for primer in main_db:
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
			sample=row['s_Sample']
			haplotype=row['h_popUID']
			for rep_number in range(1, len(main_db[primer][sample])+1):
				rep_nm_header=f'R{rep_number}_name'
				rep_cnt_header=f'R{rep_number}_ReadCnt'
				rep_name=''
				if rep_nm_header in row:
					rep_name=row[rep_nm_header][3:]
					count=row[rep_cnt_header]
#				print(primer, sample, rep_name, haplotype)
				if rep_name!='' and haplotype!='':
					if haplotype not in main_db[primer][sample][rep_name]:
						main_db[primer][sample][rep_name][haplotype]=0
					if count!='':
						main_db[primer][sample][rep_name][haplotype]+=int(count)
	return main_db

main_db=make_empty_db(sample_file)
main_db=populate_db(main_db, snakemake_folder)
yaml.dump(main_db, open(database_yaml, 'w'))
#pickle.dump(main_db, open(database_pickle, 'wb'))
