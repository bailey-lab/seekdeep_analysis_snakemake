'''
creates a database of read counts for every haplotype of every replicate of
every sample of every primer of a seekdeep run.
'''
import csv
import pickle
import gzip
import yaml

snakemake_folder=snakemake.input['seekdeep_results']
#database_pickle=snakemake.output['pickle_database']
database_yaml=snakemake.output['yaml_database']
sample_file=f'{snakemake_folder}/sampleNames.tab.txt'
zipped=snakemake.params['zipped_status']
seekdeep_subdir=snakemake.params['seekdeep_subdir']

#def get_extraction_totals:
	

def make_empty_db(sample_file):
	empty_db={}
	for line_number, line in enumerate(open(sample_file)):
		if line_number>0:
			line=line.strip().split()
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
		if zipped:
			overall_summary=f'{snakemake_folder}/{seekdeep_subdir}/popClustering/{primer}/analysis/selectedClustersInfo.tab.txt.gz'
			file_handle=gzip.open(overall_summary, mode='rt')
		else:
			overall_summary=f'{snakemake_folder}/{seekdeep_subdir}/popClustering/{primer}/analysis/selectedClustersInfo.tab.txt'
			file_handle=open(overall_summary)
		summary_reader=csv.DictReader(file_handle, delimiter='\t')
		for row in summary_reader:
			sample=row['s_Sample']
			haplotype=row['h_popUID']
			for rep_number in range(1, len(main_db[primer][sample])+1):
				rep_nm_header=f'R{rep_number}_name'
				rep_cnt_header=f'R{rep_number}_ReadCnt'
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
