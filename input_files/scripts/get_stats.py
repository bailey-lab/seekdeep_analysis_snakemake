'''
options:
a primer can fail for all samples
a sample can fail for all primers (completely fail)
a sample can fail for a subset of primers if all replicates fail for that primer
a replicate of a sample can fail for all primers (completely fail)
a replicate of a sample can fail for a subset of primers
'''

#import pickle
import yaml
#db_pickle=snakemake.input['pickle_database']
db_yaml=snakemake.input['yaml_database']
main_db=yaml.safe_load(open(db_yaml))
#main_db=pickle.load(open(db_pickle, 'rb'))


empty_primers=open(snakemake.output['empty_primers'], 'w')
empty_samples=open(snakemake.output['empty_samples'], 'w')
empty_replicates=open(snakemake.output['empty_reps'], 'w')

failed_primers=[]
failed_samples={}
failed_replicates={}

for primer in main_db:
	empty_primer=True
	for sample in main_db[primer]:
		empty_sample=True
		for replicate in main_db[primer][sample]:
			if len(main_db[primer][sample][replicate])==0:
				if sample not in failed_replicates:
					failed_replicates[sample]={}
				if replicate not in failed_replicates[sample]:
					failed_replicates[sample][replicate]=[]
				failed_replicates[sample][replicate].append(primer)
			else:
				empty_primer=False
				empty_sample=False
		if empty_sample:
			if sample not in failed_samples:
				failed_samples[sample]=[]
			failed_samples[sample].append(primer)
	if empty_primer:
		failed_primers.append(primer)

empty_primers.write('failed_primers\n')
for primer in failed_primers:
	empty_primers.write(primer+'\n')

empty_samples.write('sample\tfailed_primers\n')
for sample in sorted(list(failed_samples)):
	sample_primers=failed_samples[sample]
	if len(sample_primers)==len(main_db):
		empty_samples.write(f'{sample}\tcompletely_failed\n')
	else:
		empty_samples.write(f'{sample}\t{",".join(sample_primers)}\n')

empty_replicates.write('replicate\tfailed_primers\n')
for sample in sorted(list(failed_replicates)):
	for replicate in sorted(list(failed_replicates[sample])):
		replicate_primers=failed_replicates[sample][replicate]
		if len(replicate_primers)==len(main_db):
			empty_replicates.write(f'{replicate}\tcompletely_failed\n')
		else:
			empty_replicates.write(f'{replicate}\t{",".join(replicate_primers)}\n')
