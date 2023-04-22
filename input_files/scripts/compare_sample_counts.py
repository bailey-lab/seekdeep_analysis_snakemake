'''
Tests whether read counts from extraction table are the same as read counts from
summary file.

Compares read counts per sample replicate from a summary file(summed across
haplotypes and amplicons of a given summary) against read counts per sample
replicate overall (from the last column of the extraction table).
'''
import csv
overall_summaries=snakemake.input['overall_summaries']
extraction_tables=snakemake.input['extraction_tables']
compared_folder=snakemake.output['compared_folder']
reps=snakemake.params['reps']

read_headers=[f'R{number}_ReadCnt' for number in range(1, reps+1)]
file_name_headers=[f'R{number}_name' for number in range(1, reps+1)]
#print(read_headers)

def make_overall_dict(overall_summaries):
	overall_dict={}
	working_reps=set([])
	for overall_summary in overall_summaries:
		primer=overall_summary.split('popClustering/')[1].split('/')[0]
#		print(overall_summary)
		summary_reader=csv.DictReader(open(overall_summary), delimiter='\t')
		for row in summary_reader:
			sample=row['s_Sample']
			if sample not in overall_dict:
				overall_dict[sample]={}
			for name in summary_reader.fieldnames:
				if name in file_name_headers and name!='':
					working_reps.add(primer+'_'+row[name])
				if name in read_headers:
					if name not in overall_dict[sample]:
						overall_dict[sample][name]=0
					if row[name]=='':
						row[name]=0
					overall_dict[sample][name]+=int(row[name])
	return overall_dict, sorted(list(working_reps))

overall_dict, working_reps=make_overall_dict(overall_summaries)
#print('overall_dict is', overall_dict)
#print('working reps are', '\n'.join(working_reps))
