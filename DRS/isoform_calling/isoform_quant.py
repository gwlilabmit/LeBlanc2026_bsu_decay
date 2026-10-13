import polars as pl
# from multiprocessing import Process, Queue
import multiprocessing as mp
import sys

def split_annotations(isoform_df, chunk_size): # use chunk size of 400 for 12 threads
	res = []
	new_dict = {}
	for row in isoform_df.iter_rows(named=True):
		if len(new_dict) < chunk_size:
			new_dict[row['transcript_id']] = row
		else:
			res.append(new_dict)
			new_dict = {row['transcript_id']:row}

	res.append(new_dict)
	return res

def map_reads(q, reads_df, annotation_chunk, strand): # need to set up to work with multiprocessing

	isoform_ids = annotation_chunk.keys()
	for isoform in isoform_ids:
		info = annotation_chunk.get(isoform)
		start_list = list(map(int, info['associated_starts'].split(',')))
		end_list = list(map(int, info['associated_ends'].split(',')))

		start_min = min(start_list) # furthest left starting point associated with isoform
		start_max = max(start_list)

		end_min = min(end_list)
		end_max = max(end_list) # furthest right ending point associated with isoform

		start_left_bound = start_min - 3 # wiggle room for counting
		start_right_bound = start_max + 3
		
		end_left_bound = end_min - 3
		end_right_bound = end_max + 3

		start_range = range(start_left_bound, start_right_bound+1)
		end_range = range(end_left_bound, end_right_bound+1)

		filtered = reads_df.filter((pl.col('strand') == strand) & (pl.col('leftmost_coordinate').is_in(start_range)) & (pl.col('rightmost_coordinate').is_in(end_range)))

		count = len(filtered)
		info['count'] = count
		out = {isoform: info}
		q.put(out)

	q.put(None)
	print('chunk complete')

def mp_map_reads(annotation_dictionary_list, reads_df, strand):
	if strand == 'f':
		strand = '+'
	elif strand == 'r':
		strand = '-'

	q=mp.get_context('spawn').Queue()
	process_list = []
	for i in range(len(annotation_dictionary_list)):
		p=mp.get_context('spawn').Process(target=map_reads, args=(q, reads_df, annotation_dictionary_list[i], strand))
		p.start()
		process_list.append(p)
		print(i)
		sys.stdout.flush()

	output = {}
	count = 0
	while any([mp.Process.is_alive(p) for p in process_list]) == True:
		out=q.get()
		if type(out) == dict:
			output = output | out
		elif out == None:
			count+=1
		if count == len(annotation_dictionary_list):
			for p in process_list:
				p.join()
			print('all genes calculated')

	return output

def main():

	isoform_df = pl.read_csv('gene_matching/all_timepoint_gene_matches.txt', separator='\t')
	f_isoforms = isoform_df.filter(pl.col('strand') == 'f')
	r_isoforms = isoform_df.filter(pl.col('strand') == 'r')

	f_isoforms_list = split_annotations(f_isoforms, 400)
	r_isoforms_list = split_annotations(r_isoforms, 400)

	timepoint = sys.argv[1]

	id_dict = {'t0':'2578', 't05':'2579-all', 't2':'2581-all', 't4':'2582', 't8':'2583', 't12':'2584', 't16': '2585'}
	sample_id = id_dict[timepoint]
	infile = '../streamlined_sam_files/R25-'+sample_id+'.txt'
	timepoint_df = pl.read_csv(infile, separator='\t')

	f_output = mp_map_reads(f_isoforms_list, timepoint_df, '+')
	r_output = mp_map_reads(r_isoforms_list, timepoint_df, '-')

	annotation_cols = {
		'transcript_id':pl.Int64, 
		'strand':pl.String,
		'covered_genes':pl.String, 
		'overlapping_genes':pl.String, 
		'associated_starts':pl.String, 
		'associated_ends':pl.String, 
		'isoform_ids':pl.String,
		'timepoints':pl.String,
		'timepoint_count':pl.Int64,
		'count':pl.Int64}

	output = f_output | r_output

	out_df = pl.DataFrame(schema=annotation_cols)

	for i in output.keys():
		row = output.get(i)
		pl_row = pl.from_dict(row, schema=annotation_cols)
		out_df = pl.concat([out_df, pl_row], how='vertical')

	out_df = out_df.sort(by='transcript_id')

	out_df.write_csv('isoform_counts_by_timepoint/'+timepoint+'_iso_counts.txt', separator='\t', include_header=True)

	print(timepoint, 'done')

	#timepoint_list = ['t0', 't05', 't2', 't4', 't8', 't12', 't16']
	#id_dict = {'t0':'2578', 't05':'2579-all', 't2':'2581-all', 't4':'2582', 't8':'2583', 't12':'2584', 't16': '2585'}
	# timepoint_dic = {}
	# for i, timepoint in enumerate(timepoint_list):
	# 	print('reading ' + timepoint)
	# 	sample_id = sample_id_list[i]
	# 	infile = '../streamlined_sam_files/R25-'+sample_id+'.txt'
	# 	timepoint_df = pl.read_csv(infile, separator='\t')
	# 	timepoint_dic[timepoint] = timepoint_df

	# test_list = [f_isoforms_list[0]]
	# test_df = timepoint_dic['t0']
	# output = mp_map_reads(test_list, test_df, '+')
	# print(output)
	# sys.stdout.flush()

	# annotation_cols = {
	# 	'transcript_id':pl.Int64, 
	# 	'strand':pl.String,
	# 	'covered_genes':pl.String, 
	# 	'overlapping_genes':pl.String, 
	# 	'associated_starts':pl.String, 
	# 	'associated_ends':pl.String, 
	# 	'isoform_ids':pl.String,
	# 	'timepoints':pl.String,
	# 	'timepoint_count':pl.Int64,
	# 	'count':pl.Int64}

	# for timepoint in timepoint_dic.keys():
	# 	print(timepoint, 'running')
	# 	df = timepoint_dic.get(timepoint)
	# 	f_output = mp_map_reads(f_isoforms_list, df, '+')
	# 	r_output = mp_map_reads(r_isoforms_list, df, '-')

		# output = f_output | r_output

		# # output = f_output

		# out_df = pl.DataFrame(schema=annotation_cols)

		# for i in output.keys():
		# 	row = output.get(i)
		# 	pl_row = pl.from_dict(row, schema=annotation_cols)

		# 	out_df = pl.concat([out_df, pl_row], how='vertical')

		# out_df = out_df.sort(by='transcript_id')

		# out_df.write_csv('isoform_counts_by_timepoint/'+timepoint+'_iso_counts.txt', separator='\t', include_header=True)

		# print(timepoint, 'done')

if __name__ == '__main__':
	main()