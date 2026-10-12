import polars as pl
import pandas as pd
import sys

def read_genome(infile):
	genome={}
	with open(infile, 'r') as f:
		f.readline()
		for line in f:
			info=line.strip().split()
			start=int(info[1])
			end=int(info[2])
			strand_char=info[3]
			if strand_char == '+':
				strand = 'f'
			elif strand_char == '-':
				strand = 'r'
			name=info[5]
			type=info[6]
			genome[name]=[start, end, strand, type]
			
	return genome

def run_timepoint(annotations_df, annotation_cols, timepoint_id, timepoint_df, genome, transcript_id):

	print('starting ', timepoint_id)
	sys.stdout.flush()

	for row in timepoint_df.iter_rows(named=True):
		x = row['x']
		y = row['y']
		strand = row['strand']
		row_id = 'ID='+str(row['cluster']) + '_' + str(row['winsor']) + '_' + row['strand']

		if strand == 'f':
			strand_sign = '+'
		elif strand == 'r':
			strand_sign = '-'

		iso_range = range(x, y+1)

		genes = genome.filter(pl.col('STRAND') == strand_sign)

		## check if isoform should be added to an existing transcript ID, based on transcript that have been described so far
		## potentially problematic to add to the same df that I'm looping through, but since I'll break as soon as I add to it, I think it'll be fine?
		check = 0
		for line in annotations_df.filter(pl.col('strand') == strand).iter_rows(named=True):
			line_id = line['transcript_id']
			starts_list = line['associated_starts']
			ends_list = line['associated_ends']

			starts_distances = [abs(i - x) for i in starts_list]
			ends_distances = [abs(i - y) for i in ends_list]

			## distance of 6 or less means the two points are within each others +/- 3 windows
			check_starts = any(distance <=6 for distance in starts_distances)
			check_ends = any(distance <=6 for distance in ends_distances)

			### not super feasible to update dataframe inplace
			### should make a new row w/ new list, drop the old one and then add to the end (keep the transcript ID)
			## if you want things in order you can re-sort at the end
			if check_starts and check_ends:
				if x not in starts_list:
					starts_list.append(x)
				if y not in ends_list:
					ends_list.append(y)

				id_list = line['isoform_ids']
				id_list.append(row_id)

				timepoints_list = line['timepoints']
				timepoint_count = line['timepoint_count']
				if timepoint_id not in timepoints_list:
					timepoints_list.append(timepoint_id)
					timepoint_count+=1

				line['associated_starts'] = starts_list
				line['associated_ends'] = ends_list
				line['isoform_ids'] = id_list
				line['timepoints'] = timepoints_list
				line['timepoint_count'] = timepoint_count

				new_line = pl.DataFrame(
					{
						'transcript_id':line['transcript_id'], 
						'strand':line['strand'],
						'covered_genes':[line['covered_genes']], 
						'overlapping_genes':[line['overlapping_genes']], 
						'associated_starts':[line['associated_starts']], 
						'associated_ends':[line['associated_ends']], 
						'isoform_ids':[line['isoform_ids']],
						'timepoints':[line['timepoints']],
						'timepoint_count':line['timepoint_count']
					}, schema=annotation_cols)

				annotations_df = annotations_df.remove(pl.col('transcript_id') == line_id)
				annotations_df = pl.concat([annotations_df, new_line], how='vertical')
				check = 1
				break

		if check == 1:
			continue
		else:
			pass

		# https://nedbatchelder.com/blog/201310/range_overlap_in_two_compares.html

		if strand == 'f': ## add extra wiggle room around 5' end of read (assume that all RNAs have 10 more nt of RNA beyond the 5' end of the read)
			overlap_genes = genes.filter((iso_range.stop >= pl.col('START')) & (pl.col('END') >= (iso_range.start-10)))
			cover_genes = genes.filter((iso_range.start-10 <= pl.col('START')) & (iso_range.stop >= pl.col('END')))

		elif strand == 'r':
			overlap_genes = genes.filter((iso_range.stop+10 >= pl.col('START')) & (pl.col('END') >= iso_range.start))
			cover_genes = genes.filter((iso_range.start <= pl.col('START')) & (iso_range.stop+10 >= pl.col('END')))

		overlap_list = pl.Series(overlap_genes.select('GENE_NAME')).to_list()
		cover_list = pl.Series(cover_genes.select('GENE_NAME')).to_list()

		if len(overlap_list) == 0:
			overlap_list=['None']
		
		if len(cover_list) == 0:
			cover_list=['None']

		transcript = pl.DataFrame(
			{
				'transcript_id':transcript_id, 
				'strand':strand,
				'covered_genes':[cover_list], 
				'overlapping_genes':[overlap_list], 
				'associated_starts':[[x]], 
				'associated_ends':[[y]], 
				'isoform_ids':[[row_id]],
				'timepoints':[[timepoint_id]],
				'timepoint_count':1
			}, schema=annotation_cols)
		
		transcript_id+=1
		annotations_df = pl.concat([annotations_df, transcript], how='vertical')

	print(timepoint_id, 'done')

	return annotations_df, transcript_id

def fix_difference(list1, list2):
	diff = set(list1).symmetric_difference(set(list2))
	out = list1 + list(diff)
	return out

def second_pass_merge(annotations_df, annotation_cols): # for catching isoforms that were initially established as distinct but actually overlap with each other
	
	out_annotations = pl.DataFrame(schema=annotation_cols)
	merge_count = 0 # tracks how many merges have been done overall
	transcript_id = 0

	merged_id_list = [] # once a transcript has been merged with another one, don't consider it again

	## for each isoform
	for row in annotations_df.iter_rows(named=True):
		row_transcript_id = row['transcript_id']
		if row_transcript_id in merged_id_list:
			continue
		strand = row['strand']
		covered_genes = row['covered_genes']
		overlapping_genes = row['overlapping_genes']
		row_starts = row['associated_starts']
		row_ends = row['associated_ends']
		isoform_ids = row['isoform_ids']
		timepoints = row['timepoints']
		timepoint_count = row['timepoint_count']

		merge_check = 0 # tracks how many times this row is merged with other isoforms

		## compare to each other isoform
		for row_i in annotations_df.iter_rows(named=True):

			## don't compare isoform to itself
			transcript_id_i = row_i['transcript_id']
			if row_transcript_id == transcript_id_i:
				continue

			else:
				## get starts and ends of each other isoform
				starts_i = row_i['associated_starts']
				ends_i = row_i['associated_ends']

				## get all distances between start and end points for each isoform
				starts_distances = []
				for start in starts_i:
					distances = [abs(x - start) for x in row_starts]
					starts_distances = starts_distances + distances

				ends_distances = []
				for end in ends_i:
					distances = [abs(y - end) for y in row_ends]
					ends_distances = ends_distances + distances

				## check if there's any overlap between ends of the isoforms being compared
				check_starts = any(distance <=6 for distance in starts_distances)
				check_ends = any(distance <=6 for distance in ends_distances)

				## if the two isoforms are sufficiently close to each other on each end
				if check_starts and check_ends:
					merged_id_list.append(row_transcript_id)
					merged_id_list.append(transcript_id_i)
					# combine ends 
					row_starts = fix_difference(row_starts, starts_i)
					row_ends = fix_difference(row_ends, ends_i)

					## make sure genes linked to isoform match (may never come up)
					## concatenate if they don't already match
					if row['covered_genes'] != row_i['covered_genes']:
						covered_genes = fix_difference(row['covered_genes'], row_i['covered_genes'])
					if row['overlapping_genes'] != row_i['overlapping_genes']:
						overlapping_genes = fix_difference(row['overlapping_genes'], row_i['overlapping_genes'])

					## definitely needs to be updated, will definitely be different
					isoform_ids = fix_difference(isoform_ids, row_i['isoform_ids'])

					if row['timepoints'] != row_i['timepoints']:
						timepoints = fix_difference(row['timepoints'], row_i['timepoints'])

					timepoint_count = len(timepoints)

					## make new dataframe row

					transcript = pl.DataFrame(
						{
							'transcript_id':transcript_id, 
							'strand':strand,
							'covered_genes':[covered_genes], 
							'overlapping_genes':[overlapping_genes], 
							'associated_starts':[row_starts], 
							'associated_ends':[row_ends], 
							'isoform_ids':[isoform_ids],
							'timepoints':[timepoints],
							'timepoint_count':timepoint_count
						}, schema=annotation_cols)

					out_annotations = pl.concat([out_annotations, transcript], how='vertical')

					transcript_id+=1
					merge_check+=1
					merge_count+=1

				else: ## isoforms don't overlap
					continue

		if merge_check == 0: ## if you've looped through all the other isoforms and nothing overlapped

			## add as-is to new dataframe with updated transcript id
			transcript = pl.DataFrame(
				{
					'transcript_id':transcript_id, 
					'strand':strand,
					'covered_genes':[covered_genes], 
					'overlapping_genes':[overlapping_genes], 
					'associated_starts':[row_starts], 
					'associated_ends':[row_ends], 
					'isoform_ids':[isoform_ids],
					'timepoints':[timepoints],
					'timepoint_count':timepoint_count
				}, schema=annotation_cols)

			out_annotations = pl.concat([out_annotations, transcript], how='vertical')
			transcript_id+=1

	return out_annotations, merge_count

def main():

	print('running')
	sys.stdout.flush()
	
	timepoint_list = ['t0', 't05', 't2', 't4', 't8', 't12', 't16']
	# timepoint_list = ['t8_test', 't12_test', 't16_test']
	timepoint_dic = {}
	for timepoint in timepoint_list:
		infile = 'combined_isoforms/'+timepoint+'_combined.txt'
		# infile = 'combined_isoforms/'+timepoint+'.txt'
		df = pl.read_csv(infile, separator='\t').sort(by='x')

		timepoint_dic[timepoint] = df

	genome = read_genome('CDS_168.mochiview_v13.txt')

	genome = pl.read_csv('CDS_168.mochiview_v13.txt', separator = '\t').select(pl.col('START'), pl.col('END'), pl.col('STRAND'), pl.col('FEATURE_NAME'), pl.col('GENE_NAME'), pl.col('ALIASES'))

	annotation_cols = {
		'transcript_id':pl.Int64, 
		'strand':pl.String,
		'covered_genes':pl.List(pl.String), 
		'overlapping_genes':pl.List(pl.String), 
		'associated_starts':pl.List(pl.Int64), 
		'associated_ends':pl.List(pl.Int64), 
		'isoform_ids':pl.List(pl.String),
		'timepoints':pl.List(pl.String),
		'timepoint_count':pl.Int64}
	
	annotations = pl.DataFrame(schema=annotation_cols)

	transcript_id = 0

	for timepoint in timepoint_dic.keys():
		timepoint_df = timepoint_dic.get(timepoint)
		annotations, transcript_id = run_timepoint(annotations, annotation_cols, timepoint, timepoint_df, genome, transcript_id)

	while True:
		annotations, merge_count = second_pass_merge(annotations, annotation_cols)
		if merge_count == 0:
			break ## continue merging until you do a full pass through without merging anything

	annotations = annotations.sort(by='transcript_id')

	## conversion to pandas for saving to csv
	annotations_pd = annotations.to_pandas()

	## modifying format for more streamlined output
	annotations_pd['covered_genes'] = [','.join(map(str, l)) for l in annotations_pd['covered_genes']]
	annotations_pd['overlapping_genes'] = [','.join(map(str, l)) for l in annotations_pd['overlapping_genes']]
	annotations_pd['associated_starts'] = [','.join(map(str, l)) for l in annotations_pd['associated_starts']]
	annotations_pd['associated_ends'] = [','.join(map(str, l)) for l in annotations_pd['associated_ends']]
	annotations_pd['isoform_ids'] = [','.join(map(str, l)) for l in annotations_pd['isoform_ids']]
	annotations_pd['timepoints'] = [','.join(map(str, l)) for l in annotations_pd['timepoints']]

	annotations_pd.to_csv('gene_matching/all_timepoint_gene_matches.txt', sep='\t', index=False)
	# annotations_pd.to_csv('gene_matching/test_gene_matches.txt', sep='\t', index=False)

	# annotations.write_ndjson('test_annotations.txt')
	# annotations.write_csv('test_annotations.txt', include_header=True, separator='\t')

if __name__ == '__main__':
	main()