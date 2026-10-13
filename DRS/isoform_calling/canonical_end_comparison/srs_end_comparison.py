import polars as pl

def main():

	filtered_isoforms=pl.read_csv('../combined_isoform_quant/all_genes_filtered_by_count.txt', separator='\t')

	## primary transcript ends
	primary_ends_5=pl.read_csv('primary_5_ends.txt', separator='\t')
	primary_ends_3=pl.read_csv('primary_3_ends.txt', separator='\t')

	compare_primary_5, primary_5_count = compare_ends(filtered_isoforms, primary_ends_5, 5)
	write_matches('primary_5_matches.txt', compare_primary_5)

	compare_primary_3, primary_3_count = compare_ends(filtered_isoforms, primary_ends_3, 3)
	write_matches('primary_3_matches.txt', compare_primary_3)

	## all ends
	pyRAP_ends_5=pl.read_csv('pyRAP_5_ends.txt', separator='\t')
	pyRAP_ends_3=pl.read_csv('pyRAP_3_ends.txt', separator='\t')

	compare_pyRAP_5, all_5_count = compare_ends(filtered_isoforms, pyRAP_ends_5, 5)
	write_matches('all_5_matches.txt', compare_pyRAP_5)

	compare_pyRAP_3, all_3_count = compare_ends(filtered_isoforms, pyRAP_ends_3, 3)
	write_matches('all_3_matches.txt', compare_pyRAP_3)

	## quantify degree of matching

	with open('match_report.txt', 'w') as o:
		o.write("# of isoforms that match a published 5' end: " + str(all_5_count) + '\n')
		o.write("# of isoforms that match a published primary 5' end: " + str(primary_5_count) + '\n')
		o.write("# of isoforms that match a published 3' end: " + str(all_3_count) + '\n')
		o.write("# of isoforms that match a published primary 3' end: " + str(primary_3_count) + '\n')


	all_5_match, no_match_5 = count_matches(compare_pyRAP_5)
	all_3_match, no_match_5 = count_matches(compare_pyRAP_3)

	with open('match_report.txt', 'a') as o:
		o.write("# of published 5' ends with at least 1 isoform match: " + str(all_5_match) + '\n')
		o.write("# of published 3' ends with at least 1 isoform match: " + str(all_3_match) + '\n')

def count_matches(matches):
	match = 0
	no_match = 0
	for i in matches['f']:
		end_obj = matches['f'][i]
		if len(end_obj.iso_ids) > 0:
			match+=1
		elif len(end_obj.iso_ids) == 0:
			no_match+=1

	for i in matches['r']:
		end_obj = matches['r'][i]
		if len(end_obj.iso_ids) > 0:
			match+=1
		elif len(end_obj.iso_ids) == 0:
			no_match+=1

	return match, no_match


def compare_ends(isoform_df, srs_ends, end_type):

	ends_f=srs_ends.filter(pl.col('Strand') == 'f')
	ends_r=srs_ends.filter(pl.col('Strand') == 'r')

	verified_ends = {'f':{}, 'r':{}}

	for end in ends_f.iter_rows(named=True):
		e = end['Position']
		verified_ends['f'][e] = End(e, 'f')

	for end in ends_r.iter_rows(named=True):
		e = end['Position']
		verified_ends['r'][e] = End(e, 'r')

	match_count = 0
	for iso in isoform_df.iter_rows(named=True):
		start_list = list(map(int,iso['associated_starts'].split(',')))
		end_list = list(map(int,iso['associated_ends'].split(',')))
		strand = iso['strand']

		if end_type == 5:
			for end in verified_ends[strand]:
				e = verified_ends[strand][end].srs_end
				if strand == 'f':
					distances = [(start - e) for start in start_list]
				elif strand == 'r':
					distances = [(e - end) for end in end_list]
				dist_bool = [(0<i<20) for i in distances]
				if any(dist_bool):
					verified_ends[strand][end].add_ids(iso['isoform_id'])
					match_count+=1

		elif end_type == 3:
			for end in verified_ends[strand]:
				e = verified_ends[strand][end].srs_end
				if strand == 'f':
					distances = [abs(end - e) for end in end_list]
				elif strand == 'r':
					distances = [abs(start - e) for start in start_list]
				dist_bool = [(i<10) for i in distances]
				if any(dist_bool):
						verified_ends[strand][end].add_ids(iso['isoform_id'])
						match_count+=1

	return verified_ends, match_count

def write_matches(outfile, compared_ends):
	with open(outfile, 'w') as o:
		o.write('published_end\tstrand\tmatching_isoform_ids\n')
		for end in compared_ends['f']:
			end_obj = compared_ends['f'][end]
			o.write(str(end)+'\t')
			o.write('f\t')
			o.write(','.join(list(map(str, end_obj.iso_ids))))
			o.write('\t')
			o.write('\n')
		for end in compared_ends['r']:
			end_obj = compared_ends['r'][end]
			o.write(str(end)+'\t')
			o.write('r\t')
			o.write(','.join(list(map(str, end_obj.iso_ids))))
			o.write('\t')
			o.write('\n')

class End:
	def __init__(self, srs_end, strand):
		self.srs_end=srs_end
		self.strand=strand
		self.iso_ids=[]

	def add_ids(self, iso_id):
		self.iso_ids.append(iso_id)


if __name__ == '__main__':
	main()