import polars as pl
import os

def main():
	primary_isoforms = pl.read_csv('../canonical_end_comparison/filtered_primary_isoforms.txt', separator='\t')
	ends_5 = pl.read_csv('../canonical_end_comparison/primary_5_matches.txt', separator='\t', truncate_ragged_lines=True)
	ends_3 = pl.read_csv('../canonical_end_comparison/primary_3_matches.txt', separator='\t', truncate_ragged_lines=True)

	newdir_list = ['end_matches_5/', 'end_matches_5/f/', 'end_matches_5/r/', 'end_matches_3/', 'end_matches_3/f/', 'end_matches_3/r/']
	for i in newdir_list:
		if not os.path.exists(i):
			os.makedirs(i)

	for i in ends_5.iter_rows(named=True):
		end = i['published_end']
		strand = i['strand']
		if i['matching_isoform_ids'] != None:
			isoform_ids = list(map(int, i['matching_isoform_ids'].split(',')))
		else:
			continue

		primary_temp = primary_isoforms.filter(pl.col('isoform_id').is_in(isoform_ids)).filter(pl.col('k_coef_var')<0.3)
		if primary_temp.select(pl.len()).item() > 1:
			primary_temp.write_csv('end_matches_5/'+strand+'/'+str(end)+'.txt', separator='\t', include_header=True)

	for i in ends_3.iter_rows(named=True):
		end = i['published_end']
		strand = i['strand']
		if i['matching_isoform_ids'] != None:
			isoform_ids = list(map(int, i['matching_isoform_ids'].split(',')))
		else:
			continue

		primary_temp = primary_isoforms.filter(pl.col('isoform_id').is_in(isoform_ids)).filter(pl.col('k_coef_var')<0.3)
		if primary_temp.select(pl.len()).item() > 1:
			primary_temp.write_csv('end_matches_3/'+strand+'/'+str(end)+'.txt', separator='\t', include_header=True)

class End:
	def __init__(self, identifier, strand):
		self.identifier = identifier
		self.strand = strand
		#self.end_type = end_type
		self.associated_ends = []
		self.associated_ids = []

	def add_ends(self, ends):
		for i in ends:
			if i not in self.associated_ends:
				self.associated_ends.append(i)

	def add_ids(self, iso_id):
		self.associated_ids.append(iso_id)

if __name__=='__main__':
	main()