import polars as pl

def main():

	isoforms = pl.read_csv('../combined_isoform_quant/all_genes_filtered_by_count.txt', separator='\t')
	primary_5_matches = pl.read_csv('primary_5_matches.txt', separator='\t', truncate_ragged_lines=True)
	primary_3_matches = pl.read_csv('primary_3_matches.txt', separator='\t', truncate_ragged_lines=True)

	primary_isoform_ids = []

	for i in primary_5_matches.iter_rows(named=True):
		end_i = i['published_end']
		strand_i = i['strand']
		if i['matching_isoform_ids'] != None:
			ids_i = list(map(int, i['matching_isoform_ids'].split(',')))

		## filtering to reduce amount of looping through that you need to do
		## don't expect any isoform to be longer than 20kb
		temp_3 = primary_3_matches.filter(pl.col('strand') == strand_i).filter(abs(pl.col('published_end') - end_i) < 20000)

		for k in temp_3.iter_rows(named=True):
			end_k = k['published_end']
			if k['matching_isoform_ids'] != None:
				ids_k = list(map(int, k['matching_isoform_ids'].split(',')))

				for isoform in ids_i:
					if isoform in ids_k:
						primary_isoform_ids.append(isoform)

	primary_isoforms = isoforms.filter(pl.col('isoform_id').is_in(primary_isoform_ids))
	primary_ids = pl.Series(primary_isoforms.select(pl.col('isoform_id'))).to_list()
	primary_isoforms.write_csv('filtered_primary_isoforms.txt', separator='\t', include_header=True)

	## primary designation takes precedence over nascent designation
	nascent_isoforms = pl.read_csv('jayasinghe_net_tss_matches.txt', separator='\t')
	nascent_isoforms = nascent_isoforms.filter(pl.col('isoform_id').is_in(primary_ids).not_())
	nascent_ids = pl.Series(nascent_isoforms.select(pl.col('isoform_id'))).to_list()
	nascent_isoforms.write_csv('filtered_nascent_isoforms.txt', separator='\t', include_header=True)

	## nascent designation takes precedence over secondary designation when 5' end matches a TSS (already baked into nascent list)
	secondary_isoforms = isoforms.filter(pl.col('isoform_id').is_in(primary_isoform_ids).not_())
	secondary_isoforms = secondary_isoforms.filter(pl.col('isoform_id').is_in(nascent_ids).not_())
	secondary_isoforms.write_csv('filtered_secondary_isoforms.txt', separator='\t', include_header=True)

	print(isoforms.shape)
	print(primary_isoforms.shape)
	print(nascent_isoforms.shape)
	print(secondary_isoforms.shape)

if __name__=='__main__':
	main()