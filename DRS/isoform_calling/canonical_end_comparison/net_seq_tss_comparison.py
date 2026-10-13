import polars as pl

def main():
	
	y_matches = pl.read_csv('yakhnin_rnet_matches.txt', separator='\t', truncate_ragged_lines=True) # "pause sites"
	j_matches = pl.read_csv('jayasinghe_rnet_matches.txt', separator='\t', truncate_ragged_lines=True) # all peaks, no pause site identification
	
	tss_matches = pl.read_csv('primary_5_matches.txt', separator='\t', truncate_ragged_lines=True)

	isoforms = pl.read_csv('../combined_isoform_quant/all_genes_filtered_by_count.txt', separator='\t')

	y_list = flat_match_list(y_matches)
	j_list = flat_match_list(j_matches)
	tss_list = flat_match_list(tss_matches)

	y_match_count = 0
	y_tss_matches = []
	for i in y_list:
		if i in tss_list:
			y_match_count +=1
			y_tss_matches.append(i)

	j_match_count = 0
	j_tss_matches = []
	for i in j_list:
		if i in tss_list:
			j_match_count +=1
			j_tss_matches.append(i)

	y_all_net = isoforms.filter(pl.col('isoform_id').is_in(y_list))
	yakhnin_tss_isos = y_all_net.filter(pl.col('isoform_id').is_in(y_tss_matches))
	yakhnin_dd_isos = y_all_net.filter(pl.col('isoform_id').is_in(y_tss_matches).not_())

	yakhnin_tss_isos.write_csv('yakhnin_net_tss_matches.txt', separator='\t', include_header=True)
	yakhnin_dd_isos.write_csv('yakhnin_net_non_tss.txt', separator='\t', include_header=True)

	j_all_net = isoforms.filter(pl.col('isoform_id').is_in(j_list))
	jayasinghe_tss_isos = j_all_net.filter(pl.col('isoform_id').is_in(j_tss_matches))
	jayasinghe_dd_isos = j_all_net.filter(pl.col('isoform_id').is_in(j_tss_matches).not_())

	jayasinghe_tss_isos.write_csv('jayasinghe_net_tss_matches.txt', separator='\t', include_header=True)
	jayasinghe_dd_isos.write_csv('jayasinghe_net_non_tss.txt', separator='\t', include_header=True)

def flat_match_list(match_df):

	# get messy list of matches
	temp = pl.Series(match_df.select('matching_isoform_ids')).to_list()

	# remove Nones
	temp = [x for x in temp if x is not None]

	# split isoform ids into individual items
	temp = [x.split(',') for x in temp]

	# flatten list & convert to integers
	temp = [int(x) for xs in temp for x in xs]

	# remove duplicates
	match_list = list(set(temp))

	return match_list

if __name__ == '__main__':
	main()