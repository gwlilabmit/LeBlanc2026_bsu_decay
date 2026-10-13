import polars as pl

def main():

	gff3 = pl.read_csv('Lawaetz_2024_annotation/B_subtilis_pyRAP.gff3', separator='\t', has_header=False)
	gff3 = gff3.drop(['column_6', 'column_8']).rename({'column_1':'chr', 'column_2':'source', 'column_3':'type', 'column_4':'start', 'column_5':'end', 'column_7':'strand', 'column_9':'info'})
	
	## filter to only 5' and 3' ends of transcripts; remove other kinds of annotations

	gff3 = gff3.filter(((pl.col('source') == 'RendSeq')) & 
		(pl.col('type') != 'CDS') & 
		(pl.col('type') != 'intergenic_UTR') & 
		(pl.col('type') != 'tRNA') & 
		(pl.col('type') != 'rRNA') & 
		(pl.col('type') != 'Novel_rRNA_or_tRNA_adjacent') &
		(pl.col('type') != 'rRNA_or_tRNA_adjacent_transcript'))

	## filter to only LB exponential growth (either WT_LB_25s or WT_Bsub_LB_exp)

	gff3 = gff3.filter(pl.col('info').str.contains_any(['WT_LB_25s', 'WT_Bsub_LB_exp']))

	#### BELOW IS CLOSE, BUT NEEDS TO TAKE STRAND INTO ACCOUNT
	## list of only 5' ends deriving from peaks

	ends_5_f = pl.Series(gff3.filter(((pl.col('strand') == '+') &
			(pl.col('type') == "5'UTR") | 
			(pl.col('type') == "sRNA") | 
			(pl.col('type') == "Novel_transcript") | 
			(pl.col('type') == "transcript") | 
			(pl.col('type') == "putative_sRNA")))
	.select(pl.col('start'))).to_list()

	ends_5_r = pl.Series(gff3.filter(((pl.col('strand') == '-') &
			(pl.col('type') == "5'UTR") | 
			(pl.col('type') == "sRNA") | 
			(pl.col('type') == "Novel_transcript") | 
			(pl.col('type') == "transcript") | 
			(pl.col('type') == "putative_sRNA")))
	.select(pl.col('end'))).to_list()

	## list of only 3' ends deriving from peaks

	ends_3_f = pl.Series(gff3.filter(((pl.col('strand') == '+') &
			(pl.col('type') == "3'UTR") | 
			(pl.col('type') == "sRNA") | 
			(pl.col('type') == "Novel_transcript") | 
			(pl.col('type') == "transcript") | 
			(pl.col('type') == "putative_sRNA")))
	.select(pl.col('end'))).to_list()

	ends_3_r = pl.Series(gff3.filter(((pl.col('strand') == '-') &
			(pl.col('type') == "3'UTR") | 
			(pl.col('type') == "sRNA") | 
			(pl.col('type') == "Novel_transcript") | 
			(pl.col('type') == "transcript") | 
			(pl.col('type') == "putative_sRNA")))
	.select(pl.col('start'))).to_list()

	with open('pyRAP_5_ends.txt', 'w') as o:
		o.write('Position\tStrand\n')
		for i in ends_5_f:
			o.write(str(i) + '\tf\n')
		for i in ends_5_r:
			o.write(str(i) + '\tr\n')

	with open('pyRAP_3_ends.txt', 'w') as o:
		o.write('Position\tStrand\n')
		for i in ends_3_f:
			o.write(str(i) + '\tf\n')
		for i in ends_3_r:
			o.write(str(i) + '\tr\n')

	# -------------------------------------

	
	## read in TSS and intrinsic terminator lists from other methods
	## keep end if it's called in at least one of them

	## filter 5' end and 3' end lists to only ends that correspond to documented TSSes or intrinsic terminators
	## save these TSS and intrinsic terminator lists

	tss_forrest = pl.read_csv('Forrest_2022_cappable_seq/cappable_seq_tss_forrest.txt', separator='\t')
	tss_warman = pl.read_csv('Warman_2021_cappable_seq/cappable_seq_tss_warman.txt', separator=' ')

	it_lalanne = pl.read_csv('Lalanne_2018_terminators/3_ends_lalanne.txt', separator='\t').drop('Closest upstream gene').rename({'Position (WT)':'Position'})
	it_mandell = pl.read_csv('Mandell_2022_term_seq/term_seq_ITs_mandell.txt', separator='\t').drop('Relative Position').rename({'POT':'Position'})

	it_lalanne = it_lalanne.rename({'Strand':'Old_strand'}).with_columns(pl.when(pl.col('Old_strand') == 1).then(pl.lit('+')).otherwise(pl.lit('-')).alias('Strand')).drop('Old_strand')


	def id_primary_ends(pyRAP_f, pyRAP_r, data_1, data_2):

		f_1 = pl.Series(data_1.filter(pl.col('Strand') == '+').select(pl.col('Position'))).to_list()
		f_2 = pl.Series(data_2.filter(pl.col('Strand') == '+').select(pl.col('Position'))).to_list()

		primary_ends = []
		for end in pyRAP_f:
			if end in f_1 or end in f_2:
				primary_ends.append((end, 'f'))

			## check for close by ends that are probably the same end
			else:
				end_range = list(range(end-3, end+3+1))
				for i in end_range:
					if i in f_1 or end in f_2:
						temp = [i[0] for i in primary_ends]
						if end not in temp:
							primary_ends.append((end, 'f'))

		r_1 = pl.Series(data_1.filter(pl.col('Strand') == '-').select(pl.col('Position'))).to_list()
		r_2 = pl.Series(data_2.filter(pl.col('Strand') == '-').select(pl.col('Position'))).to_list()

		for end in pyRAP_r:
			if end in r_1 or end in r_2:
				primary_ends.append((end, 'r'))

			## check for close by ends that are probably the same end
			else:
				end_range = list(range(end-3, end+3+1))
				for i in end_range:
					if i in r_2 or end in r_2:
						temp = [i[0] for i in primary_ends]
						if end not in temp:
							primary_ends.append((end, 'r'))

		df = pl.DataFrame({'Position':[i[0] for i in primary_ends], 'Strand':[i[1] for i in primary_ends]})
		return df


	primary_5_df = id_primary_ends(ends_5_f, ends_5_r, tss_forrest, tss_warman)
	primary_5_df.write_csv('primary_5_ends.txt', separator='\t')

	primary_3_df = id_primary_ends(ends_3_f, ends_3_r, it_lalanne, it_mandell)
	primary_3_df.write_csv('primary_3_ends.txt', separator='\t')

if __name__ == '__main__':
	main()