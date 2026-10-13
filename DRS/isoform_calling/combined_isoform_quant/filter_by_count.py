import polars as pl

def main():

	counts = pl.read_csv('isoform_counts.txt', separator='\t')
	norm_vals = pl.read_csv('isoform_norm_vals.txt', separator='\t')
	fits = pl.read_csv('isoform_fit_cross_all.txt', separator='\t')

	## filter to only isoforms with detected reads in at least the first 3 time points

	filtered_counts = counts.filter((pl.col('t0') > 2) & (pl.col('t0.5') > 2) & (pl.col('t2') > 2))

	filtered_ids = pl.Series(filtered_counts.select(pl.col('isoform_id'))).to_list()
	
	filtered_norm_vals = norm_vals.filter(pl.col('isoform_id').is_in(filtered_ids))

	all_filtered = filtered_norm_vals.join(fits, on='isoform_id', how='inner')

	## compress tau above 15 to 15 minutes
	all_filtered = all_filtered.with_columns(
		k=pl.when(pl.col('k') > 15).then(pl.lit(15)).otherwise(pl.col('k')))

	all_filtered.write_csv('all_genes_filtered_by_count.txt', separator='\t')

	filtered_fits = all_filtered.filter(pl.col('k_coef_var')<0.3)

	filtered_fits.write_csv('fits_filtered_by_cv.txt', separator='\t')

	### remove points based on maximum filter – same as used for rend-seq

	max_passing_id_list = []

	for i in filtered_fits.iter_rows(named=True):
		amounts = [i['t0'], i['t0.5'], i['t2'], i['t4'], i['t8'], i['t12'], i['t16']]
		alpha = i['alpha']
		iso_id = i['isoform_id']

		amounts_rel = [(x/alpha) for x in amounts]

		maximum = max(amounts_rel)
		temp = amounts_rel.copy()
		temp.remove(maximum)
		second_max = max(temp)

		#if (maximum > 1.2 and second_max > 1.1):
		if (maximum > 2):
			continue
		else:
			max_passing_id_list.append(iso_id)

	max_filtered_fits = filtered_fits.filter(pl.col('isoform_id').is_in(max_passing_id_list))

	max_filtered_fits.write_csv('maximum_filtered_out.txt', separator='\t')

	### categorization

	## categories:
	### fully cover >1 cds and no partial overlaps
	### fully cover >1 cds, partially overlap >1 cds
	### >1 partial overlap and no fully covered CDSes
	### no CDS overlap

	class1=0 ## only fully covers cdses
	class2=0 ## full coverage and partial overlap
	class3=0 ## only partial overlap
	class4=0 ## no cds overlap

	for isoform in all_filtered.iter_rows(named=True):
		cover_list = isoform['covered_genes'].split(',')
		overlap_list = isoform['overlapping_genes'].split(',')

		if overlap_list==['None']:
			class4+=1

		elif cover_list==['None'] and overlap_list!=['None']:
			class3+=1

		elif cover_list!=['None'] and overlap_list!=['None'] and cover_list!=overlap_list:
			class2+=1

		elif cover_list!=['None'] and overlap_list!=['None'] and cover_list==overlap_list:
			class1+=1

		else:
			other+=1

	with open('isoform_categorization_report.txt', 'w') as o:
		o.write('Covers >1 annotated CDS with no partial overlaps: '+str(class1)+'\n')
		o.write('Covers >1 annotated CDS with >1 partial overlap: '+str(class2)+'\n')
		o.write('>1 partial CDS overlap with no full CDSes covered: '+str(class3)+'\n')
		o.write('No overlap with annotated CDSes: '+str(class4)+'\n')


if __name__ == '__main__':
	main()