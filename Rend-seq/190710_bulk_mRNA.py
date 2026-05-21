import rifseq.count_ends as count

def main():
	file_name_CDS='fits/annotations/30_cut_CDS_168.txt'
	c_dic, c_index=count.parse_CDS_file_by_type(file_name_CDS, 'gene')

	wig_front = '190710/aligned_wigs/190710Li_D19-'
	wig_rear = '_NA_sequence-stripped_v2_bsubtilis_v1k1.align_5_'
	wig_indices = ['7618', '7619', '7620', '7621', '7879', '7880', '7624', '7625']

	for w in wig_indices:
		c_dic, c_index = count.get_total_CDS(c_dic, c_index, w, wig_front, wig_rear)

	count.write_dic(c_dic, c_index, '190710/bulk/190710_WT_counts.txt')

	header, total=count.get_bulk_counts('190710/bulk/190710_WT_counts.txt')
		
	count.write_bulk_counts('190710_WT', header, total, '190710/bulk/bulk_counts.txt')

	spikein=count.normalize_spike_in('190710/bulk/bulk_counts.txt', '190710/norm_factors/spikein_averages.txt')
	ncRNA=count.normalize_ncRNA('190710/bulk/bulk_counts.txt', '190710/norm_factors/ncRNA_averages.txt')

	count.write_dic(spikein, header, '190710/bulk/bulk_spikein_normalized.txt')
	count.write_dic(ncRNA, header, '190710/bulk/bulk_ncRNA_normalized.txt')

if __name__ == '__main__':
	main()