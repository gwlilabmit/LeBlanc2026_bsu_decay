import rifseq.count_ends as count

def main():
	file_name_CDS='fits/annotations/CDS_ivt_spike_v3.txt'

	c_dic, c_index=count.parse_CDS_file_by_type(file_name_CDS, 'spikeIN')

	wig_front='190710/unaligned_v2k1_wigs/190710Li_D19-'
	wig_rear='_NA_sequence-stripped_v2_bsubtilis_v1k1_ivt_spikein_v2k1.align_5_'
	wig_indices=['7618', '7619', '7620', '7621', '7879', '7880', '7624', '7625']

	for w in wig_indices:
		c_dic, c_index = count.get_total_CDS(c_dic, c_index, w, wig_front, wig_rear)

	count.write_dic(c_dic, c_index, '190710/norm_factors/spikein_counts.txt')

	count_dic, count_header=count.read_count_file('190710/norm_factors/spikein_counts.txt', bin_cut_off=10, first_cut_off=100)
	count.write_dic(count_dic, count_header, '190710/norm_factors/spikein_averages.txt')

if __name__ == '__main__':
	main()