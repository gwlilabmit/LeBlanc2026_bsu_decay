import rifseq.count_ends as count

def main():

	## 190710 wild type
	file_name_CDS='fits/annotations/ncRNA_internal_168.mochiview.txt'

	c_dic, c_index=count.parse_CDS_file_by_type(file_name_CDS, 'ncRNA')

	wig_front='190710/aligned_wigs/190710Li_D19-'
	wig_rear='_NA_sequence-stripped_v2_bsubtilis_v1k1.align_5_'
	wig_indices=['7618', '7619', '7620', '7621', '7879', '7880', '7624', '7625']

	for w in wig_indices:
		c_dic, c_index = count.get_total_CDS(c_dic, c_index, w, wig_front, wig_rear)

	count.write_dic(c_dic, c_index, '190710/norm_factors/ncRNA_counts.txt')

	count_dic, count_header=count.read_count_file('190710/norm_factors/ncRNA_counts.txt', bin_cut_off=10, first_cut_off=100)
	count.write_dic(count_dic, count_header, '190710/norm_factors/ncRNA_averages.txt')

	## 211103 wild type
	file_name_CDS='fits/annotations/ncRNA_internal_168.mochiview.txt'

	c_dic, c_index=count.parse_CDS_file_by_type(file_name_CDS, 'ncRNA')

	wig_front = '211103/wigs/211103Li_D21-'
	wig_rear = '_NA_sequence-stripped_v2_bsubtilis_v1k1_5_'
	wig_indices = ['9974', '10011', '10012', '10013', '9978', '10014', '10015', '9981']

	for w in wig_indices:
		c_dic, c_index = count.get_total_CDS(c_dic, c_index, w, wig_front, wig_rear)

	count.write_dic(c_dic, c_index, '211103/norm_factors/ncRNA_counts.txt')

	count_dic, count_header=count.read_count_file('211103/norm_factors/ncRNA_counts.txt', bin_cut_off=10, first_cut_off=100)
	count.write_dic(count_dic, count_header, '211103/norm_factors/ncRNA_averages.txt')

	## RNase Y depletion
	file_name_CDS='fits/annotations/ncRNA_internal_168.mochiview.txt'

	c_dic, c_index=count.parse_CDS_file_by_type(file_name_CDS, 'ncRNA')

	wig_front = 'Y-depletion/wigs/y_dep_'
	wig_rear = '_NA_sequence-stripped_v2_bsubtilis_v1k1_5_'
	wig_indices = ['1954', '1793', '5041', '1795', '1796', '1797', '1798', '5042']

	for w in wig_indices:
		c_dic, c_index = count.get_total_CDS(c_dic, c_index, w, wig_front, wig_rear)

	count.write_dic(c_dic, c_index, 'Y-depletion/norm_factors/ncRNA_counts.txt')

	count_dic, count_header=count.read_count_file('Y-depletion/norm_factors/ncRNA_counts.txt', bin_cut_off=10, first_cut_off=100)
	count.write_dic(count_dic, count_header, 'Y-depletion/norm_factors/ncRNA_averages.txt')
	
	## RNase Y expression
	file_name_CDS='fits/annotations/ncRNA_internal_168.mochiview.txt'

	c_dic, c_index=count.parse_CDS_file_by_type(file_name_CDS, 'ncRNA')

	wig_front = 'Y-expression/wigs/y_exp_'
	wig_rear = '_NA_sequence-stripped_v2_bsubtilis_v1k1_5_'
	wig_indices = ['1154', '1785', '2121', '5038', '5039', '1789', '5040']

	for w in wig_indices:
		c_dic, c_index = count.get_total_CDS(c_dic, c_index, w, wig_front, wig_rear)

	count.write_dic(c_dic, c_index, 'Y-expression/norm_factors/ncRNA_counts.txt')

	count_dic, count_header=count.read_count_file('Y-expression/norm_factors/ncRNA_counts.txt', bin_cut_off=10, first_cut_off=100)
	count.write_dic(count_dic, count_header, 'Y-expression/norm_factors/ncRNA_averages.txt')
	
	## RNase J depletion
	file_name_CDS='fits/annotations/ncRNA_internal_168.mochiview.txt'

	c_dic, c_index=count.parse_CDS_file_by_type(file_name_CDS, 'ncRNA')

	wig_front = 'J-depletion/wigs/rnjA_dep_'
	wig_rear = '_NA_sequence-stripped_v2_bsubtilis_v1k1_5_'
	wig_indices = ['7347', '7348', '7349', '7350', 't4', '222012', '222013', '222014']

	for w in wig_indices:
		c_dic, c_index = count.get_total_CDS(c_dic, c_index, w, wig_front, wig_rear)

	count.write_dic(c_dic, c_index, 'J-depletion/norm_factors/ncRNA_counts.txt')

	count_dic, count_header=count.read_count_file('J-depletion/norm_factors/ncRNA_counts.txt', bin_cut_off=10, first_cut_off=100)
	count.write_dic(count_dic, count_header, 'J-depletion/norm_factors/ncRNA_averages.txt')

	## RNase J expression
	file_name_CDS='fits/annotations/ncRNA_internal_168.mochiview.txt'

	c_dic, c_index=count.parse_CDS_file_by_type(file_name_CDS, 'ncRNA')

	wig_front = 'J-expression/wigs/rnjA_exp_'
	wig_rear = '_NA_sequence-stripped_v2_bsubtilis_v1k1_5_'
	wig_indices = ['222015', '222016', '222017', '222018', '222019', '222020', '222021', '222022']

	for w in wig_indices:
		c_dic, c_index = count.get_total_CDS(c_dic, c_index, w, wig_front, wig_rear)

	count.write_dic(c_dic, c_index, 'J-expression/norm_factors/ncRNA_counts.txt')

	count_dic, count_header=count.read_count_file('J-expression/norm_factors/ncRNA_counts.txt', bin_cut_off=10, first_cut_off=100)
	count.write_dic(count_dic, count_header, 'J-expression/norm_factors/ncRNA_averages.txt')

if __name__ == '__main__':
	main()