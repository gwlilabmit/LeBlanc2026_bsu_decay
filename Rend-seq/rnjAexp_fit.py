import rifseq.count_ends as count
import rifseq.fit_decay as fit

def main():
	file_name_CDS='fits/annotations/30_cut_CDS_168.txt'

	c_dic, c_index=count.parse_CDS_file_by_type(file_name_CDS, 'gene')

	wig_front = 'J-expression/wigs/rnjA_exp_'
	wig_rear = '_NA_sequence-stripped_v2_bsubtilis_v1k1_5_'
	wig_indices = ['222015', '222016', '222017', '222018', '222019', '222020', '222021', '222022']

	for w in wig_indices:
		c_dic, c_index = count.get_total_CDS(c_dic, c_index, w, wig_front, wig_rear)

	count.write_dic(c_dic, c_index, 'fits/J_exp/J_exp_counts.txt')

	count_dic, count_header=count.read_count_file('fits/J_exp/J_exp_counts.txt', bin_cut_off=10, first_cut_off=100)
	count.write_dic(count_dic, count_header, 'fits/J_exp/J_exp_averages.txt')

	normed=fit.read_avg_file('fits/J_exp/J_exp_averages.txt', 'J-expression/norm_factors/ncRNA_averages.txt')
	gene_list=normed.keys()
	file_out='fits/J_exp/J_exp_norm_vals.txt'
	header=['0min', '0.5min', '1.25min', '2min', '4min', '6min', '8min', '10min']
	count.write_dic(normed, header, file_out)

	time=[0,.5,1.25,2,4,6,8,10]
	parameters_dic={}
	for i in gene_list:
		new_y, new_x=fit.remove_zeros(normed[i], time)
		parameters=fit.fit_parameters(new_y, new_x)
		errors=fit.crossvalidate(parameters, new_y, new_x, time)
		parameters_dic[i]=parameters+errors

	file_out='fits/J_exp/J_exp_fit_cross.txt'
	header=['k', 't0', 'alpha', 'B', 'plateau', 'k_coef_var', 'delay_coef_var', 'alpha_coef_var', 'B_coef_var', 'plateau_coef_var']

	count.write_dic(parameters_dic, header, file_out)

if __name__ == '__main__':
	main()