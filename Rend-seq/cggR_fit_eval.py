import rifseq.count_ends as count
import rifseq.fit_decay as fit
import rifseq.fit_decay_evaluation as evaluate

def main():
	file_name_CDS='fits/30_cut_CDS_168.txt'

	c_dic, c_index=count.parse_CDS_file_by_name(file_name_CDS, 'cggR')

	wig_front = '211103/wigs/211103Li_D21-'
	wig_rear = '_NA_sequence-stripped_v2_bsubtilis_v1k1_5_'
	wig_indices = ['9974', '10011', '10012', '10013', '9978', '10014', '10015', '9981']

	for w in wig_indices:
		c_dic, c_index = count.get_total_CDS(c_dic, c_index, w, wig_front, wig_rear)

	count.write_dic(c_dic, c_index, 'fits/cggR_fit_evaluation/cggR_counts.txt')

	count_dic, count_header=count.read_count_file('fits/cggR_fit_evaluation/cggR_counts.txt', bin_cut_off=10, first_cut_off=100)
	count.write_dic(count_dic, count_header, 'fits/cggR_fit_evaluation/cggR_averages.txt')

	normed=fit.read_avg_file('fits/cggR_fit_evaluation/cggR_averages.txt', '211103/norm_factors/ncRNA_averages.txt')
	gene_list=normed.keys()
	file_out='fits/cggR_fit_evaluation/cggR_norm_vals.txt'
	header=['0min', '0.5min', '1.25min', '2min', '4min', '8min', '12min', '16min']
	count.write_dic(normed, header, file_out)

	file_in=open('fits/cggR_fit_evaluation/cggR_norm_vals.txt', 'r')
	header=file_in.readline().strip().split()
	normed={}

	for line in file_in:
		info = line.strip().split()
		gene = info[0]
		t0 = float(info[1])
		t05 = float(info[2])
		t1 = float(info[3])
		t2 = float(info[4])
		t4 = float(info[5])
		t8 = float(info[6])
		t12 = float(info[7])
		t16 = float(info[8])

		normed[gene] = [t0, t05, t1, t2, t4, t8, t12, t16]

	time=[0,.5,1,2,4,8,12,16]
	gene_list = ['cggR']

	# simple exponential + time delay (no cross-validation)
	delay_param_dic={}
	for i in gene_list:
		new_y, new_x=evaluate.remove_zeros(normed[i], time)
		parameters=evaluate.fit_delay_parameters(new_y, new_x)
		errors=evaluate.crossvalidate_delay(parameters, new_y, new_x, time)
		delay_param_dic[i]=parameters+errors

	delay_file_out='fits/cggR_fit_evaluation/cggR_simple_delay_model.txt'
	delay_header='gene\tk\tt0\ta\tk_cv\tt0_cv\talpha_cv\n'
	evaluate.write_fit(delay_param_dic, delay_header, delay_file_out)

	# full model fitting (no cross-validation)
	full_param_dic={}
	for i in gene_list:
		new_y, new_x=evaluate.remove_zeros(normed[i], time)
		parameters=evaluate.fit_full_parameters(new_y, new_x)
		errors=evaluate.crossvalidate_full(parameters, new_y, new_x, time)
		full_param_dic[i]=parameters+errors

	full_file_out='fits/cggR_fit_evaluation/cggR_full_model.txt'
	full_header='gene\tk\tt0\talpha\tB\tplateau\tk_cv\tt0_cv\talpha_cv\tB_cv\tplateau_cv\n'
	evaluate.write_fit(full_param_dic, full_header, full_file_out)

if __name__ == '__main__':
	main()