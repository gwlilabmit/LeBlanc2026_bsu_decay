import rifseq.count_ends as count
import rifseq.fit_decay as fit

def main():
	normed=count.normalize_spike_in('190710/norm_factors/ncRNA_averages.txt', '190710/norm_factors/spikein_averages.txt')

	file_out='190710/ncRNA_normalization/3_spikein_normalized_ncRNAs.txt'
	header=['0min', '0.5min', '1min', '2min', '4min', '8min', '16min', '32min']
	count.write_dic(normed, header, file_out)

if __name__ == '__main__':
	main()