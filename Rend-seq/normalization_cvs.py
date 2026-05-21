import rifseq.coev as coev

def main():

	ncRNA_avg, ncRNA_index=coev.read_file('190710/norm_factors/ncRNA_averages.txt')
	spikein_avg, spikein_index=coev.read_file('190710/norm_factors/spikein_averages.txt')

	coev.write_plot_cv(spikein_avg, spikein_avg, '190710/ncRNA_normalization/spike_vs_spike_cv.txt')
	coev.write_plot_cv(ncRNA_avg, spikein_avg, '190710/ncRNA_normalization/ncRNA_vs_spike_cv.txt')

if __name__ == '__main__':
	main()