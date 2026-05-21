import numpy as np

def main():

	# reads values file and produces a dictionary with the gene name as the key and a list of normalized RNA amounts as the value
	def read_values(infile):
		f=open(infile, 'r')
		header=f.readline()

		gene_dictionary={}
		for line in f:
			info=line.strip().split('\t')
			gene=info[0]
			if len(info)==9:
				timepoints=[info[1], info[2], info[3], info[4], info[5], info[6], info[7], info[8]]
			elif len(info)==8:
				timepoints=[info[1], info[2], info[3], info[4], info[5], info[6], info[7]]
			gene_dictionary[gene]=timepoints

		f.close()
		return gene_dictionary

	# reads fits file and produces a dictionary of gene names and parameters
	def read_fits(infile):
		f=open(infile, 'r')
		header=f.readline()

		gene_dictionary={}
		for line in f:
			info=line.strip().split('\t')
			gene=info[0]
			# order of params list: tau, t0, alpha, B, plateau, tau coefficient of variation
			params=[info[1], info[2], info[3], info[4], info[5], info[6]]
			gene_dictionary[gene]=params

		f.close()
		return gene_dictionary


	# takes dictionaries of fits and values from read_values and read_fits and combines them into one dictionary
	def combine_data(values_dictionary, fits_dictionary):
		combined_dictionary={}
		for gene in values_dictionary.keys():
			gene_dictionary={'parameters':[], 'values':[]}

			values_list=values_dictionary.get(gene)
			gene_dictionary['values']=values_list

			parameters_list=fits_dictionary.get(gene)
			gene_dictionary['parameters']=parameters_list

			combined_dictionary[gene]=gene_dictionary

		return combined_dictionary
			

	# combines all the data into a combined dictionary using the functions above
	# then puts each gene through filtering conditions to decide whether to write its info to a final output file
	def filter_out(values_infile, fits_infile, outfile, count_outfile):
		values=read_values(values_infile)
		parameters=read_fits(fits_infile)
		combined=combine_data(values, parameters)

		cv_over_count=0
		max_over_count=0
		tau_over_count=0
		zero_alpha_count=0

		# header is currently agnostic to what the time points actually are that are being used
		file=open(outfile,'w')
		header=['gene', 'tau', 't0', 'alpha', 'B', 'plateau', 'tau_CV', 't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8']
		for i in header:
			file.write(i + '\t')
		file.write('\n')

		for gene in combined.keys():
			gene_dictionary=combined.get(gene)
			parameters=gene_dictionary.get('parameters')
			values=gene_dictionary.get('values')

			# extract values that are important for filtering and determine whether they pass the filter or not
			cv_check=0
			maximum_check=0
			zero_alpha_check=0
			## first: CV filter
			if (parameters==None):
				continue
			else:
				CV=float(parameters[5])
				if (CV>0.3):
					cv_check+=1
					cv_over_count+=1

			## next: maximum value
			## filters out genes if:
			### their maximum value is greater than 120% of alpha AND
			### their second largest value is greater than 110% of alpha
			# this seems to be the best way to capture genes that genuinely have a significant regime of increasing expression
			# without losing genes that have a one-off high value that might be from noise and doesn't seem to hinder fitting
			# but it's possible we'd want to use a more stringent threshold -- using only the 120% maximum value threshold seems to work pretty well too
			# it just loses more genes that are probably being fit well 
			if (parameters==None):
				continue

			else:
				alpha=float(parameters[2])
				beta=float(parameters[3])
				maximum=float(max(values))

				temp=[float(i) for i in values]
				temp_max=max(temp)
				temp.remove(temp_max)
				second_maximum=max(temp)

				if alpha==0:
					zero_alpha_check+=1
					zero_alpha_count+=1
				else:
					ratio=maximum/alpha
					second_ratio=second_maximum/alpha
				
					if (ratio>1.2 and second_ratio>1.1):
						maximum_check+=1
						max_over_count+=1
					elif (beta>alpha):
						maximum_check+=1
						max_over_count+=1

			# both check values should be 0 if the gene passes both of these thresholds
			# if check is greater than 0, end the loop and move on to the next gene
			check=cv_check+maximum_check
			if (check>0):
				continue

			## finally: compress tau values greater than 15 to 15
			tau=float(parameters[0])
			# if tau is greater than 15, change its value to 15 in the list of parameters
			if (tau>15):
				new_tau=15
				new_parameters=[new_tau, parameters[1], parameters[2], parameters[3], parameters[4], parameters[5]]
				combined[gene].update({'parameters':new_parameters})
				#combined[gene]['parameters'] = new_parameters
				tau_over_count+=1

			file.write(gene + '\t')

			new_parameters=combined[gene].get('parameters')
			for parameter in new_parameters:
				file.write(str(parameter) + '\t')

			for value in values:
				file.write(str(value) + '\t')

			file.write('\n')

		file.close()

		count=open(count_outfile, 'w')
		count.write('Filtered out by coefficient of variation threshold: ' + str(cv_over_count) + ' genes\n')
		count.write('Filtered out by maximum value threshold: ' + str(max_over_count) + ' genes\n')
		count.write('Tau compressed to 15: ' + str(tau_over_count) + ' genes\n')
		count.write('Compressed to 0 reads in fitting: ' + str(zero_alpha_count) + ' genes')
		count.close()

	# note: in the filtered_out_reports, there can be overlap between CV and maximum filtered out genes

	filter_out('21_WT/21_WT_norm_vals.txt', '21_WT/21_WT_fit_cross.txt', '21_WT/21_WT_filtered.txt', '21_WT/21_WT_filtered_out_report.txt')
	filter_out('J_dep/J_dep_norm_vals.txt', 'J_dep/J_dep_fit_cross.txt', 'J_dep/J_dep_filtered.txt', 'J_dep/J_dep_filtered_out_report.txt')
	filter_out('J_exp/J_exp_norm_vals.txt', 'J_exp/J_exp_fit_cross.txt', 'J_exp/J_exp_filtered.txt', 'J_exp/J_exp_filtered_out_report.txt')


if __name__ == '__main__':
	main()
