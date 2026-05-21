## this file contains functions for normalizing to normalization factors and fitting decay model to data

from scipy import optimize
import numpy as np
import math
import copy

# for normalizing average read per base quantifications by ncRNA abundance (all analyses except initial ncRNA vs. spike-in comparisons)
def read_avg_file(avg_file_name, norm_file_name):
	norm_file=open(norm_file_name, 'r')
	norm_header=norm_file.readline()
	norm_dic={}

	for line in norm_file:
		info= line.strip().split()
		length=len(info)
		gene= info[0]
		t0 = float(info[1])
		t05 = float(info[2])
		t125 = float(info[3])
		t2 = float(info[4])
		t4 = float(info[5])
		t8 = float(info[6])
		t12 = float(info[7])
		if (length==9):
			t16 = float(info[8])
			t_list = [t0,t05,t125,t2,t4,t8,t12,t16]
			norm_dic[gene] = t_list
		elif (length==8):
			t_list = [t0,t05,t125,t2,t4,t8,t12]
			norm_dic[gene] = t_list

	file_in=open(avg_file_name, 'r')
	header=file_in.readline()
	data_dic={}
	for line in file_in:
		info=line.strip().split()
		gene=info[0]
		length=len(info)
		t0=float(info[1])
		t05=float(info[2])
		t125=float(info[3])
		t2=float(info[4])
		t4=float(info[5])
		t8=float(info[6])
		t12=float(info[7])
		if (length==9):
			t16 = float(info[8])
			t_list = [t0,t05,t125,t2,t4,t8,t12,t16]
			data_dic[gene]=t_list
		elif (length==8):
			t_list=[t0,t05,t125,t2,t4,t8,t12]
			data_dic[gene]=t_list

	temp={}
	for gene in data_dic.keys():
		temp[gene]=[]
		avg_list=data_dic.get(gene)
		for i in range(len(avg_list)):
			gene_i=avg_list[i]
			ssrA_i=norm_dic['ssrA'][i]
			rnaP_i=norm_dic['rnaP'][i]
			(temp.get(gene)).append(gene_i/(ssrA_i+rnaP_i))

	norm_file.close()
	file_in.close()
	return temp

# removes time points that have 0 counts (either inherently or ones that have been zeroed out due to low coverage)
def remove_zeros(y_values, x_values):
	new_y=[]
	new_x=[]
	check=1
	for i in range(len(y_values)):
		y=y_values[i]
		if(y>0 and check>0):
			new_y.append(y)
			new_x.append(x_values[i])
		else:
			check=-1
	return new_y, new_x

# fits full decay model from given parameters
# to use with fit parameters and residuals
def fit_model(params, x_value, y_observed):
	k=params[0]
	t0=params[1]
	# need to raise 2 to the power of this value because y_observed has previously been log2 transformed
	# having this here rather than being an output of optimize.leastsquares restricts alpha to be equal to the RNA amount at t=0
	a=2**(y_observed[0])
	B=params[2]

	# separate time points into pre- and post-time delay
	x_val=list(x_value)
	x_delay=[]
	x_fit=[]
	for x in x_val:
		if(x<t0):
			x_delay.append(x)
		else:
			x_fit.append(x)
	# calculate predicted/calculated value for each time point as alpha (if pre-time delay)
	# or using the full model (if post-time delay)
	returned_vals=[]
	for x in x_delay:
		returned_vals.append(math.log2(a))
	for x in x_fit:
		returned_vals.append(math.log2(((a-B)*np.power(2,(-(x-t0)/k)))+B))
	return np.array(returned_vals)

## calculates residuals for full decay model
def residuals(parameters, x_values, y_observed):
	y_calculated=fit_model(parameters, x_values, y_observed)
	return abs(y_calculated-y_observed)

# optimizes fit parameters using least squares
# y(t=0) as the steady state amount
def fit_parameters(new_y, new_x):
	# if there are enough non-zero time points to fit (at least 3)
	# take log2 of normalized RNA amount at each time point
	if (len(new_y)>=3):
		log_y=[]
		for y in new_y:
			log_y.append(math.log2(y))
		log_y=np.array(log_y)

		# bounds within which to fit the model
		k_low=0.1
		k_high=100

		t_low=0
		t_high=9.0

		B_low=0.0
		B_high=new_y[0]

		params_init=[1.0, 0, 0.0]
		param_bounds=([k_low, t_low, B_low], [k_high, t_high, B_high])

		OptimizeResult=optimize.least_squares(residuals, params_init, args=(new_x, log_y), bounds=param_bounds)

		parameters=OptimizeResult.x
		k=parameters[0]
		t0=parameters[1]
		a=new_y[0]
		B=parameters[2]
		pl=B/a

		return [k, t0, a, B, pl]
		
	# if there aren't enough non-zero time points to fit, return 0
	else:
		return [0]


# currently not reporting out the individual values that coefficients of variation are coming from to keep the output more manageable
def crossvalidate(fit_params, y_val, x_val, time):
	points=len(y_val)
	if fit_params==[0]:
		fit_params=[-1, -1, -1, -1, -1]

	internal_timepoints=range(0,points)
	ks=[]
	delays=[]
	alphas=[]
	Bs=[]
	plateaus=[]

	for i in internal_timepoints:
		x_copy=copy.deepcopy(x_val)
		y_copy=copy.deepcopy(y_val)
		time_pop=x_copy.pop(i)
		val_pop=y_copy.pop(i)

		new_params=fit_parameters(y_copy, x_copy)
		if new_params==[0]:
			new_params=[-1, -1, -1, -1, -1]

		ks.append(new_params[0])
		delays.append(new_params[1])
		alphas.append(new_params[2])
		Bs.append(new_params[3])
		plateaus.append(new_params[4])

	CoeV_k=np.std(ks, ddof=1)/np.mean(ks)
	CoeV_delay=np.std(delays, ddof=1)/np.mean(delays)
	CoeV_alpha=np.std(alphas, ddof=1)/np.mean(alphas)
	CoeV_B=np.std(Bs, ddof=1)/np.mean(Bs)
	CoeV_plateau=np.std(plateaus, ddof=1)/np.mean(plateaus)

	CoeV_list=[CoeV_k, CoeV_delay, CoeV_alpha, CoeV_B, CoeV_plateau]

	return CoeV_list

# for calculating differential stability of overlapping isoforms
# subracts average reads per base of longer isoforms from shorter ones they overlap with
# to subtract multiple overlapping isoforms, use a list for the longer_genes argument
# appends result of subtraction directly to the average reads per base file
# should be able to handle all data sets
def subtract_isoforms(avg_file_name, shorter_gene, longer_genes):
	file_in = open(avg_file_name, 'r')
	header = file_in.readline().strip().split()
	avg_dic = {}

	for line in file_in:
		info = line.strip().split()
		length=len(info)
		gene=info[0]
		t0=float(info[1])
		t05=float(info[2])
		t1=float(info[3])
		t2=float(info[4])
		t4=float(info[5])
		t8=float(info[6])
		t12=float(info[7])
		if (length==9):
			t16=float(info[8])
			t_list = [t0,t05,t1,t2,t4,t8,t12,t16]
			avg_dic[gene] = t_list
		elif (length==8):
			t_list = [t0,t05,t1,t2,t4,t8,t12]
			avg_dic[gene] = t_list
	file_in.close()

	short_sub = []
	for i in range(len(avg_dic[shorter_gene])):
		shorter_i = avg_dic[shorter_gene][i]
		longer_i = []
		for gene in longer_genes:
			temp = avg_dic[gene][i]
			longer_i.append(temp)
		longer_i = sum(longer_i)
		sub = shorter_i - longer_i
		short_sub.append(sub)

	with open(avg_file_name, 'a') as f:
		f.write(shorter_gene + '-subtracted')
		for i in short_sub:
			f.write('\t' + str(i))
		f.write('\n')

	f.close()

# for aggregating multiple genes together – for when you have multiple coding sequences in an isoform
# example for gene list: [['cggR'], ['gapA'], ['pgk', 'tpi', 'pgm', 'eno']]
# this doesn't filter for coverage, even though it replaces the averaging step that usually does that
# because i'm assuming you'll always be doing this on a case by case basis
def aggregate_averages(count_file_name, gene_list):
	file_in = open(count_file_name, 'r')
	header = file_in.readline().strip().split()
	out_dic = {}

	in_dic = {}
	for line in file_in:
		info = line.strip().split()
		gene = info[0]
		in_dic[gene] = info[1:]

	# now you have a dictionary of all of the info in the file, with genes as keys

	for l in gene_list:
		if len(l) == 1:
			gene = l[0]
			info = in_dic.get(gene)
			start = float(info[0])
			stop = float(info[1])
			gene_length = stop-start+1.0

			temp = []
			for i in range(3, len(info)):
				entry = float(info[i])
				temp.append(entry/gene_length)
			out_dic[gene] = temp

		if len(l) > 1:
			list_length = 0
			for cds in l:
				gene = cds
				info = in_dic.get(gene)
				start = float(info[0])
				stop = float(info[1])
				gene_length = stop-start+1
				list_length = list_length + gene_length


			read_count = [0, 0, 0, 0, 0, 0, 0, 0]
			gene_names = []
			for cds in l:
				gene = cds
				info = in_dic.get(gene)
				t0 = float(info[3])
				t05 = float(info[4])
				t1 = float(info[5])
				t2 = float(info[6])
				t4 = float(info[7])
				t8 = float(info[8])
				t12 = float(info[9])
				if len(info)==11:
					t16 = float(info[10])
					t_list = [t0,t05,t1,t2,t4,t8,t12,t16]
					read_count = [read_count[i] + t_list[i] for i in range(len(t_list))]
					gene_names.append(gene)
				elif len(info)==10:
					t_list = [t0,t05,t1,t2,t4,t8,t12]
					read_count = [read_count[i] + t_list[i] for i in range(len(t_list))]
					gene_names.append(gene)
				
			gene_names = '-'.join(gene_names)

			temp = [read_count[i]/list_length for i in range(len(read_count))]
			out_dic[gene_names] = temp
	
	file_in.close()
	out_header = header[4:]
	return out_dic, out_header

