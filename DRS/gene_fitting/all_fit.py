from scipy import optimize
import numpy as np
import math
import copy

def main():

	normed=read_file('../normalization/sbo_norm_vals.txt')
	gene_list=normed.keys()

	time=[0,0.5,2,4,8,12,16]
	parameters_dic={}
	for i in gene_list:
		new_y, new_x=remove_zeros(normed[i], time)
		parameters=fit_parameters(new_y, new_x)
		errors=crossvalidate(parameters, new_y, new_x, time)
		parameters_dic[i]=parameters+errors

	file_out='ont_fit_cross_all.txt'
	header=['k', 't0', 'alpha', 'B', 'plateau', 'k_coef_var', 'delay_coef_var', 'alpha_coef_var', 'B_coef_var', 'plateau_coef_var']
	write_dic(parameters_dic, header, file_out)

### ----------------------------

def read_file(normalized_file_name):
	norm_file=open(normalized_file_name, 'r')
	norm_header=norm_file.readline()
	norm_dic={}

	for line in norm_file:
		info= line.strip().split()
		length=len(info)
		gene= info[0]
		t0 = float(info[1])
		t05 = float(info[2])
		#t125 = float(info[3])
		t2 = float(info[3])
		t4 = float(info[4])
		t8 = float(info[5])
		t12 = float(info[6])
		t16 = float(info[7])
		t_list = [t0,t05,t2,t4,t8,t12,t16]
		norm_dic[gene] = t_list

	norm_file.close()
	return norm_dic

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
# adjusted from fitting alpha as its own free parameter to just using y(t=0) as the steady state amount on 6/21/24
# adjustment applies both here and in fit_model
# compare to mzt_halflives/halflives/fit_decay_CV
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

# takes a dictionary of data and an array describing what's in the dictionary and writes it to a txt file
# this should be broadly applicable regardless of the type of data (counts, averages, fits) as long as you have a dictionary of data
# and a header in array format
def write_dic(dic, index_array, out_file_name):
	out_file=open(out_file_name, 'w')
	out_file.write('gene')
	for ind in index_array:
		out_file.write('\t'+ind)
	out_file.write('\n')

	dic_keys=dic.keys()
	for c in dic_keys:
		out_file.write(c)
		info=dic[c]
		for i in info:
			out_file.write('\t'+str(i))
		out_file.write('\n')
	out_file.close

if __name__ == '__main__':
	main()