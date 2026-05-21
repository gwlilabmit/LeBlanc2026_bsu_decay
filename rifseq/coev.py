import scipy as sp
import numpy as np
import math
import statistics as stats

# function for calculating an average
def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

# read file of counts or average reads per base
# produces a dictionary from the file as well as an index
# indices refers to the columns of the 5avg input file. the first 4 columns are gene info (name, start, stop and strand) and are thus excluded from this analysis
# columns 5-20 are the columns of average 5' end read counts from each library (timepoint/treatment)
# specifying indices here allows you to exclude some columns from consideration if necessary (e.g. indices=range(4,13) + range(15,20) excludes columns 13 and 14)
def read_file(avg_file_name, indices= range(1,9)):
	file_in= open(avg_file_name, "r")
	header= file_in.readline().strip().split("\t")
	index_output=[]
	for h in indices:
		index= header[h]
		index_output.append(index)
	out_dic={}

	for line in file_in:
		info= line.strip().split()
		gene= info[0]
		temp=[]
		for i in indices:
			entry= float(info[i])
			temp.append(entry)
		out_dic[gene]=temp

	file_in.close()
	return out_dic, index_output

####### calculating coefficient of variation
## for two data sets, normalize each value from one data set by the other
## take the mean and standard deviation of that set of ratios to determine the coefficient of variation
def cv(x,y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    ratios = []
    for i in range(n):
        x_i = x[i]
        y_i = y[i]
        ratio = x_i/y_i
        ratios.append(ratio)
    mean = average(ratios)
    sd = stats.stdev(ratios)
    coev = sd/mean
    sigfigs = '%.3f' % coev
    return sigfigs


## uses the cv function above to find coefficient of variation between each transcript in different dictionaries
def cv_two_dic(dic1, dic2):
    key1 = dic1.keys()
    key2 = dic2.keys()

    cv_dic = {}
    for k1 in key1:
        info1 = dic1[k1]
        out_cv = []
        for k2 in key2:
            info2 = dic2[k2]
            out_cv.append(cv(info1, info2))
        cv_dic[k1] = out_cv

    cv_header = 'gene'
    for k2 in key2:
        cv_header = cv_header + '\t' + k2
    cv_header = cv_header + '\n'

    return cv_dic, cv_header

## outputs a file with coefficients of variation -- not actually being used
def write_cv(out_dic, header, out_file_name):
	out_file= open(out_file_name, "w")
	out_file.write( header)

	out_keys= out_dic.keys()
	for o in out_keys:
		info= out_dic[o]
		out_file.write(o)
		for i in info:
			out_file.write("\t" + str(i))
		out_file.write("\n")
	out_file.close()

## outputs a file with coefficients of variation structured for plotting with R
def write_plot_cv(dic1, dic2, out_file_name, header=["gene1", "gene2", "cv"]):
	out_file= open(out_file_name, "w")
	header_line=""
	for h in header:
		header_line=header_line+"\t"+h
	out_file.write(header_line.strip()+ "\n")

	key1= dic1.keys()
	key2= dic2.keys()
	for k1 in key1:
		info1= dic1[k1]
		for k2 in key2:
			info2= dic2[k2]
			out_corr=str(cv(info1, info2))
			out_file.write(k1 +"\t"+ k2 +"\t"+ out_corr +"\n")
	out_file.close()