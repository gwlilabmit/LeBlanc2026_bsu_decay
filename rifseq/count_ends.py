## this file contains functions for quantifying genes: counting 5' ends within annotated regions and getting average number of reads per base
## to add: subtraction and aggregation functions (when needed)

# reads a file containing CDS annotations and parses into a dictionary for read counting
# this command does it by alias (gene, rRNA, tRNA, etc) and NOT by the name of the gene
# returns a dictionary with gene names as keys and genomic coordinate information as values, and an index that describes the info in the values
def parse_CDS_file_by_type(file_CDS, criterium):
	file_in=open(file_CDS, 'r')
	header=file_in.readline()

	CDS_dictionary={}
	CDS_index=['start', 'stop', 'strand']
	for line in file_in:
		info=line.strip().split('\t')
		g_start=info[1]
		g_stop=info[2]
		g_strand=info[3]
		gene_name=info[5]
		alias=info[6]
		if(alias in criterium):
			CDS_dictionary[gene_name]=[g_start, g_stop, g_strand]

	file_in.close()
	return CDS_dictionary, CDS_index

# same as above, but parses CDS by gene name rather than broad category
# use if you're trying to get information from a custom annotation file (i.e. if you are quantifying regions that don't map exactly to annotated coding sequences)
def parse_CDS_file_by_name(file_CDS, criterium):
	file_in=open(file_CDS, 'r')
	header=file_in.readline()

	CDS_dictionary={}
	CDS_index=['start', 'stop', 'strand']
	for line in file_in:
		info=line.strip().split('\t')
		g_start=info[1]
		g_stop=info[2]
		g_strand=info[3]
		gene_name=info[5]
		alias=info[6]
		if(gene_name in criterium):
			CDS_dictionary[gene_name]=[g_start, g_stop, g_strand]

	file_in.close()
	return CDS_dictionary, CDS_index

# counts total number of reads between a given start and stop position in a wig file
# used below in get_total_CDS to get counts for each annotated unit parsed from a CDS file
def get_total_end_count(wig_name, start, stop):
	file_in=open(wig_name, 'r')
	file_in.readline()
	file_in.readline()

	total=0
	for line in file_in:
		pos,count=line.strip().split('\t')
		pos=int(pos)
		if(pos>=start and pos<=stop):
			total+=int(count)

	file_in.close()
	return total

# uses get_total_end_count to count 5' end mapping reads for each coding sequence from a parsed CDS file
# examples of wig front and wig rear:
# - wig front: "190710Li_D19-"
# - wig rear: "_NA_sequence-stripped_v2_bsubtilis_v1k1_ivt_spikein_v2k1.align_"
def get_total_CDS(CDS_dic, CDS_index, wig_index, wig_front, wig_rear):
	CDS_keys=CDS_dic.keys()
	out_dic={}
	for c in CDS_keys:
		info=CDS_dic[c]
		start=int(info[0])
		stop=int(info[1])
		strand=info[2]
		wig_name5=wig_front+wig_index+wig_rear
		if strand=='+':
			wig_name5+='f.wig'
		else:
			wig_name5+='r.wig'

		total=get_total_end_count(wig_name5, start, stop)

		out_dic[c]=CDS_dic[c]+[str(total)]

	out_index=CDS_index+[wig_index+'_5']
	return out_dic, out_index

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

# reads in a file with read counts and calculates average number of reads per based from the read count and the start/end coordinates
# performs filtering based on the bin_cut_off and first_cut_off inputs
# bin_cut_off dictates the minimum number of reads that must be present in each time point
# first_cut_off dictates the minimum number of reads that must be present at t=0min
# if a gene is under either cut off, it gets "zeroed out" meaning that its average is rounded down to 0 in the dictionary created by this function
def read_count_file(count_file_name, bin_cut_off=0, first_cut_off=0):
	file_in=open(count_file_name, 'r')
	header=file_in.readline().strip().split()
	out_dic={}

	for line in file_in:
		info=line.strip().split()
		gene=info[0]
		start=float(info[1])
		stop=float(info[2])
		gene_length=stop-start+1.0

		temp=[]
		first=float(info[4])
		for i in range(4, len(info)):
			entry=float(info[i])
			if(entry>=bin_cut_off):
				temp.append(entry/gene_length)
			else:
				temp.append(0)
		if(first>first_cut_off):
			out_dic[gene]=temp

	file_in.close()
	out_header=header[4:]
	return out_dic, out_header

## reads in a file with read counts and adds up the counts mapping to each gene
## to return the total number of mRNA-mapping counts at each time point
def get_bulk_counts(infile):
		f=open(infile, 'r')
		header=f.readline().strip().split('\t')
		header=header[4:]

		total = [0,0,0,0,0,0,0,0]
		for line in f:
			info = line.strip().split('\t')
			t0=int(info[4])
			t1=int(info[5])
			t2=int(info[6])
			t3=int(info[7])
			t4=int(info[8])
			t5=int(info[9])
			t6=int(info[10])
			t7=int(info[11])
			l=[t0,t1,t2,t3,t4,t5,t6,t7]
			total = [x + y for x, y in zip(total, l)]

		f.close()

		return header, total


## writes output of get_bulk_counts to a text file
def write_bulk_counts(sample_name, inheader, intotal, outfile):
	f = open(outfile, 'w')
	f.write('sample_name' + '\t')
	for i in inheader:
		f.write(str(i) + '\t')
	f.write('\n')

	f.write(sample_name + '\t')
	for i in intotal:
		f.write(str(i) + '\t')

	f.close()


## normalizes bulk mRNA from a write_bulk_counts text file by ssrA and rnaP averages
def normalize_ncRNA(count_file, norm_file):
		file_in=open(count_file, 'r')
		header=file_in.readline().strip().split()
		data_dic={}

		for line in file_in:
			info=line.strip().split()
			length=len(info)
			gene=info[0]
			t0 = float(info[1])
			t1 = float(info[2])
			t2 = float(info[3])
			t3 = float(info[4])
			t4 = float(info[5])
			t5 = float(info[6])
			t6 = float(info[7])
			if (length==9):
				t7 = float(info[8])
				t_list = [t0,t1,t2,t3,t4,t5,t6,t7]
				data_dic[gene] = t_list
			elif (length==8):
				t_list = [t0,t1,t2,t3,t4,t5,t6]
				data_dic[gene] = t_list

		norm_file=open(norm_file, 'r')
		norm_header=norm_file.readline()
		norm_dic={}
		for line in norm_file:
			info= line.strip().split()
			length=len(info)
			gene= info[0]
			t0 = float(info[1])
			t1 = float(info[2])
			t2 = float(info[3])
			t3 = float(info[4])
			t4 = float(info[5])
			t5 = float(info[6])
			t6 = float(info[7])
			if (length==9):
				t7 = float(info[8])
				t_list = [t0,t1,t2,t3,t4,t5,t6,t7]
				norm_dic[gene] = t_list
			elif (length==8):
				t_list = [t0,t1,t2,t3,t4,t5,t6]
				norm_dic[gene] = t_list


		out_dic={}
		for gene in data_dic.keys():
			out_dic[gene]=[]
			bulk_list=data_dic.get(gene)
			for i in range(len(bulk_list)):
				gene_i=bulk_list[i]
				ssrA_i=norm_dic['ssrA'][i]
				rnaP_i=norm_dic['rnaP'][i]
				(out_dic.get(gene)).append(gene_i/(ssrA_i+rnaP_i))

		norm_file.close()
		file_in.close()
		return out_dic


## same as above, but normalizes to spike-ins instead
def normalize_spike_in(count_file, norm_file):
		file_in=open(count_file, 'r')
		header=file_in.readline().strip().split()
		data_dic={}

		for line in file_in:
			info=line.strip().split()
			length=len(info)
			gene=info[0]
			t0 = float(info[1])
			t1 = float(info[2])
			t2 = float(info[3])
			t3 = float(info[4])
			t4 = float(info[5])
			t5 = float(info[6])
			t6 = float(info[7])
			if (length==9):
				t7 = float(info[8])
				t_list = [t0,t1,t2,t3,t4,t5,t6,t7]
				data_dic[gene] = t_list
			elif (length==8):
				t_list = [t0,t1,t2,t3,t4,t5,t6]
				data_dic[gene] = t_list

		norm_file=open(norm_file, 'r')
		norm_header=norm_file.readline()
		norm_dic={}
		for line in norm_file:
			info= line.strip().split()
			length=len(info)
			gene= info[0]
			t0 = float(info[1])
			t1 = float(info[2])
			t2 = float(info[3])
			t3 = float(info[4])
			t4 = float(info[5])
			t5 = float(info[6])
			t6 = float(info[7])
			if (length==9):
				t7 = float(info[8])
				t_list = [t0,t1,t2,t3,t4,t5,t6,t7]
				norm_dic[gene] = t_list
			elif (length==8):
				t_list = [t0,t1,t2,t3,t4,t5,t6]
				norm_dic[gene] = t_list


		out_dic={}
		for gene in data_dic.keys():
			out_dic[gene]=[]
			bulk_list=data_dic.get(gene)
			for i in range(len(bulk_list)):
				gene_i=bulk_list[i]
				fLuc_i = norm_dic['fLuc'][i]
				pdr5_i = norm_dic['pdr5'][i]
				txb1_i = norm_dic['txb1'][i]
				(out_dic.get(gene)).append(gene_i/(fLuc_i+pdr5_i+txb1_i))

		norm_file.close()
		file_in.close()
		return out_dic
