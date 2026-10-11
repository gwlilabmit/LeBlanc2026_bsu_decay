import sys

def main():

	t0='../streamlined_sam_files/R25-2578.txt'
	t05='../streamlined_sam_files/R25-2579-all.txt'
	#t125='streamlined_sam_files/R25-2580.txt'
	t2='../streamlined_sam_files/R25-2581-all.txt'
	t4='../streamlined_sam_files/R25-2582.txt'
	t8='../streamlined_sam_files/R25-2583.txt'
	t12='../streamlined_sam_files/R25-2584.txt'
	t16='../streamlined_sam_files/R25-2585.txt'
	
	t0_counts=read_count(t0)
	print('t0')
	sys.stdout.flush()

	t05_counts=read_count(t05)
	print('t05')
	sys.stdout.flush()

	# t125_counts=read_count(t125)
	# print('t125')
	# sys.stdout.flush()

	t2_counts=read_count(t2)
	print('t2')
	sys.stdout.flush()

	t4_counts=read_count(t4)
	print('t4')
	sys.stdout.flush()

	t8_counts=read_count(t8)
	print('t8')
	sys.stdout.flush()

	t12_counts=read_count(t12)
	print('t12')
	sys.stdout.flush()

	t16_counts=read_count(t16)
	print('t16')
	sys.stdout.flush()

	#counts_list=[t0_counts, t05_counts, t125_counts, t2_counts, t4_counts, t8_counts, t12_counts, t16_counts]
	counts_list=[t0_counts, t05_counts, t2_counts, t4_counts, t8_counts, t12_counts, t16_counts]

	overlap_dictionary={'ssrA':[], 'rnpB':[]}
	subset_dictionary={'ssrA':[], 'rnpB':[]}

	for i in counts_list:
		ssrA=i[0]
		rnpB=i[1]

		overlap_dictionary['ssrA'].append(ssrA.overlap)
		overlap_dictionary['rnpB'].append(rnpB.overlap)

		subset_dictionary['ssrA'].append(ssrA.subset)
		subset_dictionary['rnpB'].append(rnpB.subset)


	print('counts compiled')

	header=['gene', '0min', '0.5min', '2min', '4min', '8min', '12min', '16min']
	write_dict(overlap_dictionary, header, 'norm_gene_overlap.txt')
	write_dict(subset_dictionary, header, 'norm_gene_subset.txt')

	print('files written')

	# normalize overlapping and subsetting counts for each gene by each normalization approach


#################################

def range_subset(read_range, gene_range): # check if one range (gene) is a subset of another range (read)
	# read start must be less than gene start
	# read end must be greater than gene end
	r1, r2=read_range.start, read_range.stop # read range
	g1, g2=gene_range.start, gene_range.stop # gene range
	return r1<=g1 and r2>=g2 # should return True if the gene is a subset of the read – tested and confirmed

def range_overlap(read_range, gene_range): # check if one range (gene) overlaps at all with another range (read)
	# to not overlap, either range1 is entirely less than range2 or range1 is entirely greater than range2
	r1, r2=read_range.start, read_range.stop
	g1, g2=gene_range.start, gene_range.stop
	# ranges do not overlap if r2<g1 (end of read is less than start of gene)
	# ranges do not overlap if g2<r1 (end of gene is less than start of read)
	# return not (r2<=g1 and g2<=r1)
	# https://nedbatchelder.com/blog/201310/range_overlap_in_two_compares.html
	return r2>=g1 and g2>=r1 # should return True if the ranges overlap – tested and confirmed

class Gene:
	def __init__(self, name, start, end):
		self.name=name
		self.start=start
		self.end=end
		self.overlap=0
		self.subset=0

def read_count(infile):
	# actual annotations from subtiwiki
	# doesn't seem to be a good representation of where the reads actually are
	# ssrA=Gene('ssrA', 3450710, 3451071)
	# rnpB=Gene('rnpB', 2331330, 2331715)

	ssrA=Gene('ssrA', 3450712, 3451021)
	rnpB=Gene('rnpB', 2331370, 2331670)

	ssrA_range=range(ssrA.start, ssrA.end)
	rnpB_range=range(rnpB.start, rnpB.end)

	with open(infile, 'r') as f:
		f.readline()
		for line in f:
			read_info=line.strip().split('\t')
			read_strand=read_info[0]
			if read_strand == '+':
				continue
			elif read_strand == '-':
				read_start=int(read_info[1])
				if read_start < 2000000:
					continue
				elif read_start > 4000000:
					continue
				else:
					read_length=int(read_info[2])
					read_end=read_start+read_length+1
					read_range=range(read_start, read_end)

					ssrA_overlap=range_overlap(read_range, ssrA_range)
					rnpB_overlap=range_overlap(read_range, rnpB_range)

					if ssrA_overlap==True:
						ssrA.overlap+=1
						ssrA_subset=range_subset(read_range, ssrA_range)
						if ssrA_subset==True:
							ssrA.subset+=1

					if rnpB_overlap==True:
						rnpB.overlap+=1
						rnpB_subset=range_subset(read_range, rnpB_range)
						if rnpB_subset==True:
							rnpB.subset+=1

	return [ssrA, rnpB]

def write_dict(dictionary, header, outfile):
	with open(outfile, 'w') as o:
		o.write('\t'.join(header) + '\n')
		for key in dictionary.keys():
			val=dictionary.get(key)
			if type(val)==list:
				val=[str(i) for i in val]
				val='\t'.join(val)
			elif type(val)==float:
				val=str(val)
			elif type(val)==int:
				val=str(val)
			o.write(str(key) + '\t' + val + '\n')


if __name__ == '__main__':
	main()