import sys
import numpy as np

def main():

	t0='../mRNA_overlap/R25-2578.txt'
	t05='../mRNA_overlap/R25-2579-all.txt'
	#t125='../mRNA_overlap/R25-2580.txt'
	t2='../mRNA_overlap/R25-2581-all.txt'
	t4='../mRNA_overlap/R25-2582.txt'
	t8='../mRNA_overlap/R25-2583.txt'
	t12='../mRNA_overlap/R25-2584.txt'
	t16='../mRNA_overlap/R25-2585.txt'

	timepoint_list=[t0, t05, t2, t4, t8, t12, t16]

	header=['gene', '0min', '0.5min', '2min', '4min', '8min', '12min', '16min']

	# initialize list of Gene objects with all genes at 0 min

	# initialize dictionary to put Gene objects in
	gene_dictionary={}
	with open(t0, 'r') as f:
		f.readline()
		for line in f:
			info=line.strip().split()
			gene_name=info[0]
			gene_start=info[1]
			gene_end=info[2]
			gene_subset=info[3]
			gene_overlap=info[4]

			# add gene object with start and end information to gene dictionary
			# add subset and overlap information for 0 min timepoint
			gene_i=Gene(gene_start, gene_end)
			gene_i.add_subset(0,gene_subset)
			gene_i.add_overlap(0,gene_overlap)

			gene_dictionary[gene_name]=gene_i

	# update dictionary with read counts for each subsequent timepoint
	# adds subsetting and overlapping read counts to appropriate array associated with each Gene object
	ind=0
	for timepoint in timepoint_list:
		with open(timepoint, 'r') as f:
			f.readline()
			for line in f:
				info=line.strip().split()
				gene_name=info[0]
				gene_subset=info[3]
				gene_overlap=info[4]

				gene_object=gene_dictionary.get(gene_name)
				if gene_object != None:
					gene_object.add_subset(ind, gene_subset)
					gene_object.add_overlap(ind, gene_overlap)

		ind+=1

	overlap_count_dictionary={}
	subset_count_dictionary={}
	for gene in gene_dictionary.keys():
		obj=gene_dictionary.get(gene)

		overlap_count_dictionary[gene]=list(obj.overlap)
		subset_count_dictionary[gene]=list(obj.subset)

	write_dict(overlap_count_dictionary, header, 'overlap_counts.txt')
	write_dict(subset_count_dictionary, header, 'subset_counts.txt')


	# import ssrA and rnpB subset and overlap read counts as arrays
	
	norm_genes_overlap={}
	with open('norm_gene_overlap.txt', 'r') as f:
		f.readline()
		for line in f:
			info=line.strip().split()
			gene_name=info[0]
			temp=info[1:]
			temp=[int(x) for x in temp]
			overlap_counts=np.array(temp)
			norm_genes_overlap[gene_name]=overlap_counts

	norm_genes_subset={}
	with open('norm_gene_subset.txt', 'r') as f:
		f.readline()
		for line in f:
			info=line.strip().split()
			gene_name=info[0]
			temp=info[1:]
			temp=[int(x) for x in temp]
			subset_counts=np.array(temp)
			norm_genes_subset[gene_name]=subset_counts

	# sum read counts from ssrA and rnpB to get final normalization factors
	overlap_norm_factors=norm_genes_overlap['ssrA'] + norm_genes_overlap['rnpB']
	subset_norm_factors=norm_genes_subset['ssrA'] + norm_genes_subset['rnpB']

	overlap_dictionary={}
	subset_dictionary={}
	sbo_dictionary={}
	for gene in gene_dictionary.keys():
		obj=gene_dictionary.get(gene)

		overlap_norm=obj.overlap / overlap_norm_factors
		overlap_dictionary[gene]=list(overlap_norm)

		subset_norm= obj.subset / subset_norm_factors
		subset_dictionary[gene]=list(subset_norm)

		subset_by_overlap_norm=obj.subset / overlap_norm_factors
		sbo_dictionary[gene]=list(subset_by_overlap_norm)

	write_dict(overlap_dictionary, header, 'overlap_norm_vals.txt')
	write_dict(subset_dictionary, header, 'subset_norm_vals.txt')
	write_dict(sbo_dictionary, header, 'sbo_norm_vals.txt')

#################################

class Gene:
	def __init__(self, start, end):
		self.start=start
		self.end=end
		self.subset=np.zeros(7)
		self.overlap=np.zeros(7)

	def add_subset(self, index, value):
		self.subset[index]=value

	def add_overlap(self, index, value):
		self.overlap[index]=value

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


#################################

if __name__ == '__main__':
	main()