import sys

# -------------------------------------------

def main():

	original_annotation = '/home/leblanch/data/genomes/NC_000964.3/CDS_168.mochiview_v13.txt'

	annotation = read_annotation(original_annotation)

	sd_buffer = 11 ## number of nucleotides used for conservative estimate of functional Shine-Dalgarno sequence
	
	new_annotation = {}
	for gene in annotation.keys():
		info = annotation.get(gene)
		strand = info[2]

		left = info[0]
		right = info[1]
		
		if strand == '+':
			new_left = left - sd_buffer
			new_right = right

		elif strand == '-':
			new_left = left
			new_right = right + sd_buffer

		feature_name=info[3]
		gene_name=info[4]
		type=info[5]

		new_annotation[gene] = [new_left, new_right, strand, feature_name, gene_name, type]

	with open('functional_CDS_168.txt', 'w') as o:
		o.write('SEQ_NAME\tSTART\tEND\tSTRAND\tFEATURE_NAME\tGENE_NAME\n')
		for gene in new_annotation.keys():
			o.write('NC_000964.3\t')
			info = new_annotation.get(gene)
			info = [str(i) for i in info]
			o.write('\t'.join(info) + '\n')

	print('functional_CDS_168.txt written')

# -------------------------------------------

def read_annotation(annotation_file):
	f=open(annotation_file, 'r')
	f.readline()
	annotations={}
	for line in f:
		info=line.strip().split('\t')
		start=int(info[1])
		end=int(info[2])
		strand=info[3]
		feature_name=info[4]
		gene_name=info[5]
		type=info[6]
		if type=='gene':
			annotations[gene_name]=[start, end, strand, feature_name, gene_name, type]

	return annotations


# -------------------------------------------

if __name__=="__main__":
    main()