from argparse import ArgumentParser
from multiprocessing import Process, Queue
import sys

# -------------------------------------------

def main():
	parser=ArgumentParser()
	parser.add_argument('--infile', help='input text file')
	parser.add_argument('--outfile', help='name of output text file')
	args = parser.parse_args()

	annotations=read_annotation('/home/leblanch/data/genomes/NC_000964.3/30_cut_CDS_168.txt')
	annotation_dictionary_list=split_dictionary(annotations, 375)

	output=mp_map_reads(annotation_dictionary_list, args.infile)

	header=['gene_name', 'start', 'end', 'subsetting_reads', 'overlapping_reads', 'fraction_full_length']
	write_dict(output, header, args.outfile)

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
		name=info[5]
		type=info[6]
		if type=='gene':
			annotations[name]=[start, end, strand]

	return annotations

# credit: https://gist.github.com/nz-angel/31890d2c6cb1c9105e677cacc83a1ffd
def split_dictionary(input_dict, chunk_size):
    res = []
    new_dict = {}
    for k, v in input_dict.items():
        if len(new_dict) < chunk_size:
            new_dict[k] = v
        else:
            res.append(new_dict)
            new_dict = {k: v}
    res.append(new_dict)
    return res

def read_infile(infile):
	reads_list=[]

	f=open(infile, 'r')
	f.readline()
	for line in f:
		read_info=line.strip().split('\t')
		read_strand=read_info[0]
		read_start=int(read_info[1]) # comes out correct
		length=int(read_info[2])
		out_info=[read_strand, read_start, length]
		reads_list.append(out_info)

	return reads_list

# before running this, need to split up the annotation dictionary into 12 chunks
# modified from initial version to work with multiprocessing
def map_reads(q, reads_list, annotation_chunk):
	#reads_list=read_infile(infile)
	annotation_keys=annotation_chunk.keys()

	for gene in annotation_keys:
		info=annotation_chunk.get(gene)
		gene_start=int(info[0])
		gene_end=int(info[1])
		gene_strand=info[2]
		gene_range=range(gene_start,gene_end)

		subset_check=0
		overlap_check=0

		for read in reads_list:
			read_strand=read[0]

			if read_strand!=gene_strand:
					continue # don't want to count any reads on the opposite strand and this should save some work

			read_start=int(read[1]) # comes out correct
			length=int(read[2])
			read_end=read_start+length # correct but will lose 1 nt when put into a range

			read_range=range(read_start, read_end+1) # should be matching the way the read aligns

			overlap=range_overlap(read_range, gene_range)
			if overlap==True:
				overlap_check+=1

				subset=range_subset(read_range, gene_range)
				if subset==True:
					subset_check+=1

		if overlap_check>0:
			fraction_full_length=subset_check / overlap_check
			out={gene:[gene_start, gene_end, subset_check, overlap_check, fraction_full_length]}
			q.put(out)

	q.put(None)
	print('chunk complete')

def mp_map_reads(annotation_dictionary_list, infile):
	q=Queue()
	process_list=[]
	reads_list=read_infile(infile)
	for i in range(len(annotation_dictionary_list)):
		p=Process(target=map_reads, args=(q, reads_list, annotation_dictionary_list[i]))
		p.start()
		process_list.append(p)
		print(i)
		sys.stdout.flush()

	output={}
	count=0
	while any([Process.is_alive(p) for p in process_list])==True:
		out=q.get()
		if type(out)==dict:
			output=output | out
		elif out==None:
			count+=1
		if count==len(annotation_dictionary_list):
			for p in process_list:
				p.join()
			print('all genes calculated')

	return output

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


# -------------------------------------------

if __name__=="__main__":
    main()