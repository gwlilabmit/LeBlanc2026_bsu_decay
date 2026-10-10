from argparse import ArgumentParser
from itertools import groupby

# -------------------------------------------

def main():
	parser=ArgumentParser()
	parser.add_argument('--infile', help='input sam file')
	parser.add_argument('--outfile', help='name of output text file')
	args = parser.parse_args()

	shrink_sam(args.infile, args.outfile)

# -------------------------------------------

def shrink_sam(infile, outfile):
	f=open(infile, 'r')
	# commented out because these files don't have headers
	# f.readline()
	# f.readline()
	# f.readline()

	o=open(outfile, 'w')
	o.write('strand\tleftmost_coordinate\tlength\trightmost_coordinate\n')

	for line in f:
		info=line.strip().split('\t')
		flag=info[1]
		if flag=='0':
			strand='+'
		elif flag=='16':
			strand='-'
		else:
			continue
		map_pos=int(info[3])
		cigar_string=info[5]
		length=query_len(cigar_string) # calculates length from cigar string based on number of matches/mismatches & deletions, ignores softclipped and inserted nts
		right_end = map_pos + length-1 # subtracting 1 should avoid off by 1 issue

		sequence = info[9]
		trimmed_map_pos, trimmed_length, trimmed_right_end = trim_A(map_pos, cigar_string, length, right_end, sequence, strand)

		o.write(strand + '\t' + str(trimmed_map_pos) + '\t' + str(trimmed_length) + '\t' + str(trimmed_right_end) + '\n')

	o.close()
	f.close()

## since approach for trimming will depend on strand, should do all of it in this function and have the right/length/left be the outputs of this function
def trim_A(map_pos, cigar_string, length, right_end, sequence, strand):
	if strand == '+':
		if cigar_string[-1] == 'S':
			## allows for soft clipping up to 9999 (expect small number of >100 clips, accounts for >1000 clips although I'm not sure there are any)
			try:
				soft_clip = int(cigar_string[-5:-1])
			except ValueError:
				try:
					soft_clip = int(cigar_string[-4:-1])
				except ValueError:
					try:
						soft_clip = int(cigar_string[-3:-1])
					except ValueError:
						soft_clip = int(cigar_string[-2:-1])
			clipped_seq = sequence[:-soft_clip]
		else:
			soft_clip = 0
			clipped_seq = sequence

		### clipped seq should now be only the actually *aligned* sequence, not anything additional (i.e. extra polyA tail)

		A_count = 0
		for i in reversed(clipped_seq):
			if i == 'A':
				A_count+=1
			elif i!='A':
				break

		## clip As from right end
		length = length - A_count
		right_end = right_end - A_count

		return map_pos, length, right_end

	elif strand == '-':
		try:
			soft_clip = int(cigar_string[:4])
			pos=4
		except ValueError:
			try:
				soft_clip = int(cigar_string[:3])
				pos=3
			except ValueError:
				try:
					soft_clip = int(cigar_string[:2])
					pos=2
				except ValueError:
					soft_clip = int(cigar_string[:1])
					pos=1
		
		## accounting for cases where the first section is not soft clipped – not sure if this is something that will come up
		if cigar_string[pos] != 'S':
			soft_clip = 0

		clipped_seq = sequence[soft_clip:]

		## As are recorded as Ts because reverse-stranded reads are reverse complemented in SAM files
		T_count = 0
		for i in clipped_seq:
			if i == 'T':
				T_count+=1
			elif i != 'T':
				break

		## clip As (Ts) from left end
		length = length - T_count
		map_pos = map_pos + T_count

		return map_pos, length, right_end

# [here](https://stackoverflow.com/questions/39710796/infer-the-length-of-a-sequence-using-the-cigar)
def query_len(cigar_string):
	real_bases=('M', 'D') # cigar string operators for bases that are actually present in the reference genome (i.e. no clipped or insert nts)
	result=0
	cig_iter=groupby(cigar_string, lambda chr: chr.isdigit())
	for _, length_digits in cig_iter:
		length = int(''.join(length_digits))
		op = next(next(cig_iter)[1])
		if op in real_bases:
			result += length

	return result
# -------------------------------------------

if __name__=="__main__":
    main()