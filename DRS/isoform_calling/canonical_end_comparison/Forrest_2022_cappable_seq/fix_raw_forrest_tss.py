
def main():

	infile = 'cappable_seq_tss_forrest_raw.txt'
	outfile = 'cappable_seq_tss_forrest.txt'

	tss_list = []
	with open(infile, 'r') as f:
		f.readline()
		for line in f:
			info = line.strip().split('_')
			position = info[1]
			strand = info[2]
			tss_list.append((position, strand))

	with open(outfile, 'w') as o:
		o.write('Position\tStrand\n')
		for i in tss_list:
			o.write(str(i[0]) + '\t' + str(i[1]) + '\n')

if __name__ == '__main__':
	main()