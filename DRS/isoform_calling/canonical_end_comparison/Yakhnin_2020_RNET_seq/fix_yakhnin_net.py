
def main():

	infile = 'yakhnin_raw.txt'
	outfile = 'rnet_yakhnin.txt'
	genome = '../B_subtilis_168_NC_000964.3.gff3'

	strands = get_strands(genome)

	## manual additions of genes with different names
	## via subtiwiki
	strands['yaaJ'] = '+'
	strands['ricT'] = '+'
	strands['ispE'] = '+'
	strands['gcaD'] = '+'
	strands['spoVC'] = '+'
	strands['yacF'] = '+'
	strands['rplGB'] = '+'
	strands['map'] = '+'
	strands['ybxA'] = '+'
	strands['pdaB'] = '+'
	strands['cdaA'] = '+'
	strands['ycbE'] = '+'
	strands['yckE'] = '+'
	strands['fpbN'] = '+'
	strands['amj'] = '+'
	strands['ndoAI'] = '+'
	strands['int'] = '-'
	strands['conG'] = '-'
	strands['ydhU/1'] = '-'
	strands['ydiG'] = '+'
	strands['ydzT/1'] = '-'
	strands['fruC'] = '+'
	strands['dgkB'] = '+'
	strands['yfnI'] = '+'
	strands['mtrA'] = '-'
	strands['fecF'] = '-'
	strands['yfjN'] = '+'
	strands['rpsNB'] = '+'
	strands['bcaP'] = '+'
	strands['bmrB'] = '+'
	strands['yhfR'] = '-'
	strands['appA/1'] = '+'
	strands['cwlQ'] = '-'
	strands['htpX'] = '+'
	strands['pfeT'] = '+'
	strands['papB'] = '+'
	strands['ylaG'] = '+'
	strands['ylbI'] = '+'
	strands['bshC'] = '+'
	strands['mraZ'] = '+'
	strands['ylyB'] = '+'
	strands['ylnD'] = '+'
	strands['ylnE'] = '+'
	strands['ylnF'] = '+'
	strands['remA'] = '+'
	strands['yloH'] = '+'
	strands['yloM'] = '+'
	strands['fapR'] = '+'
	strands['clpY'] = '+'
	strands['flgE'] = '+'
	strands['flhG'] = '+'
	strands['rasP'] = '+'
	strands['ymfK/1'] = '+'
	strands['rodZ'] = '+'
	strands['rny'] = '+'
	strands['fosB'] = '+'
	strands['plsY'] = '-'
	strands['gltB'] = '-'
	strands['azoR1'] = '-'
	strands['pghZ'] = '-'
	strands['yoyK/1'] = '+'
	strands['ypnP'] = '-'
	strands['sco'] = '-'
	strands['thyB'] = '-'
	strands['ypcP'] = '-'
	strands['bpsB'] = '-'
	strands['gpsB'] = '-'
	strands['asnS'] = '-'
	strands['bshA'] = '-'
	strands['menH'] = '-'
	strands['ypfD'] = '-'
	strands['yloI'] = '+'
	strands['odhA'] = '-'
	strands['motI'] = '-'
	strands['ribT'] = '-'
	strands['rnz'] = '-'
	strands['yidC2'] = '+'
	strands['yqiB'] = '-'
	strands['lipM'] = '-'
	strands['tapA'] = '-'
	strands['yqgP'] = '-'
	strands['ypfD'] = '-'
	strands['rpmGA'] = '-'
	strands['nfo'] = '-'
	strands['pgpH'] = '-'
	strands['sigKN'] = '+'
	strands['txpA'] = '+'
	strands['sknR'] = '+'
	strands['zinT'] = '-'
	strands['yrhJ'] = '-'
	strands['trmU'] = '-'
	strands['yrbF'] = '-'
	strands['obg'] = '-'
	strands['ysxC'] = '-'
	strands['zapA'] = '-'
	strands['ytaG'] = '-'
	strands['ytkK'] = '-'
	strands['thiI'] = '-'
	strands['sftA'] = '-'
	strands['murJ'] = '-'
	strands['floT'] = '-'
	strands['yubC'] = '+'
	strands['uppP'] = '-'
	strands['pncB'] = '-'
	strands['citG'] = '-'
	strands['catR'] = '+'
	strands['tpi'] = '-'
	strands['lutC'] = '-'
	strands['lutB'] = '-'
	strands['lutA'] = '-'
	strands['lutP'] = '+'
	strands['psdB'] = '-'
	strands['yvoE'] = '-'
	strands['lysP (yvsH)'] = '-'
	strands['cysI'] = '-'
	strands['swrAA/1'] = '-'
	strands['hpf'] = '-'
	strands['smiA'] = '-'
	strands['degV'] = '-'
	strands['tagV'] = '+'
	strands['tagT'] = '+'
	strands['spoIIT'] = '-'
	strands['maeA'] = '-'
	strands['lipL'] = '-'
	strands['hemQ'] = '+'
	strands['efeB'] = '-'
	strands['wapI'] = '-'
	strands['dra'] = '-'
	strands['walJ'] = '-'
	strands['gdpP'] = '-'
	strands['thdF'] = '-'
	strands['yycH'] = '-'
	strands['trmU'] = '-'

	peak_list = []
	with open(infile, 'r') as f:
		f.readline()
		for line in f:
			info = line.strip().split('\t')
			gene = info[0]
			rel_strand = info[-1] ## strand relative to genome annotation as reported in paper

			abs_strand = strands.get(gene) ## absolute strand (i.e. + or - strand)
			if abs_strand == '+':
				abs_strand = 'f'
			elif abs_strand == '-':
				abs_strand = 'r'
			if abs_strand == None:
				print(gene)
			peak_pos = info[1]
			peak_list.append((peak_pos, abs_strand))

	with open(outfile, 'w') as o:
		o.write('Position\tStrand\n')
		for i in peak_list:
			o.write(str(i[0]) + '\t' + str(i[1]) + '\n')

def get_strands(genome):
	strands = {}
	with open(genome, 'r') as f:
		f.readline()
		f.readline()
		f.readline()
		for line in f:
			info = line.strip().split('\t')
			biotype = info[2]
			if biotype == 'CDS':
				continue
			else:
				strand = info[6]
				tags = info[8].split(';')
				name = tags[2].split('=')[-1]
				strands[name] = strand

	return strands

if __name__ == '__main__':
	main()