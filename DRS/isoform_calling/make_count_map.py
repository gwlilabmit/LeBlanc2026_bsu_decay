import polars as pl
import numpy as np
import os
import sys

## point of this script: make "end map" file WITH:
## rRNAs, tRNAs, highly expressed ncRNAs (e.g. ssrA) removed -- i.e. no reads that overlap these annotations
## intervening regions broken into different chunks to support parallelization
## so that section breaks are in places without isoforms of interest to avoid any edge effects

def read_genome(infile):
	genome={}
	to_remove={}
	with open(infile, 'r') as f:
		f.readline()
		for line in f:
			info=line.strip().split()
			start=int(info[1])
			end=int(info[2])
			strand=info[3]
			name=info[5]
			type=info[6]
			genome[name]=[start, end, strand, type]

			if type=='rRNA':
				to_remove[name]=[start, end, strand, type]
			elif type=='tRNA':
				to_remove[name]=[start, end, strand, type]
			elif type=='tmRNA':
				to_remove[name]=[start, end, strand, type]

	return genome, to_remove

def range_overlap(read_range, gene_range): # check if one range (gene) overlaps at all with another range (read)
	# to not overlap, either range1 is entirely less than range2 or range1 is entirely greater than range2
	r1, r2=read_range.start, read_range.stop
	g1, g2=gene_range.start, gene_range.stop
	# ranges do not overlap if r2<g1 (end of read is less than start of gene)
	# ranges do not overlap if g2<r1 (end of gene is less than start of read)
	# return not (r2<=g1 and g2<=r1)
	# https://nedbatchelder.com/blog/201310/range_overlap_in_two_compares.html
	return r2>=g1 and g2>=r1 # should return True if the ranges overlap – tested and confirmed

def main():
	
	timepoint = sys.argv[1]
	in_strand = sys.argv[2] ## must be + or -

	if in_strand == 'f':
		strand='+'
	elif in_strand == 'r':
		strand='-'
	else:
		strand=in_strand

	input_path = '~/data/isoform_calling/streamlined_sam_files/R25-'
	timepoint_keys = {'t0':'2578', 't05':'2579-all', 't2':'2581-all', 't4':'2582', 't8':'2583', 't12':'2584', 't16':'2585'}

	infile = input_path + timepoint_keys[timepoint] + '.txt'

	reads = pl.read_csv(infile, separator='\t')
	reads = reads.filter(pl.col('strand') == strand)
	## removed this line because input file now already has the rightmost coordinate included
	# reads = reads.with_columns((pl.col('leftmost_coordinate') + pl.col('length')).alias('rightmost_coordinate'))

	genome_length=4215606
	genome_coord_list=list(range(1, genome_length+1))

	contact_dict=dict((x,{}) for x in genome_coord_list)

	## this is the loop from the other code that uses the duckdbs, so the conversion to rightmost has already happened
	for read in reads.iter_rows():
		leftmost=read[1]
		rightmost=read[3]
		left_key=contact_dict.get(leftmost)
		if rightmost in left_key.keys():
			val=left_key.get(rightmost)
			new_val=val+1
			contact_dict[leftmost][rightmost]=new_val
		else:
			contact_dict[leftmost][rightmost]=1

	g, to_remove=read_genome('CDS_168.mochiview_v13.txt')
	to_remove['rnpB']=[2331330, 2331715, '-', 'ncRNA']
	to_remove['scr']=[26533, 26635, '+', 'ncRNA']
	to_remove['bsrA']=[2814490, 2814682, '-', 'ncRNA']
	to_remove['bsrB']=[2095927, 2096113, '-', 'ncRNA']

	ncRNAs=[]
	for i in to_remove.keys():
		temp=to_remove.get(i)
		left=temp[0]
		right=temp[1]
		strand_i=temp[2]
		if strand_i==strand:
			ncRNAs.append((left,right))

	ncRNAs=sorted(ncRNAs)

	## defining regions between highly-expressed ncRNAs
	## only works if input list is sorted (the sorted function should sort by the first element in the tuples I think)
	temp_sections=[]
	edge=0
	for i in ncRNAs:
		left=i[0]
		right=i[1]
		potential_section=(edge, left)
		temp_sections.append(potential_section)
		edge=right

	final_section=(edge, 4215606)
	temp_sections.append(final_section)

	regions=[]
	for i in temp_sections:
		size=i[1]-i[0]
		if size<200:
			continue
		elif size>=200:
			regions.append(i)

	output_list=[]
	for left in contact_dict.keys():
		contacts=contact_dict.get(left)
		if len(contacts)>0:
			for right in contacts.keys():
				contact_range=range(left, right+1)
				## determine if the contacts overlap with any of the ncRNA annotations
				## if so add it to the output list as an ncRNA contact
				check=0 ## keeps track of if this particular contact overlaps with at least one ncRNA annotation
				for annotation in ncRNAs:
					if check==1:
						break
					l=annotation[0]
					r=annotation[1]
					anno_range=range(l,r)
					overlap_bool=range_overlap(contact_range, anno_range)
					if overlap_bool==True:
						count=contacts.get(right)
						output_list.append(['ncRNA', str(left), str(right), str(count)])
						check+=1
						break
					else:
						continue

				if check>=1: ## if an overlap was found with an ncRNA, no further consideration of this contact
					break
				if check==0: ## if no overlaps were found with ncRNAs, start checking for overlap with other regions
					for re in range(len(regions)):
						reg=regions[re]
						region_range=range(reg[0], reg[1]+1)
						overlap_bool=range_overlap(contact_range, region_range)
						if overlap_bool==True:
							count=contacts.get(right)
							output_list.append(['sec'+str(re), str(left), str(right), str(count)])
						else:
							continue

	if strand == '+':
		out_strand = 'f'
	elif strand == '-':
		out_strand = 'r'

	with open('count_maps/'+timepoint+'_'+out_strand+'.txt', 'w') as o:
		for i in output_list:
			o.write('\t'.join(i) + '\n')



if __name__ == '__main__':
	main()