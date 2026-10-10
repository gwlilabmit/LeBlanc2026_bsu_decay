from argparse import ArgumentParser

# -------------------------------------------

def main():
    
    parser=ArgumentParser()
    parser.add_argument('--infile', help='input text file')
    parser.add_argument('--outfile', help='name of output text file')
    args = parser.parse_args()
    
    ssrA=extract_annotation('../genomes/NC_000964.3/CDS_168.mochiview_v13.txt', 'ssrA')
    rnpB=['rnpB', 2331320, 2331720, '-']
    #rnpB=extract_annotation('CDS_168.mochiview_v13.txt', 'rnpB')
    
    o=open(args.outfile,'w')
    o.write('gene_name\tstart\tend\ttotal_mapping_reads\tsubsetting_reads\toverlapping_reads\tfraction_full_length\n')
    o.close()
    
    map_reads_by_gene(args.infile, ssrA, args.outfile)
    map_reads_by_gene(args.infile, rnpB, args.outfile)

# -------------------------------------------

def extract_annotation(annotation_file, target_name):
    f=open(annotation_file, 'r')
    f.readline()
    for line in f:
        info=line.strip().split('\t')
        start=int(info[1])
        end=int(info[2])
        strand=info[3]
        name=info[5]
        if name==target_name:
            annotation=[target_name, start, end, strand]
            return annotation
            break
        else:
            continue

def map_reads_by_gene(infile, gene, outfile):
    
    gene_name=gene[0]
    gene_start=gene[1]
    gene_end=gene[2]
    gene_strand=gene[3]
    
    gene_range=range(gene_start, gene_end)
    
    subset_check=0
    overlap_check=0
    
    f=open(infile, 'r')
    f.readline()
    for line in f:
        read_info=line.strip().split('\t')
        read_strand=read_info[0]

        if read_strand!=gene_strand:
            continue # don't want to count any reads on the opposite strand and this should save some work       
        
        read_start=int(read_info[1]) # double check how this works out in alignment
            
        length=int(read_info[2])
        read_end=read_start+length # double check how this works out in alignment
        
        read_range=range(read_start, read_end)
        
        overlap=range_overlap(read_range, gene_range)
        if overlap==True:
            overlap_check+=1

        subset=range_subset(read_range, gene_range)
        if subset==True:
            subset_check+=1
                
    f.close()
        
    fraction_full_length=subset_check / overlap_check

    if overlap_check>0:
        o=open(outfile,'a')
        o.write(str(gene_name)+'\t'+str(gene_start)+'\t'+str(gene_end)+'\t'+'\t'+str(subset_check)+'\t'+str(overlap_check)+'\t'+str(fraction_full_length)+'\n')
        o.close()

    print(gene_name + ' done')
            
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


# -------------------------------------------

if __name__=="__main__":
    main()