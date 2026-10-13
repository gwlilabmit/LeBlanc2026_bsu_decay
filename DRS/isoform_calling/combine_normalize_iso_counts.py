import polars as pl
import numpy as np
import sys

def main():

	timepoint_list = ['t0', 't05', 't2', 't4', 't8', 't12', 't16']
	timepoint_dic={}

	for timepoint in timepoint_list:
		print('reading ' + timepoint)
		infile = 'isoform_counts_by_timepoint/'+timepoint+'_iso_counts.txt'
		timepoint_df = pl.read_csv(infile, separator='\t')
		timepoint_dic[timepoint] = timepoint_df

	iso_dic={}
	for transcript in timepoint_dic['t0'].iter_rows(named=True):
		transcript_id = transcript['transcript_id']
		strand = transcript['strand']
		cover_genes = transcript['covered_genes'].split(',')
		overlap_genes = transcript['overlapping_genes'].split(',')
		starts = list(map(int, transcript['associated_starts'].split(',')))
		ends = list(map(int, transcript['associated_ends'].split(',')))
		called_timepoints = transcript['timepoints']

		iso_dic[transcript_id]=Isoform(transcript_id, strand, cover_genes, overlap_genes, starts, ends, called_timepoints)

	norm_genes_overlap={}
	print('getting normalization factors')
	with open('../normalization/norm_gene_overlap.txt', 'r') as f:
		f.readline()
		for line in f:
			info=line.strip().split()
			gene_name=info[0]
			temp=info[1:]
			temp=[int(x) for x in temp]
			overlap_counts=np.array(temp)
			norm_genes_overlap[gene_name]=overlap_counts

	norm_factors=norm_genes_overlap['ssrA'] + norm_genes_overlap['rnpB']

	iso_count_dic={}
	print('normalizing')
	for i in iso_dic.keys():
		isoform = iso_dic.get(i)
		for index, timepoint in enumerate(timepoint_list):
			timepoint_df = timepoint_dic.get(timepoint)
			count = timepoint_df.filter(pl.col('transcript_id') == i).select(pl.col('count')).item()
			isoform.add_count(index, count)

		normed = isoform.counts_array / norm_factors
		isoform.norm_array = normed
		iso_count_dic[i] = isoform

	with open('combined_isoform_quant/isoform_counts.txt', 'w') as o:
		print('writing counts')
		o.write('isoform_id\tstrand\tcovered_genes\toverlapping_genes\tassociated_starts\tassociated_ends\tcalled_timepoints\tt0\tt0.5\tt2\tt4\tt8\tt12\tt16\n')
		for i in iso_count_dic.keys():
			isoform = iso_count_dic.get(i)
			o.write(str(isoform.iso_id))
			o.write('\t')
			o.write(str(isoform.strand))
			o.write('\t')
			o.write(','.join(isoform.cover_genes))
			o.write('\t')
			o.write(','.join(isoform.overlap_genes))
			o.write('\t')
			o.write(','.join(list(map(str, isoform.starts))))
			o.write('\t')
			o.write(','.join(list(map(str, isoform.ends))))
			o.write('\t')
			o.write(isoform.called_timepoints)
			o.write('\t')
			o.write('\t'.join(list(map(str, isoform.counts_array))))
			o.write('\n')

	with open('combined_isoform_quant/isoform_norm_vals.txt', 'w') as o:
		print('writing norm vals')
		o.write('isoform_id\tstrand\tcovered_genes\toverlapping_genes\tassociated_starts\tassociated_ends\tcalled_timepoints\tt0\tt0.5\tt2\tt4\tt8\tt12\tt16\n')
		for i in iso_count_dic.keys():
			isoform = iso_count_dic.get(i)
			o.write(str(isoform.iso_id))
			o.write('\t')
			o.write(str(isoform.strand))
			o.write('\t')
			o.write(','.join(isoform.cover_genes))
			o.write('\t')
			o.write(','.join(isoform.overlap_genes))
			o.write('\t')
			o.write(','.join(list(map(str, isoform.starts))))
			o.write('\t')
			o.write(','.join(list(map(str, isoform.ends))))
			o.write('\t')
			o.write(isoform.called_timepoints)
			o.write('\t')
			o.write('\t'.join(list(map(str, isoform.norm_array))))
			o.write('\n')


class Isoform:
	def __init__(self, iso_id, strand, cover_genes, overlap_genes, starts, ends, called_timepoints):
		self.iso_id=iso_id
		self.strand=strand
		self.cover_genes=cover_genes
		self.overlap_genes=overlap_genes
		self.starts=starts
		self.ends=ends
		self.called_timepoints=called_timepoints
		self.counts_array=np.zeros(7, dtype=np.int64)
		self.norm_array=np.zeros(7, dtype=np.int64)

	def add_count(self, index, count):
		self.counts_array[index]=count

if __name__ == '__main__':
	main()