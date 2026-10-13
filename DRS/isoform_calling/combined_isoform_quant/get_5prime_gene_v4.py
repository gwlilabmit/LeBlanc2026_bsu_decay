import polars as pl
from Bio.Seq import Seq
import regex

def main():
	increasers = pl.read_csv('increaser_fits.txt', separator='\t')
	annotations = pl.read_csv('../CDS_168.mochiview_v13.txt', separator='\t')
	annotations_f = annotations.filter((pl.col('STRAND')=='+') & (pl.col('ALIASES')=='gene'))
	annotations_r = annotations.filter((pl.col('STRAND')=='-') & (pl.col('ALIASES')=='gene'))

	genes_list = []

	sch = increasers.schema
	gene_iso_df = pl.DataFrame(schema=sch)
	non_gene_iso_df = pl.DataFrame(schema=sch)

	for i in increasers.iter_rows(named=True):
		if i['strand'] == 'f':
			ends = [int(x) for x in i['associated_starts'].split(',')]
			end = min(ends)
			end_window = list(range(end-5, end+6))
			gene = annotations_f.filter(pl.col('START').is_in(end_window))
			if gene.shape[0]>0:
				genes_list.append(gene.select(pl.col('GENE_NAME')).item())
				gene_iso_df=pl.concat([gene_iso_df, pl.DataFrame(i)], how='vertical')
			else:
				non_gene_iso_df=pl.concat([non_gene_iso_df, pl.DataFrame(i)], how='vertical')

		elif i['strand'] == 'r':
			ends = [int(x) for x in i['associated_ends'].split(',')]
			end = max(ends)
			end_window = list(range(end-5, end+6))
			gene = annotations_r.filter(pl.col('END').is_in(end_window))
			if gene.shape[0]>0:
				genes_list.append(gene.select(pl.col('GENE_NAME')).item())
				gene_iso_df=pl.concat([gene_iso_df, pl.DataFrame(i)], how='vertical')
			else:
				non_gene_iso_df=pl.concat([non_gene_iso_df, pl.DataFrame(i)], how='vertical')

	genes_df=pl.DataFrame(genes_list)
	genes_df.write_csv('increaser_genes_list.txt', separator='\t')

	gene_iso_df.write_csv('increaser_genes_data.txt', separator='\t')
	non_gene_iso_df.write_csv('increaser_non_genes_data.txt', separator='\t')

	with open('B_subtilis_168_NC_000964.3.fa', 'r') as g:
		g.readline()
		genome = g.readline()

	gene_seqs=[]
	with open('gene_upstream_seqs.txt', 'w') as o:

		for i in gene_iso_df.iter_rows(named=True):
			iso_id = i['isoform_id']
			
			if i['strand'] == 'f':
				ends = [int(x) for x in i['associated_starts'].split(',')]
				end = min(ends)
				upstream = genome[end-21:end+19]

				o.write(upstream + '\n')
				gene_seqs.append(upstream)

			elif i['strand'] == 'r':
				ends = [int(x) for x in i['associated_ends'].split(',')]
				end = max(ends)
				upstream = genome[end-20:end+20]
				upstream=Seq(upstream)
				upstream_rc = upstream.reverse_complement()

				o.write(str(upstream_rc )+ '\n')
				gene_seqs.append(str(upstream_rc))

	non_gene_seqs=[]
	with open('non_gene_upstream_seqs.txt', 'w') as o:

		for i in non_gene_iso_df.iter_rows(named=True):
			iso_id = i['isoform_id']
			
			if i['strand'] == 'f':
				ends = [int(x) for x in i['associated_starts'].split(',')]
				end = min(ends)
				upstream = genome[end-21:end+19]

				o.write(upstream + '\n')
				gene_seqs.append(upstream)

			elif i['strand'] == 'r':
				ends = [int(x) for x in i['associated_ends'].split(',')]
				end = max(ends)
				upstream = genome[end-20:end+20]
				upstream=Seq(upstream)
				upstream_rc = upstream.reverse_complement()

				o.write(str(upstream_rc )+ '\n')
				non_gene_seqs.append(str(upstream_rc))

if __name__ == '__main__':
	main()