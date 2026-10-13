import polars as pl
from Bio.Seq import Seq
import regex

## original version set of increasers is based on actual RNA increase amount
## this one is based on having k>15

## this modified version filters first by k>15, then by increases
## uses same set of secondary isoforms used to plot the CDF in figure 4 (this is the change from the last version)

## modification from v3: removes k CV > 0.3 – probably not as important here
## but will ensure that the set i use here is the same as the set in the CDF in the last figure

def main():

	info = pl.read_csv('../canonical_end_comparison/filtered_secondary_isoforms.txt', separator='\t')

	id_list = []

	for i in info.iter_rows(named=True):
		k = i['k']
		CV = i['k_coef_var']
		iso_id = i['isoform_id']

		if k == 0:
			continue

		if k>=15 and CV<0.3:
			id_list.append(iso_id)

	stable = info.filter(pl.col('isoform_id').is_in(id_list))

	stable.write_csv('stable_fits.txt', separator='\t')

	### now from the very stable isoforms, specifically identify the increasers
	### not as unexpected that some things would just be pretty stable
	inc_id_list = []
	for i in stable.iter_rows(named=True):
		amounts = [i['t0'], i['t0.5'], i['t2'], i['t4'], i['t8'], i['t12'], i['t16']]
		alpha = i['t0']
		iso_id = i['isoform_id']

		if alpha == 0:
			continue

		amounts_rel = [(x/alpha) for x in amounts]

		if ((amounts_rel[-2] > 2) and (amounts_rel[-3] > 2)):
			inc_id_list.append(iso_id)

	increaser_fits = info.filter(pl.col('isoform_id').is_in(inc_id_list))

	increaser_fits.write_csv('increaser_fits.txt', separator='\t')

	with open('B_subtilis_168_NC_000964.3.fa', 'r') as g:
		g.readline()
		genome = g.readline()

	upstream_seqs=[]
	with open('upstream_seqs.txt', 'w') as o:
		for i in increaser_fits.iter_rows(named=True):
			iso_id = i['isoform_id']
			
			if i['strand'] == 'f':
				ends = [int(x) for x in i['associated_starts'].split(',')]
				end = min(ends)
				upstream = genome[end-21:end+19]

				o.write(upstream + '\n')
				upstream_seqs.append(upstream)

			elif i['strand'] == 'r':
				ends = [int(x) for x in i['associated_ends'].split(',')]
				end = max(ends)
				upstream = genome[end-20:end+20]
				upstream=Seq(upstream)
				upstream_rc = upstream.reverse_complement()

				o.write(str(upstream_rc )+ '\n')
				upstream_seqs.append(str(upstream_rc))

	primary_isos=pl.read_csv('../canonical_end_comparison/filtered_primary_isoforms.txt', separator='\t')
	all_isos=pl.concat([primary_isos, info], how='vertical')

	all_seqs = []
	with open('all_seqs.txt', 'w') as o:
		for i in all_isos.iter_rows(named=True):
			iso_id = i['isoform_id']
			CV = i['k_coef_var']
			
			if i['strand'] == 'f' and CV < 0.3:
				ends = [int(x) for x in i['associated_starts'].split(',')]
				end = min(ends)
				upstream = genome[end-21:end+19]

				o.write(upstream + '\n')
				all_seqs.append(upstream)

			elif i['strand'] == 'r' and CV < 0.3:
				ends = [int(x) for x in i['associated_ends'].split(',')]
				end = max(ends)
				upstream = genome[end-20:end+20]
				upstream=Seq(upstream)
				upstream_rc = upstream.reverse_complement()

				o.write(str(upstream_rc)+ '\n')
				all_seqs.append(str(upstream_rc))

	perfect_matches=0
	one_mismatch=0
	# upstream_only = [x[:21] for x in upstream_seqs]
	upstream_only = [x[:16] for x in upstream_seqs]
	for i in upstream_only:
		if 'AGGAGG' in i:
			perfect_matches+=1
		temp = regex.findall("(AGGAGG){e<=2}", i)
		if len(temp) >=1:
			one_mismatch+=1
		

	print(perfect_matches)
	print(perfect_matches/len(upstream_seqs))

	print(one_mismatch)
	print(one_mismatch/len(upstream_seqs))

	perfect_matches=0
	one_mismatch=0
	#upstream_only = [x[:21] for x in all_seqs]
	upstream_only = [x[:16] for x in all_seqs]
	for i in upstream_only:
		if 'AGGAGG' in i:
			perfect_matches+=1
		temp = regex.findall("(AGGAGG){e<=2}", i)
		if len(temp) >=1:
			one_mismatch+=1	

	print(perfect_matches)
	print(perfect_matches/len(all_seqs))

	print(one_mismatch)
	print(one_mismatch/len(all_seqs))


	### check 3' end seqs too

	downstream_seqs=[]
	with open('downstream_seqs.txt', 'w') as o:
		for i in increaser_fits.iter_rows(named=True):
			iso_id = i['isoform_id']
			
			if i['strand'] == 'f':
				ends = [int(x) for x in i['associated_ends'].split(',')]
				end = min(ends)
				upstream = genome[end-21:end+19]

				o.write(upstream + '\n')
				upstream_seqs.append(upstream)

			elif i['strand'] == 'r':
				ends = [int(x) for x in i['associated_starts'].split(',')]
				end = max(ends)
				upstream = genome[end-20:end+20]
				upstream=Seq(upstream)
				upstream_rc = upstream.reverse_complement()

				o.write(str(upstream_rc )+ '\n')
				upstream_seqs.append(str(upstream_rc))

	primary_isos=pl.read_csv('../canonical_end_comparison/filtered_primary_isoforms.txt', separator='\t')
	all_isos=pl.concat([primary_isos, info], how='vertical')

	all_down_seqs = []
	with open('all_3_seqs.txt', 'w') as o:
		for i in all_isos.iter_rows(named=True):
			iso_id = i['isoform_id']
			CV = i['k_coef_var']
			
			if i['strand'] == 'f' and CV < 0.3:
				ends = [int(x) for x in i['associated_ends'].split(',')]
				end = min(ends)
				upstream = genome[end-21:end+19]

				o.write(upstream + '\n')
				all_seqs.append(upstream)

			elif i['strand'] == 'r' and CV < 0.3:
				ends = [int(x) for x in i['associated_starts'].split(',')]
				end = max(ends)
				upstream = genome[end-20:end+20]
				upstream=Seq(upstream)
				upstream_rc = upstream.reverse_complement()

				o.write(str(upstream_rc)+ '\n')
				all_seqs.append(str(upstream_rc))


if __name__ == '__main__':
	main()