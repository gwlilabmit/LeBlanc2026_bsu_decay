import polars as pl
import os

def main():

	## same 5' end different 3' end

	ends_5 = pl.read_csv('../canonical_end_comparison/primary_5_matches.txt', separator='\t', truncate_ragged_lines=True)
	ratio_dic = {}

	used_list = []

	for e in ends_5.iter_rows():
		end=e[0]
		strand=e[1]
		if e[2]==None:
			continue

		end_file = 'end_matches_5/'+strand+'/'+str(end)+'.txt'
		if os.path.exists(end_file):
			df = pl.read_csv(end_file, separator='\t')

			## remove instances of ends that are very close to each other. mainly seems to be the case for some minor 3' ends
			if df.select(pl.len()).item() == 2:
				if strand == 'f':
					isos = df.sort(by='associated_ends', descending=False) ## this works to sort from higher start coordinate (shorter iso) to lower (longer iso), even when there are multiple associated ends
					end_list = []
					for i in isos.iter_rows(named=True):
						if type(i['associated_ends'])==str:
							ends_3 = list(map(int, i['associated_ends'].split(',')))
							end_list.append(ends_3)
						else:
							ends_3 = i['associated_ends']
							end_list.append(ends_3)

					iso_1_ends=end_list[0]
					iso_2_ends=end_list[1]
					if (type(iso_1_ends)==int) and (type(iso_2_ends)==int):
						distance = abs(iso_1_ends-iso_2_ends)
						if distance < 10:
							continue ## don't do analysis on 3' ends that are very close to each other

				elif strand == 'r':
					isos = df.sort(by='associated_starts', descending=True)
					end_list = []
					for i in isos.iter_rows(named=True):
						if type(i['associated_starts'])==str:
							ends_3 = list(map(int, i['associated_starts'].split(',')))
							end_list.append(ends_3)
						else:
							ends_3 = i['associated_starts']
							end_list.append(ends_3)

					iso_1_ends=end_list[0]
					iso_2_ends=end_list[1]
					if (type(iso_1_ends)==int) and (type(iso_2_ends)==int):
						distance = abs(iso_1_ends-iso_2_ends)
						if distance < 10:
							continue ## don't do analysis on 3' ends that are very close to each other

				shorter_k = isos[0].select(pl.col('k')).item() ## shorter
				longer_k = isos[1].select(pl.col('k')).item() ## longer

				ratio = (1/(shorter_k))/(1/(longer_k))

				shorter_id = isos[0].select(pl.col('isoform_id')).item()
				longer_id = isos[1].select(pl.col('isoform_id')).item()

				## if this pair of ends has already been calculated, skip it
				## can happen if annotation file has two ends that are very close together
				if (shorter_id, longer_id) in used_list:
					continue

				used_list.append((shorter_id, longer_id))

				info = [strand, shorter_id, longer_id, ratio]
				ratio_dic[end] = info

	outpath = 'primary_5_ratios.txt'
	with open(outpath, 'w') as o:
		o.write('end\tstrand\tiso_id1\tiso_id2\tratio\n')
		for i in ratio_dic.keys():
			o.write(str(i) + '\t')
			info = list(map(str, ratio_dic.get(i)))
			o.write('\t'.join(info))
			o.write('\n')

	## same 3' end different 5' end

	ends_3 = pl.read_csv('../canonical_end_comparison/primary_3_matches.txt', separator='\t', truncate_ragged_lines=True)
	ratio_dic = {}

	used_list = []

	for e in ends_3.iter_rows():
		end=e[0]
		strand=e[1]
		if e[2]==None:
			continue

		end_file = 'end_matches_3/'+strand+'/'+str(end)+'.txt'
		if os.path.exists(end_file):
			df = pl.read_csv(end_file, separator='\t')

			if df.select(pl.len()).item() == 2:
				if strand == 'f':
					isos = df.sort(by='associated_starts', descending=False)
					end_list = []
					for i in isos.iter_rows(named=True):
						if type(i['associated_starts'])==str:
							ends_3 = list(map(int, i['associated_starts'].split(',')))
							end_list.append(ends_3)
						else:
							ends_3 = i['associated_starts']
							end_list.append(ends_3)

					iso_1_ends=end_list[0]
					iso_2_ends=end_list[1]
					if (type(iso_1_ends)==int) and (type(iso_2_ends)==int):
						distance = abs(iso_1_ends-iso_2_ends)
						if distance < 10:
							continue ## don't do analysis on 3' ends that are very close to each other

				elif strand == 'r':
					isos = df.sort(by='associated_ends', descending=True)
					end_list = []
					for i in isos.iter_rows(named=True):
						if type(i['associated_ends'])==str:
							ends_3 = list(map(int, i['associated_ends'].split(',')))
							end_list.append(ends_3)
						else:
							ends_3 = i['associated_ends']
							end_list.append(ends_3)

					iso_1_ends=end_list[0]
					iso_2_ends=end_list[1]
					if (type(iso_1_ends)==int) and (type(iso_2_ends)==int):
						distance = abs(iso_1_ends-iso_2_ends)
						if distance < 10:
							continue ## don't do analysis on 3' ends that are very close to each other

				shorter_k = isos[0].select(pl.col('k')).item() ## shorter
				longer_k = isos[1].select(pl.col('k')).item() ## longer

				ratio = (1/(shorter_k))/(1/(longer_k))

				shorter_id = isos[0].select(pl.col('isoform_id')).item()
				longer_id = isos[1].select(pl.col('isoform_id')).item()

				if (shorter_id, longer_id) in used_list:
					continue

				used_list.append((shorter_id, longer_id))

				info = [strand, shorter_id, longer_id, ratio]
				ratio_dic[end] = info

	outpath = 'primary_3_ratios.txt'
	with open(outpath, 'w') as o:
		o.write('end\tstrand\tiso_id1\tiso_id2\tratio\n')
		for i in ratio_dic.keys():
			o.write(str(i) + '\t')
			info = list(map(str, ratio_dic.get(i)))
			o.write('\t'.join(info))
			o.write('\n')

if __name__ == '__main__':
	main()
