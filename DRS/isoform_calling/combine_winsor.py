import sys
import polars as pl

winsor_enum = pl.Enum(['100', '90'])
strand_enum = pl.Enum(['f', 'r'])

def remove_duplicates(df_100, df_90):
    out = df_90.join(df_100, how='anti', on=('x', 'y', 'strand'))
    return out

def compare_winsor(df_100, df_90, outfile):

    differences = remove_duplicates(df_100, df_90)

    for row in differences.iter_rows(named=True):
        ## if the isoform is not exactly the same as one already called, but it's very similar, it doesn't get merged
        ## may want to have different cutoffs for 5' and 3' ends, but will set those later, so come back and adjust
        ## expect that the cutoff will be >=2 for both kinds of ends, so use that for the disqualifier threshold
        disqualifiers = df_100.filter((abs(pl.col('x') - row['x']) <= 2) & (abs(pl.col('y') - row['y']) <= 2))
        if disqualifiers.shape[0] > 0:
            continue ## if there are no-winsor isoforms sufficiently similar to the 90% isoform in question, move on to the next

        ## merged isoforms must meet the criteria and also be sufficiently different from each other!
        temp_x = df_100.filter(((pl.col('x') - row['x']) > 5) & ((pl.col('x') - row['x']) < 175))
        ends_x = temp_x.filter(abs(pl.col('y') - row['y']) <= 5)
        if ends_x.shape[0] > 0:
            out = pl.DataFrame(row).cast({'strand':strand_enum, 'winsor':winsor_enum})
            with open(outfile, 'a') as o:
                out.write_csv(o, include_header=False, separator='\t')
            #df_100 = pl.concat([df_100, out], how='vertical')

        temp_y = df_100.filter((row['y'] - pl.col('y') > 5) & (row['y'] - pl.col('y') < 175))
        ends_y = temp_y.filter(abs(pl.col('x') - row['x']) <= 5)
        if ends_y.shape[0] > 0:
            out = pl.DataFrame(row).cast({'strand':strand_enum, 'winsor':winsor_enum})
            with open(outfile, 'a') as o:
                out.write_csv(o, include_header=False, separator='\t')
            #df_100 = pl.concat([df_100, out], how='vertical')

    #return df_100

def main():

    timepoint=sys.argv[1]
    out_prefix = 'combined_isoforms/'+timepoint

    infile_100 = 'isoform_clusters/'+timepoint+'_100wins_isoform_clusters.txt'
    infile_90 = 'isoform_clusters/'+timepoint+'_90wins_isoform_clusters.txt'

    in_100 = pl.read_csv(infile_100, separator='\t').cast({'strand':strand_enum}).with_columns(pl.lit('100').cast(winsor_enum).alias('winsor'))
    in_90 = pl.read_csv(infile_90, separator='\t').cast({'strand':strand_enum}).with_columns(pl.lit('90').cast(winsor_enum).alias('winsor'))

    print('files read')

    in_100.write_csv(out_prefix+'_combined.txt', separator='\t')
    compare_winsor(in_100, in_90, out_prefix+'_combined.txt')

if __name__ == '__main__':
    main()