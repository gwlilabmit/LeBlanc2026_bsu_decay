import sys
import polars as pl

# winsor_enum = pl.Enum(['100', '90'])
strand_enum = pl.Enum(['f', 'r'])

def remove_short(df):
    df = df.remove((pl.col('y') - pl.col('x')) <= 200)
    return df

def final_z(df):
    h = df.filter(pl.col('hmean') >= pl.col('vmean'))
    h = h.rename({'hz':'zscore'})
    h = h.drop('vz')

    v = df.filter(pl.col('vmean') > pl.col('hmean'))
    v = v.rename({'vz':'zscore'})
    v = v.drop('hz')

    out = pl.concat([h, v], how='vertical')
    return out

def cleanup(df):
    df = remove_short(df)
    z_clean = final_z(df)
    cleanedup = z_clean.drop(['region', 'hmean', 'vmean', 'hstd', 'vstd'])

    return cleanedup

def z_filter(df, zscore_threshold):
    df = df.filter(pl.col('zscore') >= zscore_threshold)
    return df

def process_infile(infile, strand, column_names, zscore_threshold):
    try:
        df = pl.read_csv(infile, has_header=False, separator='\t').rename(column_names)
    except pl.exceptions.ComputeError: ## addresses rare error in normalized map files
        df = pl.read_csv(infile, has_header=False, separator='\t', infer_schema=False).rename(column_names)
        df = df.remove(pl.col('x') == '1.0')
        sch = {'region':pl.String, 'x':pl.Int64, 'y':pl.Int64, 'count':pl.Int64, 'hmean':pl.Float64, 'vmean':pl.Float64, 'hstd':pl.Float64, 'vstd':pl.Float64, 'hz':pl.Float64, 'vz':pl.Float64}
        df = df.cast(sch)
    cleaned = cleanup(df).with_columns(pl.lit(strand).cast(strand_enum).alias('strand'))
    filtered = z_filter(cleaned, zscore_threshold)

    return filtered

# def remove_duplicates(df_100, df_90):
#     out = df_90.join(df_100, how='anti', on=('x', 'y', 'strand'))
#     return out

# def compare_winsor(df_100, df_90, outfile):

#     differences = remove_duplicates(df_100, df_90)

#     for row in differences.iter_rows(named=True):
#         temp_x = df_100.filter(((pl.col('x') - row['x']) > 5) & ((pl.col('x') - row['x']) < 160))
#         ends_x = temp_x.filter(abs(pl.col('y') - row['y']) <= 5)
#         if ends_x.shape[0] > 0:
#             out = pl.DataFrame(row).cast({'strand':strand_enum, 'winsor':winsor_enum})
#             with open(outfile, 'a') as o:
#                 out.write_csv(o, include_header=False, separator='\t')
#             #df_100 = pl.concat([df_100, out], how='vertical')

#         temp_y = df_100.filter((row['y'] - pl.col('y') > 5) & (row['y'] - pl.col('y') < 160))
#         ends_y = temp_y.filter(abs(pl.col('x') - row['x']) <= 5)
#         if ends_y.shape[0] > 0:
#             out = pl.DataFrame(row).cast({'strand':strand_enum, 'winsor':winsor_enum})
#             with open(outfile, 'a') as o:
#                 out.write_csv(o, include_header=False, separator='\t')
#             #df_100 = pl.concat([df_100, out], how='vertical')

#     return df_100

def main():

    timepoint=sys.argv[1]
    z_thresh=int(sys.argv[2])
    out_prefix='filtered_isoforms/'+timepoint+'_'+str(z_thresh)

    col_names={'column_1':'region', 'column_2':'x', 'column_3':'y', 'column_4':'count', 
               'column_5':'hmean', 'column_6':'vmean', 'column_7':'hstd', 'column_8':'vstd',
               'column_9':'hz', 'column_10':'vz'}

    f_infile_100='normalized_maps_0.0/'+timepoint+'_f_normalized_0.0.txt'
    r_infile_100='normalized_maps_0.0/'+timepoint+'_r_normalized_0.0.txt'

    f_100 = process_infile(f_infile_100, 'f', col_names, z_thresh)
    r_100 = process_infile(r_infile_100, 'r', col_names, z_thresh)

    all_100 = pl.concat([f_100, r_100], how='vertical')
    # all_100 = all_100.with_columns(pl.lit('100').cast(winsor_enum).alias('winsor'))

    all_100.write_csv(out_prefix+'_100wins.txt', separator='\t')

    ####

    f_infile_90='normalized_maps_0.1/'+timepoint+'_f_normalized_0.1.txt'
    r_infile_90='normalized_maps_0.1/'+timepoint+'_r_normalized_0.1.txt'

    f_90 = process_infile(f_infile_90, 'f', col_names, z_thresh)
    r_90 = process_infile(r_infile_90, 'r', col_names, z_thresh)

    all_90 = pl.concat([f_90, r_90], how='vertical')
    all_90.write_csv(out_prefix+'_90wins.txt', separator='\t')
    # all_90 = all_90.with_columns(pl.lit('90').cast(winsor_enum).alias('winsor'))

    ## only keep points from 90% winsorization that are 150 nt upstream (5'end) or 150 nt downstream (3' end)
    ## of a point called in 100% winsorization

    # compare_winsor(all_100, all_90, outfile) ## will append to the already-created output file

if __name__ == '__main__':
    main()