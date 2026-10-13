import os
import sys
import polars as pl
import numpy as np
from scipy.stats.mstats import winsorize

def winsorize_mean(arr, i, winsor_amt, res):
    if len(arr) == False:
        mean = 0.0
    else:
        padded = np.zeros((150,))
        padded[0:len(arr)] = arr ## make all arrays 150 long, with zeros making up difference
        if winsor_amt > 0.0:
            mean = np.mean(winsorize(padded, limits = [0, winsor_amt]))
            res[i] = mean
        elif winsor_amt == 0.0:
            mean = np.mean(padded)
            res[i] = mean
        else:
            print('invalid winsorization threshold')

def winsorize_std(arr, i, winsor_amt, res):
    if len(arr) == False:
        std = 0
    else:
        padded = np.zeros((150,))
        padded[0:len(arr)] = arr ## make all arrays 150 long, with zeros making up difference
        if winsor_amt > 0.0:
            std = np.std(winsorize(padded, limits = [0, winsor_amt]))
            res[i] = std
        elif winsor_amt == 0.0:
            std = np.std(padded)
            res[i] = std
        else:
            print('invalid winsorization threshold')

def hmean(df, x, y, winsor_amt):
    res = np.zeros_like(x, dtype = np.float64)
    for i in range(len(x)):
        x_i = x[i]
        y_i = y[i]
        # lower = x_i - 160
        # upper = x_i - 10

        lower = x_i + 10
        upper = x_i + 160

        background_rows = np.array(df.filter(pl.col('y') == y_i,
                                             pl.col('x').is_between(lower, upper, closed='right')).select('count')[:,0])
        winsorize_mean(background_rows, i, winsor_amt, res)

    return res

def vmean(df, x, y, winsor_amt):
    res = np.zeros_like(y, dtype = np.float64)
    for i in range(len(y)):
        x_i = x[i]
        y_i = y[i]
        lower = y_i - 160
        upper = y_i - 10

        background_rows = np.array(df.filter(pl.col('x') == x_i,
                                             pl.col('y').is_between(lower, upper, closed='left')).select('count')[:,0])
        winsorize_mean(background_rows, i, winsor_amt, res)
    return res

def hstd(df, x, y, winsor_amt):
    res = np.ones_like(x, dtype = np.float64)
    for i in range(len(x)):
        x_i = x[i]
        y_i = y[i]
        lower = x_i + 10
        upper = x_i + 160

        background_rows = np.array(df.filter(pl.col('y') == y_i,
                                             pl.col('x').is_between(lower, upper, closed='right')).select('count')[:,0])
        winsorize_std(background_rows, i, winsor_amt, res)
    return res

def vstd(df, x, y, winsor_amt):
    res = np.ones_like(y, dtype = np.float64)
    for i in range(len(y)):
        x_i = x[i]
        y_i = y[i]
        lower = y_i - 160
        upper = y_i - 10

        background_rows = np.array(df.filter(pl.col('x') == x_i,
                                             pl.col('y').is_between(lower, upper, closed='left')).select('count')[:,0])
        winsorize_std(background_rows, i, winsor_amt, res)
    return res

def zscore(val, mean, std):
    difference = np.subtract(val, mean)
    score = np.divide(difference, std)
    score = np.nan_to_num(score, nan=0, posinf=0, neginf=0)
    return score

# these two are basically the same as above, but only use non-zero values in determining mean and doesn't winsorize
def hmean_nz(df, x, y):
    res = np.zeros_like(x, dtype = np.float64)
    for i in range(len(x)):
        x_i = x[i]
        y_i = y[i]
        # lower = x_i - 160
        # upper = x_i - 10

        lower = x_i + 10
        upper = x_i + 160

        background_rows = np.array(df.filter(pl.col('y') == y_i,
                                    pl.col('x').is_between(lower, upper, closed='right'),
                                    pl.col('count')>0).select('count')[:,0])
        if len(background_rows) == 0:
            res[i] = 0
        else:
            nz_avg = np.mean(background_rows)
            res[i] = nz_avg

    return res

def vmean_nz(df, x, y):
    res = np.zeros_like(y, dtype = np.float64)
    for i in range(len(y)):
        x_i = x[i]
        y_i = y[i]
        lower = y_i - 160
        upper = y_i - 10

        background_rows = np.array(df.filter(pl.col('x') == x_i,
                                             pl.col('y').is_between(lower, upper, closed='left'),
                                             pl.col('count')>0).select('count')[:,0])
        if len(background_rows) == 0:
            res[i] = 0
        else:
            nz_avg = np.mean(background_rows)
            res[i] = nz_avg

    return res

def main():
    
    sample_name = sys.argv[1] ## e.g. t0_f
    section = sys.argv[2] ## e.g. sec0
    winsor = float(sys.argv[3]) # eg. 0.05 or 0.1
    if section == 'None':
        return ''
    
    else:
        infile = 'count_maps/'+sample_name+'.txt'
        
        df = pl.read_csv(infile, has_header=False, separator='\t')
        df = df.rename({'column_1':'region', 'column_2':'x', 'column_3':'y', 'column_4':'count'})
        df = df.filter(pl.col('region') == section)

        ## calculate horizontal and vertical means
        norm = df.with_columns(pl.struct(['x', 'y']).map_batches(
            lambda row: hmean(df, row.struct.field('x'), row.struct.field('y'), winsor
                            ), return_dtype=pl.Float64, returns_scalar=True).alias('hmean'))
        norm = norm.with_columns(pl.struct(['x', 'y']).map_batches(
            lambda row: vmean(df, row.struct.field('x'), row.struct.field('y'), winsor
                            ), return_dtype=pl.Float64, returns_scalar=True).alias('vmean'))
        
        ## calculate horizontal and vertical standard deviations
        norm = norm.with_columns(pl.struct(['x', 'y']).map_batches(
            lambda row: hstd(df, row.struct.field('x'), row.struct.field('y'), winsor
                            ), return_dtype=pl.Float64, returns_scalar=True).alias('hstd'))                                 
        norm = norm.with_columns(pl.struct(['x', 'y']).map_batches(
            lambda row: vstd(df, row.struct.field('x'), row.struct.field('y'), winsor
                            ), return_dtype=pl.Float64, returns_scalar=True).alias('vstd'))
        
        ## calculate horizontal and vertical z-scores
        norm = norm.with_columns(pl.struct(['count', 'hmean', 'hstd']).map_batches(
            lambda row: zscore(row.struct.field('count'), row.struct.field('hmean'), row.struct.field('hstd')
                            ), return_dtype=pl.Float64, returns_scalar=True).alias('hz'))
        
        norm = norm.with_columns(pl.struct(['count', 'vmean', 'vstd']).map_batches(
            lambda row: zscore(row.struct.field('count'), row.struct.field('vmean'), row.struct.field('vstd')
                            ), return_dtype=pl.Float64, returns_scalar=True).alias('vz'))
        
        # ## calculate non-zero horizontal and vertical means for downstream filtering
        # norm = norm.with_columns(pl.struct(['x', 'y']).map_batches(
        #     lambda row: hmean_nz(df, row.struct.field('x'), row.struct.field('y')
        #                          ), return_dtype=pl.Float64, returns_scalar=True).alias('hmean_nz'))
        
        # norm = norm.with_columns(pl.struct(['x', 'y']).map_batches(
        #     lambda row: vmean_nz(df, row.struct.field('x'), row.struct.field('y')
        #                          ), return_dtype=pl.Float64, returns_scalar=True).alias('vmean_nz'))


        ## done in next step instead
        # h = norm.filter(pl.col('hmean') >= pl.col('vmean'))
        # h = h.rename({'hz':'zscore'})
        # h = h.drop('vz')

        # v = norm.filter(pl.col('vmean') > pl.col('hmean'))
        # v = v.rename({'vz':'zscore'})
        # v = v.drop('hz')

        # out = pl.concat([h, v], how='vertical')

        outfile = 'normalized_maps_'+str(winsor)+'/'+sample_name+'_normalized_'+str(winsor)+'.txt'

        if os.path.exists(outfile):
            with open(outfile, 'a') as f:
                norm.write_csv(f, separator='\t', include_header = False)

        else:
            norm.write_csv(outfile, separator='\t', include_header = False)

if __name__ == '__main__':
    main()
