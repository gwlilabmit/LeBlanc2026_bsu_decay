import sys
import polars as pl
import numpy as np
import scipy
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering


def strand_clusters(df, timepoint, winsor, strand):

    strand_calls = df.filter(pl.col('strand') == strand)
    calls_arr = np.array(strand_calls.select(pl.col('x'), pl.col('y')))
    clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=4, linkage='single').fit(calls_arr)
    print('clustering complete')

    with plt.rc_context({'lines.linewidth':0.0005}):
        plot_dendrogram(clustering)
        plt.yscale('symlog')
        plt.axhline(y=4, linewidth=0.5, linestyle = 'dotted', c='k')
        plt.savefig('dendrograms/'+timepoint+'_'+str(winsor)+'_'+strand+'.pdf')
        print('dendrograms/'+timepoint+'_'+str(winsor)+'_'+strand+'.pdf saved')
        plt.close()

    strand_calls = strand_calls.with_columns(pl.Series(name='cluster', values=clustering.labels_))
    
    strand_isoforms = pl.DataFrame(schema={'x':pl.Int64, 'y':pl.Int64, 'count':pl.Int64, 'zscore':pl.Float64, 'strand':pl.String, 'cluster':pl.Int64, 'count_sum':pl.Int64,
                                            'x_min':pl.Int32, 'x_max':pl.Int32, 'y_min':pl.Int32, 'y_max':pl.Int32})
    cluster_ids = strand_calls['cluster'].explode().unique()
    for i in cluster_ids:
        cluster_i = strand_calls.filter(pl.col('cluster') == i)
        read_sum = np.sum(np.array(cluster_i.select(pl.col('count'))))
        if read_sum > 10:
            cluster_i = cluster_i.with_columns(pl.lit(read_sum).cast(pl.Int64).alias('count_sum'))

            x_sort = cluster_i.sort(by=pl.col('x'), descending=True)
            x_min = x_sort.tail(1).select(pl.col('x')).item()
            x_max = x_sort.head(1).select(pl.col('x')).item()

            y_sort = cluster_i.sort(by=pl.col('y'), descending=True)
            y_min = y_sort.tail(1).select(pl.col('y')).item()
            y_max = y_sort.head(1).select(pl.col('y')).item()

            count_sort = cluster_i.sort(by=pl.col('count'), descending=True)
            top_count = count_sort[0].with_columns(pl.lit(x_min).alias('x_min'), pl.lit(x_max).alias('x_max'),
                                           pl.lit(y_min).alias('y_min'), pl.lit(y_max).alias('y_max'))

            strand_isoforms = pl.concat([strand_isoforms, top_count], how='vertical')

    return strand_isoforms

def plot_dendrogram(model):
    ## from scikit-learn documentation example
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    
    scipy.cluster.hierarchy.dendrogram(linkage_matrix, color_threshold = 0)

def main():
    timepoint = sys.argv[1]
    z_thresh = sys.argv[2]
    winsor_level = sys.argv[3]

    infile = 'filtered_isoforms/'+timepoint+'_'+str(z_thresh)+'_'+str(winsor_level)+'wins.txt'
    calls = pl.read_csv(infile, has_header=True, separator='\t')

    f_isos = strand_clusters(calls, timepoint, winsor_level, 'f')
    r_isos = strand_clusters(calls, timepoint, winsor_level, 'r')

    isoforms = pl.concat([f_isos, r_isos], how='vertical')

    outfile = 'isoform_clusters/'+timepoint+'_'+str(winsor_level)+'wins_isoform_clusters.txt'
    isoforms.write_csv(outfile, separator='\t')

if __name__ == '__main__':
    main()