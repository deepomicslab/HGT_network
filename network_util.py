import numpy as np
import pandas as pd
from copy import deepcopy as dc
import os

# normalize by log10 and rescale to [0,1]
def log_min_max(df, log_max, log_min):
    range_interval = log_max - log_min
    df = np.array(list(map(np.log10, df)))
    df -= log_min
    df /= range_interval
    df[df == -np.inf] = 0
    return df

def simple_min_max(df):
    return (df - df.min().min())/(df.max().max() - df.min().min())


def log_rescale(df):
    np.fill_diagonal(df.values, 0)
    max_num = df.max().max()
    min_num = sorted(list(set(df.values.flatten())))[1]
    log_max = np.log10(max_num)
    log_min = np.log10(min_num)
    matrix = log_min_max(df.values, log_max, log_min)
    df = pd.DataFrame(matrix, columns=df.columns, index=df.index)
    np.fill_diagonal(df.values, 0)
    return df

def make_norm_net(data, taxonomy_df, level, odir):
    taxonomy_df = dc(taxonomy_df[~taxonomy_df[level].isna()])
    sp_list = list(set(taxonomy_df[level]))
    graph = pd.DataFrame(columns=sp_list, index=sp_list)
    graph.fillna(0, inplace=True)
    for idx in data.index:
        tid1 = data.loc[idx, 0]
        tid2 = data.loc[idx, 2]
        if not (tid1 in taxonomy_df.index and tid2 in taxonomy_df.index):
            continue
        sp1 = taxonomy_df.loc[tid1, level]
        sp2 = taxonomy_df.loc[tid2, level]
        graph.loc[sp1, sp2] += data.loc[idx, 4]
    if '{}__'.format(level) in graph.index:
        graph.drop('{}__'.format(level), axis=0, inplace=True)
        graph.drop('{}__'.format(level), axis=1, inplace=True)
    graph.to_csv(os.path.join(odir, 'raw_directed_graph_dig.tsv'), sep='\t')

    np.fill_diagonal(graph.values, 0)
    undirected = graph + graph.T
    undirected = undirected.loc[:, (undirected != 0).any(axis=0)]
    undirected = undirected.loc[(undirected != 0).any(axis=1), :]

    # normalized by degree
    col_sum = np.sqrt(undirected.sum())
    undirected = undirected.div(col_sum, axis=1)
    undirected = undirected.div(col_sum, axis=0)
    undirected.to_csv(os.path.join(odir, 'degree_norm.tsv'), sep='\t') 
