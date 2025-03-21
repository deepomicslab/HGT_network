import os 
import numpy as np
import pandas as pd
import networkx as nx 
from pyseat.SEAT import SEAT
from copy import deepcopy as dc
import matplotlib.pyplot as plt
from scipy import optimize
import powerlaw
from scipy.stats import pearsonr

# entropy 
def pe(G):
    degree = G.sum()
    m2 = G.sum().sum()
    value = 0
    for sp in degree.index:
        e = degree[sp]/m2
        if e == 0:
            continue
        value -= e * np.log2(e)
    return value

def pr(G):
    G = nx.from_pandas_adjacency(G)
    pagerank = nx.pagerank(G)
    print(pagerank.values())
    max_value = max(pagerank.values())
    max_index = list(pagerank.values()).index(max_value)
    sp = list(pagerank.keys())[max_index]
    pr_df = pd.DataFrame(index=pagerank.keys(), columns=['pr'])
    for node, pr_v in pagerank.items():
        pr_df.loc[node, 'pr']= pr_v
    pr_df.sort_values(by='pr', inplace=True, ascending=False)
    return pr_df, sp 

def community_detection(G):## G is pandas dataframe
    #community_detection(graph)
    sp_list = list(G.index)
    seat = SEAT(affinity="precomputed",
                #sparsification="precomputed",
                objective="SE",
                max_k=20,
                strategy="bottom_up")
    seat.fit_predict(G)
    labels = seat.labels_
    communities = {}
    for i, l in enumerate(labels):
        if l not in communities:
            communities[l] = []
        sp = sp_list[i]
        communities[l].append(sp)
    return communities

def connected_component(G):
    G = nx.from_pandas_adjacency(G)
    connected = nx.is_connected(G)
    ncomponent = nx.number_connected_components(G)
    components = nx.connected_components(G)
    component_dict = {}
    for i, component in enumerate(components):
        component_dict[i] = []
        for j in component:
            component_dict[i].append(j)
    return connected, ncomponent, component_dict

def louvain_community(G):
    communities = {}
    G = nx.from_pandas_adjacency(G)
    result = nx.community.louvain_communities(G)
    for i, node_set in enumerate(result):
        communities[i] = []
        for n in node_set:
            communities[i].append(n)
    return communities

def node_connectivity(G):
    G = nx.from_pandas_adjacency(G)
    connectivity = nx.node_connectivity(G)
    return connectivity

def betweenness(G):
    G = nx.from_pandas_adjacency(G)
    betweenness = nx.betweenness_centrality(G)
    result = pd.DataFrame()
    for node, centrality in betweenness.items():
        result.loc[node, 'betweenness']= centrality
    result.sort_values(by='betweenness', inplace=True, ascending=False)
    return result

def clustering_coef(G):
    G = nx.from_pandas_adjacency(G)
    avg_clustering = nx.average_clustering(G)
    clustering = nx.clustering(G)
    result = pd.DataFrame()
    for node, coefficient in clustering.items():
        result.loc[node, 'coef']= coefficient
    result.sort_values(by='coef', inplace=True, ascending=False)
    return avg_clustering, result

def distribution(data, title):
    plt.hist(data, bins=30, density=True, alpha=0.7, color='skyblue')
    plt.title(title)
    plt.xlabel("Value")
    plt.ylabel("Density")
    plt.show()

def scale_free(graph):
    G = nx.from_pandas_adjacency(graph)
    degrees = dict(G.degree())
    weights = nx.get_edge_attributes(G, 'weight')
    #correlation, p = pearsonr(degrees.values(), weights.values())
    fit_degrees = powerlaw.Fit(list(degrees.values()))
    fit_weights = powerlaw.Fit(list(weights.values()))
    p1 = fit_degrees.power_law.alpha
    p2 = fit_weights.power_law.alpha
    fig2 = fit_degrees.plot_pdf(color='b', linewidth=2)
    fit_degrees.power_law.plot_pdf(color='g', linestyle='--', ax=fig2)
    plt.show()
    is_scale_free = p1 > 2.0 and p2 > 2.0
    R, p = fit_degrees.distribution_compare('power_law', 'exponential', normalized_ratio=True)
    return is_scale_free

def transitivity(G):
    G = nx.from_pandas_adjacency(G)
    transitivity = nx.transitivity(G)
    return transitivity


def density(G):
    G = nx.from_pandas_adjacency(G)
    density = nx.density(G)
    return density

def inter_intra(level, taxonomy_df, graph, normalized=False):
    classify_dict = {}
    for sp in graph.index:
        c = taxonomy_df.loc[sp, level]
        if c not in classify_dict.keys():
            classify_dict[c] = []
        classify_dict[c].append(sp)
    
    clist = list(classify_dict.keys())
    intra_df = pd.DataFrame(index=clist, columns=['edge_n', 'nor_edge_n', 'sp_n'])
    for c, select_list in classify_dict.items():
        intra_part = graph.loc[select_list, select_list]
        np.fill_diagonal(intra_part.values, 0)
        intra_df.loc[c, 'edge_n'] = intra_part.sum().sum()/2
        intra_df.loc[c, 'sp_n'] = len(select_list)
        intra_df.loc[c, 'nor_edge_n'] = intra_df.loc[c, 'edge_n']/(intra_df.loc[c, 'sp_n']**2)
    
    inter_df = pd.DataFrame(index=clist, columns=clist)
    inter_df.fillna(0, inplace=True)
    for i in range(len(clist)):
        c1 = clist[i]
        selected_list1 = classify_dict[c1]
        n1 = intra_df.loc[c1, 'sp_n']
        if not normalized:             
            inter_df.loc[c1, c1] = intra_df.loc[c, 'edge_n']
        else:
            inter_df.loc[c1, c1] = intra_df.loc[c, 'nor_edge_n']
        for j in range(i+1, len(clist)):
            c2 = clist[j]
            selected_list2 = classify_dict[c2]
            n2 = intra_df.loc[c2, 'sp_n']
            inter_part = graph.loc[selected_list1, selected_list2]
            if not normalized:
                inter_df.loc[c1, c2] = inter_part.sum().sum()
            else:
                inter_df.loc[c1, c2] = inter_part.sum().sum()/(n1*n2)
            inter_df.loc[c2, c1] = inter_df.loc[c1, c2]

    return intra_df, inter_df


def assortativity(G):
    G = nx.from_pandas_adjacency(G)
    assortativity = nx.degree_assortativity_coefficient(G)
    return assortativity


def alge_connectivity(G):
    G = nx.from_pandas_adjacency(G)
    connectivity = nx.algebraic_connectivity(G)
    return connectivity
