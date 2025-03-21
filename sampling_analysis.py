import os 
import sys
import numpy as np
import pandas as pd
import networkx as nx 
from copy import deepcopy as dc
import matplotlib.pyplot as plt
import network_analysis as nta
import random
from scipy.stats import mannwhitneyu
import tree_util
import copy
import network_util as nutil

def decribe_tree(newick_tree):
    json_tree = tree_util.parse(newick_tree)
    largest = {'largest': 0}
    leaf_list, l = tree_util.recu_compute(json_tree, 0, largest)
    largest_level = largest['largest']
    nlayer = largest_level
    leaf_list, l = tree_util.recu_compute(json_tree, 0, largest)
    layer_leaves_dict = tree_util.make_layer_dict(nlayer)

    tree_util.recu_layer(json_tree, layer_leaves_dict)
    tree_util.to_layer_leaves(layer_leaves_dict, nlayer)
    name_dict = {}
    result = {}
    # compute leaf layer
    result['leaves_dict'] = copy.deepcopy(layer_leaves_dict)
    parent_dict = {}
    tree_util.parents(json_tree, parent_dict)
    node_leaves = {}
    for level in layer_leaves_dict.keys():
        for node, sp_list in layer_leaves_dict[level].items():
            if node in node_leaves.keys():
                continue
            node_leaves[node] = copy.deepcopy(sp_list)
    subtree_nodes = {}
    for l in leaf_list:
        parent = parent_dict[l]
        if parent not in subtree_nodes.keys():
            subtree_nodes[parent] = []
        subtree_nodes[parent].append(l)

    for node in node_leaves.keys():
        parent = parent_dict[node]
        if parent not in subtree_nodes.keys():
            subtree_nodes[parent] = []
        subtree_nodes[parent] += subtree_nodes[node]
        subtree_nodes[parent].append(node)

    for node in subtree_nodes.keys():
        subtree_nodes[node].append(node)

    direct_children_dict = {}
    for node, parent in parent_dict.items():
        if parent not in direct_children_dict:
            direct_children_dict[parent] = []
        direct_children_dict[parent].append(node)

    result['leaves_dict']['root'] = {'root': list(name_dict.values())}
    node_leaves['root'] = leaf_list
    return node_leaves, parent_dict, subtree_nodes, direct_children_dict

def report(norm_df, n_sample, newick_tree=None):
    if newick_tree:
        node_leaves, parent_dict, subtree_nodes, direct_children_dict = decribe_tree(newick_tree)
    sampled_list = random.sample(list(norm_df.index), min(norm_df.shape[0], n_sample))
    norm_df = norm_df.loc[sampled_list, sampled_list]
    graph = norm_df
    # some features
    # number of genera 
    report_dict = {}
    report_dict['edge_number'] = (norm_df > 0).sum().sum()/2
    report_dict['size'] = graph.shape[0]
    pe = nta.pe(graph)
    report_dict['pe'] = pe
    pr = nta.pr(graph)
    report_dict['largest_pr'] = pr
    btw = nta.betweenness(graph)
    #btw.to_csv(os.path.join(odir, 'betweeness.tsv'), sep='\t', header=True, index=True)
    density = nta.density(graph)
    report_dict['density'] = density
    cc = nta.clustering_coef(graph)
    report_dict['node_connect_coef'] = cc[0]
    report_dict['transitivity'] = nta.transitivity(graph)
    report_dict['assortativity'] = nta.assortativity(graph)
    report_dict['alge'] = nta.alge_connectivity(graph)
    #cc[1].to_csv(os.path.join(odir, 'node_connect_coef.tsv'), sep='\t', header=True, index=True)
    report_dict['norm_df'] = dc(norm_df)
    if newick_tree:
        rename_dict = {}
        for name in norm_df.index:
            rename_dict[name] = name.replace(' ', '-')
            rename_dict[name] = name.replace('_', '-')
        norm_df.rename(index=rename_dict, columns=rename_dict, inplace=True)
        ## log rescalse
        norm_df = nutil.log_rescale(norm_df)
        common_sp = list(set(norm_df.index).intersection(set(node_leaves['root'])))
        norm_df = norm_df.loc[common_sp, common_sp]
    return report_dict                                                                              

def check_min_shape(norm_df_dict):
    min_shape = 1000000
    for key in norm_df_dict:
        if norm_df_dict[key].shape[0] < min_shape:
            min_shape = norm_df_dict[key].shape[0]
    return min_shape

def report_dict(norm_df_dict, n_sample=None, newick_tree=None):
    result_dict = {}
    if not n_sample:
        n_sample = check_min_shape(norm_df_dict)
    for key in norm_df_dict:
        #odir = os.path.join(outer_dir, key)
        print(key)
        result_dict[key] = report(norm_df_dict[key], n_sample, newick_tree)
        
    return result_dict

def repeat_exp(repeat_time, norm_df_dict, n_sample=None, newick_tree=None):
    result_dict = {}
    for i in range(repeat_time):
        result_dict[i] = dc(report_dict(norm_df_dict, n_sample, newick_tree))
    return result_dict

def feature_vector(repeat_dict, dtype, feature):
    result = []
    for i in repeat_dict:
        result.append(repeat_dict[i][dtype][feature])
    return result

def test_feature(repeat_dict, feature, dtype_list):
    result_diff_df = pd.DataFrame(index=dtype_list, columns=dtype_list)
    result_p = pd.DataFrame(index=dtype_list, columns=dtype_list)
    for dtype1 in dtype_list:
        v1 = feature_vector(repeat_dict, dtype1, feature)
        for dtype2 in dtype_list:
            v2 = feature_vector(repeat_dict, dtype2, feature)
            t, p = mannwhitneyu(v1, v2)
            print(p)
            result_p.loc[dtype1, dtype2] = p
            diff = np.mean(v1) - np.mean(v2)
            if p < 0.05:
                result_diff_df.loc[dtype1, dtype2] = diff
    return result_diff_df, result_p
    
def fr(d_df, profile, sname, is_d):
    sp_list = list(set(profile.columns).intersection(set(d_df.index)))

    sp_d_df = d_df.loc[sp_list, sp_list]

    sp_profile = np.array(profile.loc[sname, sp_list])
    value = np.dot(sp_profile.reshape(len(sp_profile), 1),sp_profile.reshape(1, len(sp_profile)))
    width = value.shape[0]
    if is_d:
        cor_df = np.ones(shape=(width, width)) - sp_d_df.values
    else:
        cor_df = sp_d_df.values
    for i in range(width):
        cor_df[i][i] = 0
    value = np.multiply(value, cor_df)
    fr_df = pd.DataFrame(value, index=sp_list, columns=sp_list)
    return fr_df
