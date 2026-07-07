#!/usr/bin/env python
# coding: utf-8

import numpy as np
import torch
from torch_geometric.utils import to_undirected, remove_isolated_nodes
from sklearn.preprocessing import MinMaxScaler
from torch_geometric.data import Data
import pickle

from torch_geometric.data import Data
from torch_geometric.utils import to_undirected, subgraph


def to_undirected_mean(edge_index, edge_weight):

    edge_weight = edge_weight[edge_index[0] != edge_index[1]]
    edge_index = edge_index[:,edge_index[0] != edge_index[1]]

    edge_index = torch.tensor(edge_index).long()
    edge_weight = torch.tensor(edge_weight).float()

    edge_index, edge_weight = to_undirected(edge_index, edge_weight, reduce = "mean")
    
    return edge_index, edge_weight


def remove_large_index(edge_index, edge_weight, max_index):
    large_index = edge_index > max_index
    keep_index = torch.logical_not(torch.logical_or(large_index[0], large_index[1]))
    
    edge_index = edge_index[:,keep_index]
    edge_weight = edge_weight[keep_index]
    
    return edge_index, edge_weight
    

def data_process(config, chr_name):
    
    hic_dir_path = config["hic_dir_path"]
    feature_dir_path = config["feature_dir_path"]

    edge_file_name = config["edge_file_name"].replace("chr0", chr_name)
    weight_file_name = config["weight_file_name"].replace("chr0", chr_name)
    feature_file_name = config["feature_file_name"].replace("chr0", chr_name)
    
    
    x = np.load(feature_dir_path+feature_file_name, allow_pickle=True).astype(float).T
    x = np.nan_to_num(x)

    scalerx = MinMaxScaler()
    x = scalerx.fit_transform(x)
    x = torch.tensor(x).float()

    edge_hic = np.load(hic_dir_path+edge_file_name)
    weight_hic = np.load(hic_dir_path+weight_file_name)
    edge_hic, weight_hic = to_undirected_mean(edge_hic, weight_hic)
    
    edge_hic, weight_hic = remove_large_index(edge_hic, weight_hic, x.shape[0]-1)
            
    
    return x, edge_hic, weight_hic
    

def sliding_data_func(graph, size):

    sliding_data = []
    sliding_index = []

    ptr = torch.nonzero(graph.x)[-1][0].item()
    
    start = torch.nonzero(graph.x)[0][0].item()
    
    for s in range(start,ptr,int(size)):
        
        if s+size < ptr:
            sub_nodes = torch.tensor(np.arange(s,s+size))
        else:
            sub_nodes = torch.tensor(np.arange(ptr-size,ptr))

        try:
            sub_edge_index, sub_edge_attr = subgraph(sub_nodes, graph.edge_index, graph.edge_attr, relabel_nodes=True)
        except:
            continue
        
        sub_x = graph.x[sub_nodes]
            
        
        sliding_data.append(Data(x=sub_x, edge_index=sub_edge_index, edge_attr=sub_edge_attr))
        sliding_index.append(sub_nodes.tolist())
    
    return sliding_data, sliding_index


def data_process_all_chrom(config, size):

    data_list = []
    index_list = []
    
    for chr_name in config["chr_names"]:
    
        x, edge_hic, weight_hic = data_process(config, chr_name)

        graph = Data(x=x, edge_index=edge_hic, edge_attr=weight_hic)
        
        sliding_data, sliding_index = sliding_data_func(graph, size)
        
        data_list.extend(sliding_data)
        index_list.extend(sliding_index)
    
    return data_list, index_list
    
