import pickle
import numpy as np
import torch
from torch_geometric.utils import to_undirected, subgraph
from torch_geometric.data import Data
from sklearn.preprocessing import MinMaxScaler
from make_splits_val import neg_generator

def to_undirected_mean(edge_index, edge_weight):
    edge_weight = edge_weight[edge_index[0] != edge_index[1]]
    edge_index = edge_index[:, edge_index[0] != edge_index[1]]

    edge_index = torch.tensor(edge_index).long()
    edge_weight = torch.tensor(edge_weight).float()

    edge_index, edge_weight = to_undirected(edge_index, edge_weight, reduce="mean")
    
    return edge_index, edge_weight

def remove_large_index(edge_index, edge_weight, max_index):
    large_index = edge_index > max_index
    keep_index = torch.logical_not(torch.logical_or(large_index[0], large_index[1]))
    
    edge_index = edge_index[:, keep_index]
    edge_weight = edge_weight[keep_index]
    
    return edge_index, edge_weight
    
def remove_large_index_hyperedge(hyperedges, max_index):
    result = []
    for p in hyperedges:
        r = []
        for pp in p:
            if pp > max_index:
                continue
            else:
                r.append(pp)
    
        if len(set(r)) > 2:
            result.append(tuple(set(r)))
    
    return result

def subgraph_hyperedge(hyperedges, sub_nodes):
    result = []
    
    sub_nodes_set = set(sub_nodes.tolist())
    nodes_min = min(sub_nodes_set)
    
    for p in hyperedges:
        r = []
        for pp in p:
            if pp in sub_nodes_set:
                r.append(pp - nodes_min)
            else:
                continue
    
        if len(r) > 2:
            result.append(tuple(r))
    
    return result

def data_process(config, chr_name):
    hic_dir_path = config["hic_dir_path"]
    porec_dir_path = config["porec_dir_path"]
    feature_dir_path = config["feature_dir_path"]

    edge_file_name = config["edge_file_name"].replace("chr0", chr_name)
    weight_file_name = config["weight_file_name"].replace("chr0", chr_name)
    feature_file_name = config["feature_file_name"].replace("chr0", chr_name)
    hyperedge_file_name = config["hyperedge_file_name"].replace("chr0", chr_name)
    
    porec_dirs = config["porec_dirs"]
    
    x = np.load(feature_dir_path + feature_file_name, allow_pickle=True).astype(float).T
    x = np.nan_to_num(x)

    scalerx = MinMaxScaler()
    x = scalerx.fit_transform(x)
    x = torch.tensor(x).float()

    edge_hic = np.load(hic_dir_path + edge_file_name)
    weight_hic = np.load(hic_dir_path + weight_file_name)
    edge_hic, weight_hic = to_undirected_mean(edge_hic, weight_hic)
    
    porec_hypers = []
    
    for porec_dir in porec_dirs:
        with open(porec_dir + hyperedge_file_name, 'rb') as f:
            porec_hyper = pickle.load(f)        
        
        for p in porec_hyper:
            if len(p) == 2:
                continue
            porec_hypers.append(p)
            
    porec_hypers = list(set(porec_hypers))
    porec_hypers = remove_large_index_hyperedge(porec_hypers, x.shape[0] - 1)
    edge_hic, weight_hic = remove_large_index(edge_hic, weight_hic, x.shape[0] - 1)
            
    return x, edge_hic, weight_hic, porec_hypers
    
def sliding_data_func(graph, size):
    sliding_data = []
    ptr = torch.nonzero(graph.x)[-1][0].item()
    start = torch.nonzero(graph.x)[0][0].item()
    
    for s in range(start, ptr, int(size)):
        if s + size < ptr:
            sub_nodes = torch.tensor(np.arange(s, s + size))
        else:
            sub_nodes = torch.tensor(np.arange(ptr - size, ptr))
            
        try:
            sub_edge_index, sub_edge_attr = subgraph(sub_nodes, graph.edge_index, graph.edge_attr, relabel_nodes=True)
        except:
            continue
        
        sub_x = graph.x[sub_nodes]
            
        sub_hyper_edge = subgraph_hyperedge(graph.hyper_edge, sub_nodes)
        
        mns, sns, cns = neg_generator(sub_hyper_edge, len(sub_hyper_edge))
        
        sliding_data.append(Data(x=sub_x, edge_index=sub_edge_index, edge_attr=sub_edge_attr, hyper_edge=sub_hyper_edge,
                                 mns=mns, sns=sns, cns=cns))
    
    return sliding_data

def data_process_all_chrom(config, size):
    data_list = []
    for chr_name in config["chr_names"]:
        x, edge_hic, weight_hic, porec_hypers = data_process(config, chr_name)
        graph = Data(x=x, edge_index=edge_hic, edge_attr=weight_hic, hyper_edge=porec_hypers)
        
        sliding_data = sliding_data_func(graph, size)
        
        data_list.extend(sliding_data)
    
    return data_list
