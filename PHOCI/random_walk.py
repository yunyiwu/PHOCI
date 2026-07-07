#!/usr/bin/env python
# coding: utf-8

import random
import torch
from torch_cluster import random_walk
import pickle


def random_walking(data, points, device, num_experiments = 5000000):
    
    size_dist_all = {}

    with open("annotation/size_dist_mean", "rb") as f:
        size_dist_mean = pickle.load(f)

    p_dict = size_dist_mean

    data = data.to(device)

    rands = []

    for _ in range(num_experiments):
        num_samples = random.choices(list(p_dict.keys()), weights=list(p_dict.values()))[0]
        sampled_start = random.choices(points, k=1)
    
        start = torch.tensor(sampled_start).to(device)
        walk = random_walk(data.edge_index[0], data.edge_index[1], start, num_samples-1, num_nodes=data.num_nodes)
        walk = walk.cpu().tolist()[0]
    
        rands.append(list(set(walk)))
    
    return rands
