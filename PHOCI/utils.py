'''
Utilities functions for the framework.
'''
import numpy as np
import argparse
import torch
import warnings
warnings.filterwarnings('ignore')
from sklearn import metrics
from torchmetrics import AveragePrecision

     
def gen_size_dist(hyperedges):
    size_dist = {}
    for edge in hyperedges:
        leng = len(edge)
        if leng not in size_dist :
            size_dist[leng] = 0
        size_dist[leng] += 1
    if 1 in size_dist:
        del size_dist[1]
    if 2 in size_dist:
        del size_dist[2]
    total = sum(v for k, v in size_dist.items())
    for i in size_dist:
        size_dist[i] = float(size_dist[i]) / total
    return size_dist  


def measure(label, pred):
    #average_precision = AveragePrecision()
    auc_roc = metrics.roc_auc_score(np.array(label), np.array(pred))
    #ap = average_precision(torch.tensor(pred), torch.tensor(label))
    return auc_roc, auc_roc #ap
