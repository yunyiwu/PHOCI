#!/usr/bin/env python
# coding: utf-8

#####Using GM12878 cell-line as example

import time
import random
import pickle
import torch
import torch.nn as nn
import numpy as np
from sklearn import metrics

from config import config_train
import models
from aggregator import MaxminMeanAggregator
from data_load import load_val
from batch import HEBatchGenerator

print("GM12878")

def model_train(bce_loss, clip, bs, data, model, Aggregator, optim, hedges, labels, train_pred, train_label, device):
    batch_size = len(hedges) 

    model.train()
    Aggregator.train()
    optim.zero_grad()
        
    v = model(data.x, data.edge_index, data.edge_attr)
    
    preds = []
    embeds = []
    for hedge in hedges:
        embeddings = v[hedge]
        max_min_dist = torch.max(hedge) - torch.min(hedge)
        pred, embed = Aggregator(max_min_dist, embeddings)
        preds.append(pred)
        embeds.append(embed)
        train_pred.append(pred.detach())
        
    train_label.append(labels.detach())
    labels = labels.type(torch.FloatTensor).to(device)
    preds = torch.stack(preds)
    preds = preds.squeeze()
    
    loss = bce_loss(preds, labels)
      
    loss.backward()
    optim.step()

    for _, param in model.named_parameters():
        param.clamp(-clip, clip)
    for _, param in Aggregator.named_parameters():
        param.clamp(-clip, clip)
        
    return loss.item(), train_pred, train_label


def model_eval(data, test_batchloader, model, Aggregator):
    model.eval()
    Aggregator.eval()
    with torch.no_grad():
        total_pred = []
        test_label = []
        num_data = 0

        v = model(data.x, data.edge_index, data.edge_attr)
        
        while True:
            hedges, labels, is_last = test_batchloader.next()
            batch_size = len(hedges)
            num_data += batch_size
                
            for hedge in hedges:
                embeddings = v[hedge]
                max_min_dist = torch.max(hedge) - torch.min(hedge)
                pred, _ = Aggregator(max_min_dist, embeddings)
                total_pred.append(pred.detach())
            test_label.append(labels.detach())
            
            if is_last:
                break
                
        total_pred = torch.stack(total_pred)       
        test_label = torch.cat(test_label, dim=0)
        
    return total_pred.tolist(), test_label.tolist()


device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

dim_edge = 400
dim_vertex = 400

alpha_e = 0
alpha_v = 0

input_dim = 15 #features_num

# Hyperparameters
n_layers = 3
size = 1000
nv = size
memory_size = 0
lr = 0.0001
bs = 1024
clip = 0.01

sliding_data = []

for chr_name in config_train["chr_names"]:
    with open("chrom_data/GM12878/"+chr_name+"_sliding_data", "rb") as f:
        sliding_data_ = pickle.load(f)
    sliding_data.extend(sliding_data_)
        
print(len(sliding_data))

random.shuffle(sliding_data)

with open("chrom_data/GM12878/chr15_sliding_data", "rb") as f:
    test_data1 = pickle.load(f)
    
with open("chrom_data/K562/chr15_sliding_data", "rb") as f:
    test_data2 = pickle.load(f)

model_dir = False

model = models.GCNEncoder(input_dim, dim_vertex, dim_vertex, n_layers)
model.to(device)

if model_dir:
    model.load_state_dict(torch.load(model_dir+"model_epoch_4"))

cls_layers = [dim_vertex*2+2, dim_vertex, 256, 128, 32, 8, 1] 
Aggregator = MaxminMeanAggregator(dim_vertex, cls_layers)
Aggregator.to(device)

if model_dir:
    Aggregator.load_state_dict(torch.load(model_dir+"Aggregator_epoch_4"))


def evaluate(data):
    data = data.to(device)
        
    test_batchloader = load_val(data.hyper_edge, bs, device, label="pos")
    val_pred_pos, total_label_pos = model_eval(data, test_batchloader, model, Aggregator)
                
    test_batchloader = load_val(data.mns, bs, device, label="mns")
    val_pred_mns, total_label_mns = model_eval(data, test_batchloader, model, Aggregator)
    mns_auc = metrics.roc_auc_score(total_label_mns+total_label_pos, val_pred_mns+val_pred_pos)
    mns_ap = metrics.average_precision_score(total_label_mns+total_label_pos, val_pred_mns+val_pred_pos)
        
    test_batchloader = load_val(data.sns, bs, device, label="sns")
    val_pred_sns, total_label_sns = model_eval(data, test_batchloader, model, Aggregator)
    sns_auc = metrics.roc_auc_score(total_label_sns+total_label_pos, val_pred_sns+val_pred_pos)
    sns_ap = metrics.average_precision_score(total_label_sns+total_label_pos, val_pred_sns+val_pred_pos)       
        
    test_batchloader = load_val(data.cns, bs, device, label="cns")
    val_pred_cns, total_label_cns = model_eval(data, test_batchloader, model, Aggregator)
    cns_auc = metrics.roc_auc_score(total_label_cns+total_label_pos, val_pred_cns+val_pred_pos)
    cns_ap = metrics.average_precision_score(total_label_cns+total_label_pos, val_pred_cns+val_pred_pos)       
        
    print("Val mns AUC:"+str(mns_auc)+",Val mns AP:"+str(mns_ap))
    print("Val sns AUC:"+str(sns_auc)+",Val sns AP:"+str(sns_ap))
    print("Val cns AUC:"+str(cns_auc)+",Val cns AP:"+str(cns_ap))


best_roc = 0
best_epoch = 0 
optim = torch.optim.Adam(list(model.parameters())+list(Aggregator.parameters()), lr=lr)

bce_loss = nn.BCELoss()

epoch_num = 100

for epoch in range(epoch_num):
    
    s = time.time()
    
    for data in sliding_data:
        
        if len(data.hyper_edge) < bs:
            continue
        
        idcs = np.arange(len(data.hyper_edge))
        np.random.shuffle(idcs)
        
        third = int(len(data.hyper_edge)/3)
        
        mns = [data.mns[i] for i in idcs[0:third]]
        sns = [data.sns[i] for i in idcs[third:2*third]]
        cns = [data.cns[i] for i in idcs[2*third:]]
        
        train = data.hyper_edge + mns + sns + cns
        train_label = [1 for i in range(len(data.hyper_edge))] + [0 for i in range(len(data.hyper_edge))]
        
        train_batchloader = HEBatchGenerator(train, train_label, bs, device, test_generator=False) 

        data = data.to(device)

        train_pred, train_label = [], []
        loss_sum, count  = 0.0, 0.0
            
        # Train
        while True :
            s1 = time.time()
            hedges, labels, is_last = train_batchloader.next()
            loss,  train_pred, train_label = model_train(bce_loss, clip, bs, data, model, Aggregator, optim, hedges, labels, train_pred, train_label, device)
            loss_sum += loss
            count += 1
        
            e1 = time.time()
        
            if is_last:
                break
                
        train_pred = torch.stack(train_pred)
        train_pred = train_pred.squeeze()
        train_label = torch.round(torch.cat(train_label, dim=0))        
        train_auc = metrics.roc_auc_score(np.array(train_label.cpu()), np.array(train_pred.cpu()))
        train_ap = metrics.average_precision_score(np.array(train_label.cpu()), np.array(train_pred.cpu()))
        
        print("Epoch:"+str(epoch)+",Dloss:"+str(loss_sum)+",Train_AUC:"+str(train_auc)+",Train_AP:"+str(train_ap))
        
    torch.save(model.state_dict(), "models/GM12878/model_epoch_"+str(epoch))
    torch.save(Aggregator.state_dict(), "models/GM12878/Aggregator_epoch_"+str(epoch))
        
    if epoch%10 == 0:
        rand_index = random.randint(0, len(test_data1) - 1)
        print(rand_index)
        print("Same cell line:")
        evaluate(test_data1[rand_index])
        print("Cross cell line:")
        evaluate(test_data2[rand_index])        
    
    e = time.time()
    print(e-s)
