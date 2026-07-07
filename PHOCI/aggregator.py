import torch
import torch.nn as nn
import numpy as np
import torch.nn.functional as F
    
class MaxminMeanAggregator(nn.Module):
    def __init__(self, dim_vertex, layers):
        super(MaxminMeanAggregator, self).__init__()
        Layers = []
        for i in range(len(layers)-1):
            Layers.append(nn.Linear(layers[i], layers[i+1]))
            if i != len(layers)-2:
                Layers.append(nn.ReLU(True))
        self.cls = nn.Sequential(*Layers)
    
    def aggregate(self, max_min_dist, embeddings):
        
        mean = embeddings.mean(dim=0).squeeze()
        
        max_val, _ = torch.max(embeddings, dim=0)
        min_val, _ = torch.min(embeddings, dim=0)
        
        max_min = max_val - min_val
        
        v_num = torch.tensor(embeddings.shape[0]).to(embeddings.device)
        
        return torch.cat((max_min_dist.view((1)), v_num.view((1)), max_min, mean))
    
    def classify(self, embedding):
        return F.sigmoid(self.cls(embedding))
    
    def forward(self, max_min_dist, embeddings):
        embedding = self.aggregate(max_min_dist, embeddings)
        pred = self.classify(embedding)
        return pred, embedding