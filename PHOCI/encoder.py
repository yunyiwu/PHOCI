import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv

class GCNEncoder(nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels, num_layers):
        super(GCNEncoder, self).__init__()
        
        self.linear1 = nn.Linear(in_channels, hidden_channels)
        
        self.propagate = torch.nn.ModuleList()
        self.propagate.append(SAGEConv(hidden_channels, hidden_channels))
        
        for _ in range(num_layers - 2):
            self.propagate.append(SAGEConv(hidden_channels, hidden_channels))
        
        self.propagate.append(SAGEConv(hidden_channels, hidden_channels))
        
        
    def forward(self, x, edge_index, edge_attr):
        
        x = self.linear1(x)
        x = F.normalize(x,p=2,dim=1)*0.5
        
        for propagate in self.propagate:
            x = propagate(x, edge_index)
        
        return x