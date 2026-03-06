import torch
import torch.nn as nn
import torch.nn.functional as F

class Channel_mlp_CNN(nn.Module):

    def __init__(self, c_in, c_out):
        super().__init__()
        self.fc = nn.Linear(c_in, c_out)
        
    def forward(self, x):
        # Shape: [Batch, dim, n_chann, T]
        [Batch, _, n_chann, T] = x.shape
        x = x.permute(0, 1, 3, 2)
        # Shape: [Batch, 1, T, n_chann]
        x = F.relu(self.fc(x))
        x = x.permute(0, 1, 3, 2)
        return x