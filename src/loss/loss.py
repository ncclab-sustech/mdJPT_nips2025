import torch
import torch.nn as nn
from src.utils import top_k_accuracy
class SimCLRLoss(nn.Module):

    def __init__(self, temperature):
        super(SimCLRLoss, self).__init__()
        self.temperature = temperature
        self.CEL = torch.nn.CrossEntropyLoss()
        self.device = torch.device('cpu')
    
    def to(self, device):
        self.device = device
        self.CEL = self.CEL.to(device)
        return self

    def info_nce_loss(self, features):
        
        device = self.device

        # print(features.shape)
        bs = int(features.shape[0] // 2)
        labels = torch.cat([torch.arange(bs) for i in range(2)], dim=0)
        labels = (labels.unsqueeze(0) == labels.unsqueeze(1)).float()
        labels = labels.to(device)

        similarity_matrix = torch.matmul(features, features.T)

        # discard the main diagonal from both: labels and similarities matrix
        mask = torch.eye(labels.shape[0], dtype=torch.bool).to(device)
        labels = labels[~mask].view(labels.shape[0], -1)
        similarity_matrix = similarity_matrix[~mask].view(similarity_matrix.shape[0], -1)

        # select and combine multiple positives
        positives = similarity_matrix[labels.bool()].view(labels.shape[0], -1)

        # select only the negatives
        negatives = similarity_matrix[~labels.bool()].view(similarity_matrix.shape[0], -1)

        # Put the positive column at the end (when all the entries are the same, the top1 acc will be 0; while if the
        # positive column is at the start, the top1 acc might be exaggerated)
        logits = torch.cat([negatives, positives], dim=1)
        # The label means the last column contain the positive pairs
        labels = torch.ones(logits.shape[0], dtype=torch.long)*(logits.shape[1]-1)
        labels = labels.to(device)

        logits = logits / self.temperature
        return logits, labels

    def forward(self, features):
        # fea need to be normalized to 1
        self.to(features.device)
        logits, labels = self.info_nce_loss(features)
        loss = self.CEL(logits, labels)
        [acc_1, acc_5] = top_k_accuracy(logits, labels)
        return loss, logits, labels, [acc_1, acc_5]
