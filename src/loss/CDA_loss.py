import torch
import torch.nn as nn
from src.utils import top_k_accuracy
class CDALoss(nn.Module):

    def __init__(self, cfg):
        super(CDALoss, self).__init__()
        self.CEL = torch.nn.CrossEntropyLoss()
        self.device = torch.device('cpu')
        self.cfg = cfg
    def cov_mat(self, data):
            # 计算每个dim维度下的通道间协方差矩阵，使用批量矩阵乘法。
            # 参数：
            #     data (Tensor): 输入数据，形状为 [B, dim, C, T]。
            # 返回：
            #     cov_matrix (Tensor): 每个dim的协方差矩阵，形状为 [B, dim, C, C]。
            B, dim, C, T = data.shape
            
            # 计算均值，形状为 [B, dim, C, 1]
            mean_data = data.mean(dim=-1, keepdim=True)  # 计算每个通道在T维度上的均值
            data_centered = data - mean_data  # 去中心化，形状 [B, dim, C, T]
            # save_batch_images(data[0], 'cov_data')
            # 转置数据，形状变为 [B, dim, T, C]
            data_centered = data_centered.permute(0, 1, 3, 2)  # shape: [B, dim, T, C]
            # 将数据reshape成 [B * dim, C, T] 和 [B * dim, T, C]
            data_centered = data_centered.reshape(B * dim, C, T)  # shape: [B * dim, C, T]
            # 使用批量矩阵乘法计算协方差矩阵，结果形状为 [B * dim, C, C]
            cov_matrix = torch.bmm(data_centered, data_centered.permute(0, 2, 1)) / (T - 1)  # shape: [B * dim, C, C]
            # 将协方差矩阵重新reshape为 [B, dim, C, C]
            cov_matrix = cov_matrix.reshape(B, dim, C, C)  # shape: [B, dim, C, C]
            cov_matrix[:, :, range(C), range(C)] = 0
            return cov_matrix
    
    def CDA_loss(self, cov_mats):
        dis = 0
        ndim = cov_mats[0].shape[1]
        for dim in range(ndim):
            ind_cen = []
            for cov_sub_i in cov_mats:
                ind_cen.append(self.get_ind_cen(cov_sub_i[:, dim]))
            for i in range(len(ind_cen)):
                for j in range(i+1, len(ind_cen)):
                    dis = dis + self.frobenius_distance(ind_cen[i], ind_cen[j])
        loss = torch.log(dis + 1.0)
        return loss
    
    def frobenius_distance(self, matrix_a, matrix_b):
        return torch.linalg.norm(matrix_a-matrix_b, 'fro')
    
    def get_ind_cen(self, mat):
        if not self.cfg.train.loss.to_riem:
            return torch.mean(mat, dim=0)
        mat = torch.squeeze(mat)
        mat = self.frechet_mean(mat)
        return mat
    
    def forward(self, features):
        cov_mats = []
        for i, fea_i in enumerate(features):
            cov_mat_i = self.cov_mat(fea_i)
            cov_mats.append(cov_mat_i[:cov_mat_i.shape[0]//2])
            cov_mats.append(cov_mat_i[cov_mat_i.shape[0]//2:])
        align_loss = self.CDA_loss(cov_mats)
        return align_loss
