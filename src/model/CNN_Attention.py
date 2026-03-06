import torch.nn as nn
import torch.nn.functional as F
import torch
import numpy as np
from src.model.Channel_MLP import Channel_mlp_CNN

def stratified_layerNorm(out, n_samples):
    n_subs = int(out.shape[0] / n_samples)
    out_str = out.clone()
    for i in range(n_subs):
        out_oneSub = out[n_samples*i: n_samples*(i+1)]
        out_oneSub = out_oneSub.reshape(out_oneSub.shape[0], -1, out_oneSub.shape[-1]).permute(0,2,1)
        out_oneSub = out_oneSub.reshape(out_oneSub.shape[0]*out_oneSub.shape[1], -1)
        # out_oneSub[torch.isinf(out_oneSub)] = -50
        # out_oneSub[torch.isnan(out_oneSub)] = -50
        out_oneSub_str = out_oneSub.clone()
        # We don't care about the channels with very small activations
        # out_oneSub_str[:, out_oneSub.abs().sum(dim=0) > 1e-4] = (out_oneSub[:, out_oneSub.abs().sum(dim=0) > 1e-4] - out_oneSub[
        #     :, out_oneSub.abs().sum(dim=0) > 1e-4].mean(dim=0)) / (out_oneSub[:, out_oneSub.abs().sum(dim=0) > 1e-4].std(dim=0) + 1e-3)
        out_oneSub_str = (out_oneSub - out_oneSub.mean(dim=0)) / (out_oneSub.std(dim=0) + 1e-3)
        out_str[n_samples*i: n_samples*(i+1)] = out_oneSub_str.reshape(n_samples, -1, out_oneSub_str.shape[1]).permute(0,2,1).reshape(n_samples, out.shape[1], out.shape[2], -1)
        # out_str[torch.isnan(out_str)]=1
    return out_str

class cnn_PatchTST(nn.Module):
    # 配置说明 125Hz采样率基线  使用参数  dilation_array=[1,3,6,12]      seg_att = 15  avgPoolLen = 15  timeSmootherLen=3 mslen = 2,3   如果频率变化请在基线上乘以相应倍数
    def __init__(self, n_timeFilters, timeFilterLen, n_msFilters, msFilter_timeLen, n_channs=64, dilation_array=np.array([1,6,12,24]), seg_att=30, avgPoolLen = 30,
                  timeSmootherLen=6, multiFact=2, stratified=[], activ='softmax', temp=1.0, saveFea=True, has_att=True, extract_mode='me', global_att=False):
        super().__init__()
        self.stratified = stratified
        self.msFilter_timeLen = msFilter_timeLen
        self.activ = activ
        self.temp = temp
        self.dilation_array = np.array(dilation_array)   
        self.saveFea = saveFea
        self.has_att = has_att
        self.extract_mode = extract_mode
        self.global_att = global_att
        

        # time and spacial conv
        self.timeConv = nn.Conv2d(1, n_timeFilters, (1, timeFilterLen), padding=(0, (timeFilterLen-1)//2))
        self.msConv1 = nn.Conv2d(n_timeFilters, n_timeFilters*n_msFilters, (n_channs, msFilter_timeLen), groups=n_timeFilters)
        self.msConv2 = nn.Conv2d(n_timeFilters, n_timeFilters*n_msFilters, (n_channs, msFilter_timeLen), dilation=(1,self.dilation_array[1]), groups=n_timeFilters)
        self.msConv3 = nn.Conv2d(n_timeFilters, n_timeFilters*n_msFilters, (n_channs, msFilter_timeLen), dilation=(1,self.dilation_array[2]), groups=n_timeFilters)
        self.msConv4 = nn.Conv2d(n_timeFilters, n_timeFilters*n_msFilters, (n_channs, msFilter_timeLen), dilation=(1,self.dilation_array[3]), groups=n_timeFilters)

        n_msFilters_total = n_timeFilters * n_msFilters * 4

        # Attention
        self.seg_att = seg_att               #  *2 等比缩放
        self.att_conv = nn.Conv2d(n_msFilters_total, n_msFilters_total, (1, self.seg_att), groups=n_msFilters_total)
        self.att_pool = nn.AvgPool2d((1, self.seg_att), stride=1)
        self.att_pointConv = nn.Conv2d(n_msFilters_total, n_msFilters_total, (1, 1))

        # projector avepooling+timeSmooth
        self.avgpool = nn.AvgPool2d((1, avgPoolLen))
        self.timeConv1 = nn.Conv2d(n_msFilters_total, n_msFilters_total * multiFact, (1, timeSmootherLen), groups=n_msFilters_total)
        self.timeConv2 = nn.Conv2d(n_msFilters_total * multiFact, n_msFilters_total * multiFact * multiFact, (1, timeSmootherLen), groups=n_msFilters_total * multiFact)
        # # pooling  时间上的max pooling目前不需要，因为最后输出层特征会整体做个时间上的平均,时间上用ave比max更符合直觉
        # self.maxPoolLen = maxPoolLen
        # self.maxpool = nn.MaxPool2d((1, self.maxPoolLen),self.maxPoolLen)
        # # self.flatten = nn.Flatten()
    
    def forward(self, input):
        # input.shape should be [B, dim, n_channel, T]
        if 'initial' in self.stratified:
            input = stratified_layerNorm(input, int(input.shape[0]/2))
        # out = self.timeConv(input)
        out = input
        p = self.dilation_array * (self.msFilter_timeLen - 1)
        out1 = self.msConv1(F.pad(out, (int(p[0]//2), p[0]-int(p[0]//2)), "constant", 0))
        out2 = self.msConv2(F.pad(out, (int(p[1]//2), p[1]-int(p[1]//2)), "constant", 0))
        out3 = self.msConv3(F.pad(out, (int(p[2]//2), p[2]-int(p[2]//2)), "constant", 0))
        out4 = self.msConv4(F.pad(out, (int(p[3]//2), p[3]-int(p[3]//2)), "constant", 0))
        # print(out1.shape, out2.shape, out3.shape, out4.shape)
        out = torch.cat((out1, out2, out3, out4), 1) # (B, dims, 1, T)

        # Attention
        if self.has_att:
            att_w = F.relu(self.att_conv(F.pad(out, (self.seg_att-1, 0), "constant", 0)))
            if self.global_att:
                att_w = torch.mean(F.pad(att_w, (self.seg_att-1, 0), "constant", 0),-1).unsqueeze(-1) # (B, dims, 1, 1)
            else:
                att_w = self.att_pool(F.pad(att_w, (self.seg_att-1, 0), "constant", 0)) # (B, dims, 1, T)
            att_w = self.att_pointConv(att_w)
            if self.activ == 'relu':
                att_w = F.relu(att_w)
            elif self.activ == 'softmax':
                att_w = F.softmax(att_w / self.temp, dim=1)
            out = att_w * F.relu(out)          # (B, dims, 1, T)
        else:
            if self.extract_mode == 'me':
                out = F.relu(out)
        if self.saveFea:
            return out
        else:         # projecter
            if self.extract_mode == 'de':
                out = F.relu(out)
            out = self.avgpool(out)    # B*(t_dim*n_msFilters*4)*1*t_pool
            if 'middle1' in self.stratified:
                out = stratified_layerNorm(out, int(out.shape[0]/2))
            out = F.relu(self.timeConv1(out))
            out = F.relu(self.timeConv2(out))          #B*(t_dim*n_msFilters*4*multiFact*multiFact)*1*t_pool
            if 'middle2' in self.stratified:
                out = stratified_layerNorm(out, int(out.shape[0]/2))     
            proj_out = out.reshape(out.shape[0], -1)
            return F.normalize(proj_out, dim=1)

    
    def set_saveFea(self, saveFea):
        self.saveFea = saveFea

    def set_stratified(self,stratified):
        self.stratified = stratified

class Conv_att_simple_mlp(nn.Module):
    # 配置说明 125Hz采样率基线  使用参数  dilation_array=[1,3,6,12]      seg_att = 15  avgPoolLen = 15  timeSmootherLen=3 mslen = 2,3   如果频率变化请在基线上乘以相应倍数
    def __init__(self, n_timeFilters, timeFilterLen, n_msFilters, msFilter_timeLen, n_channs=64, dilation_array=np.array([1,6,12,24]), seg_att=30, avgPoolLen = 30,
                  timeSmootherLen=6, multiFact=2, stratified=[], activ='softmax', temp=1.0, saveFea=True, has_att=True, extract_mode='me', global_att=False,
                  c_mlps=None):
        super().__init__()
        self.stratified = stratified
        self.msFilter_timeLen = msFilter_timeLen
        self.activ = activ
        self.temp = temp
        self.dilation_array = np.array(dilation_array)   
        self.saveFea = saveFea
        self.has_att = has_att
        self.extract_mode = extract_mode
        self.global_att = global_att
        

        # time and spacial conv
        self.timeConv = nn.Conv2d(1, n_timeFilters, (1, timeFilterLen), padding=(0, (timeFilterLen-1)//2))
        self.msConv1 = nn.Conv2d(n_timeFilters, n_timeFilters*n_msFilters, (n_channs, msFilter_timeLen), groups=n_timeFilters)
        self.msConv2 = nn.Conv2d(n_timeFilters, n_timeFilters*n_msFilters, (n_channs, msFilter_timeLen), dilation=(1,self.dilation_array[1]), groups=n_timeFilters)
        self.msConv3 = nn.Conv2d(n_timeFilters, n_timeFilters*n_msFilters, (n_channs, msFilter_timeLen), dilation=(1,self.dilation_array[2]), groups=n_timeFilters)
        self.msConv4 = nn.Conv2d(n_timeFilters, n_timeFilters*n_msFilters, (n_channs, msFilter_timeLen), dilation=(1,self.dilation_array[3]), groups=n_timeFilters)

        n_msFilters_total = n_timeFilters * n_msFilters * 4

        # Attention
        self.seg_att = seg_att               #  *2 等比缩放
        self.att_conv = nn.Conv2d(n_msFilters_total, n_msFilters_total, (1, self.seg_att), groups=n_msFilters_total)
        self.att_pool = nn.AvgPool2d((1, self.seg_att), stride=1)
        self.att_pointConv = nn.Conv2d(n_msFilters_total, n_msFilters_total, (1, 1))

        # projector avepooling+timeSmooth
        self.avgpool = nn.AvgPool2d((1, avgPoolLen))
        self.timeConv1 = nn.Conv2d(n_msFilters_total, n_msFilters_total * multiFact, (1, timeSmootherLen), groups=n_msFilters_total)
        self.timeConv2 = nn.Conv2d(n_msFilters_total * multiFact, n_msFilters_total * multiFact * multiFact, (1, timeSmootherLen), groups=n_msFilters_total * multiFact)
        # # pooling  时间上的max pooling目前不需要，因为最后输出层特征会整体做个时间上的平均,时间上用ave比max更符合直觉
        # self.maxPoolLen = maxPoolLen
        # self.maxpool = nn.MaxPool2d((1, self.maxPoolLen),self.maxPoolLen)
        # # self.flatten = nn.Flatten()
        self.c_mlps = c_mlps
    
    def forward(self, input, dataset=0):
        # input.shape should be [B, dim, n_channel, T]
        if 'initial' in self.stratified:
            input = stratified_layerNorm(input, int(input.shape[0]/2))
        out = self.timeConv(input)
        out = self.c_mlps[dataset](out)
        out = stratified_layerNorm(out, int(out.shape[0]/2))
        p = self.dilation_array * (self.msFilter_timeLen - 1)
        out1 = self.msConv1(F.pad(out, (int(p[0]//2), p[0]-int(p[0]//2)), "constant", 0))
        out2 = self.msConv2(F.pad(out, (int(p[1]//2), p[1]-int(p[1]//2)), "constant", 0))
        out3 = self.msConv3(F.pad(out, (int(p[2]//2), p[2]-int(p[2]//2)), "constant", 0))
        out4 = self.msConv4(F.pad(out, (int(p[3]//2), p[3]-int(p[3]//2)), "constant", 0))
        # print(out1.shape, out2.shape, out3.shape, out4.shape)
        out = torch.cat((out1, out2, out3, out4), 1) # (B, dims, 1, T)

        # Attention
        if self.has_att:
            att_w = F.relu(self.att_conv(F.pad(out, (self.seg_att-1, 0), "constant", 0)))
            if self.global_att:
                att_w = torch.mean(F.pad(att_w, (self.seg_att-1, 0), "constant", 0),-1).unsqueeze(-1) # (B, dims, 1, 1)
            else:
                att_w = self.att_pool(F.pad(att_w, (self.seg_att-1, 0), "constant", 0)) # (B, dims, 1, T)
            att_w = self.att_pointConv(att_w)
            if self.activ == 'relu':
                att_w = F.relu(att_w)
            elif self.activ == 'softmax':
                att_w = F.softmax(att_w / self.temp, dim=1)
            out = att_w * F.relu(out)          # (B, dims, 1, T)
        else:
            if self.extract_mode == 'me':
                out = F.relu(out)
        if self.saveFea:
            return out
        else:         # projecter
            if self.extract_mode == 'de':
                out = F.relu(out)
            out = self.avgpool(out)    # B*(t_dim*n_msFilters*4)*1*t_pool
            if 'middle1' in self.stratified:
                out = stratified_layerNorm(out, int(out.shape[0]/2))
            out = F.relu(self.timeConv1(out))
            out = F.relu(self.timeConv2(out))          #B*(t_dim*n_msFilters*4*multiFact*multiFact)*1*t_pool
            if 'middle2' in self.stratified:
                out = stratified_layerNorm(out, int(out.shape[0]/2))     
            proj_out = out.reshape(out.shape[0], -1)
            return F.normalize(proj_out, dim=1)

    
    def set_saveFea(self, saveFea):
        self.saveFea = saveFea

    def set_stratified(self,stratified):
        self.stratified = stratified



class cnn_MLLA(nn.Module):
    # 配置说明 125Hz采样率基线  使用参数  dilation_array=[1,3,6,12]      seg_att = 15  avgPoolLen = 15  timeSmootherLen=3 mslen = 2,3   如果频率变化请在基线上乘以相应倍数
    def __init__(self, n_timeFilters, timeFilterLen, n_msFilters, msFilter_timeLen, n_channs=64, dilation_array=np.array([1,6,12,24]), seg_att=30, avgPoolLen = 30,
                  timeSmootherLen=6, multiFact=2, stratified=[], activ='softmax', temp=1.0, saveFea=True, has_att=True, extract_mode='me', global_att=False):
        super().__init__()
        self.stratified = stratified
        self.msFilter_timeLen = msFilter_timeLen
        self.activ = activ
        self.temp = temp
        self.dilation_array = np.array(dilation_array)   
        self.saveFea = saveFea
        self.has_att = has_att
        self.extract_mode = extract_mode
        self.global_att = global_att
        

        # time and spacial conv
        self.timeConv = nn.Conv2d(1, n_timeFilters, (1, timeFilterLen), padding=(0, (timeFilterLen-1)//2))
        self.msConv1 = nn.Conv2d(n_timeFilters, n_timeFilters*n_msFilters, (n_channs, msFilter_timeLen), groups=n_timeFilters)
        self.msConv2 = nn.Conv2d(n_timeFilters, n_timeFilters*n_msFilters, (n_channs, msFilter_timeLen), dilation=(1,self.dilation_array[1]), groups=n_timeFilters)
        self.msConv3 = nn.Conv2d(n_timeFilters, n_timeFilters*n_msFilters, (n_channs, msFilter_timeLen), dilation=(1,self.dilation_array[2]), groups=n_timeFilters)
        self.msConv4 = nn.Conv2d(n_timeFilters, n_timeFilters*n_msFilters, (n_channs, msFilter_timeLen), dilation=(1,self.dilation_array[3]), groups=n_timeFilters)

        n_msFilters_total = n_timeFilters * n_msFilters * 4

        # Attention
        self.seg_att = seg_att               #  *2 等比缩放
        self.att_conv = nn.Conv2d(n_msFilters_total, n_msFilters_total, (1, self.seg_att), groups=n_msFilters_total)
        self.att_pool = nn.AvgPool2d((1, self.seg_att), stride=1)
        self.att_pointConv = nn.Conv2d(n_msFilters_total, n_msFilters_total, (1, 1))

        # projector avepooling+timeSmooth
        self.avgpool = nn.AvgPool2d((1, avgPoolLen))
        self.timeConv1 = nn.Conv2d(n_msFilters_total, n_msFilters_total * multiFact, (1, timeSmootherLen), groups=n_msFilters_total)
        self.timeConv2 = nn.Conv2d(n_msFilters_total * multiFact, n_msFilters_total * multiFact * multiFact, (1, timeSmootherLen), groups=n_msFilters_total * multiFact)
        # # pooling  时间上的max pooling目前不需要，因为最后输出层特征会整体做个时间上的平均,时间上用ave比max更符合直觉
        # self.maxPoolLen = maxPoolLen
        # self.maxpool = nn.MaxPool2d((1, self.maxPoolLen),self.maxPoolLen)
        # # self.flatten = nn.Flatten()
    
    def forward(self, input):
        # input.shape should be [B, dim, n_channel, T]
        if 'initial' in self.stratified:
            input = stratified_layerNorm(input, int(input.shape[0]/2))
        # out = self.timeConv(input)
        out = input
        p = self.dilation_array * (self.msFilter_timeLen - 1)
        out1 = self.msConv1(F.pad(out, (int(p[0]//2), p[0]-int(p[0]//2)), "constant", 0))
        out2 = self.msConv2(F.pad(out, (int(p[1]//2), p[1]-int(p[1]//2)), "constant", 0))
        out3 = self.msConv3(F.pad(out, (int(p[2]//2), p[2]-int(p[2]//2)), "constant", 0))
        out4 = self.msConv4(F.pad(out, (int(p[3]//2), p[3]-int(p[3]//2)), "constant", 0))
        # print(out1.shape, out2.shape, out3.shape, out4.shape)
        out = torch.cat((out1, out2, out3, out4), 1) # (B, dims, 1, T)

        # Attention
        if self.has_att:
            att_w = F.relu(self.att_conv(F.pad(out, (self.seg_att-1, 0), "constant", 0)))
            if self.global_att:
                att_w = torch.mean(F.pad(att_w, (self.seg_att-1, 0), "constant", 0),-1).unsqueeze(-1) # (B, dims, 1, 1)
            else:
                att_w = self.att_pool(F.pad(att_w, (self.seg_att-1, 0), "constant", 0)) # (B, dims, 1, T)
            att_w = self.att_pointConv(att_w)
            if self.activ == 'relu':
                att_w = F.relu(att_w)
            elif self.activ == 'softmax':
                att_w = F.softmax(att_w / self.temp, dim=1)
            out = att_w * F.relu(out)          # (B, dims, 1, T)
        else:
            if self.extract_mode == 'me':
                out = F.relu(out)
        if self.saveFea:
            return out
        else:         # projecter
            if self.extract_mode == 'de':
                out = F.relu(out)
            out = self.avgpool(out)    # B*(t_dim*n_msFilters*4)*1*t_pool
            if 'middle1' in self.stratified:
                out = stratified_layerNorm(out, int(out.shape[0]/2))
            out = F.relu(self.timeConv1(out))
            out = F.relu(self.timeConv2(out))          #B*(t_dim*n_msFilters*4*multiFact*multiFact)*1*t_pool
            if 'middle2' in self.stratified:
                out = stratified_layerNorm(out, int(out.shape[0]/2))     
            proj_out = out.reshape(out.shape[0], -1)
            return F.normalize(proj_out, dim=1)

    
    def set_saveFea(self, saveFea):
        self.saveFea = saveFea

    def set_stratified(self,stratified):
        self.stratified = stratified