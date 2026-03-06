import torch
import pytorch_lightning as pl
import torch.nn as nn
import os
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from timm.models.layers import DropPath
from src.utils import stratified_layerNorm

class   channel_MLLA(nn.Module):
    def __init__(self, patch_size, hidden_dim, out_dim, depth, patch_stride, drop_path, n_filter, filterLen, n_heads):
        super().__init__()
        self.patch_size = patch_size
        self.patch_stride = patch_stride
        self.avgpool = nn.AdaptiveAvgPool1d(1)
        # self.encoders = nn.ModuleList([MLLA_BasicLayer(
        #                             in_dim=patch_size, hidden_dim=hidden_dim, out_dim=out_dim,
        #                             depth=1, num_heads=8) 
        #                             for _ in range(self.n_channels)])
        self.encoder = MLLA_BasicLayer(
                                    in_dim=patch_size, hidden_dim=hidden_dim, out_dim=out_dim,
                                    depth=depth, num_heads=n_heads)

    def forward(self, x):
        # print(x.shape)
        # x has shape [Batch, D1, n_channels, T]
        B, D1, n_channels, T = x.shape

        # Permute and reshape to combine batch and channel dimensions
        x = x.permute(0, 2, 1, 3)  # Shape: [Batch, n_channels, D1, T]
        x = x.reshape(B * n_channels, D1, T)  # Shape: [Batch * n_channels, D1, T]

        # Permute to match expected input shape for to_patch
        x = x.permute(0, 2, 1)  # Shape: [Batch * n_channels, T, D1]
        # print(x.shape)

        # Apply to_patch to all channels at once
        x = to_patch(x, patch_size=self.patch_size, stride=self.patch_stride)  # Shape: [Batch * n_channels, N, D1]

        # Pass through the encoder
        x = self.encoder(x)  # Shape: [Batch * n_channels, N, out_dim]

        # Reshape back to separate batch and channel dimensions
        x = x.view(B, n_channels, x.shape[1], x.shape[2])  # Shape: [Batch, n_channels, N, out_dim]

        # Permute dimensions as required
        x = x.permute(0, 1, 3, 2)  # Shape: [Batch, n_channels, out_dim, N]

        # Apply stratified_layerNorm
        x = stratified_layerNorm(x, n_samples=B // 2)

        # Permute back to original dimension order
        x = x.permute(0, 1, 3, 2)  # Shape: [Batch, n_channels, N, out_dim]

        return x

def to_patch(data, patch_size=50, stride=25):
    batchsize, timelen, dim = data.shape
    num_patches = (timelen - patch_size) // stride + 1
    patches = torch.zeros((batchsize, num_patches, patch_size)).to('cuda')
    for i in range(num_patches):
        start_idx = i * stride
        end_idx = start_idx + patch_size
        patches[:, i, :] = data[:, start_idx:end_idx, 0]
    return patches

class RoPE(nn.Module):
    def __init__(self, dim, max_seq_len=512, base=10000):
        super().__init__()
        position = torch.arange(max_seq_len, dtype=torch.float32).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, dim, 2, dtype=torch.float32) * (-torch.log(torch.tensor(base)) / dim))
        sinusoid_inp = position * div_term
        sin, cos = torch.sin(sinusoid_inp), torch.cos(sinusoid_inp)
        self.register_buffer('sin', sin)
        self.register_buffer('cos', cos)

    def forward(self, x):
        sin = self.sin[:x.size(-2), :].unsqueeze(0).unsqueeze(0)
        cos = self.cos[:x.size(-2), :].unsqueeze(0).unsqueeze(0)
        x1, x2 = x[..., ::2], x[..., 1::2]
        x = torch.stack([x1 * cos - x2 * sin, x1 * sin + x2 * cos], dim=-1)
        return x.flatten(-2)

class Mlp(nn.Module):
    def __init__(self, in_features, hidden_features=None, out_features=None, act_layer=nn.GELU, drop=0.):
        super().__init__()
        out_features = out_features or in_features
        hidden_features = hidden_features or in_features
        self.fc1 = nn.Linear(in_features, hidden_features)
        self.act = act_layer()
        self.drop1 = nn.Dropout(drop)
        self.fc2 = nn.Linear(hidden_features, out_features)
        self.drop2 = nn.Dropout(drop)

    def forward(self, x):
        x = self.fc1(x)
        x = self.act(x)
        x = self.drop1(x)
        x = self.fc2(x)
        x = self.drop2(x)
        return x

class LinearAttention(nn.Module):
    """Optimized Linear Attention with LePE and RoPE."""

    def __init__(self, dim, num_heads, qkv_bias=True):
        super().__init__()
        nn.Transformer
        self.dim = dim
        self.num_heads = num_heads
        self.head_dim = dim // num_heads
        self.scale = self.head_dim ** -0.5

        self.qkv = nn.Linear(dim, dim * 3, bias=qkv_bias)
        self.elu = nn.ELU()
        self.lepe = nn.Conv1d(dim, dim, kernel_size=3, padding=1, groups=dim)
        self.rope = RoPE(dim=self.head_dim, max_seq_len=512, base=10000)

    def forward(self, x):
        """
        Args:
            x: input features with shape of (B, N, C)
        """
        print(x.shape)
        B, N, C = x.shape
        H = self.num_heads
        D = self.head_dim

        # Compute Q, K, V
        qkv = self.qkv(x).reshape(B, N, 3, H, D).permute(2, 0, 3, 1, 4)  # (3, B, H, N, D)
        q, k, v = qkv[0], qkv[1], qkv[2]

        # Apply ELU activation and add 1
        q = self.elu(q) + 1
        k = self.elu(k) + 1

        # Apply RoPE
        q = self.rope(q)
        k = self.rope(k)

        # Scaled dot-product attention
        attn_scores = torch.matmul(q, k.transpose(-2, -1)) * self.scale  # (B, H, N, N)
        attn_weights = torch.softmax(attn_scores, dim=-1)
        attn_output = torch.matmul(attn_weights, v)  # (B, H, N, D)

        # Combine heads
        attn_output = attn_output.transpose(1, 2).reshape(B, N, C)  # (B, N, C)

        # Apply LePE
        v_conv = v.transpose(1, 2).reshape(B, N, C).permute(0, 2, 1)  # (B, C, N)
        lepe_output = self.lepe(v_conv).permute(0, 2, 1)  # (B, N, C)

        x = attn_output + lepe_output  # (B, N, C)

        return x


class vanillaMultiHeadAtt(nn.Module):
    """Optimized Linear Attention with LePE and RoPE."""

    def __init__(self, dim, num_heads, qkv_bias=True):
        super().__init__()
        self.dim = dim
        self.num_heads = num_heads
        self.head_dim = dim // num_heads
        self.scale = self.head_dim ** -0.5

        self.qkv = nn.Linear(dim, dim * 3, bias=qkv_bias)
        self.elu = nn.ELU()
        self.lepe = nn.Conv1d(dim, dim, kernel_size=3, padding=1, groups=dim)
        self.rope = RoPE(dim=self.head_dim, max_seq_len=512, base=10000)
        self.att = nn.MultiheadAttention(embed_dim=dim, num_heads=num_heads, batch_first=True)

    def forward(self, x):
        """
        Args:
            x: input features with shape of (B, N, C)
        """
        B, N, C = x.shape
        H = self.num_heads
        D = self.head_dim

        # Compute Q, K, V
        qkv = self.qkv(x).reshape(B, N, 3, C).permute(2, 0, 1, 3)
        q, k, v = qkv[0], qkv[1], qkv[2]
        attn_output, attn_output_weights = self.att(q, k, v)


        x = attn_output + x  # (B, N, C)

        return x

class MLLA_EEG_Block(nn.Module):
    def __init__(self, dim, num_heads, mlp_ratio=1., qkv_bias=True, drop=0., drop_path=0.,
                 act_layer=nn.GELU):
        super().__init__()
        self.cpe = nn.Conv1d(dim, dim, kernel_size=3, padding=1, groups=dim)
        self.norm1 = nn.LayerNorm(dim)
        self.attn = LinearAttention(dim=dim, num_heads=num_heads, qkv_bias=qkv_bias)
        # self.attn = vanillaMultiHeadAtt(dim=dim, num_heads=num_heads, qkv_bias=qkv_bias)
        self.drop_path = DropPath(drop_path) if drop_path > 0. else nn.Identity()
        self.norm2 = nn.LayerNorm(dim)
        self.mlp = Mlp(in_features=dim, hidden_features=int(dim * mlp_ratio),
                       act_layer=act_layer, drop=drop)

    def forward(self, x):
        # Convolutional Position Encoding
        x = x + self.cpe(x.permute(0, 2, 1)).permute(0, 2, 1)

        # Attention Block
        x_res = x
        x = self.norm1(x)
        x = self.attn(x)
        x = x_res + self.drop_path(x)

        # MLP Block
        x_res = x
        x = self.norm2(x)
        x = self.mlp(x)
        x = x_res + self.drop_path(x)

        return x

class MLLA_BasicLayer(nn.Module):
    """ A basic MLLA layer for one stage."""

    def __init__(self, in_dim, hidden_dim, out_dim, depth, num_heads, drop_path=0., qkv_bias=True):
        super().__init__()
        self.in_dim = in_dim
        self.hidden_dim = hidden_dim
        self.out_dim = out_dim
        self.depth = depth

        # Input and Output Projections
        self.read_in = nn.Linear(in_dim, hidden_dim)
        self.read_out = nn.Linear(hidden_dim, out_dim)

        # Build blocks
        self.blocks = nn.ModuleList([
            MLLA_EEG_Block(dim=hidden_dim, num_heads=num_heads, qkv_bias=qkv_bias, drop_path=drop_path)
            for _ in range(depth)
        ])

    def forward(self, x):
        x = F.relu(self.read_in(x))
        # print(x.shape)
        for blk in self.blocks:
            x = blk(x)
        x = F.relu(self.read_out(x))
        return x

