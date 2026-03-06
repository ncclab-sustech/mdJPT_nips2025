import torch
import pytorch_lightning as pl
import torch.nn as nn
import os
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from timm.models.layers import DropPath
from src.utils import stratified_layerNorm, save_tensor_or_ndarray, report_vram
from torch import Tensor
from typing import Callable, Optional
from torch import nn, einsum

class   channel_MLLA(nn.Module):
    def __init__(self, context_window, patch_size, hidden_dim, out_dim, depth, patch_stride, n_heads):
        super().__init__()
        self.patch_size = patch_size
        self.patch_stride = patch_stride
        self.avgpool = nn.AdaptiveAvgPool1d(1)
        # self.encoders = nn.ModuleList([MLLA_BasicLayer(
        #                             in_dim=patch_size, hidden_dim=hidden_dim, out_dim=out_dim,
        #                             depth=1, num_heads=8) 
        #                             for _ in range(self.n_channels)])
        patch_num = int((context_window - patch_size)/patch_stride + 1)
        self.encoder = MLLA_BasicLayer(q_len=patch_num,
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
        x_res = torch.clone(x)
        x = self.norm1(x)
        print(x.shape)
        x = self.attn(x)
        x = x_res + self.drop_path(x)

        # MLP Block
        x_res = torch.clone(x)
        x = self.norm2(x)
        x = self.mlp(x)
        x = x_res + self.drop_path(x)

        return x

class MLLA_BasicLayer(nn.Module):
    """ A basic MLLA layer for one stage."""

    def __init__(self, q_len, in_dim, hidden_dim, out_dim, depth, num_heads, encoder_type='MLLA'):
        super().__init__()
        self.in_dim = in_dim
        self.hidden_dim = hidden_dim
        self.out_dim = out_dim
        self.depth = depth

        # Input and Output Projections
        self.read_in = nn.Linear(in_dim, hidden_dim)
        self.read_out = nn.Linear(hidden_dim, out_dim)

        # Build blocks
        if(encoder_type == 'Transformer'):
            self.blocks = nn.ModuleList([
                TransformerEncoderLayer(q_len=q_len, d_model=hidden_dim, n_heads=num_heads,
                                        d_k=None, d_v=None, d_ff=256, norm='BatchNorm',
                                                        attn_dropout=0, dropout=0,
                                                        activation='gelu', res_attention=False,
                                                        pre_norm=False, store_attn=False)
                for _ in range(depth)
            ])
        if(encoder_type == 'MLLA'):
            self.blocks = nn.ModuleList([
                MLLAEncoderLayer(q_len=q_len, d_model=hidden_dim, n_heads=num_heads,
                                        d_k=None, d_v=None, d_ff=256, norm='BatchNorm',
                                                        attn_dropout=0, dropout=0,
                                                        activation='gelu', res_attention=False,
                                                        pre_norm=False, store_attn=False)
                for _ in range(depth)
            ])

    def forward(self, x):
        x = self.read_in(x)
        # print(x.shape)
        for blk in self.blocks:
            x = blk(x)
        x = self.read_out(x)
        return x


class TransformerEncoderLayer(nn.Module):
    def __init__(self, q_len, d_model, n_heads, d_k=None, d_v=None, d_ff=256, store_attn=False,
                 norm='BatchNorm', attn_dropout=0, dropout=0., bias=True, activation="gelu", res_attention=False, pre_norm=False):
        super().__init__()
        assert not d_model%n_heads, f"d_model ({d_model}) must be divisible by n_heads ({n_heads})"
        d_k = d_model // n_heads if d_k is None else d_k
        d_v = d_model // n_heads if d_v is None else d_v

        # Multi-Head attention
        self.res_attention = res_attention
        self.self_attn = _MultiheadAttention(d_model, n_heads, d_k, d_v, attn_dropout=attn_dropout, proj_dropout=dropout, res_attention=res_attention)

        # Add & Norm
        self.dropout_attn = nn.Dropout(dropout)
        if "batch" in norm.lower():
            self.norm_attn = nn.Sequential(Transpose(1,2), nn.BatchNorm1d(d_model), Transpose(1,2))
        else:
            self.norm_attn = nn.LayerNorm(d_model)

        # Position-wise Feed-Forward
        self.ff = nn.Sequential(nn.Linear(d_model, d_ff, bias=bias),
                                get_activation_fn(activation),
                                nn.Dropout(dropout),
                                nn.Linear(d_ff, d_model, bias=bias))

        # Add & Norm
        self.dropout_ffn = nn.Dropout(dropout)
        if "batch" in norm.lower():
            self.norm_ffn = nn.Sequential(Transpose(1,2), nn.BatchNorm1d(d_model), Transpose(1,2))
        else:
            self.norm_ffn = nn.LayerNorm(d_model)

        self.pre_norm = pre_norm
        self.store_attn = store_attn


    def forward(self, src:Tensor, prev:Optional[Tensor]=None, key_padding_mask:Optional[Tensor]=None, attn_mask:Optional[Tensor]=None) -> Tensor:

        # Multi-Head attention sublayer
        if self.pre_norm:
            src = self.norm_attn(src)
        ## Multi-Head attention
        if self.res_attention:
            src2, attn, scores = self.self_attn(src, src, src, prev, key_padding_mask=key_padding_mask, attn_mask=attn_mask)
        else:
            src2, attn = self.self_attn(src, src, src, key_padding_mask=key_padding_mask, attn_mask=attn_mask)
        if self.store_attn:
            self.attn = attn
        ## Add & Norm
        src = src + self.dropout_attn(src2) # Add: residual connection with residual dropout
        if not self.pre_norm:
            src = self.norm_attn(src)

        # Feed-forward sublayer
        if self.pre_norm:
            src = self.norm_ffn(src)
        ## Position-wise Feed-Forward
        src2 = self.ff(src)
        ## Add & Norm
        src = src + self.dropout_ffn(src2) # Add: residual connection with residual dropout
        if not self.pre_norm:
            src = self.norm_ffn(src)

        if self.res_attention:
            return src, scores
        else:
            return src


class MLLAEncoderLayer(nn.Module):
    def __init__(self, q_len, d_model, n_heads, d_k=None, d_v=None, d_ff=256, store_attn=False,
                 norm='BatchNorm', attn_dropout=0, dropout=0., bias=True, activation="gelu", res_attention=False, pre_norm=False):
        super().__init__()
        assert not d_model%n_heads, f"d_model ({d_model}) must be divisible by n_heads ({n_heads})"
        d_k = d_model // n_heads if d_k is None else d_k
        d_v = d_model // n_heads if d_v is None else d_v

        # Multi-Head attention
        self.res_attention = res_attention
        self.self_attn = _MultiheadAttention(d_model, n_heads, d_k, d_v, attn_dropout=attn_dropout, proj_dropout=dropout, res_attention=res_attention)

        self.ff = nn.Sequential(nn.Linear(d_model, d_ff, bias=bias),
                                nn.GELU(),
                                nn.Dropout(dropout),
                                nn.Linear(d_ff, d_model, bias=bias))
        self.cpe1 = nn.Conv1d(d_model, d_model, 3, padding=1, groups=d_model)
        self.cpe2 = nn.Conv1d(d_model, d_model, 3, padding=1, groups=d_model)
        self.dwc = nn.Conv1d(d_model, d_model, 3, padding=1, groups=d_model)
        norm_layer=nn.LayerNorm
        self.norm1 = norm_layer(d_model)
        self.norm2 = norm_layer(d_model)
        self.act_proj = nn.Linear(d_model, d_model)
        self.in_proj = nn.Linear(d_model, d_model)
        self.out_proj = nn.Linear(d_model, d_model)
        self.act = nn.SiLU()


    def forward(self, src:Tensor, prev:Optional[Tensor]=None, key_padding_mask:Optional[Tensor]=None, attn_mask:Optional[Tensor]=None) -> Tensor:
        src = src + self.cpe1(src.permute(0, 2, 1)).permute(0, 2, 1)
        shortcut = src
        src = self.norm1(src)
        act_res = self.act(self.act_proj(src))
        src = self.in_proj(src)
        src = self.act(self.dwc(src.permute(0, 2, 1))).permute(0, 2, 1)
        # Linear Attention
        if self.res_attention:
            src2, attn, scores = self.self_attn(src, src, src, prev, key_padding_mask=key_padding_mask, attn_mask=attn_mask)
        else:
            src2, attn = self.self_attn(src, src, src, key_padding_mask=key_padding_mask, attn_mask=attn_mask)
        src2 = self.out_proj(src2 * act_res)
        src2 = shortcut + src2
        src2 = src2 + self.cpe2(src2.permute(0, 2, 1)).permute(0, 2, 1)
        # FFN
        src2 = src2 + self.ff(src2)
        if self.res_attention:
            return src2, scores
        else:
            return src2
class _MultiheadAttention(nn.Module):
    def __init__(self, d_model, n_heads, d_k=None, d_v=None, res_attention=False, attn_dropout=0., proj_dropout=0., qkv_bias=True, lsa=False):
        """Multi Head Attention Layer
        Input shape:
            Q:       [batch_size (bs) x max_q_len x d_model]
            K, V:    [batch_size (bs) x q_len x d_model]
            mask:    [q_len x q_len]
        """
        super().__init__()
        d_k = d_model // n_heads if d_k is None else d_k
        d_v = d_model // n_heads if d_v is None else d_v

        self.n_heads, self.d_k, self.d_v = n_heads, d_k, d_v

        self.W_Q = nn.Linear(d_model, d_k * n_heads, bias=qkv_bias)
        self.W_K = nn.Linear(d_model, d_k * n_heads, bias=qkv_bias)
        self.W_V = nn.Linear(d_model, d_v * n_heads, bias=qkv_bias)

        # Scaled Dot-Product Attention (multiple heads)
        self.res_attention = res_attention
        # self.sdp_attn = _ScaledDotProductAttention(d_model, n_heads, attn_dropout=attn_dropout, res_attention=self.res_attention, lsa=lsa)
        self.lnr_attn = _LinearAttention(d_model, n_heads, attn_dropout=attn_dropout, res_attention=self.res_attention, lsa=lsa)

        # Poject output
        self.to_out = nn.Sequential(nn.Linear(n_heads * d_v, d_model), nn.Dropout(proj_dropout))


    def forward(self, Q:Tensor, K:Optional[Tensor]=None, V:Optional[Tensor]=None, prev:Optional[Tensor]=None,
                key_padding_mask:Optional[Tensor]=None, attn_mask:Optional[Tensor]=None):

        bs = Q.size(0)
        if K is None: K = Q
        if V is None: V = Q

        # Linear (+ split in multiple heads)
        q_s = self.W_Q(Q).view(bs, -1, self.n_heads, self.d_k).transpose(1,2)       # q_s    : [bs x n_heads x max_q_len x d_k]
        k_s = self.W_K(K).view(bs, -1, self.n_heads, self.d_k).permute(0,2,3,1)     # k_s    : [bs x n_heads x d_k x q_len] - transpose(1,2) + transpose(2,3)
        v_s = self.W_V(V).view(bs, -1, self.n_heads, self.d_v).transpose(1,2)       # v_s    : [bs x n_heads x q_len x d_v]

        # Apply Scaled Dot-Product Attention (multiple heads)
        if self.res_attention:
            output, attn_weights, attn_scores = self.lnr_attn(q_s, k_s, v_s, prev=prev, key_padding_mask=key_padding_mask, attn_mask=attn_mask)
        else:
            output, attn_weights = self.lnr_attn(q_s, k_s, v_s, key_padding_mask=key_padding_mask, attn_mask=attn_mask)
        # output: [bs x n_heads x q_len x d_v], attn: [bs x n_heads x q_len x q_len], scores: [bs x n_heads x max_q_len x q_len]

        # back to the original inputs dimensions
        output = output.transpose(1, 2).contiguous().view(bs, -1, self.n_heads * self.d_v) # output: [bs x q_len x n_heads * d_v]
        output = self.to_out(output)

        if self.res_attention: return output, attn_weights, attn_scores
        else: return output, attn_weights


class _ScaledDotProductAttention(nn.Module):
    r"""Scaled Dot-Product Attention module (Attention is all you need by Vaswani et al., 2017) with optional residual attention from previous layer
    (Realformer: Transformer likes residual attention by He et al, 2020) and locality self sttention (Vision Transformer for Small-Size Datasets
    by Lee et al, 2021)"""

    def __init__(self, d_model, n_heads, attn_dropout=0., res_attention=False, lsa=False):
        super().__init__()
        self.attn_dropout = nn.Dropout(attn_dropout)
        self.res_attention = res_attention
        head_dim = d_model // n_heads
        self.scale = nn.Parameter(torch.tensor(head_dim ** -0.5), requires_grad=lsa)
        self.lsa = lsa

    def forward(self, q:Tensor, k:Tensor, v:Tensor, prev:Optional[Tensor]=None, key_padding_mask:Optional[Tensor]=None, attn_mask:Optional[Tensor]=None):
        '''
        Input shape:
            q               : [bs x n_heads x max_q_len x d_k]
            k               : [bs x n_heads x d_k x seq_len]
            v               : [bs x n_heads x seq_len x d_v]
            prev            : [bs x n_heads x q_len x seq_len]
            key_padding_mask: [bs x seq_len]
            attn_mask       : [1 x seq_len x seq_len]
        Output shape:
            output:  [bs x n_heads x q_len x d_v]
            attn   : [bs x n_heads x q_len x seq_len]
            scores : [bs x n_heads x q_len x seq_len]
        '''

        # Scaled MatMul (q, k) - similarity scores for all pairs of positions in an input sequence
        attn_scores = torch.matmul(q, k) * self.scale      # attn_scores : [bs x n_heads x max_q_len x q_len]

        # Add pre-softmax attention scores from the previous layer (optional)
        if prev is not None: attn_scores = attn_scores + prev

        # Attention mask (optional)
        if attn_mask is not None:                                     # attn_mask with shape [q_len x seq_len] - only used when q_len == seq_len
            if attn_mask.dtype == torch.bool:
                attn_scores.masked_fill_(attn_mask, -np.inf)
            else:
                attn_scores += attn_mask

        # Key padding mask (optional)
        if key_padding_mask is not None:                              # mask with shape [bs x q_len] (only when max_w_len == q_len)
            attn_scores.masked_fill_(key_padding_mask.unsqueeze(1).unsqueeze(2), -np.inf)

        # normalize the attention weights
        attn_weights = F.softmax(attn_scores, dim=-1)                 # attn_weights   : [bs x n_heads x max_q_len x q_len]
        attn_weights = self.attn_dropout(attn_weights)

        # compute the new values given the attention weights
        output = torch.matmul(attn_weights, v)                        # output: [bs x n_heads x max_q_len x d_v]

        if self.res_attention: return output, attn_weights, attn_scores
        else: return output, attn_weights

class _LinearAttention(nn.Module):
    r"""Scaled Dot-Product Attention module (Attention is all you need by Vaswani et al., 2017) with optional residual attention from previous layer
    (Realformer: Transformer likes residual attention by He et al, 2020) and locality self sttention (Vision Transformer for Small-Size Datasets
    by Lee et al, 2021)"""

    def __init__(self, d_model, n_heads, attn_dropout=0., res_attention=False, lsa=False):
        super().__init__()
        self.attn_dropout = nn.Dropout(attn_dropout)
        self.res_attention = res_attention
        head_dim = d_model // n_heads
        self.scale = nn.Parameter(torch.tensor(head_dim ** -0.5), requires_grad=lsa)
        self.lsa = lsa

    def forward(self, q:Tensor, k:Tensor, v:Tensor, prev:Optional[Tensor]=None, key_padding_mask:Optional[Tensor]=None, attn_mask:Optional[Tensor]=None):
        '''
        Input shape:
            q               : [bs x n_heads x max_q_len x d_k]
            k               : [bs x n_heads x d_k x seq_len]
            v               : [bs x n_heads x seq_len x d_v]
            prev            : [bs x n_heads x q_len x seq_len]
            key_padding_mask: [bs x seq_len]
            attn_mask       : [1 x seq_len x seq_len]
        Output shape:
            output:  [bs x n_heads x q_len x d_v]
            attn   : [bs x n_heads x q_len x seq_len]
            scores : [bs x n_heads x q_len x seq_len]
        '''
        bs, n_head, L, dim = q.shape
        eps = 1e-5
        k = k.permute(0, 1, 3, 2)
    
        # 重排列为 [bs, L, n_head, dim]
        q = q.permute(0, 2, 1, 3)  # [bs, L, n_head, dim]
        k = k.permute(0, 2, 1, 3)
        v = v.permute(0, 2, 1, 3)
        
        # 特征映射 (ELU+1)
        phi_Q = F.elu(q) + 1  # [bs, L, n_head, dim]
        phi_K = F.elu(k) + 1
        
        # 计算 KV = phi_K^T V (使用 einsum 避免维度混淆)
        KV = torch.einsum('blhd,blhm->bhdm', phi_K, v)  # [bs, n_head, dim, dim]
        
        # 计算归一化因子 Z = 1 / (phi_Q * sum(phi_K))
        K_sum = phi_K.sum(dim=1, keepdim=True)  # [bs, 1, n_head, dim]
        Z = 1.0 / (torch.einsum('blhd,bkhd->blh', phi_Q, K_sum) + eps)  # [bs, L, n_head]
        
        # 计算输出 V_new = phi_Q * KV * Z
        V_new = torch.einsum('blhd,bhdm->blhm', phi_Q, KV) * Z.unsqueeze(-1)  # [bs, L, n_head, dim]
        
        # 恢复原始维度 [bs, n_head, L, dim]
        return V_new.permute(0, 2, 1, 3), None


class Transpose(nn.Module):
    def __init__(self, *dims, contiguous=False): 
        super().__init__()
        self.dims, self.contiguous = dims, contiguous
    def forward(self, x):
        if self.contiguous: return x.transpose(*self.dims).contiguous()
        else: return x.transpose(*self.dims)

    
def get_activation_fn(activation):
    if callable(activation): return activation()
    elif activation.lower() == "relu": return nn.ReLU()
    elif activation.lower() == "gelu": return nn.GELU()
    raise ValueError(f'{activation} is not available. You can use "relu", "gelu", or a callable') 
    