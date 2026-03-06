import torch
import torch.nn as nn
import pytorch_lightning as pl
from src.utils import stratified_layerNorm, save_tensor_or_ndarray, report_vram
class TemporalTransformer(pl.LightningModule):
    def __init__(self, n_chann, dim, dim_out, max_len=1024, n_heads=4, num_layers=1, dropout=0.1, lr=1e-3):
        super().__init__()
        self.save_hyperparameters()

        self.input_proj = nn.Linear(n_chann, dim)  # project each time step's n_chann -> dim

        # Learnable positional encoding: [1, max_len, dim]
        self.pos_embedding = nn.Parameter(torch.randn(1, max_len, dim))

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=dim, nhead=n_heads, dropout=dropout, batch_first=True
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)

        self.output_proj = nn.Linear(dim, dim_out)
        self.lr = lr

    def forward(self, x):
        # x: [B, 1, n_chann, T]
        B, _, n_chann, T = x.shape
        x = x.view(B, n_chann, T)     # [B, n_chann, T]
        x = x.permute(0, 2, 1)        # [B, T, n_chann]

        x = self.input_proj(x)        # [B, T, dim]
        x = x + self.pos_embedding[:, :T]  # Add position encoding

        x = self.transformer(x)       # [B, T, dim]
        x = self.output_proj(x)       # [B, T, dim]

        x = x.permute(0, 2, 1)        # [B, dim, T]
        x = x.unsqueeze(2)            # [B, dim, 1, T]
        x = x.expand(-1, -1, n_chann, -1)  # [B, dim, n_chann, T]
        x = stratified_layerNorm(x, n_samples=B // 2)
        return x