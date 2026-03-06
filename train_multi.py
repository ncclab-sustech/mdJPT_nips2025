import hydra
from omegaconf import DictConfig
import os
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "0,2"
os.environ["WORLD_SIZE"]="1"
import torch
import numpy as np
import pytorch_lightning as pl
# from pytorch_lightning.loggers import WandbLogger
# import wandb
from pytorch_lightning.callbacks import ModelCheckpoint, EarlyStopping, Callback
from src.data.multi_dataloader import MultiDataModule
from src.model.MultiModel_PL import MultiModel_PL
import logging
from omegaconf import OmegaConf
from pytorch_lightning.loggers import TensorBoardLogger
from src.utils import *

log = logging.getLogger(__name__)

@hydra.main(config_path="cfgs_multi", config_name="config_multi", version_base="1.3")
def run_pipeline(cfg: DictConfig):
    pl.seed_everything(cfg.train.seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    # run_name = f'{cfg.log.run_name}_{get_current_time_string()}'
    run_name = cfg.log.run_name
    save_dir = os.path.join('log', run_name)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    config_save_path = os.path.join(save_dir, f'{run_name}.yaml')
    OmegaConf.save(config=cfg, f=config_save_path)
    print(f"Config saved to: {config_save_path}")
    n_folds = cfg.train.n_fold
    for fold in range(n_folds):
        print(f'fold {fold}')
        logger_save_dir = os.path.join(save_dir, str(fold))
        if not os.path.exists(logger_save_dir):
            os.makedirs(logger_save_dir)
            logger = TensorBoardLogger(save_dir=logger_save_dir, name=run_name)
        cfg.data_cfg_list = [cfg.data_0, cfg.data_1, cfg.data_2, cfg.data_3, cfg.data_4]
        cfg.data_cfg_list = [cfg_i for cfg_i in cfg.data_cfg_list if cfg_i.dataset_name != 'None']
        print(f'Using {len(cfg.data_cfg_list)} datasets to pretrain')
        dm = MultiDataModule(cfg.data_cfg_list, fold, n_folds, num_workers=cfg.train.num_workers,
                            n_pairs=cfg.train.n_pairs,
        )
        dm.setup("fit")
        model = MultiModel_PL(cfg)
        total_params = sum(p.numel() for p in model.parameters())
        total_size = sum(p.numel() * p.element_size() for p in model.parameters())
        log.info(f'Total number of parameters: {total_params}')
        log.info(f'Model size: {total_size} bytes ({total_size / (1024 ** 2):.2f} MB)')
        es_monitor = "loss_total/train" if n_folds == 1 else "loss_total/val"
        cp_dir = os.path.join(save_dir, 'ckpt')
        checkpoint_callback = ModelCheckpoint(dirpath=cp_dir, 
                                            filename='{epoch:02d}', every_n_epochs=cfg.train.save_interval,
                                            save_top_k=-1)
        # if all subs are used to pretrain, no will be used in validation.
        limit_val_batches = 0.0 if n_folds == 1 else 1.0
        trainer = pl.Trainer(
            callbacks=[checkpoint_callback],
            max_epochs=cfg.train.max_epochs, min_epochs=cfg.train.min_epochs, 
            accelerator='gpu', devices=cfg.train.n_gpu_use, strategy='ddp_find_unused_parameters_true',
            limit_val_batches=limit_val_batches
        )
        trainer.fit(model, dm)


    
    pass

if __name__ == '__main__':
    run_pipeline()