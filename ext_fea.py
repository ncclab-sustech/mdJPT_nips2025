import os
os.environ["CUDA_VISIBLE_DEVICES"]="0"
os.environ["WORLD_SIZE"]="1"
import numpy as np
from src.data.io_utils import load_finetune_EEG_data, get_load_data_func, load_processed_SEEDV_NEW_data
from src.data.data_process import running_norm_onesubsession, LDS, LDS_acc, LDS_gpu
from src.data.dataset import ext_Dataset 
import torch
from torch.utils.data import DataLoader
from src.model.MultiModel_PL import MultiModel_PL
import hydra
from omegaconf import DictConfig
import pytorch_lightning as pl
from tqdm import tqdm
import mne
from src.utils import video_order_load, reorder_vids_sepVideo, reorder_vids_back

def normTrain(data,data_train):
    temp = np.transpose(data_train,(0,1,3,2))
    temp = temp.reshape(-1,temp.shape[-1])
    data_mean = np.mean(temp, axis=0)
    data_var = np.var(temp, axis=0)
    data_normed = (data - data_mean.reshape(-1,1)) / np.sqrt(data_var + 1e-5).reshape(-1,1)
    return data_normed

def cal_fea(data,mode):
    if mode == 'de':
        fea = 0.5*np.log(2*np.pi*np.exp(1)*(np.var(data, 3)) + 1.0).squeeze()
        fea[fea<-40] = -40
    elif mode == 'me':
        fea = np.mean(data, axis=3).squeeze()
    return fea

@hydra.main(config_path="cfgs_multi", config_name="config_multi", version_base="1.3")
def ext_fea(cfg: DictConfig) -> None:
    pl.seed_everything(cfg.train.seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    load_dir = os.path.join(cfg.data_val.data_dir,'processed_data')
    print('data loading...')
    data, onesub_label, n_samples_onesub, n_samples_sessions = load_finetune_EEG_data(load_dir, cfg.data_val)
    print('data loaded')
    print(f'data ori shape:{data.shape}')
    print(f'n_samples2_onesub shape:{n_samples_onesub.shape}')
    print(f'n_samples2_sessions shape:{n_samples_sessions.shape}')
    # DATA SHAPE: [n_subs,session*vid*n_samples, n_chans, n_pionts]
    data = data.reshape(cfg.data_val.n_subs, -1, data.shape[-2], data.shape[-1])
    save_dir = os.path.join(cfg.data_val.data_dir,'ext_fea')
    if not os.path.exists(save_dir):
        os.makedirs(save_dir) 
    np.save(save_dir+'/onesub_label.npy',onesub_label)
    if cfg.val.extractor.normTrain:
        val_subs_all = cfg.data_val.val_subs_all
        if cfg.val.n_fold == "loo":
            val_subs_all = [[i] for i in range(cfg.data_val.n_subs)]
            n_folds = len(val_subs_all)
        else:
            n_folds = len(val_subs_all)
    else:
        n_folds = 1
    if cfg.val.extractor.use_pretrain:
        print('Use pretrain model:')
        cp_path = os.path.join('log', cfg.log.run_name, 'ckpt', f'epoch={(cfg.val.extractor.ckpt_epoch-1):02d}.ckpt')
        print(f'checkpoint load from: {cp_path}')
        cfg.data_cfg_list = [cfg.data_0, cfg.data_1, cfg.data_2, cfg.data_3, cfg.data_4, cfg.data_val]
        cfg.data_cfg_list = [cfg_i for cfg_i in cfg.data_cfg_list if cfg_i.dataset_name != 'None']
        Extractor = MultiModel_PL.load_from_checkpoint(checkpoint_path=cp_path, cfg=cfg, strict=False)
        Extractor.save_fea = True
        Extractor.cnn_encoder.set_saveFea(True)
        trainer = pl.Trainer(accelerator='gpu', devices=1)
    for fold in tqdm(range(n_folds), desc='Extracting feature......'):
        if cfg.val.extractor.normTrain:
            val_subs = val_subs_all[fold]
            train_subs = list(set(range(cfg.data_val.n_subs)) - set(val_subs))
        else:
            val_subs = list(range(cfg.data_val.n_subs))
            train_subs = list(range(cfg.data_val.n_subs))
        if cfg.val.extractor.reverse:
            train_subs, val_subs = val_subs, train_subs
        print(f'train_subs:{train_subs}')
        print(f'val_subs:{val_subs}' )
        data_train = data[train_subs]
        if cfg.val.extractor.normTrain:
            print('normTraining')
            data_fold = normTrain(data,data_train)
            print('normDone')
        else:
            data_fold = data
        if cfg.val.extractor.use_pretrain:
            print('Use pretrain model:')
            data_fold = data_fold.reshape(-1, data_fold.shape[-2], data_fold.shape[-1])
            label_fold = np.tile(onesub_label, cfg.data_val.n_subs)
            print(data_fold.shape)
            foldset = ext_Dataset(data_fold, label_fold)
            del data_fold, label_fold
            fold_loader = DataLoader(foldset, batch_size=cfg.val.extractor.batch_size, shuffle=False, num_workers=cfg.train.num_workers)
            pred = trainer.predict(Extractor, fold_loader)
            mllaout = [pred_i[2] for pred_i in pred]
            pred = [pred_i[0] for pred_i in pred]
            pred = torch.cat(pred, dim=0).cpu().numpy()
            # Pred = [n_fea, dim, 1, T]
            fea = cal_fea(pred,cfg.val.extractor.fea_mode)
            fea = fea.reshape(cfg.data_val.n_subs,-1,fea.shape[-1])
        else:
            n_subs, n_samples, n_chans, sfreqs = data_fold.shape
            freqs = [[1,4], [4,8], [8,14], [14,30], [30,47]]
            de_data = np.zeros((n_subs, n_samples, n_chans, len(freqs)))
            n_samples_onesub_cum = np.concatenate((np.array([0]), np.cumsum(n_samples_onesub)))
            for idx, band in enumerate(freqs):
                for sub in tqdm(range(n_subs)):
                    print(f'sub:{sub}')
                    for vid in tqdm(range(len(n_samples_onesub)), desc=f'Direct DE Processing sub: {sub}', leave=False):
                        data_onevid = data_fold[sub,n_samples_onesub_cum[vid]:n_samples_onesub_cum[vid+1]]
                        data_onevid = data_onevid.transpose(1,0,2)
                        data_onevid = data_onevid.reshape(data_onevid.shape[0],-1)
                        data_video_filt = mne.filter.filter_data(data_onevid, sfreqs, l_freq=band[0], h_freq=band[1])
                        data_video_filt = data_video_filt.reshape(n_chans, -1, sfreqs)
                        de_onevid = 0.5*np.log(2*np.pi*np.exp(1)*(np.var(data_video_filt, 2))).T
                        de_data[sub,  n_samples_onesub_cum[vid]:n_samples_onesub_cum[vid+1], :, idx] = de_onevid
            fea = de_data.reshape(n_subs, n_samples, -1)
        fea_train = fea[train_subs]
        data_mean = np.mean(np.mean(fea_train, axis=1),axis=0)
        data_var = np.mean(np.var(fea_train, axis=1),axis=0)
        if np.isinf(fea).any():
            print("There are inf values in the array")
        else:
            print('No inf')
        if np.isnan(fea).any():
            print("There are nan values in the array")
        else:
            print('No nan')

        # reorder
        if cfg.data_val.dataset_name == 'FACED':
            vid_order = video_order_load(cfg.data_val.n_vids)
            if cfg.data_val.n_class == 3:
                n_vids = 28
            elif cfg.data_val.n_class == 9:
                n_vids = 28
            vid_inds = np.arange(n_vids)
            fea, vid_play_order_new = reorder_vids_sepVideo(fea, vid_order, vid_inds, n_vids)
        
        n_sample_sum_sessions = np.sum(n_samples_sessions,1)
        n_sample_sum_sessions_cum = np.concatenate((np.array([0]), np.cumsum(n_sample_sum_sessions)))
        # print(f'before norm:{fea.shape}')
        for sub in tqdm(range(cfg.data_val.n_subs), desc='Running norm......'):
            for s in tqdm(range(len(n_sample_sum_sessions)), desc=f'running norm sub: {sub}', leave=False):
                fea[sub,n_sample_sum_sessions_cum[s]:n_sample_sum_sessions_cum[s+1]] = running_norm_onesubsession(
                                            fea[sub,n_sample_sum_sessions_cum[s]:n_sample_sum_sessions_cum[s+1]],
                                            data_mean,data_var,cfg.val.extractor.rn_decay)
        # print(f'before LDS:{fea.shape}')
        if np.isinf(fea).any():
            print("There are inf values in the array")
        else:
            print('no inf')
        if np.isnan(fea).any():
            print("There are nan values in the array 2")
        else:
            print('no nan')
        # order back
        if cfg.data_val.dataset_name == 'FACED':
            fea = reorder_vids_back(fea, len(vid_inds), vid_play_order_new)
        
        n_samples_onesub_cum = np.concatenate((np.array([0]), np.cumsum(n_samples_onesub)))
        # LDS
        if(cfg.val.extractor.LDS):
            for sub in tqdm(range(cfg.data_val.n_subs), desc='LDS......'):
                for vid in tqdm(range(len(n_samples_onesub)), desc=f'LDS Processing sub: {sub}', leave=False):
                    fea[sub,n_samples_onesub_cum[vid]:n_samples_onesub_cum[vid+1]] = LDS_gpu(fea[sub,n_samples_onesub_cum[vid]:n_samples_onesub_cum[vid+1]])
                # print('LDS:',fea[sub,0])
        fea = fea.reshape(-1,fea.shape[-1])
        if np.isinf(fea).any():
            print("There are inf values in the array")
        else:
            print('no inf')
        if np.isnan(fea).any():
            print("There are nan values in the array")
        else:
            print('no nan')
        if not cfg.val.extractor.normTrain:
            save_path = os.path.join(save_dir,cfg.log.run_name+f'_all_fea_{f'epoch={(cfg.val.extractor.ckpt_epoch-1):02d}.ckpt' if cfg.val.extractor.use_pretrain else ""}{cfg.val.extractor.fea_mode if cfg.val.extractor.use_pretrain else cfg.val.extractor.fea_mode}.npy')
        else:
            save_path = os.path.join(save_dir,cfg.log.run_name+f'_f{fold}_fea_{f'epoch={(cfg.val.extractor.ckpt_epoch-1):02d}.ckpt' if cfg.val.extractor.use_pretrain else ""}{cfg.val.extractor.fea_mode if cfg.val.extractor.use_pretrain else cfg.val.extractor.fea_mode}.npy')
        np.save(save_path,fea)
        print(f'fea saved to {save_path}')
        if not cfg.val.extractor.normTrain:
            break


if __name__ == '__main__':
    ext_fea()