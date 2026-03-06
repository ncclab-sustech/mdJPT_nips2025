from pytorch_lightning.utilities.types import TRAIN_DATALOADERS
import torch
import pytorch_lightning as pl
import torch.nn as nn
import os
import numpy as np
from src.data.dataset import EEG_Dataset
from torch.utils.data import Dataset, DataLoader
import random
from itertools import cycle


def get_train_subs(n_subs, fold, n_folds):
    n_per = round(n_subs / n_folds)
    if n_folds == 1:
        val_subs = []
    elif fold < n_folds - 1:
        val_subs = range(n_per * fold, n_per * (fold + 1))
    else:
        val_subs = range(n_per * fold, n_subs)
    train_subs = list(set(range(n_subs)) - set(val_subs))
    return [train_subs, val_subs]


class EEGSampler:
    def __init__(self, datasets, n_pairs):
        self.n_pairs = n_pairs
        self.datasets = datasets  # List of datasets
        self.num_datasets = len(datasets)
        self.subs_list = []
        self.n_subs_list = []
        self.n_sessions_list = []
        self.pairs_list = []
        self.n_pairs_list = []
        self.max_n_pairs = 0

        # For preloading data
        self.save_dirs = []
        self.sliced_data_dirs = []
        self.n_samples_session_list = []
        self.n_vids_list = []
        self.batch_sizes = []
        self.n_sessions = []
        self.n_per_session_list = []
        self.n_per_session_cum_list = []
        self.n_samples_per_trial_list = []
        self.n_samples_cum_session_list = []

        # Initialize datasets
        for dataset in self.datasets:
            subs = dataset.train_subs if dataset.train_subs is not None else dataset.val_subs
            self.subs_list.append(subs)
            self.n_subs_list.append(len(subs))
            self.n_sessions_list.append(dataset.n_session)

        # Create pairs for each dataset
        for idx, dataset in enumerate(self.datasets):
            n_subs = self.n_subs_list[idx]
            n_sessions = self.n_sessions_list[idx]
            print(f'dataset {idx}: n_subs={n_subs} n_sessions={n_sessions}')
            pairs = []

            for i in range(n_subs * n_sessions):
                for j in range(i + n_sessions, n_subs * n_sessions, n_sessions):
                        pairs.append((i, j))

            random.shuffle(pairs)
            self.pairs_list.append(pairs)
            self.n_pairs_list.append(len(pairs))
            self.max_n_pairs = max(self.max_n_pairs, len(pairs))

            # Preload data for the dataset
            save_dir = os.path.join(dataset.cfg.data_dir, 'sliced_data')
            sliced_data_dir = os.path.join(
                save_dir, f'sliced_len{dataset.cfg.timeLen}_step{dataset.cfg.timeStep}')
            n_samples_session = np.load(
                os.path.join(sliced_data_dir, 'metadata', 'n_samples_sessions.npy'))
            n_vid = dataset.cfg.n_vids
            batch_size = n_vid
            n_session = dataset.cfg.n_session
            n_per_session = np.sum(n_samples_session, 1).astype(int)
            n_per_session_cum = np.concatenate((np.array([0]), np.cumsum(n_per_session)))
            n_samples_per_trial = int(n_vid / n_samples_session.shape[1])
            n_samples_cum_session = np.concatenate(
                (np.zeros((n_session, 1)), np.cumsum(n_samples_session, 1)), 1)

            self.save_dirs.append(save_dir)
            self.sliced_data_dirs.append(sliced_data_dir)
            self.n_samples_session_list.append(n_samples_session)
            self.n_vids_list.append(n_vid)
            self.batch_sizes.append(batch_size)
            self.n_sessions.append(n_session)
            self.n_per_session_list.append(n_per_session)
            self.n_per_session_cum_list.append(n_per_session_cum)
            self.n_samples_per_trial_list.append(n_samples_per_trial)
            self.n_samples_cum_session_list.append(n_samples_cum_session)

        print('Dataloader Lengths:')
        print(self.n_pairs_list)

    def get_sample(self, dataset_idx, subsession_pair):
        n_per_session = self.n_per_session_list[dataset_idx]
        n_per_session_cum = self.n_per_session_cum_list[dataset_idx]
        n_samples_cum_session = self.n_samples_cum_session_list[dataset_idx]
        n_session = self.n_sessions[dataset_idx]
        n_samples_per_trial = self.n_samples_per_trial_list[dataset_idx]
        n_sub = self.n_subs_list[dataset_idx]
        batch_size = self.batch_sizes[dataset_idx]

        subsession1, subsession2 = subsession_pair

        # Ensure both subsessions are from the same session
        cur_session = int(subsession1 % n_session)
        assert cur_session == int(subsession2 % n_session), "Subsessions are from different sessions"

        cur_sub1 = int(subsession1 // n_session)
        cur_sub2 = int(subsession2 // n_session)

        ind_abs_list = []

        n_trials = len(n_samples_cum_session[cur_session]) - 2
        for i in range(n_trials):
            start = int(n_samples_cum_session[cur_session][i])
            end = int(n_samples_cum_session[cur_session][i + 1])
            ind_one = np.random.choice(np.arange(start, end), n_samples_per_trial, replace=False)
            ind_abs_list.append(ind_one)

        # For the remaining samples
        i = n_trials
        start = int(n_samples_cum_session[cur_session][i])
        end = int(n_samples_cum_session[cur_session][i + 1])
        remaining_samples = int(batch_size - n_samples_per_trial * n_trials)
        ind_one = np.random.choice(np.arange(start, end), remaining_samples, replace=False)
        ind_abs_list.append(ind_one)

        ind_abs = np.concatenate(ind_abs_list)

        # Compute indices for both subsessions
        ind_this1 = ind_abs + np.sum(n_per_session) * cur_sub1 + n_per_session_cum[cur_session]
        ind_this2 = ind_abs + np.sum(n_per_session) * cur_sub2 + n_per_session_cum[cur_session]

        return ind_this1, ind_this2

    def __len__(self):
        return self.n_pairs

    def __iter__(self):
        for i in range(self.n_pairs):
            index = random.randint(0, self.max_n_pairs - 1)
            idx_list = []
            for idx in range(self.num_datasets):
                pairs = self.pairs_list[idx]
                n_pairs = self.n_pairs_list[idx]
                pair = pairs[index % n_pairs]
                idx_1, idx_2 = self.get_sample(dataset_idx=idx, subsession_pair=pair)
                idx_combined = np.concatenate((idx_1, idx_2))
                idx_list.append(list(idx_combined.astype(int)))
            yield idx_list


class MultiDataset(Dataset):
    def __init__(self, datasets):
        self.datasets = datasets
        self.lens = [len(dataset) for dataset in datasets]
        print('Dataset lengths:', self.lens)

    def __len__(self):
        return max(self.lens)

    def __getitem__(self, idx_list):
        data_list = []
        label_list = []
        for idxs, dataset in zip(idx_list, self.datasets):
            data = [dataset[i][0] for i in idxs]
            label = [dataset[i][1] for i in idxs]
            data = torch.stack(data)
            label = torch.stack(label)
            data_list.append(data)
            label_list.append(label)
        return data_list, label_list


class MultiDataModule(pl.LightningDataModule):
    def __init__(self, data_list, fold, n_folds, n_pairs=256, num_workers=8, device='cpu', sub_list_pre = None):
        super().__init__()
        self.device = device
        self.n_pairs = n_pairs
        self.num_workers = num_workers
        self.data_list = data_list  # List of data configs
        self.train_subs_list = []
        self.val_subs_list = []
        self.trainsets = []
        self.valsets = []
        self.datasets = []
        self.fold = fold
        self.n_folds = n_folds
        self.sub_list_pre = [None] * len(data_list) if sub_list_pre is None else sub_list_pre

    def prepare_data(self) -> None:
        pass

    def setup(self, stage=None):
        for index, data_cfg in enumerate(self.data_list):
            [train_subs, val_subs] = get_train_subs(data_cfg.n_subs, self.fold, self.n_folds) if self.sub_list_pre[index] is None else [self.sub_list_pre[index][0], self.sub_list_pre[index][1]]
            self.train_subs_list.append(train_subs)
            self.val_subs_list.append(val_subs)
            print(f'Dataset {index} \nTrain subs: {train_subs} \nVal subs: {val_subs}\n')

        if stage == 'fit' or stage is None:
            self.trainsets = []
            self.valsets = []
            for idx, data_cfg in enumerate(self.data_list):
                trainset = EEG_Dataset(data_cfg, train_subs=self.train_subs_list[idx], mods='train', sliced=False)
                valset = EEG_Dataset(data_cfg, val_subs=self.val_subs_list[idx], mods='val', sliced=False)
                self.trainsets.append(trainset)
                self.valsets.append(valset)
            self.trainset = MultiDataset(self.trainsets)
            self.valset = MultiDataset(self.valsets)

        if stage == 'validate':
            self.trainsets = []
            self.valsets = []
            for idx, data_cfg in enumerate(self.data_list):
                valset = EEG_Dataset(data_cfg, val_subs=self.val_subs_list[idx], mods='val', sliced=False)
                self.valsets.append(valset)
            self.valset = MultiDataset(self.valsets)

    def train_dataloader(self):
        self.trainsampler = EEGSampler(datasets=self.trainsets, n_pairs=self.n_pairs)
        self.trainloader = DataLoader(self.trainset, sampler=self.trainsampler,
                                      pin_memory=True, num_workers=self.num_workers,
                                      prefetch_factor=10)
        return self.trainloader

    def val_dataloader(self):
        self.valsampler = EEGSampler(datasets=self.valsets, n_pairs=int(self.n_pairs // 4))
        self.valloader = DataLoader(self.valset, sampler=self.valsampler,
                                    pin_memory=True, num_workers=self.num_workers)
        return self.valloader
