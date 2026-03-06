import hydra
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "1"
os.environ["WORLD_SIZE"] = "1"
from omegaconf import DictConfig
from src.model.valMLP import simpleNN3, MLPModel
import numpy as np
from src.data.dataset import PDataset
import pytorch_lightning as pl
from pytorch_lightning.callbacks import ModelCheckpoint
from torch.utils.data import DataLoader
import torch
import logging
from sklearn.preprocessing import label_binarize
from sklearn.metrics import precision_score, recall_score, f1_score, roc_auc_score, average_precision_score

torch.set_float32_matmul_precision('high')  # optional performance tweak

@hydra.main(config_path="cfgs_multi", config_name="config_multi", version_base="1.3")
def train_mlp(cfg: DictConfig) -> None:
    pl.seed_everything(cfg.val.mlp.seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    # prepare cross-validation folds
    val_subs_all = cfg.data_val.val_subs_all
    if cfg.val.n_fold == "loo":
        val_subs_all = [[i] for i in range(cfg.data_val.n_subs)]
    n_folds = len(val_subs_all)

    # storage for metrics
    accs, precisions, recalls, f1s, aurocs, auprcs = [], [], [], [], [], []

    for fold in range(n_folds):
        print(f"=== Fold {fold} ===")
        # checkpoint callback
        cp_dir = os.path.join(cfg.log.mlp_cp_dir, cfg.log.run_name)
        os.makedirs(cp_dir, exist_ok=True)
        checkpoint_callback = ModelCheckpoint(
            monitor="mlp/val/acc", verbose=True, mode="max",
            dirpath=cp_dir,
            filename=f'{cfg.data_val.dataset_name}_mlp_f{fold}_wd={cfg.val.mlp.wd}_{{epoch}}',
            save_top_k=1,
        )

        # split subjects
        val_subs = val_subs_all[fold]
        train_subs = list(set(range(cfg.data_val.n_subs)) - set(val_subs))
        if cfg.val.extractor.reverse:
            train_subs, val_subs = val_subs, train_subs
        print(f"Finetune subjects: {train_subs}")
        print(f"Test subjects:   {val_subs}")

        # load features
        save_dir = os.path.join(cfg.data_val.data_dir, 'ext_fea')
        save_path = os.path.join(
            save_dir,
            f"{cfg.log.run_name}_{f'{fold}' if cfg.val.extractor.normTrain else 'all'}_fea_"
            + (f"epoch={(cfg.val.extractor.ckpt_epoch-1):02d}.ckpt" if cfg.val.extractor.use_pretrain else "")
            + f"{cfg.val.extractor.fea_mode}.npy"
        )
        print(f'Using feature path: {save_path}')
        data = np.load(save_path)
        data = np.nan_to_num(data)
        data = data.reshape(cfg.data_val.n_subs, -1, data.shape[-1])

        # labels
        onesub_label = np.load(os.path.join(save_dir, 'onesub_label.npy'))
        labels = np.tile(onesub_label, cfg.data_val.n_subs)
        train_labels = np.tile(onesub_label, len(train_subs))
        val_labels = np.tile(onesub_label, len(val_subs))

        # datasets & loaders
        trainset = PDataset(data[train_subs].reshape(-1, data.shape[-1]), train_labels)
        valset   = PDataset(data[val_subs].reshape(-1, data.shape[-1]), val_labels)
        trainLoader = DataLoader(trainset, batch_size=cfg.val.mlp.batch_size, shuffle=True)
        valLoader   = DataLoader(valset,   batch_size=cfg.val.mlp.batch_size, shuffle=False)

        # model & trainer
        fea_dim = data.shape[-1]
        base_model = simpleNN3(fea_dim, cfg.val.mlp.hidden_dim, cfg.val.mlp.out_dim, 0.1)
        lightning_module = MLPModel(base_model, cfg.val.mlp)
        trainer = pl.Trainer(
            callbacks=[checkpoint_callback],
            max_epochs=cfg.val.mlp.max_epochs,
            min_epochs=cfg.val.mlp.min_epochs,
            accelerator='gpu', devices=1,
            limit_val_batches=1.0
        )

        # fit
        trainer.fit(lightning_module, trainLoader, valLoader)

        # load best checkpoint for metric computation
        best_ckpt = checkpoint_callback.best_model_path
        print(f"Loading best model from: {best_ckpt}")
        best_model = MLPModel.load_from_checkpoint(
            best_ckpt, model=base_model, cfg=cfg.val.mlp
        )
        best_model.eval()
        best_model.freeze()

        # gather predictions
        y_true, y_pred, y_prob = [], [], []
        for x_batch, y_batch in valLoader:
            logits = best_model.model(x_batch.to(best_model.device))
            probs  = torch.softmax(logits, dim=1).detach().cpu().numpy()
            preds  = probs.argmax(axis=1)
            y_true.append(y_batch.numpy())
            y_pred.append(preds)
            y_prob.append(probs)
        y_true = np.concatenate(y_true)
        y_pred = np.concatenate(y_pred)
        y_prob = np.concatenate(y_prob)
        print(np.unique(y_true))
        print(y_prob.shape)

        # compute metrics
        acc = (y_pred == y_true).mean()
        precision = precision_score(y_true, y_pred, average='macro', zero_division=0)
        recall    = recall_score(y_true, y_pred, average='macro', zero_division=0)
        f1        = f1_score(y_true, y_pred, average='macro', zero_division=0)
        n_classes = len(np.unique(y_true))
        if n_classes == 2:
            # 只提取正类的概率（假设正类是 label=1）
            y_prob_pos = y_prob[:, 1]
            auroc = roc_auc_score(y_true, y_prob_pos)
            auprc = average_precision_score(y_true, y_prob_pos)
        else:
            # 多分类情况，使用多输出版本
            y_true_bin = label_binarize(y_true, classes=np.arange(n_classes))
            auroc = roc_auc_score(y_true_bin, y_prob, multi_class='ovo', average='macro')
            auprc = average_precision_score(y_true_bin, y_prob, average='macro')

        # store
        accs.append(acc)
        precisions.append(precision)
        recalls.append(recall)
        f1s.append(f1)
        aurocs.append(auroc)
        auprcs.append(auprc)

        # print fold results
        print(f"Fold {fold} results:")
        print(f"  Acc   = {acc:.2f}")
        print(f"  Prec  = {precision:.2f}")
        print(f"  Rec   = {recall:.2f}")
        print(f"  F1    = {f1:.2f}")
        print(f"  AUROC = {auroc:.2f}")
        print(f"  AUPRC = {auprc:.2f}")

    # summary
    print("\n=== Summary across folds ===")
    metrics = {
        'Acc':      accs,
        'Precision':precisions,
        'Recall':   recalls,
        'F1':       f1s,
        'AUROC':    aurocs,
        'AUPRC':    auprcs
    }
    latex = ''
    for name, vals in metrics.items():
        mean = np.mean(vals)
        std  = np.std(vals)
        print(f"{name}: {mean*100:.2f} ± {std*100:.2f}")
        latex = latex + f"{mean*100:.2f} ± {std*100:.2f} & "
    print(latex)

if __name__ == '__main__':
    train_mlp()
