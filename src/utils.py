from datetime import datetime
import numpy as np
import torch
import matplotlib.pyplot as plt
import os
import scipy
import logging
import socket
from datetime import datetime, timedelta

def get_current_time_string():
    # 获取当前时间
    now = datetime.now()
    # 格式化为形如20250221123405的字符串
    time_string = now.strftime("%Y%m%d%H%M%S")
    return time_string

import numpy as np
import torch
import matplotlib.pyplot as plt
import os

def stratified_layerNorm(out, n_samples):
    n_samples = int(n_samples)
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
        # out_str[n_samples*i: n_samples*(i+1)] = out_oneSub_str
        # out_str[torch.isnan(out_str)]=1
    return out_str

def LDS(sequence):
    print(sequence.shape)
    # shape: (timeSample, n_dims)  timesample 为一个vid的采样数
    # print(sequence.shape) # (26, 256)   

    sequence_new = torch.zeros_like(sequence) # (26, 256)
    ave = torch.mean(sequence, axis=0)

    for te in range(sequence.shape[1]):
        X = sequence[:, te].reshape(1, -1) # (26,) to (1, 26)
        u0 = ave[te]

        V0 = 0.01
        A = 1
        T = 0.0001
        C = 1
        sigma = 1
        givenAll = 1

        [m, n] = X.shape # (1, 26)
        P = torch.zeros((m, m, n)) # (1, 1, 26)
        u = torch.zeros((m, n)) # (1, 26)
        V = torch.zeros((m, m, n)) # (1, 1, 26)
        K = torch.zeros((m, m, n)) # (1, 1, 26)

        K[:, :, 0] = V0*C / (C*V0*C + sigma)
        u[:, 0] = u0 + K[:, :, 0] * (X[:, 0] - C*u0)
        V[:, :, 0] = (torch.eye(m) - K[:, :, 0] * C) * V0

        for i in range(1, n):
            P[:, :, i-1] = A * V[:, :, i-1] * A + T
            K[:, :, i] = P[:, :, i-1] * C / (C * P[:, :, i-1] * C + sigma)
            u[:, i] = A * u[:, i-1] + K[:, :, i] * (X[:, i] - C*A*u[:, i-1])
            V[:, :, i] = (torch.eye(m) - K[:, :, i] * C) * P[:, :, i-1]

        if givenAll == 1:
            uAll = torch.zeros((m, n))
            VAll = torch.zeros((m, m, n))
            J = torch.zeros((m, m, n))
            uAll[:, n-1] = u[:, n-1]
            VAll[:, :, n-1] = V[:, :, n-1]

            for ir in range(n-1):
                i = n-2 - ir
                # print(i)
                J[:, :, i] = V[:, :, i] * A / P[:, :, i]
                uAll[:, i] = u[:, i] + J[:, :, i] * \
                    (uAll[:, i+1] - A * u[:, i])
                VAll[:, :, i] = V[:, :, i] + J[:, :, i] * \
                    (VAll[:, :, i+1] - P[:, :, i]) * J[:, :, i]

            X = uAll

        else:
            X = u

        sequence_new[:, te] = X
    return sequence_new

def LDS_new(sequence):
    # sequence: (B, n, n_dims)
    B, n, n_dims = sequence.shape

    # Compute the mean over the time axis
    ave = torch.mean(sequence, dim=1)  # (B, n_dims)

    # Permute sequence to shape (B, n_dims, n)
    X = sequence.permute(0, 2, 1)  # (B, n_dims, n)

    # Initial state mean
    u0 = ave  # (B, n_dims)

    # Define constants as tensors
    V0 = torch.tensor(0.01, dtype=sequence.dtype, device=sequence.device)
    A = torch.tensor(1.0, dtype=sequence.dtype, device=sequence.device)
    T = torch.tensor(0.0001, dtype=sequence.dtype, device=sequence.device)
    C = torch.tensor(1.0, dtype=sequence.dtype, device=sequence.device)
    sigma = torch.tensor(1.0, dtype=sequence.dtype, device=sequence.device)
    givenAll = 1

    # Initialize lists to collect values
    u_list = []
    V_list = []
    P_list = []
    K_list = []

    # Initial Kalman gain
    K_init = V0 * C / (C * V0 * C + sigma)  # Scalar tensor
    K0 = K_init.expand(B, n_dims)  # (B, n_dims)

    # Initial estimates
    u_prev = u0 + K0 * (X[:, :, 0] - C * u0)  # (B, n_dims)
    V_prev = (1 - K0 * C) * V0  # (B, n_dims)
    u_list.append(u_prev)
    V_list.append(V_prev)
    K_list.append(K0)

    # Forward pass (Kalman Filter)
    for i in range(1, n):
        P_prev = A * V_prev * A + T  # (B, n_dims)
        Denominator = C * P_prev * C + sigma  # (B, n_dims)
        K_curr = P_prev * C / Denominator  # (B, n_dims)
        u_curr = A * u_prev + K_curr * (X[:, :, i] - C * A * u_prev)  # (B, n_dims)
        V_curr = (1 - K_curr * C) * P_prev  # (B, n_dims)

        # Append current values to lists
        u_list.append(u_curr)
        V_list.append(V_curr)
        K_list.append(K_curr)
        P_list.append(P_prev)

        # Update previous values
        u_prev = u_curr
        V_prev = V_curr

    # Stack lists to form tensors of shape (B, n_dims, n)
    u = torch.stack(u_list, dim=2)  # (B, n_dims, n)
    V = torch.stack(V_list, dim=2)  # (B, n_dims, n)
    K = torch.stack(K_list, dim=2)  # (B, n_dims, n)

    # Backward pass (Rauch–Tung–Striebel smoother)
    if givenAll == 1:
        uAll_list = []
        VAll_list = []

        # Initialize with the last element
        uAll_prev = u[:, :, -1]  # (B, n_dims)
        VAll_prev = V[:, :, -1]  # (B, n_dims)

        # Insert at the beginning of the list
        uAll_list.insert(0, uAll_prev)
        VAll_list.insert(0, VAll_prev)

        for ir in range(n - 1):
            i = n - 2 - ir
            V_i = V[:, :, i]  # (B, n_dims)
            P_i = A * V_i * A + T  # (B, n_dims)

            J_i = V_i * A / P_i  # (B, n_dims)
            u_i = u[:, :, i]
            uAll_curr = u_i + J_i * (uAll_prev - A * u_i)
            VAll_curr = V_i + J_i * (VAll_prev - P_i) * J_i

            # Insert at the beginning of the list
            uAll_list.insert(0, uAll_curr)
            VAll_list.insert(0, VAll_curr)

            # Update previous values
            uAll_prev = uAll_curr
            VAll_prev = VAll_curr

        # Stack lists to form tensors
        uAll = torch.stack(uAll_list, dim=2)  # (B, n_dims, n)

        X = uAll
    else:
        X = u

    # Permute back to original shape
    sequence_new = X.permute(0, 2, 1)  # (B, n, n_dims)
    return sequence_new

def save_img(data, filename='image_with_colorbar.png', cmap='viridis'):
    # return
    print('called')
    if isinstance(data, torch.Tensor):
        data = data.detach().cpu().numpy()
    fig, ax = plt.subplots()
    cax = ax.imshow(data, cmap=cmap)
    fig.colorbar(cax)
    plt.savefig(filename, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)

def save_batch_images(data, folder='heatmaps', cmap='viridis'):
    # 创建保存热力图的文件夹
    if not os.path.exists(folder):
        os.makedirs(folder)

    if isinstance(data, torch.Tensor):
        data = data.detach().cpu().numpy()

    # 假设数据的形状为 [B, W, D]
    B, W, D = data.shape

    # 循环遍历每个批次，生成热力图
    for i in range(B):
        filename = os.path.join(folder, f'heatmap_{i}.png')
        fig, ax = plt.subplots()
        cax = ax.imshow(data[i], cmap=cmap)
        fig.colorbar(cax)
        plt.savefig(filename, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig)

def top_k_accuracy(logits, labels, ks=[1, 5]):
    """
    计算给定 logits 和 labels 的 top-k 准确率。
    
    参数：
    - logits: Tensor of shape (N, num_classes), 模型的预测输出
    - labels: Tensor of shape (N,), 包含实际类别的下标格式的标签
    - ks: List of integers, 指定要计算的 k 值，如 [1, 5] 表示计算 top-1 和 top-5 准确率
    
    返回：
    - acc_list: List of accuracies corresponding to each k in ks
    """
    max_k = max(ks)  # 计算需要的最大 k 值
    batch_size = labels.size(0)

    # 获取 top-k 预测，返回的 top_preds 形状为 (N, max_k)
    _, top_preds = logits.topk(max_k, dim=1, largest=True, sorted=True)

    # 扩展 labels 的形状以便与 top_preds 比较，形状为 (N, max_k)
    top_preds = top_preds.t()  # 转置为 (max_k, N)
    correct = top_preds.eq(labels.view(1, -1).expand_as(top_preds))  # 检查 top-k 内是否有正确的标签

    # 计算每个 k 对应的准确率
    acc_list = []
    for k in ks:
        correct_k = correct[:k].reshape(-1).float().sum(0, keepdim=True)  # 获取前 k 的正确预测数量
        accuracy = correct_k.mul_(100.0 / batch_size).item()  # 计算准确率百分比
        acc_list.append(accuracy)

    return acc_list

def video_order_load(n_vids=28):
    datapath = '/mnt/dataset0/qingzhu/AutoICA_Processed_EEG/Faced/Processed_data_filter_epoch_0.50_47_Auto_ICA_def_Threshold/After_remarks'
    filesPath = os.listdir(datapath)
    filesPath.sort()
    vid_orders = np.zeros((len(filesPath), n_vids),dtype=int)
    for idx, file in enumerate(filesPath):
        # Here don't forget to arange the subjects arrangement
        # print(file)
        remark_file = os.path.join(datapath,file,'After_remarks.mat')
        subject_remark = scipy.io.loadmat(remark_file)['After_remark']
        vid_orders[idx, :] = [np.squeeze(subject_remark[vid][0][2]) for vid in range(0,n_vids)]
    # print('vid_order shape: ', vid_orders.shape)
    return vid_orders

def reorder_vids_sepVideo(data, vid_play_order, sel_vid_inds, n_vids_all):
    # data: (n_subs, n_points, n_feas)
    n_vids = len(sel_vid_inds)
    n_subs = data.shape[0]
    # print('n_subs:',n_subs)
    vid_play_order_copy = vid_play_order.copy()
    vid_play_order_new = np.zeros((n_subs, len(sel_vid_inds))).astype(np.int32)
    data_reorder = np.zeros_like(data)
    if n_vids_all == 24:
        for sub in range(n_subs):
            tmp = vid_play_order_copy[sub,:]
            tmp = tmp[(tmp<13)|(tmp>16)]
            tmp[tmp>=17] = tmp[tmp>=17] - 4
            tmp = tmp - 1

            tmp_new = []
            for i in range(len(tmp)):
                if tmp[i] in sel_vid_inds:
                    tmp_new.append(np.where(sel_vid_inds==tmp[i])[0][0])
            tmp_new = np.array(tmp_new)

            vid_play_order_new[sub, :] = tmp_new

            data_sub = data[sub, :, :]
            data_sub = data_sub.reshape(n_vids, -1, data_sub.shape[-1])
            data_sub = data_sub[tmp_new, :, :]
            data_reorder[sub, :, :] = data_sub.reshape(-1, data_sub.shape[-1])
    elif n_vids_all == 28:
        for sub in range(n_subs):
            tmp = vid_play_order_copy[sub,:]
            tmp = tmp - 1

            tmp_new = []
            for i in range(len(tmp)):
                if tmp[i] in sel_vid_inds:
                    tmp_new.append(np.where(sel_vid_inds==tmp[i])[0][0])
            tmp_new = np.array(tmp_new)

            vid_play_order_new[sub, :] = tmp_new

            data_sub = data[sub, :, :]
            data_sub = data_sub.reshape(n_vids, -1, data_sub.shape[-1])
            data_sub = data_sub[tmp_new, :, :]
            data_reorder[sub, :, :] = data_sub.reshape(-1, data_sub.shape[-1])
        
    return data_reorder, vid_play_order_new


def reorder_vids_back(data, n_vids, vid_play_order_new):
    # data: (n_subs, n_points, n_feas)
    # return:(n_subs, n_points, n_feas)
    n_subs = data.shape[0]
    n_samples = data.shape[1]//n_vids


    data_back = np.zeros((n_subs, n_vids, n_samples, data.shape[-1]))

    for sub in range(n_subs):
        data_sub = data[sub, :, :].reshape(n_vids, n_samples, data.shape[-1])
        data_back[sub, vid_play_order_new[sub, :], :, :] = data_sub
    data_back = data_back.reshape(n_subs, n_vids*n_samples, data.shape[-1])
    return data_back

def save_tensor_or_ndarray(data, save_name, save_path="./saved_data"):
    """
    保存 Tensor 或 NumPy 数组到本地文件。

    参数:
        data: 输入的 Tensor 或 NumPy 数组。
        save_name: 保存的文件名（不带后缀）。
        save_path: 保存的路径，默认为 "./saved_data"。

    返回:
        None
    """
    # 如果路径不存在，创建路径
    os.makedirs(save_path, exist_ok=True)

    # 检查输入类型
    if isinstance(data, torch.Tensor):
        # 如果 Tensor 在 GPU 上，移动到 CPU
        if data.is_cuda:
            data = data.cpu().detach()
        # 转换为 NumPy 数组
        data = data.numpy()
    elif not isinstance(data, np.ndarray):
        raise ValueError("输入必须是 PyTorch Tensor 或 NumPy 数组")

    # 构建完整文件路径
    file_path = os.path.join(save_path, f"{save_name}.npy")

    # 保存为 .npy 文件
    np.save(file_path, data)
    print(f"数据已保存到: {file_path}")

import torch

def report_vram(tensor: torch.Tensor):
    print(tensor.shape)
    """返回张量占用的 VRAM 大小（以字节为单位）
    
    Args:
        tensor (torch.Tensor): 输入的 PyTorch 张量
        
    Returns:
        int: 张量占用的 VRAM 大小（字节）
    """
    if not tensor.is_cuda:
        print("Warning: Tensor is not on GPU. Returning CPU memory usage instead.")
    
    # 计算张量占用的总字节数 = 元素数量 × 每个元素的字节大小
    bytes = tensor.numel() * tensor.element_size()
    report = f"Tensor VRAM usage: {bytes / 1024 ** 2:.2f} MB"
    return report

TIME_FORMAT_STR: str = "%b_%d_%H_%M_%S"
def trace_handler(prof: torch.profiler.profile):
   # Prefix for file names.
   host_name = socket.gethostname()
   timestamp = datetime.now().strftime(TIME_FORMAT_STR)
   file_prefix = f"{host_name}_{timestamp}"

   # Construct the trace file.
#    prof.export_chrome_trace(f"{file_prefix}.json.gz")

   # Construct the memory timeline file.
   prof.export_memory_timeline(f"{file_prefix}.html", device="cuda:0")

def get_vram_profiler():
    prof = torch.profiler.profile(
       activities=[
            torch.profiler.ProfilerActivity.CPU,
            torch.profiler.ProfilerActivity.CUDA,
        ],
        schedule=torch.profiler.schedule(wait=0, warmup=0, active=6, repeat=1),
        record_shapes=True,
        profile_memory=True,
        with_stack=True,
        on_trace_ready=trace_handler,
        )
    return prof

def read_npy(dir):
    data = np.load(dir)
    return data, data.shape