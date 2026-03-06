import numpy as np
import torch

def LDS(sequence):
    # shape: (timeSample, n_dims)  timesample 为一个vid的采样数
    # print(sequence.shape) # (26, 256)   

    sequence_new = np.zeros_like(sequence) # (26, 256)
    ave = np.mean(sequence, axis=0)

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
        P = np.zeros((m, m, n)) # (1, 1, 26)
        u = np.zeros((m, n)) # (1, 26)
        V = np.zeros((m, m, n)) # (1, 1, 26)
        K = np.zeros((m, m, n)) # (1, 1, 26)

        K[:, :, 0] = V0*C / (C*V0*C + sigma)
        u[:, 0] = u0 + K[:, :, 0] * (X[:, 0] - C*u0)
        V[:, :, 0] = (np.eye(m) - K[:, :, 0] * C) * V0

        for i in range(1, n):
            P[:, :, i-1] = A * V[:, :, i-1] * A + T
            K[:, :, i] = P[:, :, i-1] * C / (C * P[:, :, i-1] * C + sigma)
            u[:, i] = A * u[:, i-1] + K[:, :, i] * (X[:, i] - C*A*u[:, i-1])
            V[:, :, i] = (np.eye(m) - K[:, :, i] * C) * P[:, :, i-1]

        if givenAll == 1:
            uAll = np.zeros((m, n))
            VAll = np.zeros((m, m, n))
            J = np.zeros((m, m, n))
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

def LDS_acc(self,sequence):
    # print(sequence.shape) # (28*30, 256)
    # sequence_new = np.zeros_like(sequence) # (26, 256)
    ave = np.mean(sequence, axis=0) # [256,]
    u0 = ave
    X = sequence.T # [256, 26]
    V0 = 0.01
    A = 1
    T = 0.0001
    C = 1
    sigma = 1
    [m, n] = X.shape # (256,26)
    # P = np.zeros((m, n)) # (256,26) dia
    # u = np.zeros((m, n)) # (256,26)
    # V = np.zeros((m, n)) # (256,26) dia
    # K = np.zeros((m, n)) # (256,26)  todo: m,n diagonal of K
    # K[:, 0] = (V0 * C / (C * V0 * C + sigma)) * np.ones((m,))
    # u[:, 0] = u0 + K[:, 0] * (X[:, 0] - C * u0)
    # V[:, 0] = (np.ones((m,)) - K[:, 0] * C) * V0
    K = (V0*C / (C*V0*C + sigma))*np.ones((m,)).reshape(m,1)
    u = (u0 + K[:, 0] * (X[:, 0] - C*u0)).reshape(m,1)
    V = ((np.ones((m,)) - K[:, 0] * C) * V0).reshape(m,1)
    for i in range(1, n):
        # P[:, i - 1] = A * V[:, i - 1] * A + T
        # K[:, i] = P[:, i - 1] * C / (C * P[:, i - 1] * C + sigma)
        # u[:, i] = A * u[:, i - 1] + K[:, i] * (X[:, i] - C * A * u[:, i - 1])
        # V[:, i] = (np.ones((m,)) - K[:, i] * C) * P[:, i - 1]
        P_t = (A * V[ :, i-1] * A + T).reshape(m,1)
        if i == 1:
            P = P_t
        else:
            P = np.concatenate((P,P_t),-1)
        K_t = (P[:, i-1] * C / (C * P[:, i-1] * C + sigma)).reshape(m,1)
        K = np.concatenate((K,K_t),-1)
        u_t = (A * u[:, i-1] + K[:, i] * (X[:, i] - C*A*u[:, i-1])).reshape(m,1)
        u = np.concatenate((u,u_t),-1)
        V_t = ((np.ones((m,)) - K[:, i] * C) * P[:, i-1]).reshape(m,1)
        V = np.concatenate((V,V_t),-1)
    X = u
    return X.T

def LDS_gpu(sequence):
    """
    GPU-accelerated version of the LDS function using PyTorch.

    Args:
        sequence (torch.Tensor): Input tensor of shape (timeSample, n_dims).

    Returns:
        torch.Tensor: Smoothed sequence of the same shape as the input.
    """
    # Ensure sequence is a PyTorch tensor
    if not isinstance(sequence, torch.Tensor):
        sequence = torch.tensor(sequence)
    
    # Move sequence to GPU if available
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    sequence = sequence.to(device)
    sequence = sequence.float()  # Ensure floating point operations

    timeSample, n_dims = sequence.shape

    # Initialize output tensor
    sequence_new = torch.zeros_like(sequence)

    # Compute the mean over the time axis
    ave = torch.mean(sequence, dim=0)  # Shape: (n_dims,)

    # Transpose sequence to shape (n_dims, timeSample)
    X = sequence.transpose(0, 1)  # Shape: (n_dims, timeSample)

    # Initial state mean
    u0 = ave  # Shape: (n_dims,)

    # Define constants as tensors for GPU compatibility
    V0 = torch.tensor(0.01, device=device, dtype=sequence.dtype)
    A = torch.tensor(1.0, device=device, dtype=sequence.dtype)
    T = torch.tensor(0.0001, device=device, dtype=sequence.dtype)
    C = torch.tensor(1.0, device=device, dtype=sequence.dtype)
    sigma = torch.tensor(1.0, device=device, dtype=sequence.dtype)
    givenAll = 1  # Flag for backward smoothing

    n = timeSample  # Number of time steps

    # Initialize tensors for Kalman filter variables
    P = torch.zeros((n_dims, n), device=device, dtype=sequence.dtype)
    u = torch.zeros((n_dims, n), device=device, dtype=sequence.dtype)
    V = torch.zeros((n_dims, n), device=device, dtype=sequence.dtype)
    K = torch.zeros((n_dims, n), device=device, dtype=sequence.dtype)

    # Initial Kalman gain
    K_init = V0 * C / (C * V0 * C + sigma)  # Scalar tensor
    K[:, 0] = K_init  # Broadcast to all dimensions

    # Initial estimates
    u[:, 0] = u0 + K[:, 0] * (X[:, 0] - C * u0)
    V[:, 0] = (1 - K[:, 0] * C) * V0

    # Forward pass (Kalman Filter)
    for i in range(1, n):
        P_prev = A * V[:, i - 1] * A + T  # Shape: (n_dims,)
        Denominator = C * P_prev * C + sigma  # Shape: (n_dims,)
        K[:, i] = P_prev * C / Denominator  # Shape: (n_dims,)
        u[:, i] = A * u[:, i - 1] + K[:, i] * (X[:, i] - C * A * u[:, i - 1])
        V[:, i] = (1 - K[:, i] * C) * P_prev  # Shape: (n_dims,)

    # Backward pass (Rauch–Tung–Striebel smoother)
    if givenAll == 1:
        uAll = torch.zeros_like(u)
        VAll = torch.zeros_like(V)
        J = torch.zeros_like(K)

        # Initialize with the last element
        uAll[:, -1] = u[:, -1]
        VAll[:, -1] = V[:, -1]

        for ir in range(n - 1):
            i = n - 2 - ir
            P_i = A * V[:, i] * A + T  # Shape: (n_dims,)
            J[:, i] = V[:, i] * A / P_i  # Shape: (n_dims,)
            uAll[:, i] = u[:, i] + J[:, i] * (uAll[:, i + 1] - A * u[:, i])
            VAll[:, i] = V[:, i] + J[:, i] * (VAll[:, i + 1] - P_i) * J[:, i]

        X_smooth = uAll
    else:
        X_smooth = u

    # Transpose back to original shape (timeSample, n_dims)
    sequence_new = X_smooth.transpose(0, 1)
    sequence_new = sequence_new.cpu().numpy()

    return sequence_new


def running_norm(data,data_mean,data_var,decay_rate):
    # data  (subs,n_points,dim,...)
    # output data_norm:(subs,n_points,dim,...)

    data_norm = np.zeros_like(data)
    for sub in range(data.shape[0]):
        running_sum = np.zeros(data.shape[-1])
        running_square = np.zeros(data.shape[-1])
        decay_factor = 1
        for counter in range(data.shape[1]):
            data_one = data[sub, counter]
            running_sum = running_sum + data_one
            running_mean = running_sum / (counter+1)
            running_square = running_square + data_one**2
            running_var = (running_square - 2 * running_mean * running_sum) / (counter+1) + running_mean**2

            curr_mean = decay_factor*data_mean + (1-decay_factor)*running_mean
            curr_var = decay_factor*data_var + (1-decay_factor)*running_var
            decay_factor = decay_factor*decay_rate

            data_one = (data_one - curr_mean) / np.sqrt(curr_var + 1e-50)
            data_norm[sub, counter, :] = data_one
    return data_norm

def running_norm_onesubsession(data,data_mean,data_var,decay_rate):
    # data  (n_points,dim,...)
    # output data_norm:(n_points,dim,...)
    # one session represent one sub

    data_norm = np.zeros_like(data)

    running_sum = np.zeros(data.shape[-1])
    running_square = np.zeros(data.shape[-1])
    decay_factor = 1

    for counter in range(data.shape[0]):
        data_one = data[counter]
        running_sum = running_sum + data_one
        running_mean = running_sum / (counter+1)
        running_square = running_square + data_one**2
        running_var = (running_square - 2 * running_mean * running_sum) / (counter+1) + running_mean**2

        curr_mean = decay_factor*data_mean + (1-decay_factor)*running_mean
        curr_var = decay_factor*data_var + (1-decay_factor)*np.maximum(running_var,0)
        decay_factor = decay_factor*decay_rate
        data_one = (data_one - curr_mean) / np.sqrt(curr_var + 1e-10)
        
        data_norm[counter, :] = data_one
    return data_norm


