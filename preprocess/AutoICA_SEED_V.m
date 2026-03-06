% This is the script for SEED_V EEG dataset preprocess, including:
% (1) Downsample
% (2) Band pass Filter
% (3) Divide into trials, specifically, data matrix for 1 subject 1 vedio
% (4) Detect and Interpolate bad channels
% (5) Auto ICA, remove artifacts, and line/eye/muscle noise components, based on matlab package EEGLab and Discover-eeg
% (6) Rereference
% (7) Reorder trials
% (8) Save data to .mat files

clear;
close all;
clc;

bpfreq = [0.5, 47];

data_dir = "....../EEG_raw";
files = dir(fullfile(data_dir, '*.cnt'));
n_sub = numel(files);

% n_sub = 1;


bad_prop_1 = 0.4;
thresh_1 = 3;
bad_prop_2 = 0.01;
thresh_2 = 30;

for sub_id = 1:n_sub
    %% Load Data
    file_name = fullfile(data_dir, files(sub_id).name);
    file_name = char(file_name);
    cfg = [];
    cfg.dataset = file_name;
    cfg.channel = {'all', '-M1', '-M2', '-VEO', '-HEO', '-CB1', '-CB2'};
    data = ft_preprocessing(cfg);

    %% Downsample
    resamplefs = 125;
    cfg = [];
    cfg.resamplefs = resamplefs;
    cfg.detrend = 'no';
    data_downsampled = ft_resampledata(cfg, data);
    pre_info.resamplefs = cfg.resamplefs;

    %% Filter
    cfg=[];
    cfg.bpfilter = 'yes';
    bpfreq = [0.5, 47];
    cfg.bpfreq = bpfreq;
%     cfg.continuous = 'yes';
    % The default bpfiltord is 4. Consider changing it to 3 if error occurs.
%     cfg.bpfiltord = 3; 
    pre_info.bpfreq = cfg.bpfreq;
    data_filted = ft_preprocessing(cfg, data_downsampled);

    %% Divide trial
    cfg = [];
    % the begin and end time in second of trials in all 3 sessions, modify accordingly.
    beg_s1 = [30, 132, 287, 555, 773, 982, 1271, 1628, 1730, 2025, 2227, 2435, 2667, 2932, 3204];
    end_s1 = [102, 228, 524, 742, 920, 1240, 1568, 1697, 1994, 2166, 2401, 2607, 2901, 3172, 3359];
    beg_s2 = [30, 299, 548, 646, 836, 1000, 1091, 1392, 1657, 1809, 1966, 2186, 2333, 2490, 2741];
    end_s2 = [267, 488, 614, 773, 967, 1059, 1331, 1622, 1777, 1908, 2153, 2302, 2428, 2709, 2817];
    beg_s3 = [30, 353, 478, 674, 825, 908, 1200, 1346, 1451, 1711, 2055, 2307, 2457, 2726, 2888];
    end_s3 = [321, 418, 643, 764, 877, 1147, 1284, 1418, 1679, 1996, 2275, 2425, 2664, 2857, 3066];

    session = regexp(files(sub_id).name, '_([^_]+)_', 'tokens', 'once');
    session = session{1};
    if session == '1'
        trl = [beg_s1;end_s1] * resamplefs;
    elseif session == '2'
        trl = [beg_s2;end_s2] * resamplefs;
    elseif session == '3'
        trl = [beg_s3;end_s3] * resamplefs;
    end
    fprintf("Session %s\n", session)
%     cfg.begsample = beg_s1;
%     cfg.endsample = end_s1;
% 
%     cfg.begsample = cfg.begsample * resamplefs;
%     cfg.endsample = cfg.endsample * resamplefs;
    trl = transpose(trl);
    trl = padarray(trl, [0,1], "post");
    cfg.trl = trl;
    data_trialed = ft_redefinetrial(cfg, data_filted);

    %% 1st Interpolate bad channels

    load chn_coords
    load('chn_coords')
    keep_channels = [1, 3, 6, 14];

    bad_channel_type_1 = cell(1, 15);
    bad_channel_type_2 = cell(1, 15);
    for i = 1: length(data_trialed.trial)
        visual_check_data(data_trialed.trial{i}, data_trialed.hdr.label, 50, 0, 'before interpolation', resamplefs);
%         pause(40);
        x = transpose(data_trialed.trial{i}); % change to (n_points, n_channs)
        [iBad_1, toGood_1] = nt_find_bad_channels_custom(x, bad_prop_1, thresh_1);
        [iBad_2, toGood_2] = nt_find_bad_channels_custom(x, bad_prop_2, thresh_2);
        bad_channel_type_1{i} = union(bad_channel_type_1{i}, setdiff(iBad_1, keep_channels));
        bad_channel_type_2{i} = union(bad_channel_type_2{i}, setdiff(iBad_2, keep_channels));
        iBad = union(iBad_1, iBad_2);
        % keep [1, 3, 6, 14] (Fp1, Fp2, F7, F8) at first interpolation
        iBad = setdiff(iBad, keep_channels);
        if iBad
            pre_info.iBad{i} = iBad;
            fprintf('trial %d, iBad:', i)
            disp(iBad)
            disp(data_trialed.hdr.label(iBad))
%             visual_check_data(data.trial{i}, hdr.label, 50, 0, ['trial ', num2str(i)], resamplefs)
        else
            fprintf('1st bad_channel inter - Trial %d No bad channels identified\n\n', i)
        end
        SEED_coords_matrix = load('SEED_coords_matrix.mat').coords_matrix;
        if iBad
            [toGood,fromGood] = nt_interpolate_bad_channels_custom(x, iBad, SEED_coords_matrix);
            x = x * (toGood * fromGood);
        end
        data_trialed.trial{i} = transpose(x);
    end

    %% Auto ICA
    EEGtemp = [];
    EEGtemp.filepath = sprintf('C:\Emotion Classification\prep_code_clPaper\Processed\subject_%d\trial_%d', sub_id, i);
    EEGtemp.filename = sprintf('subject_%d_trial_%d', sub_id, i);
    EEGtemp.data = cat(2, data_trialed.trial{:});
    %           visual_check_data(EEGtemp.data, hdr.label, 50, 0, 'Before-ICA', resamplefs);
    EEGtemp.etc = [];
    EEGtemp.setname = [];
    EEGtemp.icawinv = [];
    EEGtemp.icaweights = [];
    EEGtemp.icasphere = [];
    EEGtemp.nbchan = size(EEGtemp.data, 1);
    EEGtemp.pnts = size(EEGtemp.data, 2);
    EEGtemp.trials = 1;
    EEGtemp.srate = resamplefs;
    EEGtemp.xmin = 0;
    EEGtemp.xmax = (EEGtemp.pnts - 1) / EEGtemp.srate;
    EEGtemp.chanlocs = readlocs('C:\Emotion Classification\AutoICA\SEED_10_20_standard.ced');
    %               visual_check_data(EEGtemp.data, hdr.label, 50, 0, ['trial ', num2str(i)], resamplefs)
    EEGtemp = pop_runica(EEGtemp,'icatype','runica','concatcond','off');
    EEGtemp = pop_iclabel(EEGtemp,'default');
    debug = 1;
    if(debug)
        pop_eegplot(EEGtemp, 0);
        pop_viewprops(EEGtemp, 0);
    end
    %             IClabel: Brain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other.
    params_ICLabel_cus = [0 0;0.5 1; 0.7 1; 0.7 1; 0.7 1; 0.7 1; 0 0];
    params_ICLabel_def = [0 0;0.8 1; 0.8 1; 0 0; 0 0; 0 0; 0 0];
    % which threshold to use
%     params_IClabel = params_ICLabel_def;
%     thresh = 'Def';
    params_IClabel = params_ICLabel_def;
    thresh = 'Def';
    EEGtemp = pop_icflag(EEGtemp, params_IClabel);
    
    classifications = EEGtemp.etc.ic_classification.ICLabel.classifications; % Keep classifications before component substraction
    
    disp('Classified as noise components:')
    disp(find(EEGtemp.reject.gcompreject == 1));
    
    disp(EEGtemp.etc.ic_classification.ICLabel.classifications)
    EEGtemp = pop_subcomp(EEGtemp,[], debug); % Subtract artifactual independent components
    EEGtemp.etc.ic_classification.ICLabel.orig_classifications = classifications;
    
    % binary classification
    binary_matrix = zeros(size(classifications));
    for i = 1:size(classifications, 1)
        [max_values, max_indices] = max(classifications(i, :));
        if(max_values > params_IClabel(max_indices))
            binary_matrix(i, max_indices) = 1;
        end
    end
    classification_count = sum(binary_matrix, 1);

    %% Divide trial and 2nd bad channel interpolate

    %  Divide
    data_ica = data_trialed;
    start_time = 1;
    end_time = 0;
    for vid_id = 1:size(trl, 1)
        end_time = start_time + (trl(vid_id, 2) - trl(vid_id, 1));
        disp(start_time);
        disp(end_time);
        data_ica.trial{vid_id} = EEGtemp.data(:, start_time:end_time);
        start_time = end_time + 1;
    end

    %  bad channel interpolate

    for i = 1: length(data_ica.trial)
%         visual_check_data(data_trialed.trial{i}, data_trialed.hdr.label, 50, 0, 'before interpolation', resamplefs);
%         pause(40);
        x = transpose(data_ica.trial{i}); % change to (n_points, n_channs)
        [iBad_1, toGood_1] = nt_find_bad_channels_custom(x, bad_prop_1, thresh_1);
        [iBad_2, toGood_2] = nt_find_bad_channels_custom(x, bad_prop_2, thresh_2);
        bad_channel_type_1{i} = union(bad_channel_type_1{i}, iBad_1);
        bad_channel_type_2{i} = union(bad_channel_type_2{i}, iBad_2);
        iBad = union(iBad_1, iBad_2);
        if iBad
            pre_info.iBad{i} = iBad;
            fprintf('trial %d, iBad:', i)
            disp(iBad)
            disp(data_ica.hdr.label(iBad))
%             visual_check_data(data.trial{i}, hdr.label, 50, 0, ['trial ', num2str(i)], resamplefs)
        else
            fprintf('2st bad_channel inter - Trial %d No bad channels identified\n\n', i)
        end
        SEED_coords_matrix = load('SEED_coords_matrix.mat').coords_matrix;
        if iBad
            [toGood,fromGood] = nt_interpolate_bad_channels_custom(x, iBad, SEED_coords_matrix);
            x = x * (toGood * fromGood);
        end
        data_ica.trial{i} = transpose(x);
    end
    if(0)
        visual_check_data(data_ica.trial{i}, data_ica.hdr.label, 50, 0, 'After 2nd interpolation', resamplefs);
%         pause(40);
    end

    %% Rereference
    nomas_ind = 1:60;
    for i = 1: length(data_ica.trial)
        data_ica.trial{i} = data_ica.trial{i}(nomas_ind, :) - repmat(mean(data_ica.trial{i}(nomas_ind, :)), 60, 1);
    end

    %% Reorder trials
%     trial_reorder = data.trial;
%     trial_reorder(num_reorder) = data.trial;
%     for i = 1:28
%         n_samples(sub,i) = size(trial_reorder{1,i}, 2) / resamplefs;
%     end
%     disp(n_samples(sub,:))
%     pre_info.num_reorder = num_reorder;
    pre_info.classification = classification_count;
    pre_info.noise_sum = sum(pre_info.classification(2:6));
    pre_info.bad_channel_1 = bad_channel_type_1;
    pre_info.bad_channel_2 = bad_channel_type_2;
    if(debug)
        visual_check_data(data_ica.trial{5}, data_ica.hdr.label(nomas_ind), 50, 0, 'after process', resamplefs);
    end

    %% Save data
    
    n_samples_one = zeros(1, 15);
    data_all_cleaned = zeros(60, 0);
    for i = 1: length(data_ica.trial)
        data_all_cleaned = cat(2, data_all_cleaned, data_ica.trial{i});
        n_samples_one(1, i) = size(data_ica.trial{i}, 2) / resamplefs;
    end
    save_folder_name = sprintf('Processed_data_filter_%.2f_%.2f_AutoICA_%s_Threshold', bpfreq(1), bpfreq(2), thresh);
    save_dir = fullfile('C:\Emotion Classification\Processed\SEED_V', save_folder_name, 'data');
    if ~exist(save_dir)
        mkdir(save_dir);
    end
    if(files(sub_id).name(2) == '_')
        subject_id = files(sub_id).name(1);
    end
    if(files(sub_id).name(3) == '_')
        subject_id = files(sub_id).name([1, 2]);
    end
    save(fullfile(save_dir, sprintf('Sub_%s_Session_%s', subject_id, session)), 'data_all_cleaned', 'pre_info', 'n_samples_one');
end

