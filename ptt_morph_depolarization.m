% Preprocess data

clc; clear; close all;
addpath(genpath(pwd));

%% Load data
output_dir = "Dataset/Depolarization_Data/raw_data";
list_dir = dir(fullfile(output_dir, '*.mat'));
data_all = cell(1, numel(list_dir));
scene_all = strings(1, numel(list_dir));
sid_all = zeros(1, numel(list_dir));

for i_data = 1:numel(list_dir)
    file_path = fullfile(output_dir, list_dir(i_data).name);
    temp = load(file_path);
    data_all{i_data} = temp.data;

    % Parse scene and sid from filename
    split_dir = strsplit(list_dir(i_data).name(1:end-4), '_');
    scene_all(i_data) = split_dir{1};
    sid_all(i_data) = data_all{i_data}.sid;
end
%%
debug_mode = false;
if debug_mode
    NIR_G = struct;
    NIR_G.SBP = {};
    NIR_G.DBP = {};
    NIR_G.MBP = {};
    
    NIR_G.sid = {};
    NIR_G.rid = {};
    NIR_G.SNR_NIR = {};
    NIR_G.SNR_G = {};
    NIR_G.mean_NIR_cycles = {};
    NIR_G.mean_G_cycles = {};
    NIR_G.n_valid_NIR_cycles = {};
    NIR_G.n_valid_G_cycles = {};
    
    nir_fit_deg = 11;
    g_fit_deg = 7;
    
    NIR_G.NIR_coef = {};
    NIR_G.G_coef = {};
    NIR_G.G11_coef = {};
    NIR_G.G_coef_pad = {};
    NIR_G.NIR_coef_norm = {};
    NIR_G.G_coef_norm = {};
    NIR_G.G11_coef_norm = {};
    NIR_G.G_coef_pad_norm = {};
    
    
    norm_cycle_len = 100;
    fs = 125;
    % param for ssa
    lf = 0.5;
    hf = 6.0;
    [B, A] = butter(4, [lf*2/125, hf*2/fs], 'bandpass');
    
    ssa_deg = 1:4;
    ssa_seg_len = 10;
    ssa_seg_step = 2;

    for i_data=2
    % for i_data=1:length(data_all)
        nir_data = data_all{i_data}.trace1(2,:);
        g_data   = data_all{i_data}.trace3(2,:);
        fs = data_all{i_data}.fs;
        ts = data_all{i_data}.ts_bp;
        sbp = data_all{i_data}.sbp;
        dbp = data_all{i_data}.dbp;
        mbp = data_all{i_data}.mbp;
        sid = data_all{i_data}.sid;
    
        nir_data = -(nir_data/mean(nir_data)-1);
        g_data = -(g_data/mean(g_data)-1);
        
    
        det_param = [400, 4];
    
        % nir_data = ssa_overlap_merge(nir_data, fs, ssa_seg_len, ssa_seg_step, ssa_deg, det_param);
        % g_data = ssa_overlap_merge(g_data, fs, ssa_seg_len, ssa_seg_step, ssa_deg, det_param);
        [nir_data, g_data] = mssa_overlap_merge(nir_data, g_data, fs, ssa_seg_len, ssa_seg_step, ssa_deg, det_param);
    
        % Upsample (not sure if needed)
        nir_data = interp_resample(nir_data, fs, 100);
        g_data = interp_resample(g_data, fs, 100);
        ts = int32(ts/fs*100);
        fs = 100;
    
        [nir_peaks, nir_onsets] = findpeaks_pyampd(nir_data, fs);
        [g_peaks, g_onsets]     = findpeaks_pyampd(g_data, fs);
        g_median_cycle_len = median(diff(g_onsets)); % Used for cycle length outlier removal
    end

    figure;
    plot(nir_data); hold on; plot(nir_peaks, nir_data(nir_peaks), 'ro');
    nir_filt = filtfilt(B, A, nir_data);
    g_filt   = filtfilt(B, A ,g_data);
    figure;
    subplot(2,1,1);
    plot(nir_filt); hold on; plot(g_filt);
    subplot(2,1,2);
    plot(nir_data); hold on; plot(g_data);

    % 1. Envelope Outlier Removal
    [anomaly_idx, env_diff, upper_env, lower_env, nir_good_peaks, nir_good_onsets] = detect_motion_artifacts(nir_data, fs, 1.75, nir_peaks, nir_onsets);
    [~, ~, ~, ~, g_good_peaks, g_good_onsets] = detect_motion_artifacts(g_data, fs, 1.75, g_peaks, g_onsets);
    
    % 2. Cycle length and peak location outlier removal
    nir_triplets = peak_onset_tripletize(nir_good_peaks, nir_good_onsets, g_median_cycle_len, 0.25, true);
    g_triplets = peak_onset_tripletize(g_good_peaks, g_good_onsets, g_median_cycle_len, 0.25, true);
    
    % 3. Correlation (Template matching) Outlier Removal
    nir_cycles = get_trip_norm_cycles(nir_data, nir_triplets, norm_cycle_len);
    g_cycles = get_trip_norm_cycles(g_data, g_triplets, norm_cycle_len);
    
    [nir_median_cycle, selected_nir_ind] = high_corr_median_cycle(nir_cycles);
    [g_median_cycle, selected_g_ind] = high_corr_median_cycle(g_cycles); 
    
    nir_triplets = nir_triplets(selected_nir_ind, :);
    g_triplets = g_triplets(selected_g_ind, :);


    figure;
    subplot(2,1,1);
    plot(nir_cycles');
    title('NIR');
    subplot(2,1,2);
    plot(g_cycles');
    title('G');
end
% %%
% n_ts = length(ts);
%     bp_measurement_range = int32(15*fs);
%     for i_ts = 15
%         start_point = max(ts(i_ts) - bp_measurement_range, 1);
%         end_point = min(ts(i_ts) + bp_measurement_range, length(nir_data));
%         nir_trip_seg = nir_triplets(nir_triplets(:,1)>start_point & nir_triplets(:,2)<end_point,:);
%         g_trip_seg = g_triplets(g_triplets(:,1)>start_point & g_triplets(:,2)<end_point,:);
% 
%         if isempty(nir_trip_seg) | isempty(g_trip_seg)
%             fprintf("Empty triplets for sid %d and rid %d\n", sid, i_ts);
%             continue;
%         end
% 
%         % Normalized cycles
%         nir_cycles_seg = get_trip_norm_cycles(nir_data, nir_trip_seg, norm_cycle_len);
%         g_cycles_seg = get_trip_norm_cycles(g_data, g_trip_seg, norm_cycle_len);
% 
%         % k=1;
%         % if size(nir_cycles_seg,1) > k
%         %     nir_cycles_seg = pca_recon(nir_cycles_seg, k, false);
%         % end
%         % if size(g_cycles_seg,1) > k
%         %     g_cycles_seg   = pca_recon(g_cycles_seg, k, false);
%         % end
% 
%         NIR_G.n_valid_NIR_cycles{end+1} = size(nir_cycles_seg, 1);
%         NIR_G.n_valid_G_cycles{end+1} = size(g_cycles_seg, 1);
% 
%         x = linspace(0, 1, norm_cycle_len);
%         all_p_NIR = [];
%         all_p_NIR_norm = [];
%         all_p_G = [];
%         all_p_G11 = [];
%         all_p_G_norm = [];
%         all_p_G11_norm = [];
%         for i=1:size(nir_cycles_seg,1)
%             p_NIR = polyfit(x, nir_cycles_seg(i,:), nir_fit_deg);
%             p_NIR_norm = p_NIR/norm(p_NIR);
%             all_p_NIR = [all_p_NIR; p_NIR];
%             all_p_NIR_norm = [all_p_NIR_norm; p_NIR_norm];
%         end
%         for i=1:size(g_cycles_seg,1)
%             p_G = polyfit(x, g_cycles_seg(i,:), g_fit_deg);
%             p_G11= polyfit(x, g_cycles_seg(i,:), nir_fit_deg);
% 
%             p_G_norm = p_G/norm(p_G);
%             p_G11_norm = p_G11/norm(p_G11);
% 
%             all_p_G = [all_p_G; p_G];
%             all_p_G_norm = [all_p_G_norm; p_G_norm];
%             all_p_G11 = [all_p_G11; p_G11];
%             all_p_G11_norm = [all_p_G11_norm; p_G11_norm];
%         end
% 
% 
%         NIR_G.NIR_coef{end+1}        = median(all_p_NIR, 1);
%         NIR_G.G_coef{end+1}          = median(all_p_G, 1);
%         NIR_G.G11_coef{end+1}        = median(all_p_G11, 1);
%         NIR_G.NIR_coef_norm{end+1}   = median(all_p_NIR_norm, 1);
%         NIR_G.G_coef_norm{end+1}     = median(all_p_G_norm, 1);
%         NIR_G.G11_coef_norm{end+1}   = median(all_p_G11_norm, 1);
% 
%         NIR_G.G_coef_pad{end+1}      = [zeros(1, nir_fit_deg - g_fit_deg), median(all_p_G, 1)];
%         NIR_G.G_coef_pad_norm{end+1} = [zeros(1, nir_fit_deg - g_fit_deg), median(all_p_G_norm, 1)];
% 
%         % Save cycles
%         NIR_G.sid{end+1} = sid;
%         NIR_G.rid{end+1} = i_ts;
%         NIR_G.SBP{end+1} = sbp(i_ts);
%         NIR_G.DBP{end+1} = dbp(i_ts);
%         NIR_G.MBP{end+1} = mbp(i_ts);
%     end
% 
% %%
% figure;
% subplot(2,1,1);
% plot(nir_cycles_seg');
% subplot(2,1,2);
% plot(g_cycles_seg');

%% Preprocess
colors = {'red','green','blue'};
all_NIR_G = {};
for color_ch = 1
    % for i_data = 1
    NIR_G = struct;
    NIR_G.SBP = {};
    NIR_G.DBP = {};
    NIR_G.MBP = {};
    
    NIR_G.sid = {};
    NIR_G.rid = {};
    % NIR_G.SNR_NIR = {};
    % NIR_G.SNR_G = {};
    % NIR_G.mean_NIR_cycles = {};
    % NIR_G.mean_G_cycles = {};
    NIR_G.n_valid_NIR_cycles = {};
    NIR_G.n_valid_G_cycles = {};
    NIR_G.NIR_coef = {};
    NIR_G.G_coef = {};
    NIR_G.G11_coef = {};
    NIR_G.G_coef_pad = {};
    NIR_G.NIR_coef_norm = {};
    NIR_G.G_coef_norm = {};
    NIR_G.G11_coef_norm = {};
    NIR_G.G_coef_pad_norm = {};
    NIR_G.ptt = {};
    for i_data=1:length(data_all)
        nir_data = data_all{i_data}.trace1(1,:);
        g_data   = data_all{i_data}.trace1(2,:);
        sid = data_all{i_data}.sid;
        if sid<11
            continue;
        end
        NIR_G = get_coef(data_all, nir_data, g_data, i_data, NIR_G);
    end
    all_NIR_G{end+1} = NIR_G;
end

function NIR_G = get_coef(data_all, nir_data, g_data, i_data, NIR_G)

    norm_cycle_len = 100;
    
    % param for ssa
    ssa_deg = 1:4;
    ssa_seg_len = 10;
    ssa_seg_step = 2;
    nir_fit_deg = 11;
    g_fit_deg = 7;

    fs = data_all{i_data}.fs;
    ts = data_all{i_data}.ts_bp;
    sbp = data_all{i_data}.sbp;
    dbp = data_all{i_data}.dbp;
    ptt = data_all{i_data}.ptt;
    % sbp = sbp - mean(sbp);
    % dbp = dbp - mean(dbp);
    % mbp = data_all{i_data}.mbp;
    sid = data_all{i_data}.sid;
    scene = data_all{i_data}.scene;
    fprintf("Working on subject %d\n", sid);

    nir_data = -(nir_data/mean(nir_data)-1);
    g_data = -(g_data/mean(g_data)-1);
    det_param = [400, 4];

    nir_data = ssa_overlap_merge(nir_data, fs, ssa_seg_len, ssa_seg_step, ssa_deg, det_param);
    g_data = ssa_overlap_merge(g_data, fs, ssa_seg_len, ssa_seg_step, ssa_deg, det_param);
    % [nir_data, g_data] = mssa_overlap_merge(nir_data, g_data, fs, ssa_seg_len, ssa_seg_step, ssa_deg, det_param);

    % Upsample (not sure if needed)
    nir_data = interp_resample(nir_data, fs, 100);
    g_data = interp_resample(g_data, fs, 100);
    ts = int32(ts/fs*100);
    fs = 100;

    [nir_peaks, nir_onsets] = findpeaks_pyampd(nir_data, fs);
    [g_peaks, g_onsets]     = findpeaks_pyampd(g_data, fs);
    g_median_cycle_len = median(diff(g_onsets)); % Used for cycle length outlier removal

    % 1. Envelope Outlier Removal
    [anomaly_idx, env_diff, upper_env, lower_env, nir_good_peaks, nir_good_onsets] = detect_motion_artifacts(nir_data, fs, 1.75, nir_peaks, nir_onsets);
    [~, ~, ~, ~, g_good_peaks, g_good_onsets] = detect_motion_artifacts(g_data, fs, 1.75, g_peaks, g_onsets);

    % 2. Cycle length and peak location outlier removal
    nir_triplets = peak_onset_tripletize(nir_good_peaks, nir_good_onsets, g_median_cycle_len, 0.25, true);
    g_triplets = peak_onset_tripletize(g_good_peaks, g_good_onsets, g_median_cycle_len, 0.25, true);

    % 3. Correlation (Template matching) Outlier Removal
    nir_cycles = get_trip_norm_cycles(nir_data, nir_triplets, norm_cycle_len);
    g_cycles = get_trip_norm_cycles(g_data, g_triplets, norm_cycle_len);

    [nir_median_cycle, selected_nir_ind] = high_corr_median_cycle(nir_cycles);
    [g_median_cycle, selected_g_ind] = high_corr_median_cycle(g_cycles); 

    nir_triplets = nir_triplets(selected_nir_ind, :);
    g_triplets = g_triplets(selected_g_ind, :);

    n_ts = length(ts);
    n_ptt = size(ptt,2);
    n_ts = min(n_ts, n_ptt);
    bp_measurement_range = int32(10*fs);
    if strcmp(scene, 'icu')
        ts_selection = 1:30:n_ts;
    else
        ts_selection = 1:n_ts;
    end
    for i_ts = ts_selection
        start_point = max(ts(i_ts) - bp_measurement_range, 1);
        end_point = min(ts(i_ts) + bp_measurement_range, length(nir_data));
        start_point_sec = start_point/fs;
        end_point_sec = end_point/fs;
        nir_trip_seg = nir_triplets(nir_triplets(:,1)>start_point & nir_triplets(:,2)<end_point,:);
        g_trip_seg = g_triplets(g_triplets(:,1)>start_point & g_triplets(:,2)<end_point,:);

        if isempty(nir_trip_seg) | isempty(g_trip_seg)
            fprintf("Empty triplets for sid %d and rid %d\n", sid, i_ts);
            continue;
        end

        % Normalized cycles
        nir_cycles_seg = get_trip_norm_cycles(nir_data, nir_trip_seg, norm_cycle_len);
        g_cycles_seg = get_trip_norm_cycles(g_data, g_trip_seg, norm_cycle_len);

        % k=1;
        % if size(nir_cycles_seg,1) > k
        %     nir_cycles_seg = pca_recon(nir_cycles_seg, k, false);
        % end
        % if size(g_cycles_seg,1) > k
        %     g_cycles_seg   = pca_recon(g_cycles_seg, k, false);
        % end

        NIR_G.n_valid_NIR_cycles{end+1} = size(nir_cycles_seg, 1);
        NIR_G.n_valid_G_cycles{end+1} = size(g_cycles_seg, 1);

        x = linspace(0, 1, norm_cycle_len);
        all_p_NIR = [];
        all_p_NIR_norm = [];
        all_p_G = [];
        all_p_G11 = [];
        all_p_G_norm = [];
        all_p_G11_norm = [];
        for i=1:size(nir_cycles_seg,1)
            p_NIR = polyfit(x, nir_cycles_seg(i,:), nir_fit_deg);
            p_NIR_norm = p_NIR/norm(p_NIR);
            all_p_NIR = [all_p_NIR; p_NIR];
            all_p_NIR_norm = [all_p_NIR_norm; p_NIR_norm];
        end
        for i=1:size(g_cycles_seg,1)
            p_G = polyfit(x, g_cycles_seg(i,:), g_fit_deg);
            p_G11= polyfit(x, g_cycles_seg(i,:), nir_fit_deg);
        
            p_G_norm = p_G/norm(p_G);
            p_G11_norm = p_G11/norm(p_G11);

            all_p_G = [all_p_G; p_G];
            all_p_G_norm = [all_p_G_norm; p_G_norm];
            all_p_G11 = [all_p_G11; p_G11];
            all_p_G11_norm = [all_p_G11_norm; p_G11_norm];
        end


        NIR_G.NIR_coef{end+1}        = median(all_p_NIR, 1);
        NIR_G.G_coef{end+1}          = median(all_p_G, 1);
        NIR_G.G11_coef{end+1}        = median(all_p_G11, 1);
        NIR_G.NIR_coef_norm{end+1}   = median(all_p_NIR_norm, 1);
        NIR_G.G_coef_norm{end+1}     = median(all_p_G_norm, 1);
        NIR_G.G11_coef_norm{end+1}   = median(all_p_G11_norm, 1);
        
        NIR_G.G_coef_pad{end+1}      = [zeros(1, nir_fit_deg - g_fit_deg), median(all_p_G, 1)];
        NIR_G.G_coef_pad_norm{end+1} = [zeros(1, nir_fit_deg - g_fit_deg), median(all_p_G_norm, 1)];

        % Save cycles
        NIR_G.sid{end+1} = sid;
        NIR_G.rid{end+1} = i_ts;
        NIR_G.ptt{end+1} = mean(ptt(:,max(i_ts-15,1):min(i_ts+15,n_ts)),2);
        if strcmp(scene, 'icu')
            sbp_temp = sbp(max(1,i_ts-15):min(i_ts+15,n_ts));
            NIR_G.SBP{end+1} = median(sbp_temp);
            dbp_temp = dbp(max(1,i_ts-15):min(i_ts+15,n_ts));
            NIR_G.DBP{end+1} = median(dbp_temp);
        else
            NIR_G.SBP{end+1} = sbp(i_ts);
            NIR_G.DBP{end+1} = dbp(i_ts);
        end
    end
end
%%
% save('all_NIR_G_with_ptt_raw_bp.mat', 'all_NIR_G');
save('all_NIR_G_red_green.mat', 'all_NIR_G');
% load('all_NIR_G2.mat', 'all_NIR_G');
% %% test PCA for waveform extraction
% recon_g_cycles_seg = pca_recon(g_cycles_seg, 1, false);
%% Find the common trial and fusion of all RGB channels
% --- Step 1: Build key map ---
keys_all = cell(1, 3);
for c = 1:3
    sid = all_NIR_G{c}.sid;
    rid = all_NIR_G{c}.rid;
    keys_all{c} = strcat(string(sid), "_", string(rid));  % e.g., "S01_15"
end

% --- Step 2: Find common keys ---
common_keys = keys_all{1};
for c = 2:3
    common_keys = intersect(common_keys, keys_all{c}, 'stable');
end

% --- Step 3: Concatenate features across channels for all properties ---
fields_to_merge = fieldnames(all_NIR_G{1});
skip_fields = {'sid', 'rid', 'SBP', 'DBP', 'MBP', 'SNR_NIR', 'SNR_G', 'mean_NIR_cycles' ,'mean_G_cycles'};  % Don't merge these

concat_NIR_G = struct;
concat_NIR_G.sid = {};
concat_NIR_G.rid = {};
concat_NIR_G.SBP = {};
concat_NIR_G.DBP = {};

for f = 1:length(fields_to_merge)
    field = fields_to_merge{f};
    if any(strcmp(field, skip_fields)), continue; end
    concat_NIR_G.(field) = {};
end

for i = 1:length(common_keys)
    key = common_keys{i};
    sid_rid = split(key, "_");
    sid = sid_rid{1};
    rid = str2double(sid_rid{2});
    concat_NIR_G.sid{end+1} = str2double(sid);
    concat_NIR_G.rid{end+1} = rid;

    % Use SBP/DBP from channel 1
    idx1 = find(strcmp(keys_all{1}, key));
    concat_NIR_G.SBP{end+1} = all_NIR_G{1}.SBP{idx1};
    concat_NIR_G.DBP{end+1} = all_NIR_G{1}.DBP{idx1};

    for f = 1:length(fields_to_merge)
        field = fields_to_merge{f};
        if any(strcmp(field, skip_fields)), continue; end

        merged_feat = [];
        for c = 1:3
            keys = keys_all{c};
            idx = find(strcmp(keys, key));
            value = all_NIR_G{c}.(field){idx};
            merged_feat = [merged_feat, value];  % concat across channels
        end
        concat_NIR_G.(field){end+1} = merged_feat;
    end
end
%%
reg_bp(all_NIR_G{1}, 1, 0, {'SVR'}, 'Red only', false);
%%
reg_bp(all_NIR_G{1}, 1, 0, {'SVR'}, 'Red w/ PTT', true);
%%
reg_bp(all_NIR_G{2}, 1, 0, {'SVR'}, 'Green only', false);
%%
reg_bp(all_NIR_G{2}, 1, 0, {'SVR'}, 'Green w/ PTT', true);
%%
reg_bp(all_NIR_G{3}, 1, 0, {'SVR'}, 'Blue only', false);
%%
reg_bp(all_NIR_G{3}, 1, 0, {'SVR'}, 'Blue w/ PTT', true);
%%
reg_bp(concat_NIR_G, 1, 0, {'RF'}, 'concated RGB');
%% Cross Merge Red and Green channel (R0 + G90 OR G0 + R90)
% --- Step 1: Build key map ---
keys_all = cell(1, 2);
for c = 1:2
    sid = all_NIR_G{c}.sid;
    rid = all_NIR_G{c}.rid;
    keys_all{c} = strcat(string(sid), "_", string(rid));  % e.g., "S01_15"
end

% --- Step 2: Find common keys ---
common_keys = keys_all{1};
for c = 2
    common_keys = intersect(common_keys, keys_all{c}, 'stable');
end

% --- Step 3: Concatenate features across R0 and G90 ---
fields_to_merge = fieldnames(all_NIR_G{1});
skip_fields = {'sid', 'rid', 'SBP', 'DBP', 'MBP', 'SNR_NIR', 'SNR_G', 'mean_NIR_cycles' ,'mean_G_cycles'};  % Don't merge these

R0_G90 = struct;
R0_G90.sid = {};
R0_G90.rid = {};
R0_G90.SBP = {};
R0_G90.DBP = {};

for f = 1:length(fields_to_merge)
    field = fields_to_merge{f};
    if any(strcmp(field, skip_fields)), continue; end
    R0_G90.(field) = {};
end

for i = 1:length(common_keys)
    key = common_keys{i};
    sid_rid = split(key, "_");
    sid = sid_rid{1};
    rid = str2double(sid_rid{2});
    R0_G90.sid{end+1} = str2double(sid);
    R0_G90.rid{end+1} = rid;

    % Use SBP/DBP from channel 1
    idx1 = find(strcmp(keys_all{1}, key));
    R0_G90.SBP{end+1} = all_NIR_G{1}.SBP{idx1};
    R0_G90.DBP{end+1} = all_NIR_G{1}.DBP{idx1};

    for f = 1:length(fields_to_merge)
        field = fields_to_merge{f};
        if any(strcmp(field, skip_fields)), continue; end
        if contains(field, 'NIR')
            idx_R = find(strcmp(keys_all{1}, key));
            R0_G90.(field){end+1} = all_NIR_G{1}.(field){idx_R}; 
        else
            idx_G = find(strcmp(keys_all{2}, key));
            R0_G90.(field){end+1} = all_NIR_G{2}.(field){idx_G};
        end
    end
end

% --- Step 4: Apply the same operation to concate R90 and G0 ---
fields_to_merge = fieldnames(all_NIR_G{1});
skip_fields = {'sid', 'rid', 'SBP', 'DBP', 'MBP', 'SNR_NIR', 'SNR_G', 'mean_NIR_cycles' ,'mean_G_cycles'};  % Don't merge these

R90_G0 = struct;
R90_G0.sid = {};
R90_G0.rid = {};
R90_G0.SBP = {};
R90_G0.DBP = {};

for f = 1:length(fields_to_merge)
    field = fields_to_merge{f};
    if any(strcmp(field, skip_fields)), continue; end
    R90_G0.(field) = {};
end

for i = 1:length(common_keys)
    key = common_keys{i};
    sid_rid = split(key, "_");
    sid = sid_rid{1};
    rid = str2double(sid_rid{2});
    R90_G0.sid{end+1} = str2double(sid);
    R90_G0.rid{end+1} = rid;

    % Use SBP/DBP from channel 1
    idx1 = find(strcmp(keys_all{1}, key));
    R90_G0.SBP{end+1} = all_NIR_G{1}.SBP{idx1};
    R90_G0.DBP{end+1} = all_NIR_G{1}.DBP{idx1};

    for f = 1:length(fields_to_merge)
        field = fields_to_merge{f};
        if any(strcmp(field, skip_fields)), continue; end
        if contains(field, 'NIR')
            idx_G = find(strcmp(keys_all{2}, key));
            R90_G0.(field){end+1} = all_NIR_G{2}.(field){idx_G}; 
        else
            idx_R = find(strcmp(keys_all{1}, key));
            R90_G0.(field){end+1} = all_NIR_G{1}.(field){idx_R};
        end
    end
end
%%
reg_bp(R0_G90, 1, 0, {'SVR'}, 'concated R0 G90');
%%
reg_bp(R90_G0, 1, 0, {'SVR'}, 'concated R90 G0');
%%
% select_sid = [1:10, 12, 14, 15, 16, 17, 18, 21, 23, 27, 29];
colors = {'red','green','blue'};
load('all_NIR_G2.mat');
mdl = 'RF';
select_sid = 11:30;
for i_ch = 1
    NIR_G = all_NIR_G{i_ch};
    all_sids = cell2mat(NIR_G.sid);
    good_sid_idx = ismember(all_sids, select_sid); 
    NIR_G.SBP = NIR_G.SBP(good_sid_idx);
    NIR_G.DBP = NIR_G.DBP(good_sid_idx);
    NIR_G.sid = NIR_G.sid(good_sid_idx);
    NIR_G.rid = NIR_G.rid(good_sid_idx);
    NIR_G.n_valid_NIR_cycles = NIR_G.n_valid_NIR_cycles(good_sid_idx);
    NIR_G.n_valid_G_cycles = NIR_G.n_valid_G_cycles(good_sid_idx);
    NIR_G.NIR_coef = NIR_G.NIR_coef(good_sid_idx);
    NIR_G.G_coef = NIR_G.G_coef(good_sid_idx);
    NIR_G.G11_coef = NIR_G.G11_coef(good_sid_idx);
    NIR_G.G_coef_pad = NIR_G.G_coef_pad(good_sid_idx);
    NIR_G.NIR_coef_norm = NIR_G.NIR_coef_norm(good_sid_idx);
    NIR_G.G_coef_norm = NIR_G.G_coef_norm(good_sid_idx);
    NIR_G.G11_coef_norm = NIR_G.G11_coef_norm(good_sid_idx);
    NIR_G.G_coef_pad_norm = NIR_G.G_coef_pad_norm(good_sid_idx);
    
    [pd_bp, gt_bp] = reg_bp(NIR_G, 1, 0, {mdl}, [colors{i_ch}, ' ', mdl]);
end


%%
figure;
plot(pd_bp{2,3}); hold on;
plot(gt_bp{2,3});
%% Polyfit coefficients seem noisy, maybe need PCA
k=1;
if size(nir_cycles_seg,1) > k
    nir_cycles_seg_pca = pca_recon(nir_cycles_seg, k, false);
end
if size(g_cycles_seg,1) > k
    g_cycles_seg_pca   = pca_recon(g_cycles_seg, k, false);
end
%%
figure; 
subplot(2,2,1); plot(nir_cycles_seg');
title('trace1')
subplot(2,2,2); plot(g_cycles_seg');
title('trace3')
subplot(2,2,3); plot(nir_cycles_seg_pca');
ylabel('PCA reduction');
subplot(2,2,4); plot(g_cycles_seg_pca');
%%
function X_recon = pca_recon(cycles, k, plot_flag)
if nargin<2
    k=min(4, size(cycles,1));
end
if nargin<3
    plot_flag=false;
end
[coeff, score, ~, ~, ~] = pca(cycles);
X_recon = score(:, 1:k) * coeff(:, 1:k)';  % Low-rank approximation
X_recon = X_recon + mean(cycles, 1);      % Add back the mean

% Plot some original and reconstructed cycles side-by-side
if plot_flag
    n_plot = min(5, size(cycles, 1));
    figure;
    for i = 1:n_plot
        subplot(n_plot, 1, i);
        plot(cycles(i,:), 'k-', 'DisplayName', 'Original'); hold on;
        plot(X_recon(i,:), 'r--', 'DisplayName', sprintf('Reconstructed (k=%d)', k));
        legend;
        xlabel('Time'); ylabel('Amplitude');
        title(sprintf('Cycle %d', i));
        grid on;
    end
    sgtitle(sprintf('Original vs Reconstructed Cycles Using Top %d PCs', k));
end
end
%%
% function clean_cycle = get_clean_cycle_by_pca(cycles, k)
% % cycles: [n x L] matrix of normalized cycles (each row is a cycle)
% % clean_cycle: [1 x L] - principal component direction projected back
% 
% % Center the data (PCA expects zero-mean rows)
% X = cycles - mean(cycles, 1);
% 
% % Perform PCA using SVD
% [U, S, V] = svd(X, 'econ');  % X = U S V^T
% 
% % First principal component direction (in original feature space)
% pc1 = V(:,1:k)';           % [1 x L]
% coeffs = U(:,1:k) * S(1:k,1:k);  % projections on PC1
% 
% % Project mean cycle onto PC1 and scale back
% clean_cycle = mean(cycles, 1) + coeffs(2,:) * pc1;
% end
% 
% k_list = [1, 2, 3, 5, 10];  % Try different numbers of PCs
% colors = lines(length(k_list));  % Distinct colors for plotting
% 
% % Plot all cycles for context
% figure;
% plot(nir_cycles_seg(2,:)', 'Color', [0.8, 0.8, 0.8]); hold on;
% 
% % Plot mean for comparison
% plot(mean(nir_cycles_seg, 1), 'b-', 'LineWidth', 2, 'DisplayName', 'Mean');
% 
% % Plot PCA-reconstructed cycles with different k
% for i = 1:length(k_list)
%     k = k_list(i);
%     clean_cycle = get_clean_cycle_by_pca(nir_cycles_seg, k);
%     plot(clean_cycle, '--', 'Color', colors(i,:), 'LineWidth', 2, ...
%         'DisplayName', sprintf('PCA-%d', k));
% end
% 
% legend('show');
% xlabel('Time (normalized)');
% ylabel('Amplitude');
% title('Clean Cycle Reconstruction with Varying Number of PCs');
% grid on;

%% M-SSA
i_data = 30;
nir_data = data_all{i_data}.nir_raw;
g_data   = data_all{i_data}.g_raw;
fs = data_all{i_data}.fs;
sid = data_all{i_data}.sid;
sbp = data_all{i_data}.sbp;
dbp = data_all{i_data}.dbp;
L = fs*10;


nir_data = -(nir_data/mean(nir_data)-1);
g_data = -(g_data/mean(g_data)-1);

if sid < 21
    det_param = [100, 1];
    % det_param = 100;
else
    det_param = [200, 2];
    % det_param = 200;
end

% Origianl SSA as baseline
ssa_nir_data = ssa_overlap_merge(nir_data, fs, ssa_seg_len, ssa_seg_step, ssa_deg, det_param);
ssa_g_data = ssa_overlap_merge(g_data, fs, ssa_seg_len, ssa_seg_step, ssa_deg, det_param);

%%
nir_data_seg = nir_data(1:L);
g_data_seg   = g_data(1:L);

A1 = tar_detrend(L, det_param(1)); % Remove low-freq only
A2 = tar_detrend(L, det_param(2)); % Remove low-freq only
A = A1 - A2;

% nir_data_seg_lp = nir_data_seg*A1;
% g_data_seg_lp   = g_data_seg*A1;
% comb_data_seg_lp = [nir_data_seg_lp', g_data_seg_lp'];
% Y_lp = mssa(comb_data_seg_lp, 1*fs);
nir_data_seg = nir_data_seg*A;
g_data_seg   = g_data_seg*A;
comb_data_seg = [nir_data_seg', g_data_seg'];

[Y, A, EV, EW] = mssa(comb_data_seg, 1*fs);


%%
figure;
for i_plot = 1:4
    subplot(4,2,(i_plot-1)*2+2);
    plot((1:size(Y,1))/fs, sum(Y(:,1,1:2*i_plot), 3)); hold on;
    plot((1:size(Y,1))/fs, sum(Y(:,2,1:2*i_plot), 3));
    if i_plot==1
        legend({'NIR', 'G'});
        title('MSSA recon')
    end
    ylabel(sprintf('RC %d', i_plot));

    subplot(4,2,(i_plot-1)*2+1);
    ssa_nir_data_seg = ssa(nir_data_seg, fs, 1:2*i_plot, false);
    ssa_g_data_seg   = ssa(g_data_seg, fs, 1:2*i_plot, false);
    plot((1:size(Y,1))/fs, ssa_nir_data_seg); hold on;
    plot((1:size(Y,1))/fs, ssa_g_data_seg);

    if i_plot==1
        legend({'NIR', 'G'});
        title('SSA recon')
    end
end
sgtitle(sprintf('Sub%d, SBP%.f, DBP%.f', sid, mean(sbp), mean(dbp)))

%%
% plot(nir_data(1:512)); hold on;
% plot(g_data(1:512))
% periodogram(nir_data(1:512), rectwin(512), 512, fs); hold on
% periodogram(g_data(1:512), rectwin(512), 512, fs)
% legend({'NIR data', 'G data'})
% title('Periodogram')
%%
function [pd_bp, gt_bp]=reg_bp(NIR_G, channel_norms, same_segs, reg_mdls, title_name, use_ptt)
    results = [];   % will become a MATLAB table
    for channel_norm = channel_norms
        for same_deg = same_segs
            if same_deg
                if channel_norm
                    G_mat = cell2mat(NIR_G.G11_coef_norm'); 
                    NIR_mat   = cell2mat(NIR_G.NIR_coef_norm');           % NIR-only, length 12
                else
                    G_mat = cell2mat(NIR_G.G11_coef');
                    NIR_mat   = cell2mat(NIR_G.NIR_coef');           % NIR-only, length 12
                end
                Concat_mat= [G_mat, NIR_mat];                     % length 24
                Diff_mat  = NIR_mat - G_mat;                      % length 12
            else
                if channel_norm
                    G_mat = cell2mat(NIR_G.G_coef_norm');
                    G_mat_pad = cell2mat(NIR_G.G_coef_pad_norm');
                    NIR_mat   = cell2mat(NIR_G.NIR_coef_norm');           % NIR-only, length 12
                else
                    G_mat = cell2mat(NIR_G.G_coef');
                    G_mat_pad = cell2mat(NIR_G.G_coef_pad');
                    NIR_mat   = cell2mat(NIR_G.NIR_coef');           % NIR-only, length 12
                end
                Concat_mat= [NIR_mat, G_mat];                     % length 24
                Diff_mat  = NIR_mat - G_mat_pad;                      % length 12
            end
            if use_ptt
                ptt_mat = cell2mat(NIR_G.ptt);
                ptt_mat = ptt_mat';
        
                G_mat = [G_mat, ptt_mat];
                NIR_mat = [NIR_mat, ptt_mat];
                Diff_mat = [Diff_mat, ptt_mat];
                Concat_mat = [Concat_mat, ptt_mat];
            end
                FeatSet   = {G_mat, NIR_mat, Diff_mat, Concat_mat};
            FeatName  = {'G only', 'NIR only', 'NIR - G', 'NIR ⊕ G'};
            
            SBP_vec   = cell2mat(NIR_G.SBP)';                 % n_seg × 1
            DBP_vec   = cell2mat(NIR_G.DBP)';                 % n_seg × 1
            sid_vec   = cell2mat(NIR_G.sid)';                 % subject IDs
            
            % reg_mdls = {'RF'};
            % reg_mdls = {'SVR'};
            % (2)  Prepare storage for metrics  .....................................
            nT   = numel(FeatSet);
            nM   = numel(reg_mdls);           % 2 models: SVR & RF
            MAE  = zeros(2,nT,nM);            % rows: SBP/DBP, cols: feat, pages: model
            RMSE = zeros(2,nT,nM);
            Rval = zeros(2,nT,nM);
            gt_bp = cell(2,nT,nM);
            pd_bp = cell(2,nT,nM);
            
            % (3)  Leave-One-Subject-Out loop ......................................
            unique_sids = unique(sid_vec);
            clrMap = lines(nM);               % one color per model
            mrkMap = {'o','^'};               % marker per model
            % clrMap_sid = lines(numel(unique_sids));  % unique color per subject
            n_subj = numel(unique_sids);  % total number of unique subjects
            clrMap_sid = parula(n_subj); % or jet(n_subj), turbo(n_subj), etc.
            sid_color_map = containers.Map(num2cell(unique_sids), ...
                                        num2cell(clrMap_sid, 2));
            
            for m = 1:nM                      % model loop (1 = SVR, 2 = RF)
                mdlName = reg_mdls{m};
                figure('Visible', true, 'Name', sprintf('SBP & DBP Regression – %s', mdlName), ...
               'Position', [100 100 2000 1200]);
                t = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
                for f = 1:nT                  % feature type loop
                    X = FeatSet{f};
            
                    % containers for aggregated predictions
                    pred_SBP = [];  true_SBP = [];
                    pred_DBP = [];  true_DBP = [];
            
                    for s = unique_sids'      % LOSO loop
                        testIdx  = sid_vec == s;   trainIdx = ~testIdx;
                        Xtrain   = X(trainIdx,:);  Xtest = X(testIdx,:);
                        ytrain_S = SBP_vec(trainIdx);  ytest_S = SBP_vec(testIdx);
                        ytrain_D = DBP_vec(trainIdx);  ytest_D = DBP_vec(testIdx);
            
                        % ----- choose regressor -----
                        switch mdlName
                            case 'SVR'
                                mdl_S = fitrsvm(Xtrain, ytrain_S, ...
                                                 'KernelFunction','rbf', ...
                                                 'Standardize', true, ...
                                                 'BoxConstraint', 10, ...
                                                 'Epsilon', 2);
                                mdl_D = fitrsvm(Xtrain, ytrain_D, ...
                                                 'KernelFunction','rbf', ...
                                                 'Standardize', true, ...
                                                 'BoxConstraint', 10, ...
                                                 'Epsilon', 2);
                            case 'RF'
                                mdl_S = fitrensemble(Xtrain, ytrain_S, 'Method','Bag');
                                mdl_D = fitrensemble(Xtrain, ytrain_D, 'Method','Bag');
                        end
            
                        pred_SBP = [pred_SBP; predict(mdl_S, Xtest)];
                        true_SBP = [true_SBP; ytest_S];
                        pred_DBP = [pred_DBP; predict(mdl_D, Xtest)];
                        true_DBP = [true_DBP; ytest_D];
                    end
            
                    % -------- compute & store metrics ----------
                    MAE (1,f,m) = mean(abs(pred_SBP-true_SBP));
                    RMSE(1,f,m) = sqrt(mean((pred_SBP-true_SBP).^2));
                    Rval(1,f,m) = corr(pred_SBP,true_SBP);
            
                    MAE (2,f,m) = mean(abs(pred_DBP-true_DBP));
                    RMSE(2,f,m) = sqrt(mean((pred_DBP-true_DBP).^2));
                    Rval(2,f,m) = corr(pred_DBP,true_DBP);
            
                    pd_bp{1,f,m} = pred_SBP;
                    pd_bp{2,f,m} = pred_DBP;
                    gt_bp{1,f,m} = true_SBP;
                    gt_bp{2,f,m} = true_DBP;
                    
    
                    newRow = table( ...
                    channel_norm   , ...                                 % Norm flag
                    same_deg       , ...                                 % Same-deg flag
                    {FeatName{f}}    , ...                                 % 'G only', ...
                    {reg_mdls{m}}    , ...                                 % 'SVR' or 'RF'
                    MAE (1,f,m)    , RMSE(1,f,m) , Rval(1,f,m) , ...     % SBP metrics
                    MAE (2,f,m)    , RMSE(2,f,m) , Rval(2,f,m) );        % DBP metrics
        
                    results = [results ; newRow];
                    % ----------- plotting ----------------------
                    % SBP subplot (row 1)
                    subplot(2,4,f);
                    for s = unique_sids'
                        sid_col = sid_color_map(s);
                        scatter(true_SBP(sid_vec==s), pred_SBP(sid_vec==s), 30, 'filled', ...
                                'Marker', mrkMap{m}, ...
                                'MarkerEdgeColor', sid_col, ...
                                'MarkerFaceColor', sid_col); hold on;
                    end
                    plot(xlim,xlim,'k--'); grid on
                    xlabel('True'); ylabel('Pred.');
                    title(sprintf('SBP – %s', FeatName{f}));      
                    if f == 2
                        legend_labels = arrayfun(@(s) sprintf('SID %d', s), unique_sids, 'UniformOutput', false);
                        legend(legend_labels, 'Location', 'east');
                    end
    
                    % DBP subplot (row 2)
                    subplot(2,4,4+f);
                    for s=unique_sids'
                        sid_col = sid_color_map(s);
                        scatter(true_DBP(sid_vec==s), pred_DBP(sid_vec==s), 30, 'filled', ...
                                'Marker', mrkMap{m}, ...
                                'MarkerEdgeColor', sid_col, ...
                                'MarkerFaceColor', sid_col); hold on;
                    end
                    plot(xlim,xlim,'k--'); grid on
                    xlabel('True'); ylabel('Pred.');
                    title(sprintf('DBP – %s', FeatName{f}));
                end
                % (4)  Add metric text boxes (once per subplot) .........................
                for f = 1:nT
                    % ---------- SBP ----------
                    subplot(2,4,f);
                    txt = {sprintf('%s  MAE %.2f  RMSE %.2f  R %.2f', ...
                                   mdlName, MAE(1,f,m), RMSE(1,f,m), Rval(1,f,m))};
                    text(min(xlim)+0.04*range(xlim), max(ylim)-0.08*range(ylim), txt, ...
                         'BackgroundColor','w','EdgeColor','k','FontSize',8);
                
                    % ---------- DBP ----------
                    subplot(2,4,4+f);
                    txt = {sprintf('%s  MAE %.2f  RMSE %.2f  R %.2f', ...
                                   mdlName, MAE(2,f,m), RMSE(2,f,m), Rval(2,f,m))};
                    text(min(xlim)+0.04*range(xlim), max(ylim)-0.08*range(ylim), txt, ...
                         'BackgroundColor','w','EdgeColor','k','FontSize',8);
                end
                sgtitle(sprintf("LOSO {%s} Norm%d SameDeg%d %s", mdlName, channel_norm, same_deg, title_name));
                % pause(5); % The figure need some time to adjust
                % fig_save_fn = sprintf("LOSO_%s_Norm%d_SameDeg%d_SSA%d-%d.png", mdlName, channel_norm, same_deg, ssa_deg(1), ssa_deg(end));
                % saveas(gcf, fullfile(fig_save_path, fig_save_fn));
            end
            
        end
    end
end
