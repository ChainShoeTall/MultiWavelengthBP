%%% Visualize Depolarization Mean Cycle Every 30 seconds
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
debug_mode = true;
if debug_mode   
    norm_cycle_len = 100;
    fs = 125;
    % param for ssa
    lf = 0.5;
    hf = 6.0;
    [B, A] = butter(4, [lf*2/125, hf*2/fs], 'bandpass');
    
    ssa_deg = 1:4;
    ssa_seg_len = 10;
    ssa_seg_step = 2;

    % for i_data=41
    for i_data=1:length(data_all)
        nir_data = data_all{i_data}.trace1(2,:);
        g_data   = data_all{i_data}.trace3(2,:);
        fs = data_all{i_data}.fs;
        ts = data_all{i_data}.ts_bp;
        sbp = data_all{i_data}.sbp;
        dbp = data_all{i_data}.dbp;
        sid = data_all{i_data}.sid;
        if sid > 10
            continue;
        end
    
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
    
    % 1. Envelope Outlier Removal
    [anomaly_idx, env_diff, upper_env, lower_env, nir_good_peaks, nir_good_onsets] = detect_motion_artifacts(nir_data, fs, 1.75, nir_peaks, nir_onsets);
    [~, ~, ~, ~, g_good_peaks, g_good_onsets] = detect_motion_artifacts(g_data, fs, 1.75, g_peaks, g_onsets);
    
    % 2. Cycle length and peak location outlier removal
    nir_triplets = peak_onset_tripletize(nir_good_peaks, nir_good_onsets, g_median_cycle_len, 0.25, true);
    g_triplets = peak_onset_tripletize(g_good_peaks, g_good_onsets, g_median_cycle_len, 0.25, true);
    

    ts_selection = 15*fs:30*fs:size(nir_data,2)-15*fs;
    bp_measurement_range = int32(15*fs);
    all_nir_median = {};
    all_g_median = {};
    for i_ts = ts_selection
        start_point = max(i_ts - bp_measurement_range, 1);
        end_point = min(i_ts + bp_measurement_range, length(nir_data));
        nir_trip_seg = nir_triplets(nir_triplets(:,1)>start_point & nir_triplets(:,2)<end_point,:);
        g_trip_seg = g_triplets(g_triplets(:,1)>start_point & g_triplets(:,2)<end_point,:);

        nir_cycles_seg = get_trip_norm_cycles(nir_data, nir_trip_seg, norm_cycle_len);
        g_cycles_seg = get_trip_norm_cycles(g_data, g_trip_seg, norm_cycle_len);
        [nir_median_cycle, selected_nir_ind] = high_corr_median_cycle(nir_cycles_seg,'mean');
        [g_median_cycle, selected_g_ind] = high_corr_median_cycle(g_cycles_seg, 'mean'); 

        all_nir_median{end+1} = nir_median_cycle;
        all_g_median{end+1} = g_median_cycle;
    end
    
    % ==== Visualization of Median Cycles ====
    figure; 
    n_plot = length(ts_selection);
    time_axis = linspace(0, 1, norm_cycle_len);  % normalized cycle x-axis
    
    for i = 1:n_plot
        subplot(ceil(n_plot/4), 4, i);  % 4 columns layout
        hold on;
    
        % Set transparency by RGBA (4th value is alpha), range 0 (transparent) to 1 (opaque)
        nir_color = [1, 0, 0, 0.6];  % red with 60% opacity
        g_color   = [0, 0.6, 0, 0.6]; % green with 60% opacity
        
        if i <= length(all_nir_median) && ~isempty(all_nir_median{i})
            plot(time_axis, all_nir_median{i}, 'LineWidth', 2.5, ...
                'Color', nir_color);
        end
        if i <= length(all_g_median) && ~isempty(all_g_median{i})
            plot(time_axis, all_g_median{i}, 'LineWidth', 2.5, ...
                'Color', g_color);
        end
    
        % Optional: add timestamp info
        ts_sec = double(ts_selection(i)) / fs;
        title(sprintf('%.1f sec', ts_sec), 'FontSize', 8);
        ylim([-0.5 1.5]);  % adjust based on your signal range
        xlim([0 1]);
        grid on;
    end
    
    sgtitle(sprintf('Depolarization Mean Cycles (Subject %d, Scene %s)', sid, scene_all(i_data)), 'FontWeight', 'bold');
    legend({'G0', 'G90'}, 'Location', 'best');
    
    % ==== Visualization: NIR and G Median Cycles Separately ====
    n_plot = length(ts_selection);
    time_axis = linspace(0, 1, norm_cycle_len);  % normalized x-axis for each cycle
    
    % Generate color map
    cmap = turbo(n_plot);  % Use visually distinct color map (e.g., turbo > jet)
    
    figure;
    tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % === NIR Plot ===
    nexttile; hold on;
    for i = 1:n_plot
        if i <= length(all_nir_median) && ~isempty(all_nir_median{i})
            plot(time_axis, all_nir_median{i}, 'Color', cmap(i,:), 'LineWidth', 2);
        end
    end
    title('G0 Mean Cycles');
    xlabel('Normalized Time'); ylabel('Amplitude');
    xlim([0, 1]); grid on; box on;
    
    % === G Plot ===
    nexttile; hold on;
    for i = 1:n_plot
        if i <= length(all_g_median) && ~isempty(all_g_median{i})
            plot(time_axis, all_g_median{i}, 'Color', cmap(i,:), 'LineWidth', 2);
        end
    end
    title('G90 Mean Cycles');
    xlabel('Normalized Time'); ylabel('Amplitude');
    xlim([0, 1]); grid on; box on;
    
    % === Add a single shared colorbar ===
    cb = colorbar('Position', [0.93 0.15 0.02 0.7]); % adjust position if needed
    colormap(turbo(n_plot));
    cb.Ticks = linspace(0, 1, 5);
    cb.TickLabels = arrayfun(@(x) sprintf('%.0f s', x), ...
        linspace(ts_selection(1), ts_selection(end), 5)/fs, 'UniformOutput', false);
    ylabel(cb, 'Timestamp (sec)');
    
    sgtitle(sprintf('Depolarization Mean Cycles Over Time (Subject %d, Scene %s)', sid, scene_all(i_data)), ...
        'FontWeight', 'bold');
    end
end
