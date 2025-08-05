clc; clear; close all;
addpath(genpath(pwd));
%%
load("/Users/cst/Documents/Research/sustech/Multi-wavelength-rPPG-for-BP-main/Dataset/Depolarization_Data/ICU/1P/1/update_trace/updated_dh1_trace1.mat");
load("/Users/cst/Documents/Research/sustech/Multi-wavelength-rPPG-for-BP-main/Dataset/Depolarization_Data/ICU/1P/1/update_trace/updated_dh1_trace3.mat");
% load("/Users/cst/Documents/Research/sustech/Multi-wavelength-rPPG-for-BP-main/Dataset/Depolarization_Data/LAB/cold/1/updated_Albert_trace1.mat");
% load("/Users/cst/Documents/Research/sustech/Multi-wavelength-rPPG-for-BP-main/Dataset/Depolarization_Data/LAB/cold/1/updated_Albert_trace3.mat");
%% Parameters
fs = 125;                      % Sampling frequency
bp_range = [0.7 4];            % Bandpass range in Hz (42–240 BPM)
[b, a] = butter(3, bp_range / (fs/2), 'bandpass');  % 3rd-order Butterworth

%% Prepare signals
signals = {updated_Albert_trace1, updated_Albert_trace3};
titles = {'RGB: 0 degrees', 'RGB: 90 degrees'};
colors = {'r','g','b'};
N = length(signals{1});
t = (0:N-1)/fs;

ssa_deg = 1:4;
ssa_seg_len = 10;
ssa_seg_step = 2;
det_param = [400, 4];

figure;
for i = 1:2
    rgb = double(signals{i});       % Ensure double precision
    rgb_dc = -(rgb ./ mean(rgb, 2)-1);   % DC-normalization (divided by mean)
    

    % Bandpass filter each channel
    for color_ch = 1:1
        % rgb_filt(color_ch,:) = filtfilt(b, a, rgb_dc(color_ch,:));
        rgb_filt(color_ch,:) = ssa_overlap_merge(rgb_dc(color_ch,:), fs, ssa_seg_len, ssa_seg_step, ssa_deg, det_param);
    end

    % Plot
    hold on;
    for color_ch = 1:1
        plot(t, rgb_filt(color_ch,:), colors{color_ch-1+i}, 'DisplayName', ['Channel ', upper(colors{color_ch})]);
    end
    title(['Filtered & DC-normalized RGB - ', titles{i}]);
    xlabel('Time (s)');
    ylabel('Normalized Intensity');
    legend('show');
    grid on;
end
%% Visualize for ICU data
% plot(zmff1_2_ppg1');

%%
parent_path = 'Dataset/Depolarization_Data';
bp_stim = 'cold';
scene = 'icu';
% scene = 'lab';


output_dir = fullfile("Dataset/Depolarization_Data/raw_data");
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
fs = 125;

if strcmp(scene, 'lab')
    T = readtable('/Users/cst/Documents/Research/sustech/Multi-wavelength-rPPG-for-BP-main/Dataset/Depolarization_Data/LAB/BP - LAB.xlsx', 'ReadVariableNames', false);
    % SBP and DBP start from row 2 (index 2) and alternate
    sbp_rows = 1:2:size(T,1);
    dbp_rows = 2:2:size(T,1);
    
    SBP = table2array(T(sbp_rows, 3:5)); % Assuming columns 3–5 are SBP data
    DBP = table2array(T(dbp_rows, 3:5)); % Columns 3–5 again for DBP
    subjects = T{sbp_rows, 6};           % Column 6 holds subject names
    sub_ids = 1:10;
    data_path = fullfile(parent_path, upper(scene), bp_stim);
    for enum_id = 1:length(sub_ids)
        sub_id = sub_ids(enum_id);
        data_fn = dir(fullfile(data_path, int2str(sub_id), '*.mat'));
    
        for i_fn = 1:length(data_fn)
            fn = data_fn(i_fn).name;
            if strfind(data_fn(i_fn).name, 'time_segments')
                splitted_str = strsplit(fn, '_');
                sub_name = splitted_str{1};
                ts = load(fullfile(data_fn(i_fn).folder, data_fn(i_fn).name));
                ts = ts.(data_fn(i_fn).name(1:end-4));
            end
        
            if strfind(data_fn(i_fn).name, 'trace1')
                trace1 = load(fullfile(data_fn(i_fn).folder, data_fn(i_fn).name));
                trace1 = trace1.(data_fn(i_fn).name(1:end-4));
            end
            if strfind(data_fn(i_fn).name, 'trace3')
                trace3 = load(fullfile(data_fn(i_fn).folder, data_fn(i_fn).name));
                trace3 = trace3.(data_fn(i_fn).name(1:end-4));
            end
        end
        % SBP_data = NaN;
        % MBP_data = NaN;
        % DBP_data = NaN;
    
        data = struct;
        data.trace1 = trace1;
        data.trace3 = trace3;
        data.ts_bp = [fs*120, fs*180, fs*540]; % Measured at 2min, 3min and 6min
        idx1 = find(strcmpi(subjects, sub_name));
        disp(idx1);
        data.sbp = SBP(idx1,:);
        data.dbp = DBP(idx1,:);
        data.sid = sub_id;
        data.sname = sub_name;
        data.scene = scene;
        data.fs = 125;
    
        fn = [scene, '_', num2str(sub_id), '.mat'];
        save(fullfile(output_dir, fn), 'data');
    end
else
    sub_ids = 1:20;
    data_path = fullfile(parent_path, upper(scene));
    for enum_id = 1:length(sub_ids)
        sub_id = sub_ids(enum_id);
        for i_rec = 1:2
            data_fn = dir(fullfile(data_path, sprintf('%dP',sub_id),int2str(i_rec), 'update_trace/*.mat'));
            for i_fn = 1:length(data_fn)
                fn = data_fn(i_fn).name;
                if strfind(data_fn(i_fn).name, 'time_segments')
                    splitted_str = strsplit(fn, '_');
                    sub_name = splitted_str{1};
                    ts = load(fullfile(data_fn(i_fn).folder, data_fn(i_fn).name));
                    ts = ts.(data_fn(i_fn).name(1:end-4));
                end
            
                if strfind(data_fn(i_fn).name, 'trace1')
                    trace1 = load(fullfile(data_fn(i_fn).folder, data_fn(i_fn).name));
                    trace1 = trace1.(data_fn(i_fn).name(1:end-4));
                end
                if strfind(data_fn(i_fn).name, 'trace3')
                    trace3 = load(fullfile(data_fn(i_fn).folder, data_fn(i_fn).name));
                    trace3 = trace3.(data_fn(i_fn).name(1:end-4));
                end
            end
            BP_path = fullfile(data_path, sprintf('%dP',sub_id), int2str(i_rec), sprintf('%s.mat', sub_name));
            BP_data = load(BP_path);
            SBP_data = BP_data.SBP;
            % MBP_data = BP_data.MBP;
            DBP_data = BP_data.DBP;
            
            data = struct;
            data.trace1 = trace1;
            data.trace3 = trace3;
            data.ts_bp = (30:length(SBP_data)-30)*125;
            data.sbp = SBP_data;
            data.dbp = DBP_data;
            data.sid = sub_id+10;
            data.sname = sub_name;
            data.scene = scene;
            data.fs = 125;

            load(fullfile(data_path, sprintf('%dP',sub_id),int2str(i_rec), '/PTT1_PTT2.mat'));
            data.ptt = PTT1;
            fn = [scene, '_', num2str(sub_id), '_', int2str(i_rec), '.mat'];
            save(fullfile(output_dir, fn), 'data');
        end
    end
end
%%
list_dir = dir(fullfile(output_dir, '*.mat'));
data_all = cell(1, numel(list_dir));
scene_all = strings(1, numel(list_dir));
sid_all = zeros(1, numel(list_dir));

for i = 1:numel(list_dir)
    file_path = fullfile(output_dir, list_dir(i).name);
    temp = load(file_path);
    data_all{i} = temp.data;

    % Parse scene and sid from filename
    % tokens = regexp(list_dir(i).name, '(lab|icu)_(\d+|\d+_\d+)\.mat', 'tokens');
    % scene_all(i) = tokens{1}{1};
    % sid_all(i) = str2double(tokens{1}{2});
    split_dir = strsplit(list_dir(i).name(1:end-4), '_');
    scene_all(i) = split_dir{1};
    sid_all(i) = str2double(split_dir{2});
end

%% 第4步：画出信号的语谱图，其中X表示时间，Y表示频率（次/min），图中第一条90附近的线为基波，第2条180附近的线为一次谐波，线中颜色强度最大（红色）表示幅值最大的地方
dataFile = "/Users/cst/Documents/Research/sustech/Multi-wavelength-rPPG-for-BP-main/Dataset/Depolarization_Data/ICU/2P/1/ACDC_ppg[0.5-5]/yjz1_1_ppg1.mat";
load(dataFile)
param.fps = 125;

trace_raw = mjs2_1_ppg3(:,1:7500*5);
plot_spec_fig = true;

% 滤波后PPG信号R、G、B
ppg_1 = trace_raw(1,:);
[ppg_spc_1, ppg_rate_1] = spec(ppg_1, param.fps);
ppg_2 = trace_raw(2,:);
[ppg_spc_2, ppg_rate_2] = spec(ppg_2, param.fps);
ppg_3 = trace_raw(3,:);
[ppg_spc_3, ppg_rate_3] = spec(ppg_3, param.fps);
if plot_spec_fig
% plot 滤波后RGB频谱
figure, set(gcf, 'Position', [100 100 350 550]);
colormap jet; %
d = [0.12, 0.08];
subplot(6,1,1);imagesc(ppg_spc_1(1:220,:)); set(gca,'YDir','normal');ylim([0, 220]); % imagesc(ppg_spc_1(1:220,:))
title('mjs2滤波后RGB-PPG语谱图');
subplot(6,1,2);plot(ppg_rate_1,'LineWidth',2); ylim([0, 220]); xlim([0, size(ppg_rate_1,2)])
subplot(6,1,3);imagesc(ppg_spc_2(1:220,:)); set(gca,'YDir','normal');ylim([0, 220]); 
subplot(6,1,4);plot(ppg_rate_2,'LineWidth',2); ylim([0, 220]); xlim([0, size(ppg_rate_2,2)])
subplot(6,1,5);imagesc(ppg_spc_3(1:220,:)); set(gca,'YDir','normal');ylim([0, 220]); 
subplot(6,1,6);plot(ppg_rate_3,'LineWidth',2); ylim([0, 220]); xlim([0, size(ppg_rate_3,2)])
end
%%
figure;
plot(ppg_2)