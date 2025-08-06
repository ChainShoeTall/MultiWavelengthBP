function [anomaly_idx, env_diff, upper_env, lower_env, good_peaks, good_troughs] = detect_motion_artifacts(ppg, fs, threshold_multiplier, pk_locs, tr_locs)
% Motion artifact detection based on envelope difference and IQR
% Inputs:
%   ppg                - raw PPG signal (1D array)
%   fs                 - sampling frequency (Hz)
%   threshold_multiplier - IQR scaling factor (default = 2)
% Outputs:
%   anomaly_idx        - indices of motion artifacts
%   env_diff           - envelope difference signal
%   upper_env, lower_env - upper and lower envelopes

    if nargin < 3
        threshold_multiplier = 1.75;
    end
    if nargin < 5
        [pk_locs, tr_locs] = findpeaks_pyampd(ppg, fs);
    end

    %% Step 1: Bandpass filter
    pks = ppg(pk_locs);
    trs = ppg(tr_locs);
    % [pks, pk_locs] = findpeaks(ppg, 'MinPeakDistance', peak_dist);
    % [trs, tr_locs] = findpeaks(-ppg, 'MinPeakDistance', peak_dist);
    % trs = -trs;

    %% Step 3: Interpolate envelopes
    upper_env = interp1(pk_locs, pks, 1:length(ppg), 'linear', 'extrap');
    lower_env = interp1(tr_locs, trs, 1:length(ppg), 'linear', 'extrap');

    %% Step 4: Envelope difference
    env_diff = abs(upper_env - lower_env);

    %% Step 5: Slope change detection
    slope_idx = [];
    for i = 2:length(env_diff)-1
        diff_prev = env_diff(i) - env_diff(i-1);
        diff_next = env_diff(i+1) - env_diff(i);
        if (diff_prev >= 0 && diff_next < 0) || (diff_prev < 0 && diff_next >= 0)
            slope_idx(end+1) = i;
        end
    end

    %% Step 6: IQR thresholding
    Q1 = prctile(env_diff, 25);
    Q3 = prctile(env_diff, 75);
    IQR = Q3 - Q1;
    upper_thresh = Q3 + threshold_multiplier * IQR;
    lower_thresh = Q1 - threshold_multiplier * IQR;

    %% Step 7: Detect anomalies at slope changes
    anomaly_idx = slope_idx(env_diff(slope_idx) > upper_thresh | ...
                            env_diff(slope_idx) < lower_thresh);

    % Filter peaks and troughs by removing any that fall within the expanded anomaly indices
    pk_mask = true(size(pk_locs));
    % Trough mask should match the number of trough locations rather than peaks
    tr_mask = true(size(tr_locs));
    for i=1:length(pk_locs)
        if ismember(pk_locs(i), anomaly_idx)
            pk_mask(max(i-2,1):min(i+2,length(pk_locs))) = false;
        end
    end
    good_peaks = pk_locs(pk_mask);
    for i=1:length(tr_locs)
        if ismember(tr_locs(i), anomaly_idx)
            tr_mask(max(i-2,1):min(i+2,length(tr_locs))) = false;
        end
    end
    good_troughs = tr_locs(tr_mask);
end
