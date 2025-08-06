% Get triplets given peaks, onsets, and cycle length.
%
% Matches the dominant peak between each pair of consecutive onsets and
% optionally removes outliers based on the relative peak position.
%
% Inputs
%   peaks  - Sorted vector of peak indices.
%   onsets - Sorted vector of onset indices.
%   cycle_len - Expected duration (in samples) of each cycle.
%   thresh_multiplier - Allowed relative deviation from cycle_len. Default: 0.25
%   peak_ratio_denoise - Enable additional ratio based denoising. Default: false
%
% Output
%   triplets - Nx3 array where each row is [onset_i, onset_{i+1}, peak_i].

function triplets = peak_onset_tripletize(peaks, onsets, cycle_len, ...
    thresh_multiplier, peak_ratio_denoise)

    if nargin < 4 || isempty(thresh_multiplier)
        thresh_multiplier = 0.25;
    end

    if nargin < 5
        peak_ratio_denoise = false;
    end

    n_onsets = numel(onsets);
    % Preallocate maximum possible size to avoid dynamic resizing
    triplets = zeros(max(n_onsets - 1, 0), 3);
    peak_ratios = zeros(1, max(n_onsets - 1, 0));
    t_idx = 0;

    for i = 1:n_onsets - 1
        start_idx = onsets(i);
        end_idx   = onsets(i + 1);
        cycle_idx = end_idx - start_idx;

        % Skip cycles that deviate too much from expected length
        if abs(cycle_idx - cycle_len) > thresh_multiplier * cycle_len
            continue;
        end

        % Find peaks in this interval
        in_range = peaks(peaks > start_idx & peaks < end_idx);

        % Require exactly one peak per cycle
        if isempty(in_range) || numel(in_range) > 1
            continue;
        end

        t_idx = t_idx + 1;
        triplets(t_idx, :) = [start_idx, end_idx, in_range];

        if peak_ratio_denoise
            peak_ratios(t_idx) = (in_range - start_idx) / cycle_idx;
        end
    end

    % Remove unused preallocated rows
    triplets = triplets(1:t_idx, :);
    peak_ratios = peak_ratios(1:t_idx);

    if peak_ratio_denoise && ~isempty(peak_ratios)
        prct25 = prctile(peak_ratios, 25);
        prct75 = prctile(peak_ratios, 75);
        iqr = prct75 - prct25;
        mask = peak_ratios > prct25 - iqr * 1.75 & ...
               peak_ratios < prct75 + iqr * 1.75;
        triplets = triplets(mask, :);
    end
end
