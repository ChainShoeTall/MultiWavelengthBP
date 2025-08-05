% Get triplets given peaks, onsets, and cycle_len
function triplets = peak_onset_tripletize(peaks, onsets, cycle_len, thresh_multiplier, peak_ratio_denoise)
% Matches peak between each onset pair.
%   Inputs:
%     onsets     - vector of onset indices (must be sorted)
%     peaks      - vector of peak indices (must be sorted)
%     cycle_len  - minimum duration (in samples) between onsets 
%     thresh_multiplier - length of each cycle
%
%   Output:
%     triplets   - Nx3 matrix, each row is [onset_i, onset_i+1, peak]
    if nargin < 4
        thresh_multiplier = 0.25;
    end

    if nargin < 5
        peak_ratio_denoise = false;
    end

    triplets = [];
    peak_ratios = [];

    n_onsets = length(onsets);
    for i = 1:n_onsets - 1
        start_idx = onsets(i);
        end_idx = onsets(i+1);
        cycle_idx = onsets(i+1) - onsets(i);

        if abs(cycle_idx - cycle_len) > thresh_multiplier*cycle_len
            continue;  % skip too short cycles
        end

        % Find peaks in this interval
        in_range = peaks(peaks > start_idx & peaks < end_idx);

        if isempty(in_range) | length(in_range) > 1
            continue;  % no peak found between this pair of onsets
        end
        
        if peak_ratio_denoise
            peak_ratios(end+1) = (in_range - start_idx)/cycle_idx;
        end
        % Append to triplet list
        triplets(end+1, :) = [start_idx, end_idx, in_range];
    end

    if peak_ratio_denoise
        prct25 = prctile(peak_ratios, 25);
        prct75 = prctile(peak_ratios, 75);
        iqr = prct75 - prct25;
        peak_ratio_mask = find(peak_ratios > prct25-iqr*1.75 & ...
            peak_ratios < prct75+iqr*1.75);
        triplets = triplets(peak_ratio_mask, :);
    end
end