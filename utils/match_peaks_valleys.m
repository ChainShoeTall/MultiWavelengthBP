function valid_triplets = match_peaks_valleys(NIR_onsets, NIR_peaks, G_onsets, G_peaks, cycle_len, match_threshold, plot_flag)
%MATCH_PEAKS_VALLEYS Match peaks and valleys (onsets) between NIR and G signals.
%
%   valid_triplets = MATCH_PEAKS_VALLEYS(NIR_onsets, NIR_peaks, G_onsets, G_peaks, cycle_len, match_threshold, plot_flag)
%
%   This function finds valid matched segments ("triplets") between the onsets (valleys)
%   and peaks of NIR and G signals. Nearest neighbor matching within a threshold is used
%   to identify corresponding onsets and peaks.
%
%   INPUTS:
%       NIR_onsets      - array of onset indices (valleys) in the NIR signal.
%       NIR_peaks       - array of peak indices in the NIR signal.
%       G_onsets        - array of onset indices (valleys) in the G signal.
%       G_peaks         - array of peak indices in the G signal.
%       cycle_len       - estimated cycle length (in samples), typically median diff of G_onsets.
%       match_threshold - (optional) matching threshold in samples (default: 10).
%       plot_flag       - (optional) true to plot the result (default: false).
%
%   OUTPUTS:
%       valid_triplets  - cell array of structures. Each structure includes:
%           - nir_start  : start index of the segment in NIR signal
%           - g_start    : start index of the segment in G signal
%           - nir_end    : end index of the segment in NIR signal
%           - g_end      : end index of the segment in G signal
%           - nir_len    : length of the segment in NIR signal
%           - g_len      : length of the segment in G signal
%           - nir_peak   : peak index within the segment in NIR signal
%           - g_peak     : peak index within the segment in G signal
%
%   EXAMPLE:
%       triplets = match_peaks_valleys(NIR_onsets, NIR_peaks, G_onsets, G_peaks, 100);
%
%   Author: Shutao Chen
%   Date: 11-Jun-2025
% ------------------------------------------------------------------------------
if nargin < 3
    cycle_len = median(diff(G_onsets));
end

if nargin < 6
    match_threshold = 10;
end
if nargin < 7
    plot_flag = false;
end


matched_onsets = [];
for i = 1:length(G_onsets)-1
    [min_dist, idx] = min(abs(NIR_onsets - G_onsets(i)));
    if min_dist <= match_threshold
        matched_onsets = [matched_onsets; NIR_onsets(idx), G_onsets(i)];
    end
end

matched_peaks = [];
for i = 1:length(G_peaks)
    [min_dist, idx] = min(abs(NIR_peaks - G_peaks(i)));
    if min_dist <= match_threshold
        matched_peaks = [matched_peaks; NIR_peaks(idx), G_peaks(i)];
    end
end

% Initialize to store valid segments
valid_triplets = {};

for i = 1:(length(matched_onsets)-1)
    % Determine the segment length
    segment_len = matched_onsets(i+1,:) - matched_onsets(i,:);
    
    % Check condition 1: segment length is within ±10 of cycle_len
    if (abs(segment_len(1) - cycle_len) <= 2*match_threshold) && ...
        (abs(segment_len(2) - cycle_len) <= 2*match_threshold)
        % Check condition 2: there is at least one matched peak in this segment
        peaks_in_segment = matched_peaks( ... 
            (matched_peaks(:, 2) > matched_onsets(i, 2)) & ...
            (matched_peaks(:, 2) < matched_onsets(i+1, 2)) & ...
            (matched_peaks(:, 1) > matched_onsets(i, 1)) & ...
            (matched_peaks(:, 1) < matched_onsets(i+1, 1)), :);
        if ~isempty(peaks_in_segment)
            % Save this segment
            valid_triplet = struct;
            valid_triplet.nir_start = matched_onsets(i,1);
            valid_triplet.g_start = matched_onsets(i,2);
            valid_triplet.nir_end = matched_onsets(i+1,1);
            valid_triplet.g_end = matched_onsets(i+1,2);
            valid_triplet.nir_len = segment_len(1)+1;
            valid_triplet.g_len = segment_len(2)+1;
            valid_triplet.nir_peak = peaks_in_segment(1,1);
            valid_triplet.g_peak = peaks_in_segment(1,2);
            valid_triplets{end+1} = valid_triplet;
        end
    end
end


if plot_flag
    figure;
    hold on;

    % Initialize plot handles for legend entries
    segment_lines = [];
    for i = 1:length(valid_triplets)
        h = plot([valid_triplets{i}.nir_start, valid_triplets{i}.nir_peak, valid_triplets{i}.nir_end], ...
             [0, 1, 0], 'b-', 'LineWidth', 1.5);
        if i==1
            segment_lines = h;  % Save only the first for the legend
        end
    end

    % Plot onsets
    h1 = plot(NIR_onsets, zeros(size(NIR_onsets)), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    h2 = plot(G_onsets, zeros(size(G_onsets)), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');

    % Plot peaks
    h3 = plot(NIR_peaks, ones(size(NIR_peaks)), 'r*', 'MarkerSize', 8, 'LineWidth', 1.5);
    h4 = plot(G_peaks, ones(size(G_peaks)), 'g*', 'MarkerSize', 8, 'LineWidth', 1.5);

    xlabel('Sample Index');
    ylabel('Marker Position');
    ylim([-1, 2]);
    title('Onsets and Peaks Alignment');
    
    % Legend only for the 5 entries
    legend([segment_lines, h1, h2, h3, h4], ...
           {'Matched Segments', 'NIR Onsets', 'G Onsets', 'NIR Peaks', 'G Peaks'}, ...
           'Location', 'best');

    grid on;
end