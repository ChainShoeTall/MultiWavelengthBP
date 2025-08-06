function [mean_cycle, selected_ind] = high_corr_median_cycle(cycles, method, ref_cycle)
%HIGH_CORR_MEAN_CYCLE Compute a robust mean cycle using high-correlation filtering
%
%   [mean_cycle, selected_ind] = HIGH_CORR_MEAN_CYCLE(cycles) computes a
%   representative mean cycle from a matrix of cycle waveforms by:
%     1. Calculating the median cycle across all input cycles.
%     2. Computing the Pearson correlation coefficient between each cycle and
%        the median cycle.
%     3. Selecting cycles whose correlation exceeds the 50th percentile
%        (median correlation).
%     4. Recomputing the median cycle from only these high-correlation cycles.
%
%   Inputs:
%     cycles         - An N-by-T matrix, where each row is a cycle of length T.
%
%   Outputs:
%     mean_cycle     - A 1-by-T row vector representing the refined mean cycle,
%                      computed from high-correlation cycles only.
%     selected_ind   - A numeric vector containing the indices of cycles
%                      used in the final mean calculation (those above
%                      median correlation).
%
%   Example:
%     [mean_cycle, idx] = high_corr_mean_cycle(all_cycles);
%     plot(mean_cycle); hold on;
%     plot(all_cycles(idx,:), '--');
%
%   Notes:
%     - Pearson correlation is computed using corrcoef().
%     - The function uses the median rather than the mean for robustness
%       against outliers.
%
%   See also: corrcoef, prctile, median

% Initial median cycle
if nargin<2
    method = 'median';
end

if nargin<3
    ref_cycle = [];
end

switch(method)
    case('median')
        if isempty(ref_cycle)
            ref_cycle = median(cycles, 1);
        end
        % Compute correlation with median cycle
        coefs = zeros(1, size(cycles, 1));
        for i_cyc = 1:size(cycles, 1)
            corr = corrcoef(cycles(i_cyc,:), ref_cycle);
            coefs(i_cyc) = corr(1,2);
        end
        
        % Select cycles above median correlation
        [sort_coefs, sort_ind] = sort(coefs, 'descend');
        selected_ind = sort_ind(1:floor(length(sort_ind)/2));
        
        % Recompute mean cycle from selected cycles
        mean_cycle = median(cycles(selected_ind,:), 1);
    case('mean')
        if isempty(ref_cycle)
            ref_cycle = mean(cycles, 1);
        end
        
        % Compute correlation with median cycle
        coefs = zeros(1, size(cycles, 1));
        for i_cyc = 1:size(cycles, 1)
            corr = corrcoef(cycles(i_cyc,:), ref_cycle);
            coefs(i_cyc) = corr(1,2);
        end
        
        % Select cycles above median correlation
        [sort_coefs, sort_ind] = sort(coefs, 'descend');
        selected_ind = sort_ind(1:floor(length(sort_ind)/2));
        
        % Recompute mean cycle from selected cycles
        mean_cycle = mean(cycles(selected_ind,:), 1);
end
