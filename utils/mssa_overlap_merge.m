function [r1, r2] = mssa_overlap_merge(s1, s2, fps, win_sec, step_sec, K, det_param)
% ssa_overlap_hann_merge  Apply SSA with overlap and Hanning-weighted merge
%
% Inputs:
%   s         - Input 1D signal
%   fps       - Sampling rate (Hz)
%   win_sec   - SSA window length in seconds (e.g., 1.5)
%   step_sec  - Step size (hop length) in seconds (e.g., 0.5)
%
% Output:
%   r         - Reconstructed full signal after SSA + weighted merge

    if nargin<5
        K = 2:6;
    end

    if nargin<6
        det_param = [100, 1];
    end

    L = round(win_sec * fps);         % window length in samples
    step = round(step_sec * fps);     % step/hop size in samples
    N = length(s1);                    % signal length

    r1 = zeros(1, N);                  % result signal
    w = zeros(1, N);                  % weight accumulator
    r2 = zeros(1, N);                  % result signal

    win = hann(L)';                   % row vector Hanning window

    if length(det_param) == 2
        A1 = tar_detrend(L, det_param(1)); % Remove low-freq
        A2 = tar_detrend(L, det_param(2)); % Remove high-freq
        A = A1 - A2; % keep between low-freq and high-freq
    else
        A = tar_detrend(L, det_param); % Remove low-freq only
    end

    for start_idx = 1:step:(N - L + 1)
        seg1 = s1(start_idx : start_idx + L - 1);
        seg2 = s2(start_idx : start_idx + L - 1);
        % seg = seg - mean(seg);             % zero-mean
        seg1 = seg1 * A;
        seg2 = seg2 * A;
        comb = [seg1', seg2'];
        [Y, ~, ~, ~] = mssa(comb, 1*fps); % Y in shape [N, D, D*M];
        seg1_denoised = sum(Y(:,1,K),3)' .* win;  %ˇ apply window
        seg2_denoised = sum(Y(:,2,K),3)' .* win;

        % overlap-add with weighting
        r1(start_idx : start_idx + L - 1) = r1(start_idx : start_idx + L - 1) + seg1_denoised;
        w(start_idx : start_idx + L - 1) = w(start_idx : start_idx + L - 1) + win;

        r2(start_idx : start_idx + L - 1) = r2(start_idx : start_idx + L - 1) + seg2_denoised;
    end

    % normalize
    w(w == 0) = 1;        % to avoid divide-by-zero
    r1 = r1 ./ w;
    r2 = r2 ./ w;
end

