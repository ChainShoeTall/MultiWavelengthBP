function r = ssa_overlap_merge(s, fps, win_sec, step_sec, K, det_param)
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
    N = length(s);                    % signal length

    r = zeros(1, N);                  % result signal
    w = zeros(1, N);                  % weight accumulator

    win = hann(L)';                   % row vector Hanning window

    if length(det_param) == 2
        A1 = tar_detrend(L, det_param(1)); % Remove low-freq
        A2 = tar_detrend(L, det_param(2)); % Remove high-freq
        A = A1 - A2; % keep between low-freq and high-freq
    else
        A = tar_detrend(L, det_param); % Remove low-freq only
    end

    for start_idx = 1:step:(N - L + 1)
        seg = s(start_idx : start_idx + L - 1);
        % seg = seg - mean(seg);             % zero-mean
        seg = seg * A;
        seg_denoised = ssa(seg, fps, K, false);      % apply SSA
        seg_denoised = seg_denoised .* win;  %ˇ apply window

        % overlap-add with weighting
        r(start_idx : start_idx + L - 1) = r(start_idx : start_idx + L - 1) + seg_denoised;
        w(start_idx : start_idx + L - 1) = w(start_idx : start_idx + L - 1) + win;
    end

    % normalize
    w(w == 0) = 1;        % to avoid divide-by-zero
    r = r ./ w;
end

