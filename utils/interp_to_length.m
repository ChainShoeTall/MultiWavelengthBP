function y_resampled = interp_to_length(y, target_len, method)
%INTERP_TO_LENGTH Interpolates the input signal to any desired length
%
% USAGE:
%   y_resampled = interp_to_length(y, target_len)
%   y_resampled = interp_to_length(y, target_len, method)
%
% INPUTS:
%   y           - Original 1D array (row or column vector)
%   target_len  - Desired length after resampling
%   method      - Interpolation method (default: 'linear')
%
% OUTPUT:
%   y_resampled - Resampled signal to target length

    % Default to linear interpolation if not specified
    if nargin < 3
        method = 'pchip';
    end

    % Create original index
    orig_idx = linspace(1, length(y), length(y));
    
    % Create target index
    target_idx = linspace(1, length(y), target_len);

    % Interpolate
    y_resampled = interp1(orig_idx, y, target_idx, method);
end