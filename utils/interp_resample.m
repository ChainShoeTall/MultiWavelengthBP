function [X_ups, ts, new_ts] = interp_resample(X, original_fs, new_fs)
%INTERP_RESAMPLE Resample a 1D signal to a new sampling frequency using interpolation
%
%   [X_ups, ts, new_ts] = interp_resample(X, original_fs, new_fs)
%
%   This function interpolates a 1D signal `X` from its original sampling
%   frequency `original_fs` to a new target frequency `new_fs`, using
%   piecewise cubic Hermite interpolation (PCHIP).
%
%   INPUTS:
%     X           - A 1D input signal (vector) to be resampled.
%     original_fs - (Optional) Original sampling frequency (Hz). Default = 30.
%     new_fs      - (Optional) Desired sampling frequency (Hz). Default = 100.
%
%   OUTPUTS:
%     X_ups   - The resampled signal at the new frequency.
%     ts      - Original time vector corresponding to X.
%     new_ts  - New time vector corresponding to X_ups.
%
%   EXAMPLE USAGE:
%     signal = sin(2*pi*1*(0:1/30:5));       % 1 Hz sine wave at 30 Hz
%     [resampled, ts, new_ts] = interp_resample(signal, 30, 100);
%
%   See also INTERP1

arguments (Input)
    X
    original_fs (1,1) double = 30
    new_fs (1,1) double = 100
end

arguments (Output)
    X_ups
    ts
    new_ts
end
ts = (0:length(X)-1)/original_fs;
new_ts = ts(1):1/new_fs:ts(end);
X_ups = interp1(ts, X, new_ts, "pchip");
end