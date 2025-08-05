function [output] = butter_bandpass_filt(X, lb, ub, fs)
[b, a]= butter(4, [lb*2/fs, ub*2/fs], 'bandpass');
output = filtfilt(b, a, X);
end