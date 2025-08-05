function [pxx, wxx] = get_periodogram(X, fs)
[pxx, wxx] = periodogram(X, hamming(length(X)), length(X), fs);
end