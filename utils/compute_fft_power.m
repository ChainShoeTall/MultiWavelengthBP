function [fxx, pxx, fft_result, N] = compute_fft_power(signal, fs)
% Compute the frequency and power spectrum using FFT
% Inputs:
%   signal - time-domain signal
%   fs     - sampling frequency
% Outputs:
%   fxx    - frequency vector (Hz)
%   pxx    - power spectrum

    N = length(signal);
    fft_result = fft(signal);
    pxx = abs(fft_result / N).^2;
    pxx = pxx(1:floor(N/2)+1);
    pxx(2:end-1) = 2 * pxx(2:end-1);
    fxx = fs * (0:floor(N/2)) / N;
end