function reconstructed_signal = reconstruct_signal_from_fft(fft_coefficients, N)
% Reconstruct the time-domain signal from FFT coefficients
% Inputs:
%   fft_coefficients - complex FFT coefficients
%   N                - original signal length
% Output:
%   reconstructed_signal - time-domain signal

    reconstructed_signal = real(ifft(fft_coefficients, N));
end