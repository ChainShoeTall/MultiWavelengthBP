function [f, Y] = simple_fft(signal, fs, plot_flag)
    % Inputs:
    %   signal - 1D time-domain signal
    %   fs     - sampling frequency (Hz)
    % Outputs:
    %   f      - frequency axis (Hz)
    %   Y      - magnitude of FFT

    if nargin < 3
        plot_flag = false;
    end
    % Ensure signal is a row vector
    signal = signal(:)';
    N = length(signal);

    % Remove mean to avoid DC offset
    signal = signal - mean(signal);

    % Compute FFT
    Y_complex = fft(signal);
    Y = abs(Y_complex / N);     % Normalize magnitude
    Y = Y(1:floor(N/2)+1);      % One-sided spectrum
    Y(2:end-1) = 2*Y(2:end-1);  % Compensate for symmetry

    % Frequency axis
    f = fs * (0:floor(N/2)) / N;

    % Plot
    if plot_flag
        figure;
        plot(f, Y, 'b-', 'LineWidth', 1.5);
        xlabel('Frequency (Hz)');
        ylabel('Magnitude');
        title('FFT Spectrum');
        grid on;
    end
end