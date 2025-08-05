function [peaks, valleys] = findpeaks_pyampd(sig, fs)
sig_struct = struct;
sig_struct.v = sig;
sig_struct.fs = fs;
[peaks, valleys, ] = detect_ppg_beats(sig_struct, 'AMPD');
end