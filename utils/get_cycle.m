function cycles = get_cycle(signal, peaks, len_before, len_after)
% Remove peaks at start and end if their cycles are incomplete
peaks = peaks((peaks+(double(len_after))<length(signal)) &(peaks-double(len_before)>1));
cycles = [];
for i_peak = 1:length(peaks)
    cycles = [cycles; signal(peaks(i_peak)-len_before:peaks(i_peak)+len_after)];
end
end