%GET_TRIP_NORM_CYCLES Extract normalized cycles from a PPG signal.
%
%   cycles = GET_TRIP_NORM_CYCLES(ppg, triplets, cycle_len) returns a matrix
%   where each row is a single pulse cycle extracted from the vector ppg.
%   The cycles are defined by the onset/offset indices in triplets and are
%   interpolated to cycle_len samples with baseline removal and amplitude
%   normalization.
%
%   Inputs
%       ppg       - Vector containing the PPG signal.
%       triplets  - Nx3 array where each row is [onset, offset, peak].
%       cycle_len - Number of samples to interpolate each cycle to.
%
%   Output
%       cycles    - Matrix of size N-by-cycle_len containing normalized cycles.

function cycles = get_trip_norm_cycles(ppg, triplets, cycle_len)

    n_cycles = size(triplets, 1);
    cycles = zeros(n_cycles, cycle_len);

    for i = 1:n_cycles
        % Extract the cycle and interpolate to the desired length
        cyc = ppg(triplets(i, 1):triplets(i, 2));
        cyc = interp_to_length(cyc, cycle_len);

        % Remove linear trend and normalise by maximum amplitude
        baseline = linspace(cyc(1), cyc(end), cycle_len);
        cyc = cyc - baseline;
        cycles(i, :) = cyc / max(cyc);
    end
end

