function [all_NIR_cycles, all_G_cycles] = get_norm_cycles(global_trips, NIR_ups, G_ups, cycle_len)
n_cycle = length(global_trips);
all_NIR_cycles = [];
all_G_cycles = [];
for i_cyc = 1:n_cycle
    trip = global_trips{i_cyc};
    NIR_cyc = NIR_ups(trip.nir_start:trip.nir_end);
    G_cyc = G_ups(trip.g_start:trip.g_end);
    NIR_peak_ind = trip.nir_peak - trip.nir_start+1;
    G_peak_ind = trip.g_peak - trip.g_start+1;

    NIR_resampled = interp_to_length(NIR_cyc, cycle_len);
    G_resampled = interp_to_length(G_cyc, cycle_len);
    NIR_debase = NIR_resampled - linspace(NIR_resampled(1), NIR_resampled(end), cycle_len);
    G_debase = G_resampled - linspace(G_resampled(1), G_resampled(end), cycle_len);
    NIR_cyc_norm = NIR_debase/max(NIR_debase);
    G_cyc_norm = G_debase/max(G_debase);
    
    if any(G_cyc_norm<0) || any(NIR_cyc_norm<0)
        continue;
    end
    
    all_NIR_cycles = [all_NIR_cycles; NIR_cyc_norm];
    all_G_cycles = [all_G_cycles; G_cyc_norm];
end
end


% figure();
% subplot(2,1,1);
% h_nir_cyc = plot(all_NIR_cycles'); hold on;
% h_nir_cyc(1).Color(4) = 0.3;
% plot(mean(all_NIR_cycles, 1), 'r-', 'LineWidth', 2);
% 
% subplot(2,1,2)
% h_g_cyc = plot(all_G_cycles'); hold on;
% h_g_cyc(1).Color(4) = 0.3;
% plot(mean(all_G_cycles, 1), 'r', 'LineWidth', 2);
