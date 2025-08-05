function all_cycles = get_trip_norm_cycles(ppg, triplets, cycle_len)
n_cycle = size(triplets,1);
all_cycles = zeros(n_cycle, cycle_len);
for i_cyc = 1:n_cycle
    cyc = ppg(triplets(i_cyc,1):triplets(i_cyc,2));
    cyc = interp_to_length(cyc, cycle_len);
    cyc_debase = cyc - linspace(cyc(1), cyc(end), cycle_len);
    all_cycles(i_cyc, :) = cyc_debase/max(cyc_debase);
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
