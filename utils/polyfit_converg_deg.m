function [NIR_deg, NIR_p, G_deg, G_p] = polyfit_converg_deg(mean_NIR_cycle, mean_G_cycle, stop_thresh)
% Get degree that polyfit converges, and corresponding coefficients
if nargin<3
    stop_thresh = 0.99;
end

fit_deg = 1;
x = linspace(0, 1, length(mean_NIR_cycle));
% x = 1:length(mean_NIR_cycle);
G_stop = false;
NIR_stop = false;
while ~G_stop || ~NIR_stop
    if ~NIR_stop
        [p_NIR, S_NIR, mu_NIR] = polyfit(x, mean_NIR_cycle, fit_deg);
        if S_NIR.rsquared > stop_thresh
            NIR_stop = true;
            NIR_deg = fit_deg;
            NIR_p = p_NIR;
        end
    end

    if ~G_stop
        [p_G, S_G, mu_G] = polyfit(x, mean_G_cycle, fit_deg);
        if S_G.rsquared > stop_thresh
            G_stop = true;
            G_deg = fit_deg;
            G_p = p_G;
        end
    end
    fit_deg = fit_deg+1;

    if fit_deg > 20
        G_deg = -1;
        fprintf('Complex curve that cannot be fitted within 20 deg\n')
        break
    end
end

% fprintf("NIR fit stop at deg: %d;\n " + ...
%     "G fit stop at deg: %d,\n", NIR_deg, G_deg);
end