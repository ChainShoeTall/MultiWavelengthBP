function [x_norm] = minmax_norm(x)
x_norm = (x-min(x))/(max(x)-min(x));
end