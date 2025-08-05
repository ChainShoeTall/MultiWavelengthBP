function D = tar_detrend(T, lambda)
% High-pass with freq ~~ fs/(2*pi*lambda)
% 1
% ~15.9 Hz (very weak smoothing)
% 10
% ~1.59 Hz
% 20
% ~0.80 Hz
% 100
% ~0.16 Hz (very strong smoothing)

if nargin<1
    lambda = 10;
end
I = speye(T); 
D2 = spdiags(ones(T-2,1)*[1 -2 1],0:2,T-2, T);
D = (I-inv(I+lambda^2*(D2'*D2)));