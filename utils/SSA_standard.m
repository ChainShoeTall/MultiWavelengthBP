function RC = SSA_standard(X, M)
% Note: This implementation is slightly difference from the reference, which
% stacks the embeding in each row; Also it computes the covariance matrix
% using the correlation and toeplitz function as an alternative.
X = X-mean(X);
X = X/std(X);
N = size(X,2);
% Compute covariance using TOEPLITZ function
% covX = xcorr(X, M-1, 'unbiased'); % 2*M-1 square mat
% Ctoep = toeplitz(covX(M:end)); % M components
% C = Ctoep;
Y=zeros(N-M+1,M);
for m=1:M    
  Y(:,m) = X((1:N-M+1)+m-1);
end
Cemb=Y'*Y / (N-M+1);
C =Cemb;
[RHO, LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);
[~, ind] = sort(LAMBDA, 'descend');
RHO = RHO(:,ind);

PC = Y*RHO;
RC = zeros(N, M);
for m=1:M
    buf = PC(:,m)*RHO(:,m)';
    buf = buf(end:-1:1,:);
    for n=1:N % anti-diagonal averaging
        RC(n,m) = mean(diag(buf, -(N-M+1)+n));
    end
end
end