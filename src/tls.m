function [Xhat] = tls(A,B)

n = size(A,2);
if size(A,1)~=size(B,1); error('Matrices are not conformant.'); end
R1 = [A B];
[~,~,V] = svd(R1,0);
r = size(A,2);
R = rq(V(:,r+1:end));Gamma = R(n+1:end,n-r+1:end);
Z = R(1:n,n-r+1:end);
Xhat = -Z/Gamma;

end