% Performs RQ decomposition
function [R,Q,varargout] = rq(A)

n = size(A,1);
if nargout<3
    [Q,R] = qr(flipud(A)',0);
else
    [Q,R,P1] = qr(flipud(A)',0);
    P(n+1-P1) = n:-1:1; % arrange permutation in right way
    varargout{1} = P;
end
R = rot90(R',2);
Q = flipud(Q');

[n,m]=size(A);

if n>m
    R = [zeros(n,n-m), R];
    Q = [zeros(n-m,m); Q];
end  

end