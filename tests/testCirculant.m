addpath('../src');

n = 500;
m = 1e2;
c = randn(1,n);
A = toeplitz([c(1) fliplr(c(2:end))], c);

X = randn(n,m);
Y = A*X;

Acirc = piDMD(X,Y,'circulant');
AcircTLS = piDMD(X,Y,'circulantTLS');

errCirc = norm(A - Acirc,'fro')/norm(A,'fro');
errCircTLS = norm(A - AcircTLS,'fro')/norm(A,'fro');

disp(['circulant-piDMD error:    ' num2str(errCirc)])
disp(['circulantTLS-piDMD error: ' num2str(errCircTLS)])