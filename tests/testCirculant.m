addpath('../src');

n = 5; m = 1e2;             % Number of states and samples
% Construct circulant matrix
c = randn(1,n);
A = toeplitz([c(1) fliplr(c(2:end))], c);

% Form data matrices
X = randn(n,m);
Y = A*X;

Acirc =    piDMD(X,Y,'circulant');
AcircTLS = piDMD(X,Y,'circulantTLS');
I = eye(n);

errCirc    = norm(A - Acirc(I),'fro')/norm(A,'fro');
errCircTLS = norm(A - AcircTLS(I),'fro')/norm(A,'fro');

disp(['circulant-piDMD error:    ' num2str(errCirc)])
disp(['circulantTLS-piDMD error: ' num2str(errCircTLS)])