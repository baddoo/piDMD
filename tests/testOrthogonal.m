nx = 4; nt = 1e3;
A0 = randn(nx);
[A,~,~] = svd(A0); % Make orthogonal
X0 = randn(nx,nt);
X = randn(nx,nt);
Y = A*X;

piA = piDMD(X,Y,'orthogonal',4);

