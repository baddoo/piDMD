%% Physics-informed dynamic mode decompositions
%
% Computes a dynamic mode decomposition when the solution matrix is
% constrained to lie in a matrix manifold. 
% The options available for the "method" so far are
%
% - "exact", "exactSVDS"
% - "orthogonal"
% - "uppertriangular, "lowertriangular"
% - "diagonal", "diagonalpinv", "diagonaltls", "symtridiagonal"
% - "circulant", "circulantTLS", "circulantunitary", "circulantsymmetric",
% "circulantskewsymmetric"
% - "BCCB", "BCCBtls", "BCCBskewsymmetric",
% "BCCBunitary", "hankel", "toeplitz"
% - "symmetric", "skewsymmetric"
%
function [A, varargout] = piDMD(X,Y,method,varargin)

[nx, nt] = size(X);

if strcmp(method,'exact') || strcmp(method,'exactSVDS')
        
    if nargin>3
        r = varargin{1};
    else
        r = min(nx,nt);
    end

    if strcmp(method,'exact')
        [Ux,Sx,Vx] = svd(X,0);
        Ux = Ux(:,1:r); Sx = Sx(1:r,1:r); Vx = Vx(:,1:r);
    elseif strcmp(method,'exactSVDS')
        [Ux,Sx,Vx] = svds(X,r);
    end

    Atilde = (Ux'*Y)*Vx*pinv(Sx);
    A = @(v) Ux*(Atilde*(Ux'*v));

    if nargout==2
    varargout{1} = eig(Atilde);
    elseif nargout>2
    [eVecs,eVals] = eig(Atilde);
    eVals = diag(eVals);
    eVecs = Y*Vx*pinv(Sx)*eVecs./eVals.';
    varargout{1} = eVals;
    varargout{2} = eVecs;
    end

elseif strcmp(method,'orthogonal')
    
    if nargin>3
        r = varargin{1}; 
    else 
        r = min(nx,nt);
    end
    [Ux,~,~] = svd(X,0); Ux = Ux(:,1:r);
    Yproj = Ux'*Y; Xproj = Ux'*X; % Project X and Y onto principal components
    [Uyx, ~, Vyx] = svd(Yproj*Xproj',0);
    Aproj = Uyx*Vyx';    
    A = @(x) Ux*(Aproj*(Ux'*x));
    
    if nargout==2
        eVals = eig(Aproj);
        varargout{1} = eVals;
    elseif nargout>2
        [eVecs,eVals] = eig(Aproj);
        varargout{1} = diag(eVals);
        varargout{2} = Ux*eVecs;
        if nargout > 3; varargout{3} = Aproj; end
    end
    
elseif strcmp(method,'uppertriangular')
    
    [R,Q] = rq(X); % Q*Q' = I
    Ut = triu(Y*Q');
    A = Ut/R;

elseif strcmp(method,'lowertriangular')
    
    A = rot90(piDMD(flipud(X),flipud(Y),'uppertriangular'),2);
    
% The codes allows for matrices of variable banded width. The fourth input,
% a 2xn matrix called d, specifies the upper and lower bounds of the
% indices of the non-zero elements. The first column corresponds to the width of
% the band below the diagonal and the second column is the width of the
% band above. For example, a diagonal matrix would have d = ones(nx,2) and 
% a tridiagonal matrix would have d = [2 2]+zeros(nx,2). If you only specify
% d as a scalar then the algorithm converts the input to obtain a banded 
% diagonal matrix of width d. 
elseif startsWith(method,'diagonal') 
    if nargin>3
        d = varargin{1}; % arrange d into an nx-by-2 matrix
        if numel(d) == 1
            d = d*ones(nx,2);
        elseif numel(d) == nx
             d = repmat(d,[1,2]);
        elseif any(size(d)~=[nx,2])
            error('Diagonal number is not in an allowable format.')
        end
    else 
        d = ones(nx,2); % default is for a diagonal matrix
    end
    % Allocate cells to build sparse matrix
    Icell = cell(1,nx); Jcell = cell(1,nx); Rcell = cell(1,nx);
    for j = 1:nx
    l1 = max(j-(d(j,1)-1),1); l2 = min(j+(d(j,2)-1),nx);
    C = X(l1:l2,:); b = Y(j,:); % preparing to solve min||Cx-b|| along each row
    if strcmp(method,'diagonal')
            sol = b/C;
    elseif strcmp(method,'diagonalpinv')
            sol = b*pinv(C);
    elseif strcmp(method,'diagonaltls')
            sol = tls(C.',b.').';
    end
    Icell{j} = j*ones(1,1+l2-l1); Jcell{j} = l1:l2; Rcell{j} = sol;
    end
    Imat = cell2mat(Icell); Jmat = cell2mat(Jcell); Rmat = cell2mat(Rcell);
    Asparse = sparse(Imat,Jmat,Rmat,nx,nx);
    A = @(v) Asparse*v;

    if nargout==2
        eVals = eigs(Asparse,nx);
        varargout{1} = eVals;
    elseif nargout>2
        [eVecs, eVals] = eigs(Asparse,nx);
        varargout{1} = diag(eVals); varargout{2} = eVecs;
    end

elseif strcmp(method,'symmetric') || strcmp(method,'skewsymmetric')
    
[Ux,S,V] = svd(X,0);
C = Ux'*Y*V;
C1 = C;
if nargin>3; r = varargin{1}; else; r = rank(X); end
Ux = Ux(:,1:r);
Yf = zeros(r);
    if strcmp(method,'symmetric') 
    for i = 1:r
        Yf(i,i) = real(C1(i,i))/S(i,i);
        for j = i+1:r
            Yf(i,j) = (S(i,i)*conj(C1(j,i)) + S(j,j)*C1(i,j)) / (S(i,i)^2 + S(j,j)^2);
        end
    end
    Yf = Yf + Yf' - diag(diag(real(Yf)));
    elseif strcmp(method,'skewsymmetric')
    for i = 1:r
        Yf(i,i) = 1i*imag(C1(i,i))/S(i,i);
        for j = i+1:nx
            Yf(i,j) = (-S(i,i)*conj(C1(j,i)) + S(j,j)*(C1(i,j))) / (S(i,i)^2 + S(j,j)^2);
        end
    end
    Yf = Yf - Yf' - 1i*diag(diag(imag(Yf)));
    end

    A = @(v) Ux*Yf*(Ux'*v);

    if nargout==2
        varargout{1} = eig(Yf);
    elseif nargout>2
        [eVecs,eVals] = eig(Yf);
        eVals = diag(eVals);
        eVecs = Ux*eVecs;
        varargout{1} = eVals;
        varargout{2} = eVecs;
    end

    
elseif strcmp(method,'toeplitz') || strcmp(method,'hankel')
   if  strcmp(method,'toeplitz'); J = eye(nx); 
   elseif strcmp(method,'hankel'); J = fliplr(eye(nx)); end
    Am = fft([eye(nx) zeros(nx)].',[],1)'/sqrt(2*nx); % Define the left matrix
    B = fft([(J*X)' zeros(nt,nx)].',[],1)'/sqrt(2*nx); % Define the right matrix
    BtB = B'*B; 
    AAt = ifft(fft([eye(nx) zeros(nx); zeros(nx,2*nx)]).').'; % Fast computation of A*A'
    y = diag(Am'*conj(Y)*B)'; % Construct the RHS of the linear system
    L = (AAt.*BtB.')'; % Construct the matrix for the linear system
    d = [y(1:end-1)/L(1:end-1,1:end-1) 0]; % Solve the linear system
    newA = ifft(fft(diag(d)).').'; % Convert the eigenvalues into the circulant matrix
    A = newA(1:nx,1:nx)*J; % Extract the Toeplitz matrix from the circulant matrix
 
elseif startsWith(method,'circulant')

 fX = fft(X); fY = fft(conj(Y));

 d = zeros(nx,1);
    if endsWith(method,'TLS') % Solve in the total least squares sense     
        for j = 1:nx
        d(j) = tls(fX(j,:)',fY(j,:)');
        end
    elseif ~endsWith(method,'TLS') % Solve the other cases
    d = diag(fX*fY')./vecnorm(fX,2,2).^2;
        if endsWith(method,'unitary'); d = exp(1i*angle(d));
        elseif endsWith(method,'symmetric'); d = real(d);
        elseif endsWith(method,'skewsymmetric'); d = 1i*imag(d);
        end
    end
  eVals = d; % These are the eigenvalues
  eVecs = fft(eye(nx)); % These are the eigenvectors
 if nargin>3
    r = varargin{1}; % Rank constraint
    res = diag(abs(fX*fY'))./vecnorm(fX')'; % Identify least important eigenvalues
    [~,idx] = mink(res,nx-r); % Remove least important eigenvalues
    d(idx) = 0; eVals(idx) = []; eVecs(:,idx) = [];
 end

 if nargout>1; varargout{1} = eVals; end
 if nargout>2; varargout{2} = eVecs; end

 A = @(v) fft(d.*ifft(v)); % Reconstruct the operator in terms of FFTs

elseif strcmp(method,'BCCB') || strcmp(method,'BCCBtls') || strcmp(method,'BCCBskewsymmetric') || strcmp(method,'BCCBunitary')
    
    if isempty(varargin); error('Need to specify size of blocks.'); end
    s = varargin{1}; p = prod(s);
    % Equivalent to applying the block-DFT matrix F 
    % defined by F = kron(dftmtx(M),dftmtx(N)) to the 
    % matrix X
    aF =  @(x) reshape(     fft2(reshape(x ,[s,size(x,2)])) ,[p,size(x,2)])/sqrt(p);
    aFt = @(x) conj(aF(conj(x)));
    fX = aF(conj(X)); fY = aF(conj(Y));
    d = zeros(p,1);
    
    if strcmp(method,'BCCB') 
    for j = 1:p; d(j) = conj(fX(j,:)*fY(j,:)')/norm(fX(j,:)').^2; end
    elseif strcmp(method,'BCCBtls')
    for j = 1:p; d(j) = tls(fX(j,:)',fY(j,:)')'; end
    elseif strcmp(method,'BCCBskewsymmetric')
    for j = 1:p; d(j) = 1i*imag(fY(j,:)/fX(j,:)); end
    elseif strcmp(method,'BCCBsymmetric')
    for j = 1:p; d(j) = real(fY(j,:)/fX(j,:)); end
    elseif strcmp(method,'BCCBunitary')
    for j = 1:p; d(j) = exp(1i*angle(fY(j,:)/fX(j,:))); end
    end

    % Returns a function handle that applies A
     if nargin>4
        r = varargin{2};
        res = diag(abs(fX*fY'))./vecnorm(fX')';
        [~,idx] = mink(res,nx-r);
        d(idx) = 0;
    end
    A = @(x) aF((conj(d).*aFt(x)));
    varargout{1} = d;
    % Eigenvalues are given by d

elseif strcmp(method,'BC') || strcmp(method,'BCtri') || strcmp(method,'BCtls')
    
    s = varargin{1}; p = prod(s);
        M = s(2); N = s(1);
    if isempty(s); error('Need to specify size of blocks.'); end
    % Equivalent to applying the block-DFT matrix F 
    % defined by F = kron(dftmtx(M),eye(N)) to the 
    % matrix X
aF  =  @(x) reshape(fft(reshape(x,[s,size(x,2)]),[],2) ,[p,size(x,2)])/sqrt(M);
aFt =  @(x) conj(aF(conj(x)));

fX = aF(X); fY = aF(Y);
    d = cell(M,1);

for j = 1:M
    ls = (j-1)*N + (1:N);
    if strcmp(method,'BC')
        d{j} = fY(ls,:)/fX(ls,:);
    elseif strcmp(method,'BCtri')
        d{j} = piDMD(fX(ls,:),fY(ls,:),'diagonal',2);
    elseif strcmp(method,'BCtls')
        d{j} = tls(fX(ls,:)',fY(ls,:)')';
    end
end 

    BD = blkdiag(d{:});
    A = @(v) aFt(BD*aF(v));        
   
elseif strcmp(method,'symtridiagonal')
    
    T1e = vecnorm(X,2,2).^2; % Compute the entries of the first block
    T1 = spdiags(T1e,0,nx,nx); % Form the leading block
    T2e = dot(X(2:end,:),X(1:end-1,:),2); % Compute the entries of the second block
    T2 = spdiags([T2e T2e],-1:0,nx,nx-1); % Form the second and third blocks
    T3e = [0; dot(X(3:end,:),X(1:end-2,:),2)]; % Compute the entries of the final block
    T3 = spdiags(T1e(1:end-1) + T1e(2:end),0,nx-1,nx-1) ...
         + spdiags(T3e,1,nx-1,nx-1) + spdiags(T3e,1,nx-1,nx-1)'; % Form the final block
    T = [T1 T2; T2' T3]; % Form the block tridiagonal matrix
    d = [dot(X,Y,2); dot(X(1:end-1,:),Y(2:end,:),2) + dot(X(2:end,:),Y(1:end-1,:),2)]; % Compute the RHS vector
    c = real(T)\real(d); % Take real parts then solve linear system
    % Form the solution matrix
    A = spdiags(c(1:nx),0,nx,nx) + spdiags([0;c(nx+1:end)],1,nx,nx) + spdiags([c(nx+1:end); 0],-1,nx,nx);
else
    error('The selected method doesn''t exist.');

end
