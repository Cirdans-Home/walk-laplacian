function [T,AVT,ERRT] = xnystraceexp(A,m,type,t,nt,varargin)
%XNYSTRACE-EXP randomized average return probability estimate
%   A adjacency matrix of the graph
%   m number of matrix vector product to be used in the estimate
%   type type of dynamic to be studied, admissible types are:
%       - 'LAPLACIAN' Standard Laplacian dynamic y'(t) = L^T y(t)
%   t = [t0,t1] time interval for the estimate, if longer uses the value
%   nt number of time points in the [t0,t1] interval, if t is longer than 2
%       takes the length of t as value

% Allocate result vectors
if length(t) == 2
    T = linspace(t(1),t(2),nt)';
elseif length(t) > 2
    T = t(:);
    nt = length(t);
end
    
AVT = zeros(nt,1);
ERRT = zeros(nt,1);

% Build randomized rectangular matrix for sampling
n = size(A,1);
tic;
[Om,improved,vtype] = generate_test_matrix(n,m,'improved');
timerandom = toc;
fprintf('Generated m = %d random %s vectors (%1.2e s).\n',...
    m,vtype,timerandom);

% We build the Laplacian matrix of interest, depends on the type entry
% variable

switch upper(type)
    case 'LAPLACIAN'
        % Standard Laplacian
        e = ones(n,1);
        d = A*e;
        L = spdiags(d,0,n,n) - A;
    case 'KATZNBT'
        % Laplacian based on the Katz matrix
        alpha = varargin{1};
        L = buildkatzbasedlaplacian(A,alpha);
        AB = pencil(L,n,[]);
    otherwise
        error('Not known or not implemented dynamic on the network: %s',type);
end

% We use here the Rational-Krylov Toolbox to generate a Krylov subspace to
% approximate the matrix exponential
tic;
F = @(x) exp(-x);
if isnumeric(L)
    lmax = eigs(L,1,'largestabs');
else
    Leig = @(v) L(false,v);
    lmax = eigs(Leig,n,1,'largestabs');
end
lmax = abs(lmax);
Z = linspace(0,t(2)*lmax,max(round(t(2)*lmax),50));
[~,xi] = util_aaa(F,Z);
xi = xi.';
timerat = toc;
fprintf('Built rational approximation of e^{-x} on [%1.2f,%1.2f] with %d poles (%1.2e s)\n',...
    0,t(2)*lmax,length(xi),timerat);
tic;
if isnumeric(L)
    [V,K,H,out] = rat_krylov(L,Om,xi);
else
    [V,K,H,out] = rat_krylov(AB,Om,xi);
end
Kthin = K(:,out.column_deflation);
Hthin = H(:, out.column_deflation);
HK = Hthin/Kthin;
Om0 = V'*Om;
timekryl = toc;
ndefl = find(out.column_deflation == 0);
ndefl = length(ndefl);
fprintf('Rational Krylov application in %1.2e s (Deflated %d columns)\n',...
    timekryl,ndefl)
% Loop over the time step for Trace Estimation
tic;
for i=1:nt
    if T(i) == 0
        AVT(1) = 1;
        ERRT(1) = 0;
    else
        E = expm(-T(i)*HK)/n;
        Y = real(V*(E*Om0));
        [AVT(i),ERRT(i)] = xnystrace_helper(Y, Om, improved);
    end
end
timestimate = toc;
fprintf('Time of the estimation phase: %1.2e s.\n',timestimate);

end

function AB = pencil(L,n,P)
% Returns the pencil structure for the Rational Krylov method

matprod = @(rho,eta,x) rho*L(false,x) - eta*x;
solve = @(nu,mu,x) lin_solve(matprod,n,nu,mu,x,P);

AB.multiply = matprod;
AB.solve = solve;
AB.isreal = true;

end

function y = lin_solve(matprod,n,nu,mu,x,P)
    y = zeros(size(x));
    
    if ~isempty(P)
       [L,U] = ilu(nu*P - mu*speye(size(P)));    
       for i=1:size(x,2)
           [v,flag,relres,iter] = gmres( @(v) matprod(nu,mu,v),x(:,i),[],1e-6,min(500,n),L,U);
           y(:,i) = v;
	   fprintf('Flag %d relres %1.2e iter %d %d\n',flag,relres,iter);
       end
    else
       for i=1:size(x,2)
           [v,flag,relres,iter] = gmres( @(v) matprod(nu,mu,v),x(:,i),[],1e-6,min(500,n));
           y(:,i) = v;
	   fprintf('Flag %d relres %1.2e iter %d %d\n',flag,relres,iter);
       end
    end

end

function [Om,improved,type] = generate_test_matrix(n,m,default,varargin)
%GENERATE_TEST_MATRIX generates the random sampling test matrix for the
%XNYSTRACE trace estimator from the repository:
%   https://github.com/eepperly/XTrace
type = default;
if ~isempty(varargin)
    for i = 1:length(varargin)
        if ischar(varargin{i}) || isstring(varargin{i})
            type = varargin{i};
            break
        end
    end
end

improved = false;
if strcmp(type, 'improved')
    improved = true;
    Om = sqrt(n) * cnormc(randn(n,m));
elseif strcmp(type, 'cimproved')
    improved = true;
    Om = sqrt(n) * cnormc(randn(n,m) + 1i*randn(n,m));
elseif strcmp(type, 'rademacher') || strcmp(type, 'signs')
    Om = -3 + 2*randi(2,n,m);
elseif strcmp(type, 'steinhaus') || strcmp(type, 'csigns') ...
        || strcmp(type, 'phases')
    Om = exp(2*pi*1i*rand(n,m));
elseif strcmp(type, 'gaussian')
    Om = randn(n,m);
elseif strcmp(type, 'cgaussian')
    Om = 1/sqrt(2) * randn(n,m) + 1i/sqrt(2)*randn(n,m);
elseif strcmp(type, 'unif') || strcmp(type, 'sphere')
    Om = sqrt(n) * normc(randn(n,m));
elseif strcmp(type, 'cunif') || strcmp(type, 'csphere')
    Om = sqrt(n) * normc(randn(n,m)+1i*randn(n,m));
elseif strcmp(type, 'orth')
    Om = sqrt(n)*orth(randn(n,m));
elseif strcmp(type, 'corth')
    Om = sqrt(n)*orth(randn(n,m) + 1i*randn(n,m));
else
    error('"%s" not recognized as matrix type', type)
end
end


function [t,err] = xnystrace_helper(Y, Om, improved)
%XNYSTRACE trace estimator from the repository:
%   https://github.com/eepperly/XTrace
[n,m] = size(Y);
nu = eps*norm(Y,'fro')/sqrt(n);
Y = Y + nu*Om; % Shift for numerical stability
[Q,R] = qr(Y,0);
H = Om'*Y; C = chol((H+H')/2);
B = R/C; % Nystrom approx is Q*B*B'*Q'

%% Normalization
if improved
    [QQ,RR] = qr(Om,0);
    WW = QQ'*Om;
    warning('off','MATLAB:nearlySingularMatrix');
    SS = cnormc(inv(RR)');
    warning('on','MATLAB:nearlySingularMatrix');
    scale = (n - m + 1) ./ (n - vecnorm(WW).^2 ...
        + abs(diag_prod(SS,WW)' .* vecnorm(SS)).^2);
else
    scale = ones(1,m);
end

% Trace estimate
warning('off','MATLAB:nearlySingularMatrix');
W = Q'*Om; S = (B/C') .* (diag(inv(H))').^(-0.5);
warning('on','MATLAB:nearlySingularMatrix');
dSW = diag_prod(S, W).';
ests = norm(B,'fro')^2 - vecnorm(S).^2 + abs(dSW).^2 .* scale - nu*n;
t = mean(ests);
err = std(ests)/sqrt(m);
end

function d = diag_prod(A,B)
d = sum(conj(A).*B,1).'; % computes diag(A'*B)
end

function M = cnormc(M)
M = M ./ vecnorm(M,2,1);
end
