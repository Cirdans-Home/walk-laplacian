function L = kwalklaplacian(A,f,varargin)
%KWALKLAPLACIAN Compute the k-walk Laplacian of a graph for a given function f. Specifically, it returns a function handle
% which computes the matrix-vector product of the k-walk Laplacian of a graph with a given vector. The returned routines takes
% two inputs: the vector to be multiplied and an optional arguments that specifies if you need the product with L or with L'.

% Parse optional arguments
% Argument to parse are the following:
% - 'kmax': maximum size of the Krylov subspace
% - 'method': method to compute the k-walk Laplacian. It can be either 'arnoldi' or 'lanczos'
% - 'reorth': boolean that specifies if the basis should be reorthogonalized
% - 'nstep': number of steps after which the basis should be reorthogonalized
% - 'tol': tolerance for the computation of the matrix functions
p = inputParser;
addParameter(p,'kmax',100,@(x) x > 0);
addParameter(p,'method','arnoldi',@(x) ismember(x,{'arnoldi','lanczos'}));
addParameter(p,'reorth',false,@islogical);
addParameter(p,'nstep',50,@(x) x > 0);
addParameter(p,'tol',1e-9,@(x) x > 0);
parse(p,varargin{:});
kmax = p.Results.kmax;
method = p.Results.method;
reorth = p.Results.reorth;
nstep = p.Results.nstep;
tol = p.Results.tol;

% Get the size of the graph
n = size(A,1);
e = ones(n,1);

% Check that f is a string representing a viable matrix function
switch upper(f)
    case "EXP"
        % This is the matrix exponential and is treated separately using the expmv routine from Matlab
        D = spdiags(expmv(A,e),0,n,n);
        L = @(x,transp) expmv_lap_apply(D,A,x,transp);
        return
    case "COS"
        fun = @(x) funm(x,@cos);
    case "SIN"
        fun = @(x) funm(x,@sin);
    case "SQRT"
        fun = @(x) sqrtm(x);
    case "LOG"
        fun = @(x) funm(x,@log);
    case "SINH"
        fun = @(x) funm(x,@sinh);
    case "COSH"
        fun = @(x) funm(x,@cosh);
    otherwise
        error('f must be a string representing a viable matrix function')
end
% Compute the diagonal matrix D as diag(f(A)e) employing Krylov method
d = run_krylov(fun,A,e,kmax,method,reorth,nstep,tol,'notransp');

% Return the function handle
L = @(x,transp) d.*x - run_krylov_block(fun,A,x,kmax,method,reorth,nstep,tol,transp);



end

function y = expmv_lap_apply(D,A,x,transp)
% Apply the matrix-vector product of the k-walk Laplacian

y = zeros(size(x));
switch transp
    case 'notransp'
        for j = 1:size(x,2)
            y(:,j) = D*x(:,j) - expmv(A,x(:,j));
        end
    case 'transp'
        for j = 1:size(x,2)
            y(:,j) = D*x(:,j) - expmv(A',x(:,j));
        end
    case 'real'
        if isreal(A)
            y = 0;
        else
            y = 1;
        end
    otherwise
        error('transp must be either notransp or transp');
end

end


function f = run_krylov_block(fun,A,b,kmax,method,reorth,nstep,tol,transp)
    if size(b,2) > 1
        f = zeros(size(b));
        for j=1:size(b,2)
            f(:,j) = run_krylov(fun,A,b(:,j),kmax,method,reorth,nstep,tol,transp);
        end
    else
        f = run_krylov(fun,A,b,kmax,method,reorth,nstep,tol,transp);
    end
end

function f = run_krylov(fun,A,b,kmax,method,reorth,nstep,tol,transp)
% Create a Krylov subspace object
switch transp
    case 'notransp'
        krylov = KrylovSubspace(A,b,kmax);
    case 'transp'
        if strcmpi(method,'lanczos')
            % If the method is Lanczos, we don't need to use the transpose of A, since A is symmetric
            krylov = KrylovSubspace(A,b,kmax);
        else
            krylov = KrylovSubspace(A',b,kmax);
        end
    case 'real'
        if isreal(A)
            y = 0;
        else
            y = 1;
        end
    otherwise
        error('transp must be either notransp or transp')
end

for i=1:kmax
    % Extend the Krylov subspace by one vector using the Arnoldi process
    switch method
        case 'arnoldi'
            krylov = arnoldi(krylov);
        case 'lanczos'
            krylov = lanczos(krylov);
        otherwise
    end
    % Reorthogonalize the basis
    if reorth && mod(i,nstep) == 0
        krylov = reorth(krylov);
    end
    % Approximate the matrix function
    fvec = fun(krylov.H(1:krylov.k-1,1:krylov.k-1));
    % Estimate the error
    err_estimate = abs(krylov.nb*fvec(i,1)*krylov.H(i+1,i));
    if err_estimate < tol
        break
    end
end
% Compute the matrix-vector product
f = krylov.nb*krylov.V(:,1:size(fvec,1))*(fvec(:,1));

end