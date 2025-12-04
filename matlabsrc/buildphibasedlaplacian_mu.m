function L = buildphibasedlaplacian_mu(A,mu)
%%BUILDPHIBASEDLAPLACIAN takes as input the adjacency matrix A and gives
%back a function handle for the computation of the matrix vector product
%with the phi-function based Laplacian downweighting the backtracking
%walsk with parameter \theta = 1 - \mu
%
% This routine needs acces to the 919 Algorithm.

n = size(A,1);

d = A*ones(n,1);
o = zeros(2*n,1);
e = [d;A*d - mu*d];
D = spdiags(d,0,n,n);
O = sparse(n,n);
I = speye(n,n);
Z = [O,I;mu*(mu*I-D),A];

% phi-based non backtracking centralities
[DL,stats] = phipm(1.0, Z, [o,e], 1e-12, false, 20);
DL = DL(1:n);

L = @(t,x) phimatvec_mu(DL,A,D,Z,n,o,x,t,mu);

fprintf('Number of substeps             : %d\n',stats(1));
fprintf('Number of rejected steps       : %d\n',stats(2));
fprintf('Number of Krylov steps         : %d\n',stats(3));
fprintf('Number of matrix exponentials  : %d\n',stats(4));

end