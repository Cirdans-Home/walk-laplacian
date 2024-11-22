function L = buildkatzbasedlaplacian(A,alpha)
%%BUILDKATZBASEDLAPLACIAN takes as input the adjacency matrix A and gives
%back a function handle for the computation of the matrix vector product
%with the Katz based walk Laplacian
%

n = size(A,1);

e = ones(n,1);
d = A*e;
D = spdiags(d,0,n,n);
I = speye(n,n);


% phi-based non backtracking centralities
matvec = @(x) x -alpha*(A*x) - alpha^2*(x - D*x);
Li = ichol(I - alpha*A - alpha^2*(I - D));

[DL,~,~,iter] = pcg(matvec,e,1e-9,n,Li,Li');

L = @(t,x) katzmatvec(DL,A,D,n,Li,alpha,x,t);

fprintf('Number of Krylov steps         : %d\n',iter);

end