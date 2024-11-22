function [alphavec,rhoz] = producealphas(A,k,check)
%PRODUCEALPHAS Prints the value of alpha to be used in the Fortran test
%for evaluating the linear system solution routine.

if ~exist('check','var')
    check = false;
end

if check
    % If it is required restricts the graph to the largest connected
    % component
    G = graph(A);
    [bin,binsize] = conncomp(G);
    idx = binsize(bin) == max(binsize);
    SG = subgraph(G, idx);
    A = SG.adjacency();
end


% Building the 2x2 block matrix
n = size(A,1);
D = spdiags(sum(A,2),0,n,n);
% Eigenvalue computation
Z = @(x) zmatvec(A,D,x,false);
rhoz = eigs(Z,2*n,1,"largestabs");
logvec = logspace(-3,-0.5,k);
alphavec = logvec/rhoz;

end

