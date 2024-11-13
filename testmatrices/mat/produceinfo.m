%% Table with Network info

clear; close all; clc;

[~,files] = system('ls -v1 *.mat');
files = split(files);
files(end) = [];

for j=1:length(files)
    load(files{j});
    % Reduce to largest connected component
    G = graph(Problem.A);
    [bin,binsize] = conncomp(G);
    idx = binsize(bin) == max(binsize);
    SG = subgraph(G, idx);
    A = SG.adjacency();
    % Building the 2x2 block matrix
    n = size(A,1);
    O = sparse(n,n);
    I = speye(n,n);
    D = spdiags(sum(A,2),0,n,n);
    Z = [O I;I-D A];
    % Eigenvalue computation
    rhoz = eigs(Z,1,"largestabs");

    fprintf(" %s & %d & %1.2f \n", ...
        Problem.name,size(A,1),rhoz);

end