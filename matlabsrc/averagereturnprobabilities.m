%% Average Return Probabilities

% Clear the environment and call 919
clear; clc; close all;
addpath('919/Matlab/Src/')

t = linspace(0,50,100);
t = [t(1),logspace(-5,log10(t(2)),50),t(3:end)];
avgret = @(l,t) sum(exp(-l*t))./length(l);

matrixfiles = {'../testmatrices/mat/USpowerGrid.mat', ...
	       '../testmatrices/mat/ca-CondMat.mat', ...
	       '../testmatrices/mat/dictionary28.mat', ...
	       '../testmatrices/mat/usroads-48.mat', ...
	       '../testmatrices/mat/human_gene1.mat'}

for II = 1:length(matrixfiles)

    filename = sprintf("batch-results-%d.mat",II);

    %% Load matrix
    load(matrixfiles{II});
    %% Who are we running on?
    fprintf("Matrix: %s\n",Problem.name);
    name = Problem.name;
    save(filename,"t","name")

    fprintf("Transforming into a connected graph\n");
    G = graph(Problem.A,"omitselfloops");
    [bin,binsize] = conncomp(G);
    idx = binsize(bin) == max(binsize);
    G = subgraph(G, idx);
    A = G.adjacency();
    e = ones(size(A,1),1);
    fprintf("Staring computations\n");

    %% Standard Laplacian
    fprintf('Standard Laplacian start\n')
    L = laplacian(G);
    l = eigs(L,size(L,1));
    stdavg = arrayfun(@(t) avgret(l,t),t);
    fprintf('Standard Laplacian end\n')
    save(filename,"l","stdavg","-append")

    %% Katz walk-based Laplacian
    fprintf('Katz (NBT) walk-based Laplacian start\n')
    [~,rhoz] = producealphas(G.adjacency,1,false);
    alpha = 1/(rhoz+1);
    Lkatz = buildkatzbasedlaplacian(G.adjacency,alpha);
    lk = eigs(@(x) Lkatz(false,x),G.numnodes,G.numnodes);
    lk(end) = 0;
    katzavg = arrayfun(@(t) avgret(lk,t),t);
    fprintf('Katz walk-based Laplacian end\n')
    save(filename,"lk","katzavg","-append")

    %% Exponential walk-based Laplacian
    fprintf('Exponential (NBT) walk-based Laplacian start\n')
    Lexp = buildphibasedlaplacian(G.adjacency);
    lexp = eigs(@(x) Lexp(false,x),G.numnodes,G.numnodes);
    lexp(end) = 0;
    expavg = arrayfun(@(t) avgret(lexp,t),t);
    fprintf('Exponential walk-based Laplacian end\n')
    save(filename,"lexp","expavg","-append")

    %% Walk-EXP Laplacian
    fprintf('Exponential walk-based Laplacian start\n')
    expA = expm(full(A));
    expL = diag(expA*e) - expA;
    eig_exp_L = eig(expL);
    eig_exp_avg = arrayfun(@(t) avgret(eig_exp_L,t),t);
    save(filename,"eig_exp_L","eig_exp_avg","-append")
    fprintf('Exponential walk-based Laplacian end\n')

    %% Katz Laplacian
    fprintf('Exponential walk-based Laplacian start\n')
    rhoA = abs(eigs(A,1,'largestabs'));
    alpha = 1/(2*rhoA);
    I = speye(size(A));
    resolvA = (I - alpha*A)\I;
    resolvL = diag(resolvA*e) - resolvA;
    eig_resolv_L = eig(full(resolvL));
    eig_resolv_avg = arrayfun(@(t) avgret(eig_resolv_L,t),t);
    save(filename,"eig_resolv_L","eig_resolv_avg","-append")
    fprintf('Exponential walk-based Laplacian end\n')
    
    %% Path Laplacian (Exp)
    fprintf('Exponential Path Laplacian start\n')
    Lpath = pathlap(G,'exp',1);
    eigLpath = eig(full(Lpath));
    eig_lpath_avg_exp = arrayfun(@(t) avgret(eigLpath,t),t); 
    save(filename,"eig_lpath_avg_exp","eigLpath","-append")
    fprintf('Exponential Path Laplacian end\n')
    
    %% Path Laplacian (Fractional)
    fprintf('Power Path Laplacian start\n')
    Lpath = pathlap(G,'power',1);
    eigLpath_pow = eig(full(Lpath));
    eig_lpath_avg_power = arrayfun(@(t) avgret(eigLpath_pow,t),t); 
    save(filename,"eig_lpath_avg_power","eigLpath_pow","-append")
    fprintf('Power Path Laplacian end\n')
    
end
