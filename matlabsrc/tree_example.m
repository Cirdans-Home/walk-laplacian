%% Example of Diffusion on a randomly generated Tree
% On your system:
% python3 -m venv <path>
% source <path>/bin/activate
% pip install networkx numpy scipy
% Then on MATLAB, you run:
% pyenv('Version','<path>/bin/python3','ExecutionMode','OutOfProcess')

clear; clc; close all;

% Generate a tree with a given branching factor
n = 100;
b = 3;
rng(b); % For reproducibility (on the MATLAB side)

% Use the Python MATLAB interface to generate a random tree using the NetworkX library
[A, pos] = pyrunfile("generate_tree.py",["A","pos"],n_nodes=int32(n),seed=int8(b));
A = sparse(double(A));

posmat = zeros(n,2);
for i=0:n-1
    posmat(i+1,:) = double(pos.get(i));
end


% Convert the Python object to a graph
G = graph(A);

% Plot the tree
figure(1);
plot(G,'XData',posmat(:,1),'YData',posmat(:,2),'NodeLabel',{});

%% Build the Laplacian(s)
X0 = zeros(n,1);    % Initial state
X0(1) = 1;
%% We start from the classical Laplacian
I = speye(n,n);
L = G.laplacian();
D = spdiags(G.degree,0,n,n);
P = I - D\L;
mc_classical = dtmc(P);

numSteps1 = 20;
numSteps2 = 40;
numSteps3 = 80;
X20 = simulate(mc_classical,numSteps1,'X0',X0);
X40 = simulate(mc_classical,numSteps2,'X0',X0);
X80 = simulate(mc_classical,numSteps3,'X0',X0);
%% Visualize
figure(2)
% 20 Steps
subplot(4,3,1)
h20 = plot(G,'XData',posmat(:,1),'YData',posmat(:,2),'NodeLabel',{},'MarkerSize',3);
highlight(h20,X20(1),'NodeColor','blue','Marker','square');
highlight(h20,X20(2:end-1),'NodeColor',[0.9290 0.6940 0.1250],'Marker','*','MarkerSize',6);
highlight(h20,X20(end),'NodeColor','red','Marker','*','MarkerSize',6);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
ylabel('$L$','Interpreter','latex','FontSize',14)
xlabel(sprintf('Visited %d nodes',length(unique(X20))));
title('20 Steps')
% 40 Steps
subplot(4,3,2)
h40 = plot(G,'XData',posmat(:,1),'YData',posmat(:,2),'NodeLabel',{},'MarkerSize',3);
highlight(h40,X40(1),'NodeColor','blue','Marker','square');
highlight(h40,X40(2:end-1),'NodeColor',[0.9290 0.6940 0.1250],'Marker','*','MarkerSize',6);
highlight(h40,X40(end),'NodeColor','red','Marker','*','MarkerSize',6);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xlabel(sprintf('Visited %d nodes',length(unique(X40))));
title('40 Steps')
% 80 Steps
subplot(4,3,3)
h80 = plot(G,'XData',posmat(:,1),'YData',posmat(:,2),'NodeLabel',{},'MarkerSize',3);
highlight(h80,X80(1),'NodeColor','blue','Marker','square');
highlight(h80,X80(2:end-1),'NodeColor',[0.9290 0.6940 0.1250],'Marker','*','MarkerSize',6);
highlight(h80,X80(end),'NodeColor','red','Marker','*','MarkerSize',6);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xlabel(sprintf('Visited %d nodes',length(unique(X80))));
title('80 Steps')
%% The ne use the walk-exp-Laplacian
beta = 1;
EA = expm(beta*A);
%beta = 1/(abs(eigs(A,1,'largestabs')+1));
%EA = inv(I - beta*A);
EA(EA < 0) = 0;
DE = spdiags(sum(EA,2),0,n,n);
LE = DE - EA;
PE = I - DE\LE;

mc_walk_exp = dtmc(PE);

numSteps1 = 20;
numSteps2 = 40;
numSteps3 = 80;
X20E = simulate(mc_walk_exp,numSteps1,'X0',X0);
X40E = simulate(mc_walk_exp,numSteps2,'X0',X0);
X80E = simulate(mc_walk_exp,numSteps3,'X0',X0);

%% Visualize
figure(2)
% 20 Steps
subplot(4,3,4)
h20E = plot(G,'XData',posmat(:,1),'YData',posmat(:,2),'NodeLabel',{},'MarkerSize',3);
highlight(h20E,X20E(1),'NodeColor','blue','Marker','square');
highlight(h20E,X20E(2:end-1),'NodeColor',[0.9290 0.6940 0.1250],'Marker','*','MarkerSize',6);
highlight(h20E,X20E(end),'NodeColor','red','Marker','*','MarkerSize',6);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
ylabel('$\bf{L}(\exp)$','Interpreter','latex','FontSize',14)
xlabel(sprintf('Visited %d nodes',length(unique(X20E))));
% 40 Steps
subplot(4,3,5)
h40E = plot(G,'XData',posmat(:,1),'YData',posmat(:,2),'NodeLabel',{},'MarkerSize',3);
highlight(h40E,X40E(1),'NodeColor','blue','Marker','square');
highlight(h40E,X40E(2:end-1),'NodeColor',[0.9290 0.6940 0.1250],'Marker','*','MarkerSize',6);
highlight(h40E,X40E(end),'NodeColor','red','Marker','*','MarkerSize',6);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xlabel(sprintf('Visited %d nodes',length(unique(X40E))));
% 80 Steps
subplot(4,3,6)
h80E = plot(G,'XData',posmat(:,1),'YData',posmat(:,2),'NodeLabel',{},'MarkerSize',3);
highlight(h80E,X80E(1),'NodeColor','blue','Marker','square');
highlight(h80E,X80E(2:end-1),'NodeColor',[0.9290 0.6940 0.1250],'Marker','*','MarkerSize',6);
highlight(h80E,X80E(end),'NodeColor','red','Marker','*','MarkerSize',6);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xlabel(sprintf('Visited %d nodes',length(unique(X80E))));

%% PATH Laplacian

Dist = G.distances();
diam = max(Dist(:));
Lpath = sparse(n,n);
%Dpath = ones(n,1);
alpha = 1;
for i=1:diam
    [p,q] = find(Dist == i);
    e = ones(size(p));
    Gi = graph(sparse(p,q,e,n,n));
    Li = Gi.laplacian;
    Di = Gi.degree;
    %Dpath = Dpath + 1/(i^alpha)*Di;
    Lpath = Lpath + exp(-alpha*i)*Li;
end
Dpath = spdiags(diag(Lpath),0,n,n);
Ppath = I - Dpath\Lpath;


mc_walk_path = dtmc(Ppath);

numSteps1 = 20;
numSteps2 = 40;
numSteps3 = 80;
X20pp = simulate(mc_walk_path,numSteps1,'X0',X0);
X40pp = simulate(mc_walk_path,numSteps2,'X0',X0);
X80pp = simulate(mc_walk_path,numSteps3,'X0',X0);

%% Visualize
figure(2)
% 20 Steps
subplot(4,3,7)
h20pp = plot(G,'XData',posmat(:,1),'YData',posmat(:,2),'NodeLabel',{},'MarkerSize',3);
highlight(h20pp,X20pp(1),'NodeColor','blue','Marker','square');
highlight(h20pp,X20pp(2:end-1),'NodeColor',[0.9290 0.6940 0.1250],'Marker','*','MarkerSize',6);
highlight(h20pp,X20pp(end),'NodeColor','red','Marker','*','MarkerSize',6);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
ylabel('$\mathcal{L}$','Interpreter','latex','FontSize',14)
xlabel(sprintf('Visited %d nodes',length(unique(X20pp))));
% 40 Steps
subplot(4,3,8)
h40pp = plot(G,'XData',posmat(:,1),'YData',posmat(:,2),'NodeLabel',{},'MarkerSize',3);
highlight(h40pp,X40pp(1),'NodeColor','blue','Marker','square');
highlight(h40pp,X40pp(2:end-1),'NodeColor',[0.9290 0.6940 0.1250],'Marker','*','MarkerSize',6);
highlight(h40pp,X40pp(end),'NodeColor','red','Marker','*','MarkerSize',6);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xlabel(sprintf('Visited %d nodes',length(unique(X40pp))));
% 80 Steps
subplot(4,3,9)
h80pp = plot(G,'XData',posmat(:,1),'YData',posmat(:,2),'NodeLabel',{},'MarkerSize',3);
highlight(h80pp,X80pp(1),'NodeColor','blue','Marker','square');
highlight(h80pp,X80pp(2:end-1),'NodeColor',[0.9290 0.6940 0.1250],'Marker','*','MarkerSize',6);
highlight(h80pp,X80pp(end),'NodeColor','red','Marker','*','MarkerSize',6);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xlabel(sprintf('Visited %d nodes',length(unique(X80pp))));

%% NBT exp Laplacian

Lnbt_fun = buildphibasedlaplacian(G.adjacency);
I = eye(G.numnodes,G.numnodes);
Lnbt = Lnbt_fun('notransp',I);
Dnbt = diag(diag(Lnbt));
Pnbt = abs(I - Dnbt\Lnbt);

mc_nbt_exp = dtmc(Pnbt);

numSteps1 = 20;
numSteps2 = 40;
numSteps3 = 80;
X20nbt = simulate(mc_nbt_exp,numSteps1,'X0',X0);
X40nbt = simulate(mc_nbt_exp,numSteps2,'X0',X0);
X80nbt = simulate(mc_nbt_exp,numSteps3,'X0',X0);

%% Visualize
figure(2)
% 20 Steps
subplot(4,3,10)
h20pp = plot(G,'XData',posmat(:,1),'YData',posmat(:,2),'NodeLabel',{},'MarkerSize',3);
highlight(h20pp,X20nbt(1),'NodeColor','blue','Marker','square');
highlight(h20pp,X20nbt(2:end-1),'NodeColor',[0.9290 0.6940 0.1250],'Marker','*','MarkerSize',6);
highlight(h20pp,X20nbt(end),'NodeColor','red','Marker','*','MarkerSize',6);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
ylabel('$\bf{L}^{NBT}(\exp)$','Interpreter','latex','FontSize',14)
xlabel(sprintf('Visited %d nodes',length(unique(X20nbt))));
% 40 Steps
subplot(4,3,11)
h40pp = plot(G,'XData',posmat(:,1),'YData',posmat(:,2),'NodeLabel',{},'MarkerSize',3);
highlight(h40pp,X40nbt(1),'NodeColor','blue','Marker','square');
highlight(h40pp,X40nbt(2:end-1),'NodeColor',[0.9290 0.6940 0.1250],'Marker','*','MarkerSize',6);
highlight(h40pp,X40nbt(end),'NodeColor','red','Marker','*','MarkerSize',6);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xlabel(sprintf('Visited %d nodes',length(unique(X40nbt))));
% 80 Steps
subplot(4,3,12)
h80pp = plot(G,'XData',posmat(:,1),'YData',posmat(:,2),'NodeLabel',{},'MarkerSize',3);
highlight(h80pp,X80nbt(1),'NodeColor','blue','Marker','square');
highlight(h80pp,X80nbt(2:end-1),'NodeColor',[0.9290 0.6940 0.1250],'Marker','*','MarkerSize',6);
highlight(h80pp,X80nbt(end),'NodeColor','red','Marker','*','MarkerSize',6);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xlabel(sprintf('Visited %d nodes',length(unique(X80nbt))));


%% Spectral gap
figure(3)
subplot(1,3,1);
eigplot(mc_walk_path)
title('k-path Laplacian')
subplot(1,3,2)
eigplot(mc_walk_exp)
title('walk-based Laplacian')
subplot(1,3,3)
eigplot(mc_nbt_exp)
title('NBT Laplacian')

