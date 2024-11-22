%% Example of a PATH graph with a trap

clear; clc; close all;

l = 5; % Length of the path (if it is even we add a node to have a center)
m = 8; % Number of edges on the central trap

%% Build and visualize the trap
if mod(l,2) == 0
    l = l + 1;
end

xpos = zeros(l+m,1);
ypos = zeros(l+m,1);
G = graph();
for i=1:l-1
    G = G.addedge(i,i+1,1);   % Build the path
    xpos(i) = i;
    ypos(i) = 0;
end
xpos(l) = l;
ypos(l) = 0;
halfnode = (l+1)/2;
for i=1:m
    G = G.addedge(halfnode,l+i,1);
    xpos(l+i) = halfnode + 0.5*cos(2*pi*i/(m+1));
    ypos(l+i) = 0 + 0.5*sin(2*pi*i/(m+1));
end
figure(1)
plot(G,'XData',xpos,'YData',ypos,'NodeLabel',{});

%% All the simulation start from the X0 node
X0    = zeros(l+m,1);
X0(1) = 1;
numsteps = 100*(l+m);

%% Build the Laplacian
L = G.laplacian;
D = spdiags(G.degree,0,G.numnodes,G.numnodes);
I = speye(G.numnodes,G.numnodes);
P = I - D\L;

mc_laplacian  = dtmc(P);
Xlap = simulate(mc_laplacian,numsteps,'X0',X0);
firstend = find(Xlap == l,1,'first');

figure(3)
subplot(1,4,1)
simplot(mc_laplacian,Xlap(1:firstend))
yticks([1,firstend])
yticklabels([1,firstend])

[pilap_two,lampi] = eigs(P',2,'largestabs');
ind = find(diag(lampi)>0);
pilap = pilap_two(:,ind);
pilap = pilap/sum(pilap);

%% Path-Laplacian
Dist = G.distances();
diam = max(Dist(:));
Lpath = sparse(G.numnodes,G.numnodes);
alpha = 1;
for i=1:diam
    [p,q] = find(Dist == i);
    e = ones(size(p));
    Gi = graph(sparse(p,q,e,G.numnodes,G.numnodes));
    Li = Gi.laplacian;
    Di = Gi.degree;
    Lpath = Lpath + exp(-alpha*i)*Li;
end
Dpath = spdiags(diag(Lpath),0,G.numnodes,G.numnodes);
Ppath = I - Dpath\Lpath;

mc_walk_path = dtmc(Ppath);

[piPath,lam] = eigs(Ppath',1,'largestabs');
piPath = piPath/sum(piPath);

Xpath = simulate(mc_walk_path,numsteps,'X0',X0);
firstend = find(Xpath == l,1,'first');
if ~isempty(firstend)
    figure(3)
    subplot(1,4,2)
    simplot(mc_walk_path,Xpath(1:firstend))
    yticks([1,firstend])
    yticklabels([1,firstend])
end
%% Walk based Laplacian
beta = 1/(2*abs(eigs(G.adjacency,1,'largestabs')));
EA = inv(I - beta*G.adjacency);
EA(EA < 0) = 0;
DE = spdiags(sum(EA,2),0,G.numnodes,G.numnodes);
LE = DE - EA;
PE = I - DE\LE;

mc_walk_exp = dtmc(PE);

Xwalk = simulate(mc_walk_exp,numsteps,'X0',X0);
firstend = find(Xwalk == l,1,'first');
if ~isempty(firstend)
    figure(3)
    subplot(1,4,4)
    simplot(mc_walk_exp,Xwalk(1:firstend))
    yticks([1,firstend])
    yticklabels([1,firstend])
end
[piwalk,lamwalk] = eigs(PE',1,'largestabs');
piwalk = piwalk/sum(piwalk);

%% Walk based Laplacian exp
beta = 1;
EA = expm(beta*full(G.adjacency()));
EA(EA < 0) = 0;
DE = spdiags(sum(EA,2),0,G.numnodes,G.numnodes);
LE = DE - EA;
PE = I - DE\LE;

mc_walk_exp = dtmc(PE);

Xwalk = simulate(mc_walk_exp,numsteps,'X0',X0);
firstend = find(Xwalk == l,1,'first');
if ~isempty(firstend)
    figure(3)
    subplot(1,3,3)
    simplot(mc_walk_exp,Xwalk(1:firstend))
    yticks([1,firstend])
    yticklabels([1,firstend])
end

[piwalkexp,lamwalkexp] = eigs(PE',1,'largestabs');
piwalkexp = piwalkexp/sum(piwalkexp);

%% EXP NBT Laplacian

Lnbt_fun = buildphibasedlaplacian(G.adjacency);
I = eye(G.numnodes,G.numnodes);
Lnbt = Lnbt_fun('notransp',I);
Dnbt = diag(diag(Lnbt));
Pnbt = abs(I - Dnbt\Lnbt);

mc_nbt_exp = dtmc(Pnbt);

Xnbt = simulate(mc_nbt_exp,numsteps,'X0',X0);
firstend = find(Xnbt == l,1,'first');
if ~isempty(firstend)
    figure(3)
    subplot(1,3,3)
    simplot(mc_nbt_exp,Xnbt(1:firstend))
    yticks([1,firstend])
    yticklabels([1,firstend])
end

[pinbtexp,lamnbtexp] = eigs(Pnbt',1,'largestabs');
pinbtexp = pinbtexp/sum(pinbtexp);

%% Visualize Steady States
figure(4)
subplot(5,2,[1;3;5;7;9])
labels = ["L","Lkatz","Lexp","Lpath","NBT (exp)"];
h = bar(labels,[pilap(l)/max(pilap),piwalk(l)/max(piwalk), ...
    piwalkexp(l)/max(piwalkexp),piPath(l)/max(piPath),...
    pinbtexp(l)/max(pinbtexp)],'FaceColor','flat');
h.CData(1,:) = [0.8500 0.3250 0.0980];
h.CData(2,:) = [0.4940 0.1840 0.5560];
h.CData(3,:) = [0.4660 0.6740 0.1880];
h.CData(4,:) = [0.6350 0.0780 0.1840];
h.CData(5,:) = [0.9529 0.2549 1.0000];
subplot(5,2,2)
hlap = bar(pilap,'FaceColor','flat');
ylim([0 0.4])
hlap.CData(l,:) = [0.8500 0.3250 0.0980];
title('Laplacian')
subplot(5,2,4)
hwalk = bar(piwalk,'FaceColor','flat');
ylim([0 0.4])
hwalk.CData(l,:) = [0.4940 0.1840 0.5560];
title('Walk Based Laplacian (Katz)')
subplot(5,2,6)
hwalk = bar(piwalkexp,'FaceColor','flat');
ylim([0 0.4])
hwalk.CData(l,:) = [0.4660 0.6740 0.1880];
title('Walk Based Laplacian (Exponential)')
subplot(5,2,8)
hpath = bar(piPath,'FaceColor','flat');
ylim([0 0.4])
hpath.CData(l,:) = [0.6350 0.0780 0.1840];
title('Transformed Path Laplacian (exp, \alpha = 1)')
subplot(5,2,10)
hnbt = bar(pinbtexp,'FaceColor','flat');
ylim([0 0.4])
hnbt.CData(l,:) = [0.9529 0.2549 1.0000];
title('NBT-exp Laplacian')

%% Visualize Steady States
figure(5)
subplot(5,2,[1;3;5;7;9])
labels = ["L","Lkatz","Lexp","Lpath","NBT (exp)"];
h = bar(labels,[pilap(l),piwalk(l), ...
    piwalkexp(l),piPath(l),...
    pinbtexp(l)],'FaceColor','flat');
h.CData(1,:) = [0.8500 0.3250 0.0980];
h.CData(2,:) = [0.4940 0.1840 0.5560];
h.CData(3,:) = [0.4660 0.6740 0.1880];
h.CData(4,:) = [0.6350 0.0780 0.1840];
h.CData(5,:) = [0.9529 0.2549 1.0000];
subplot(5,2,2)
hlap = bar(pilap,'FaceColor','flat');
hlap.CData(l,:) = [0.8500 0.3250 0.0980];
title('Laplacian')
subplot(5,2,4)
hwalk = bar(piwalk,'FaceColor','flat');
hwalk.CData(l,:) = [0.4940 0.1840 0.5560];
title('Walk Based Laplacian (Katz)')
subplot(5,2,6)
hwalk = bar(piwalkexp,'FaceColor','flat');
hwalk.CData(l,:) = [0.4660 0.6740 0.1880];
title('Walk Based Laplacian (Exponential)')
subplot(5,2,8)
hpath = bar(piPath,'FaceColor','flat');
hpath.CData(l,:) = [0.6350 0.0780 0.1840];
title('Transformed Path Laplacian (exp, \alpha = 1)')
subplot(5,2,10)
hpath = bar(pinbtexp,'FaceColor','flat');
hpath.CData(l,:) = [0.9529 0.2549 1.0000];
title('NBT-exp Laplacian')