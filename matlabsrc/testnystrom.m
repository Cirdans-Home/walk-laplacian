%% Test Nystrom Average Return Probability Estimate

clear; clc; close all;
load('../testmatrices/mat/karate.mat');

A = Problem.A;
G = graph(A,"omitselfloops");
A = G.adjacency("weighted");

%% Stochastic Estimator
avgret = @(l,t) sum(exp(-l*t))./length(l);
%% Standard Laplacian
t = [0,10];
nt = 30;
m = 4;
type = 'Laplacian';

[T,AVT,ERRT] = xnystraceexp(A,m,type,t,nt);
% Check for the error
n = size(A,1);
e = ones(n,1);
L = spdiags(A*e,0,n,n) - A;
ev = eig(full(L));
avret = arrayfun(@(T) avgret(ev,T),T);

% Visualization

figure(1)
subplot(1,2,1)
plot(T,AVT,'r-',...
    T,AVT+ERRT,'r--',...
    T,AVT-ERRT,'r--',...
    'LineWidth',2)
legend('XNysTrace-exp')
subplot(1,2,2)
semilogy(T,AVT,'r-',T,avret,'k-','LineWidth',2);
legend('XNysTrace-exp','True value')

%% Katz based Laplacian
t = [0,10];
nt = 30;
m = 3;
[~,rhoz] = producealphas(G.adjacency,1,false);
alpha = 1/(rhoz+1);
type = 'Katznbt';


[Tkatznbt,AVTkatznbt,ERRTkatznbt] = xnystraceexp(A,m,type,t,nt,alpha);

% Check for the error
n = size(A,1);
e = ones(n,1);
Lkatz = buildkatzbasedlaplacian(G.adjacency,alpha);
lk = eigs(@(x) Lkatz(false,x),G.numnodes,G.numnodes);
lk(end) = 0;
katzavg = arrayfun(@(t) avgret(lk,t),Tkatznbt);

% Visualization

figure(2)
subplot(1,2,1)
plot(Tkatznbt,AVTkatznbt,'r-',...
    Tkatznbt,AVTkatznbt+ERRT,'r--',...
    Tkatznbt,AVTkatznbt-ERRT,'r--',...
    'LineWidth',2)
legend('XNysTrace-exp')
subplot(1,2,2)
semilogy(Tkatznbt,AVTkatznbt,'r-',T,katzavg,'k-','LineWidth',2);
legend('XNysTrace-exp','True value')