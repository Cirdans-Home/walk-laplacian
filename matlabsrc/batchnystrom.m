%% Test of Randomized Approximation of Average Return Probabilities

clear; clc; close all;
addpath('rktoolbox/');

load('../testmatrices/mat/roadNet-CA.mat');

A = Problem.A;
G = graph(A,"omitselfloops");
A = G.adjacency("weighted");

t = [0,100];
nt = 90;
m = 30;
type = 'laplacian';
[~,rhoz] = producealphas(G.adjacency,1,false);
alpha = 1/(rhoz+1);
nrepetition = 10;
T = cell(nrepetition,1);
AVT = cell(nrepetition,1);
ERRT = cell(nrepetition,1);
parfor i = 1:nrepetition
    fprintf('|---------------------------------------------------------------------|\n');
    fprintf('| Executing Repetition %d                                             |\n\n',i);
    fprintf('|---------------------------------------------------------------------|\n');
    [T{i},AVT{i},ERRT{i}] = xnystraceexp(A,m,type,t,nt,alpha);
end
save('batchnystrom-laplace-CA.mat','T','AVT','ERRT');


t = [0,100];
nt = 90;
m = 30;
type = 'katznbt';
[~,rhoz] = producealphas(G.adjacency,1,false);
alpha = 1/(rhoz+1);
nrepetition = 10;
T = cell(nrepetition,1);
AVT = cell(nrepetition,1);
ERRT = cell(nrepetition,1);
parfor i = 1:nrepetition
    fprintf('|---------------------------------------------------------------------|\n');
    fprintf('| Executing Repetition %d                                             |\n\n',i);
    fprintf('|---------------------------------------------------------------------|\n');
    [T{i},AVT{i},ERRT{i}] = xnystraceexp(A,m,type,t,nt,alpha);
end
save('batchnystrom-laplace-katz-CA.mat','T','AVT','ERRT');
