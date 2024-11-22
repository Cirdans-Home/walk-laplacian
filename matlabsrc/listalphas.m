%% Produces the list of alphas and matrices for running the Fortran experiment

clear; clc; close all;

[~,files] = system('ls -v1 ../testmatrices/mm/*.mtx');
files = split(files);
files(end) = [];

kval = 4;       % Number of alphas to produce
check = false;  % Reduction to largest connected component has already 
                % been performed

j = 0;
for i=1:length(files)
    A = mmread(files{i});
    [alphavec,rhoz] = producealphas(A,kval,check);
    for k = 1:kval
        fprintf('A[%d,%d]=%1.8f\n',j,k-1,alphavec(k))
    end
    j = j + 1;
end