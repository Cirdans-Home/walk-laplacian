%% Convert matrices from SuiteSparse to Matrix-Market

clear; clc; close all hidden;

[~,files] = system('ls -v1 ../mat/*.mat');
files = split(files);
files(end) = [];


% h = waitbar(0,"Converting matrix");
for j=1:length(files)
    load(files{j});
    % Reduce to largest connected component
    G = graph(Problem.A);
    [bin,binsize] = conncomp(G);
    idx = binsize(bin) == max(binsize);
    SG = subgraph(G, idx);
    A = SG.adjacency();
    % Determine filename
    filename = sscanf(files{j},'../mat/%s.mat');
    filename = filename(1:end-4);
    filename = [filename,'.mtx'];
    % Write on file
    mmwrite(filename,A,Problem.name);
    % Update waitbar
    % waitbar(j/length(files),h,"Converting matrix");
end
