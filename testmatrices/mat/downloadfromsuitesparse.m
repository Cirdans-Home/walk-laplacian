%% This script download the test-matrices from Suite-Sparse

clear; clc; close all;

url = 'https://suitesparse-collection-website.herokuapp.com/mat/Pajek/dictionary28.mat';
filename = 'dictionary28.mat';
websave(filename,url);

url = 'https://suitesparse-collection-website.herokuapp.com/mat/Pajek/USpowerGrid.mat';
filename = 'USpowerGrid.mat';
websave(filename,url);

url = 'https://suitesparse-collection-website.herokuapp.com/mat/Belcastro/human_gene1.mat';
filename = 'human_gene1.mat';
websave(filename,url);

url = 'https://suitesparse-collection-website.herokuapp.com/mat/Belcastro/human_gene2.mat';
filename = 'human_gene2.mat';
websave(filename,url);

url = 'https://suitesparse-collection-website.herokuapp.com/mat/SNAP/ca-CondMat.mat';
filename = 'ca-CondMat.mat';
websave(filename,url);

url = 'https://suitesparse-collection-website.herokuapp.com/mat/SNAP/ca-HepPh.mat';
filename = 'ca-HepPh.mat';
websave(filename,url);

url = 'https://suitesparse-collection-website.herokuapp.com/mat/SNAP/roadNet-CA.mat';
filename = 'roadNet-CA.mat';
websave(filename,url);

url = 'https://suitesparse-collection-website.herokuapp.com/mat/SNAP/roadNet-PA.mat';
filename = 'roadNet-PA.mat';
websave(filename,url);

url = 'https://suitesparse-collection-website.herokuapp.com/mat/SNAP/roadNet-TX.mat';
filename = 'roadNet-TX.mat';
websave(filename,url);

url = 'https://suitesparse-collection-website.herokuapp.com/mat/Gleich/usroads-48.mat';
filename = 'usroads-48.mat';
websave(filename,url);

url = 'https://suitesparse-collection-website.herokuapp.com/mat/Arenas/PGPgiantcompo.mat';
filename = 'PGPgiantcompo.mat';
websave(filename,url);

url = 'https://suitesparse-collection-website.herokuapp.com/mat/DIMACS10/belgium_osm.mat';
filename = 'belgium_osm.mat';
websave(filename,url);

% Small networks
url = 'https://suitesparse-collection-website.herokuapp.com/mat/Newman/karate.mat';
filename = 'karate.mat';
websave(filename,url);