%% Visualize KatzResults

clear; clc; close all;

data_cpu = parseKatzSolutionFiles("katzsolution","*cpu*");
data_gpu = parseKatzSolutionFiles("katzsolution","*gpu*");

% Prepare a table with columns the matrix name, the value of alpha and the
% time to solve the system
alphavec = repmat(["\alpha_1";"\alpha_2";"\alpha_3";"\alpha_4"],[length(unique({data_cpu.matrix}')),1]);
resultsTableCPU = table({data_cpu.matrix}', [data_cpu.alpha]', [data_cpu.timetosolve]', alphavec, ...
    'VariableNames', {'MatrixName', 'Alpha', 'SolveTime', 'alphavec'});
resultsTableGPU = table({data_gpu.matrix}', [data_gpu.alpha]', [data_gpu.timetosolve]', alphavec, ...
    'VariableNames', {'MatrixName', 'Alpha', 'SolveTime', 'alphavec'});

% Make two heatmaps for the two dataset
figure(1);
subplot(1, 3, 1);
heatmap(resultsTableCPU,"alphavec","MatrixName",'ColorVariable',"SolveTime");
title('CPU (CSR) Solve Times (s)');
ylabel('Matrix Name');
xlabel('\alpha');
set(gca,"FontSize",16)
subplot(1, 3, 2);
heatmap(resultsTableGPU,"alphavec","MatrixName",'ColorVariable',"SolveTime");
title('GPU (CSRG) Solve Times (s)');
ylabel('');
xlabel('\alpha');
set(gca,"FontSize",16)
subplot(1,3,3); % Compute the speed-up
% Create a table with MatrixName Alpha and the speedup
speedUp = resultsTableCPU.SolveTime ./ resultsTableGPU.SolveTime;
resultsTableSpeedUp = table(resultsTableCPU.MatrixName, resultsTableCPU.Alpha, speedUp, alphavec,...
    'VariableNames', {'MatrixName', 'Alpha', 'SpeedUp','alphavec'});
heatmap(resultsTableSpeedUp,"alphavec","MatrixName",'ColorVariable',"SpeedUp");
title('Speed-Up Comparison');
ylabel('');
xlabel('\alpha');
set(gca,"FontSize",16)