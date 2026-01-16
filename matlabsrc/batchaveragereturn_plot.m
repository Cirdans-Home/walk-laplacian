%% Plot Average Return Probabilities Results
% This script has to be used in conjunction with the batchaveragereturn.m.
% The first one runs the computation on a node without X11 support and
% saves the results on 'batch-result-X.mat' files. This script loads such
% results and produces the figures that are used in the paper.
%
% To save the figures the code uses matlab2tikz. If it is not available the
% script will only plot them on screen.

close all; clear all; clc;

i = 4; % Change it for the test you want to plot
% Load the data from file
load(sprintf('batch-results-%d.mat',i));
n = length(l);  % Get the size of the matrix from the number of evs.
% Recompute avg to show all the dynamic: you may have to adapt it to the time
% interval you want to visualize.
t = linspace(0,1e6,300);
tlog = [t(1),logspace(-4,log10(t(2)),30)];
t = [tlog,t(3:end)];
avgret = @(l,t) sum(exp(-l*t))./length(l);
stdavg = arrayfun(@(t) avgret(l,t),t);
katzavg = arrayfun(@(t) avgret(lk,t),t);
expavg = arrayfun(@(t) avgret(lexp,t),t);
eig_exp_avg = arrayfun(@(t) avgret(eig_exp_L,t),t);
eig_resolv_avg = arrayfun(@(t) avgret(eig_resolv_L,t),t);
eig_lpath_avg_exp = arrayfun(@(t) avgret(eigLpath,t),t);
eig_lpath_avg_power = arrayfun(@(t) avgret(eigLpath_pow,t),t);
% Make the plots
if i == 1
    figure(i)
    loglog(t,real(stdavg),'k--',...
        t,real(katzavg),...
        t,real(expavg),'.-',...
        t,real(eig_resolv_avg),...
        t,real(eig_exp_avg),'.-',...
        t,real(eig_lpath_avg_power),'--',... %t,1/n*ones(size(t)),'k--',...
        t,real(eig_lpath_avg_exp),'--',...
        'LineWidth',2)
    legend({'$L$',...
        '$\mathbb{L}_1\left((1-\alpha x)^{-1}\right) \; \alpha = \nicefrac{1}{1 + \rho(Z)}$',...
        '$\mathbb{L}_1(\exp)$',...
        '$\mathbb{L}\left((1-\alpha x)^{-1}\right) \; \alpha = \nicefrac{1}{2\rho(A)}$',...
        '$\mathbb{L}(\exp)$',...
        '$\mathcal{L}, \; t_k = \nicefrac{1}{k}$',...
        '$\mathcal{L}, \; t_k = \exp(-k)$'},'Location','eastoutside','NumColumns',1);
    % ylabel(sprintf('%s',name)); % de-comment to put name on y-axis: now is on
    % the subcaption in the latex
    ylim([1/(1.1*n),1])
    xlim([t(2),t(end)]);
    % text(t(end),1/n,'1/|V|'); % de-comment to put limit back
    % Try using matlab2tikz to save the figure on file
elseif i == 2
    figure(i)
    loglog(t,real(stdavg),'k--',...
        t,real(katzavg),...
        t,real(expavg),'.-',...
        t,real(eig_resolv_avg),...
        t,[real(eig_exp_avg(1:21)),real(eig_exp_avg(21))*ones(1,length(t(22:end)))],'.-',...
        t,real(eig_lpath_avg_power),'--',... %t,1/n*ones(size(t)),'k--',...
        t,real(eig_lpath_avg_exp),'--',...
        'LineWidth',2)
    ylim([1/(1.1*n),1])
    xlim([t(2),t(end)]);
    % Try using matlab2tikz to save the figure on file
elseif i == 3 || i == 4
        figure(i)
        loglog(t,real(stdavg),'k--',...
            t,real(katzavg),...
            t,real(expavg),'.-',...
            t,real(eig_resolv_avg),...
            t,real(eig_exp_avg),'.-',...
            t,real(eig_lpath_avg_power),'--',... %t,1/n*ones(size(t)),'k--',...
            t,real(eig_lpath_avg_exp),'--',...
            'LineWidth',2)
        % ylabel(sprintf('%s',name)); % de-comment to put name on y-axis: now is on
        % the subcaption in the latex
        ylim([1/(1.1*n),1])
        xlim([t(2),t(end)]);
        % text(t(end),1/n,'1/|V|'); % de-comment to put limit back
        % Try using matlab2tikz to save the figure on file
end
try
    savefilename = [SafeName(name),'.tex'];
    matlab2tikz('filename',savefilename,'width','0.4\columnwidth',...
        'height','1.254in','extraAxisOptions','legend columns=1,','parseStrings',false);
catch
    warning('matlab2tikz is not responsive, do you have it?');
end


%% Make a filename safe
function name = SafeName( name )
name = matlab.lang.makeValidName(name,'ReplacementStyle','delete','Prefix','d');
end
