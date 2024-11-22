%% Plot Average Return Probabilities Results
% This script has to be used in conjunction with the batchaveragereturn.m.
% The first one runs the computation on a node without X11 support and
% saves the results on 'batch-result-X.mat' files. This script loads such
% results and produces the figures that are used in the paper.
%
% To save the figures the code uses matlab2tikz. If it is not available the
% script will only plot them on screen.


i = 6; % Change it for the test you want to plot
% Load the data from file
load(sprintf('batch-results-%d.mat',i));
n = length(l);  % Get the size of the matrix from the number of evs.
% Recompute avg to show all the dynamic: you may have to adapt it to the time
% interval you want to visualize.
t = linspace(0,1e6,300);
tlog = [t(1),logspace(-8,log10(t(2)),30)];
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
figure(i)
loglog(t,real(stdavg),...
t,real(katzavg),...
t,real(expavg),'.-',...
t,real(eig_exp_avg),'.-',...
t,real(eig_resolv_avg),...
t,real(eig_lpath_avg_exp),'--',...
t,real(eig_lpath_avg_power),'--',...
t,1/n*ones(size(t)),'k--',...
'LineWidth',2)
legend({'Laplacian','(NBT) Katz Laplacian','(NBT) Exp Laplacian',...
'Exp Laplacian','Katz Laplacian','Path Laplacian (exp)',...
'Path Laplacian (power)'},'Location','northoutside','NumColumns',2);
ylabel(sprintf('%s',name));
ylim([1/(1.1*n),1])
xlim([t(2),t(end)]);
text(t(end),1/n,'1/|V|');
% Try using matlab2tikz to save the figure on file
try
savefilename = [SafeName(name),'.tex'];
matlab2tikz('filename',savefilename,'width','0.6\columnwidth',...
    'height','1.5in','extraAxisOptions','legend columns=2,');
catch
warning('matlab2tikz is not responsive, do you have it?');
end

