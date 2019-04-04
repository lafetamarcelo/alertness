%% Using KSS/PVT on MOLI

%load the data-set
clear; close; clc;
load('regres_01.mat');
addpath('./functions/');

figure(1); hold on;
for i = 1 : length(dte.t)
    plot(dte.t{i}, dte.y{i}, 'k-', 'LineWidth', 1.4);
    plot(dte.t{i}, dte.y{i}, 'k.', 'LineWidth', 1.6);
end
hold off;
%% Determine the model structure 

[struc,dt] = struc_select('trivial',dte);

%% Determine the model parameters

parameters = est_regr(dt,struc,'16','ls');

%% Plotting some results

time_sam.init = [dt.t{1}(1), dt.init(2:end)];
time_sam.final = dt.final;
initial = dt.y{1}(1);

% simulate the alertness level for the estimate data
dts = sim_system(parameters,time_sam,initial);

figure(1);
for i = 1 : 1%length(dts.y)
    subplot(3,1,1); hold on;
    plot(dts.td{i},dts.yd{i},'r--','LineWidth',1.6);
    subplot(3,1,2); hold on;
    plot(dts.td{i},dts.yd{i} - dts.dhom{i},'r--','LineWidth',1.6);
    subplot(3,1,3); hold on;
    plot(dts.td{i},dts.dhom{i},'r--','LineWidth',1.6);
   
    if i ~= length(dts.y)
        subplot(3,1,1); hold on;
        plot(dts.tn{i},dts.yn{i},'r--','LineWidth',1.6);
        subplot(3,1,2); hold on;
        plot(dts.tn{i},dts.yn{i} - dts.nhom{i},'r--','LineWidth',1.6);
        subplot(3,1,3); hold on;
        plot(dts.tn{i},dts.nhom{i},'r--','LineWidth',1.6);    
    end
end