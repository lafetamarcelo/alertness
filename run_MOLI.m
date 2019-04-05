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

%% Run the nonlinear least squares approach

%par = est_nlLeast(dte, 'genetic algorithm');

%% Determine the model structure 

[struc,dt] = struc_select('sGolay',dte);

%% Determine the model parameters

parameters = est_regr(dt,struc,'16','ls');

%% Plotting some results

time_sam.init = [dt.t{1}(1), dt.init(2:end)];
time_sam.final = dt.final;
initial = dt.y{1}(1);

% simulate the alertness level for the estimate data
dts = sim_system(parameters,time_sam,initial);

figure(2); hold on;
for i = 1 : length(dts.y)
    subplot(3,1,1); hold on;
    plot(dts.td{i},dts.yd{i},'k.-','LineWidth',1.6);
    subplot(3,1,2); hold on;
    plot(dts.td{i},dts.yd{i} - dts.dhom{i},'k.-','LineWidth',1.6);
    subplot(3,1,3); hold on;
    plot(dts.td{i},dts.dhom{i},'k.-','LineWidth',1.6);
   
    if i ~= length(dts.y)
        subplot(3,1,1); hold on;
        plot(dts.tn{i},dts.yn{i},'b--','LineWidth',1.6);
        subplot(3,1,2); hold on;
        plot(dts.tn{i},dts.yn{i} - dts.nhom{i},'b--','LineWidth',1.6);
        subplot(3,1,3); hold on;
        plot(dts.tn{i},dts.nhom{i},'b--','LineWidth',1.6);    
    end
end

figure(3); hold on;
scatter(cell2mat(dte.t), cell2mat(dte.y), 'r', 'LineWidth', 1.2);
for i = 1 : length(dts.y)
    plot(dts.td{i},dts.yd{i},'k.-','LineWidth',1.6);
    if i ~= length(dts.y)
        plot(dts.tn{i},dts.yn{i},'b--','LineWidth',1.6);
    end
    plot([dt.init(i) dt.init(i)], [2, 16],...
                                    'Color', [.6 .6 .6], 'LineWidth', 1.2);
    plot([dt.final(i) dt.final(i)], [2, 16],...
                                    'Color', [.6 .6 .6], 'LineWidth', 1.2); 
end


%% Plotting some one step ahead predictions

dts = pred_system(parameters,dte);

figure(4); hold on;
for i = 1 : length(dts.t)
     scatter(dts.t{i}, dts.y{i}, 'r', 'LineWidth', 1.2);
     scatter(dte.t{i}, dte.y{i}, 'k', 'LineWidth', 1.4);
end
