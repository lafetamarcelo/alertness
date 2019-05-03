%load the data-set
clear; close all; clc;
addpath('./functions/');
addpath('./dataPacks/');

%% Using simulated KSS/PVT

W = autoGen(3);
[dte,dtn,A] = alertness_gen(2, W, 'random', 'not'); 

dtv.time = [[dte.t{1}(1); dte.init(2:end)'], dte.final'];
dtv.initial = dte.y{1}(1);

figure(1);
scatter(cell2mat(dte.t), cell2mat(dte.y), 'ko', 'LineWidth', 1.6);

% Include noise on the data


%% Using real KSS/PVT

load('regres_02.mat');

figure(1); hold on;
for i = 1 : length(dte.t)
   scatter(dte.t{i}, dte.y{i}, 'k', 'LineWidth', 1.4);
end
hold off;

dtv.time = [[dte.t{1}(1); dte.init(2:end)'], dte.final'];
dtv.initial = dte.y{1}(1);
%% Run the nonlinear least squares approach

preprocessing = 'None';
structure = 'Folkard'; 
algorithm = 'Force Compute'; 
identification = 'Simulated Annealing'; 
showResults = 'Sleepness'; % 'None', 'Sleepness', 'Validation'

model = alertness(preprocessing, structure, algorithm, identification, showResults);
model.fit(dte, dtv);

%% Evalute the model performance
dts = model.dts;
figure(2); hold on;
for i = 1 : length(dte.y)
    scatter(dte.t{i}, dte.y{i}, 'k', 'LineWidth', 1.4);
    plot(dts.dayTime{i}, dts.dayOut{i}, 'r-', 'LineWidth', 1.4);
    if i~= length(dte.y)
       plot(dts.nightTime{i}, dts.nightOut{i}, 'b-', 'LineWidth', 1.4); 
    end
end

figure(3); hold on;
for i = 1 : length(dte.y)
   subplot(2,1,1); hold on;
   plot(dts.dayTime{i}, dts.dayHom{i}, 'r-', 'LineWidth', 1.4);
   subplot(2,1,2); hold on;
   circ = dts.dayOut{i} - dts.dayHom{i};
   plot(dts.dayTime{i}, circ, 'r-', 'LineWidth', 1.4);
   if i ~= length(dte.y)
      subplot(2,1,1); hold on;
      plot(dts.nightTime{i}, dts.nightHom{i}, 'b-', 'LineWidth', 1.4);
      subplot(2,1,2); hold on;
      circ = dts.nightOut{i} - dts.nightHom{i};
      plot(dts.nightTime{i}, circ, 'b-', 'LineWidth', 1.4);
   end
end

%% 
%%
%%
%% Past algorithm


%% Determine the model structure 

[struc,dt] = struc_select('sGolay',dte);

%% Determine the model parameters

parameters = est_regr(dt,struc,'16','ls');

%% Plotting some results

time_sam.init = [dt.t{1}(1), dt.init(2:end)];
time_sam.final = dt.final;
initial = dt.y{1}(1);

%simulate the alertness level for the estimate data
dts = sim_system(parameters,time_sam,initial);

figure(3); hold on;
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

figure(4); hold on;
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
