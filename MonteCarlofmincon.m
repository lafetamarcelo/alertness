%% Monte Carlo of GoHorse Optimization - fmincon(.)
%% Default config
%close all;
clear;
clc;
set(0,'defaultfigurecolor',[1 1 1]);
%% Initialize variables and set the Monte Carlo configurations
MCruns = 200;

par = cell(MCruns,1);
J = cell(MCruns,1);
e_flag = cell(MCruns,1);
cruze_par = cell(MCruns,1);

%% Generate the data or load pre-saved ones
%% Define the number of information data
ppH = 2;
resolution = 10000;
space = 'random';
noise = 'colored';
SNR = 40; %dB

load('monte_data.mat');

% W = autoGen(35);
% 
% ind = cell(length(W(:,1)),1);
% for i = 1 : length(W(:,1))
% 
% if strcmp(space,'random')
%     ppW = round(W(i,2)*ppH,0);
%     cruze = .5*rand(ppW,1)./ppH - 1/ppH;
%     i_cruze = round(cruze*resolution/W(i,2),0);
%     
%     ind{i} = unique(sort(abs(round(linspace(1,resolution,ppW),0) + i_cruze')));
% elseif strcmp(space,'linear')
%     ind{i} = unique(round(linspace(1,resolution,ppW),0)); 
% end
% 
% if ind{i}(end) > resolution; ind{i}(end) = resolution; end
% if ind{i}(1) <= 0; ind{i}(1) = 1; end
% end

%% Optimization setup
options = optimoptions('fmincon','Algorithm','sqp','Display','off','MaxIter',50,...
    'UseParallel',true);
prop = .2;

%% Introducing Constrains
lb = [.1, .1, .1, -2*pi, .1, .1, .1];
ub = [pi, inf, inf, 2*pi, inf, inf, inf];

%% Real parameters
rPar(1) = pi/12; %omega
rPar(2) = 1/0.0353; %tau
rPar(3) = 2.52; %M
rPar(4) = -16.835*pi/12; %cPhase
rPar(5) = 2.4; %offset
rPar(6) = 14.3; %y_oo
rPar(7) = 2.6247; %tau_e

TexNames{1} = '$$\omega_n$$'; 
TexNames{2} = '$$\tau$$';
TexNames{3} = '$$M$$'; 
TexNames{4} = '$$\phi_{c}$$'; 
TexNames{5} = '$$A_{DC}$$'; 
TexNames{6} = '$$y_{\infty}$$';
TexNames{7} = '$$\tau_{e}$$';

%% Initialize the Monte Carlo Simulations
for mc = 1 : MCruns
    
    rng('shuffle');
    
    [dtd,dta] = alertness_sim(W,noise,SNR,resolution,ind);
    
    SNR_(mc,:) = mean(dta.SNR(:,1:4));
    
    est_d = 7;
    
    dte.y = dtd.y(1:est_d); dte.t = dtd.t(1:est_d);
    dte.init = dtd.initial(1:est_d); dte.final = dtd.final(1:est_d); 
    dte.valid.y = dtd.y(est_d+1:end); dte.valid.t = dtd.t(est_d+1:end);
    dte.valid.init = dtd.initial(est_d+1:end); dte.valid.final = dtd.final(est_d+1:end); 
    cruze_par{mc} = prop.*rPar.*(rand(1,length(rPar))-.5);
    
    init = rPar + cruze_par{mc};
   
    [par{mc},J{mc},e_flag{mc}] = fmincon(@ghCost,init,[],[],[],[],lb,ub,[],options,dte);
    
    disp(['Monte Carlo run ',num2str(mc),...
                             ' Quadratic Error: ',num2str(J{mc})]);
end

%% Monte Carlo figures
%% Plot the optimization output flags propotions and the cost function results 
labels = {'Maximum constrains violation','No feasible solution','Stopped by output funtion',...
    'Iterations exceeded','Conveged to selected options','Step Tolerance stopped',...
    'Cost function variation tolerance achieved', 'Step Tolerance to small',...
    'Derivative direction magnitude achieved'};
flags = cell2mat(e_flag)'+4;
e_flags = [count(num2str(flags),'1'), count(num2str(flags),'2'), count(num2str(flags),'3'),...
           count(num2str(flags),'4'), count(num2str(flags),'5'), count(num2str(flags),'6'),...
           count(num2str(flags),'7'), count(num2str(flags),'8'), count(num2str(flags),'9')];

figure(2); hold on;
subplot(2,1,1); 
histogram(cell2mat(J),20);
subplot(2,1,2); 
pie(e_flags,[0, 0, 0, 0, 1, 0, 0, 0, 0],labels); hold off;

%% Generate the visualization from for the cruze signal variations
figure(5); hold on;
%title('Parameters initialization cruzing','Interpreter','latex');
aux = cell2mat(cruze_par);
for i = 1 : 7
   subplot(4,2,i); hold on;
   histogram(aux(:,i)+rPar(i),20); hold on;
   title(['$$\overline{x} =',num2str(round(mean(aux(:,i)+rPar(i)),3)),...
       '\pm',num2str(round(std(aux(:,i)+rPar(i)),3)),'$$'],...
       'Interpreter','latex');
   leg = TexNames{i};
   legend(leg,'Interpreter','latex'); hold off;
end

%% Plot the proportional error for each parameter at each run
figure(12); hold on;
aux_i = cell2mat(par);
for i = 1 : 7
   subplot(4,2,i); hold on;
   e_prop = 100.*(rPar(i) - aux_i(:,i))/rPar(i); 
   histogram(e_prop,20); hold on;
   title(['$$\overline{x} =',num2str(round(mean(e_prop),3)),...
       '\pm',num2str(round(std(e_prop),3)),'$$'],'Interpreter','latex');
   leg = TexNames{i};
   legend(leg,'Interpreter','latex'); hold off; 
end

