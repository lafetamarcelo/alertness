%% preTest 6

clear;
close all;
clc;
set(0,'defaultfigurecolor',[1 1 1]);
rng('shuffle');

%%
%  Generate with GenAlertness a linear step spaced data, with high resolution
%  for (_i.e._ 80 points per hour), as a reference guide line for further 
%  analysis

WL = autoGen(21);

[dtd,dtn,A] = alertness_gen(1,WL,'random','not');
%load('data.mat')
%% Include noise to data

noise_nat = 'none';
SNR = 20;
for i = 1 : length(dtd.y)
    if strcmp(noise_nat,'white')
        dtd.y{i} = awgn(dtd.y{i},SNR);
    elseif strcmp(noise_nat,'colored')
        e = filter(1,[1 -.9],randn(length(dtd.y{i}),1));
        v = e*std(cell2mat(dtd.y))*10^(-SNR/20)/std(e);
        e = v;
        dtd.y{i} = dtd.y{i} + e; 
    end
end

figure(1)
plot(A.t,A.A,'-','Color',[.6 .6 .6]); hold on;
scatter(cell2mat(dtd.t),cell2mat(dtd.y),'r','LineWidth',1.6);
scatter(cell2mat(dtn.t),cell2mat(dtn.y),'bx','LineWidth',1.6);
legend([{'Random Ref.'},{'Day'},{'Nigth'}])
hold off;

%% Determine the estimation data
%
%
est_d = 14;

dte.y = dtd.y(1:est_d); dte.t = dtd.t(1:est_d);
dte.init = dtd.initial(1:est_d); dte.final = dtd.final(1:est_d); 
dte.valid.y = dtd.y(est_d+1:end); dte.valid.t = dtd.t(est_d+1:end);
dte.valid.init = dtd.initial(est_d+1:end); dte.valid.final = dtd.final(est_d+1:end);

A.valid.Ad = A.Ad(est_d+1:end); A.valid.td = A.td(est_d+1:end);
A.valid.An = A.An(est_d+1:end); A.valid.tn = A.tn(est_d+1:end);

%% Select the model structure  
%
%
%

struc = struc_select('trivial');

%% 

parameters = est_dcinc(dte,struc,'16');

[dts,j_c] = sim_alert_dc_inc(Parameters,dte,validation);














