%% Monte Carlo MOLI runs

close all;
clear;
clc;
set(0,'defaultfigurecolor',[1 1 1]);

%%

MCruns = 100;

par = cell(MCruns,1);
struc = cell(MCruns,1);
J = cell(MCruns,1);

%% Generate the data or load pre--saved ones

WL = autoGen(14);
[dtd,dtn,A] = alertness_gen(1,WL,'random','not');
%load('monte_data.mat');

est_d = 7;

dte.y = dtd.y(1:est_d); dte.t = dtd.t(1:est_d);
dte.init = dtd.initial(1:est_d); dte.final = dtd.final(1:est_d); 
dte.valid.y = dtd.y(est_d+1:end); dte.valid.t = dtd.t(est_d+1:end);
dte.valid.init = dtd.initial(est_d+1:end); dte.valid.final = dtd.final(est_d+1:end);

A.valid.Ad = A.Ad(est_d+1:end); A.valid.td = A.td(est_d+1:end);
A.valid.An = A.An(est_d+1:end); A.valid.tn = A.tn(est_d+1:end);

%% Create the data for filtering (including noise)
y_c = dte.y;

%%

n_s = 1;
for mc = 1 : MCruns
    
    rng('shuffle');

    noise_nat = 'white';
    SNR = 40;
    for i = 1 : length(dte.y)
        if strcmp(noise_nat,'white')
            P_s = sum(y_c{i}.*y_c{i})/length(y_c{i});
            dte.y{i} = awgn(y_c{i},SNR,'measured');
            noise = dte.y{i} - y_c{i};
            P_n = sum(noise.*noise)/length(noise);
            SNR_check(n_s) = 10*log10(P_s/P_n);
            n_s = n_s + 1;
        elseif strcmp(noise_nat,'colored')
            e = filter(1,[1 -.9],randn(length(dte.y{i}),1));
            v = e*std(cell2mat(dte.y))*10^(-SNR/20)/std(e);
            e = v;
            dte.y{i} = dte.y{i} + e; 
        end
    end
    
    struc{mc} = struc_select('barycenter',dte);
    
    par{mc} = est_regr(dte,struc{mc},'15');

    time.init = [dte.valid.t{1}(1), dte.valid.init(2:end)];
    time.final = dte.valid.final;
    initial = dte.valid.y{1}(1);
    
    dts = sim_system(par{mc},time,initial);
    
    fit = e_quad(dts,dte.valid);
     
    J{mc} = fit.quadratic;
    
    disp(['Monte Carlo run ',num2str(mc),...
                             ' Quadratic Error: ',num2str(fit.quadratic)]);
    
end

%% Signal to Noise Ratio Histogram
figure(8); hold on;
histogram(SNR_check,20);

%% Parameters error analisys

E = zeros(MCruns,9);

for mc = 1 : MCruns
   
   E(mc,1) = (par{mc}.est.omega - par{mc}.real.omega)/par{mc}.real.omega; 
   E(mc,2) = (par{mc}.est.tau - par{mc}.real.tau)/par{mc}.real.tau; 
   E(mc,3) = (par{mc}.est.cphase - ...
                                  par{mc}.real.cphase)/par{mc}.real.cphase;
   E(mc,4) = (par{mc}.est.M - par{mc}.real.M)/par{mc}.real.M;
   E(mc,6) = (par{mc}.est.offset - par{mc}.real.offset)/par{mc}.real.offset;
   E(mc,7) = (par{mc}.est.y_oo - par{mc}.real.y_oo)/par{mc}.real.y_oo;
   E(mc,8) = (par{mc}.est.tau_e - par{mc}.real.tau_e)/par{mc}.real.tau_e;
   
   for k = 1 : length(par{mc}.est.h0) 
       E(mc,5) = E(mc,5) + (par{mc}.est.h0(k) - par{mc}.real.h0(k))^2;
   end
   
   E(mc,5) = E(mc,5)^.5/length(par{mc}.est.h0);
   
   for k = 1 : length(par{mc}.est.h0) 
       E(mc,9) = E(mc,9) + (par{mc}.est.h0(k) + par{mc}.est.offset ...
                    - par{mc}.real.h0(k) - par{mc}.real.offset)^2;
   end
   
   E(mc,9) = E(mc,9)^.5/length(par{mc}.est.h0);    
end

E = 100.*E;

%% Plotting parameters error results

figure(9); hold on;
histogram(cell2mat(J),20);
title(['$$\overline{x} =',num2str(round(mean(cell2mat(J)),3)),...
     '\pm',num2str(round(std(cell2mat(J)),3)),'$$'],'Interpreter','latex');

figure(10); hold on;
for k = 1 : 8
    subplot(2,4,k); hold on;
    histogram(E(:,k),20); hold on;
    title(['$$\overline{x} =',num2str(round(mean(E(:,k)),3)),...
          '\pm',num2str(round(std(E(:,k)),3)),'$$'],'Interpreter','latex');
    leg = par{1}.TexNames{k};
    %xlim([-100,100]);
    legend([{leg}],'Interpreter','latex'); hold off;
end

figure(11); hold on;
histogram(E(:,9),20); hold on;
title(['$$\overline{x} =',num2str(round(mean(E(:,9)),3)),...
         '\pm',num2str(round(std(E(:,9)),3)),'$$'],'Interpreter','latex');
leg = 'DC + h(0)';
%xlim([-100,100]);
legend([{leg}],'Interpreter','latex'); hold off;
