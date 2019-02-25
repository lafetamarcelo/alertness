%% Monte Carlo of Go Horse Opt.

clear;
clc;
set(0,'defaultfigurecolor',[1 1 1]);

%%

MCruns = 100;

par = cell(MCruns,1);
J = cell(MCruns,1);
e_flag = cell(MCruns,1);
cruze = cell(MCruns,1);

%% Generate the data or load pre--saved ones

WL = autoGen(14);
[dtd,dtn,A] = alertness_gen(2,WL,'random','not');
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

%% Optimization setup

options = optimset('Display','off','MaxIter',100);

prop = .2;
%% Real parameters

rPar(1) = pi/12; %omega
rPar(2) = 1/0.0353; %tau
rPar(3) = 2.52; %M
rPar(4) = -16.835*pi/12; %cPhase
rPar(5) = 2.4; %offset
rPar(6) = 14.3; %y_oo
rPar(7) = 2.6247; %tau_e

%%

for mc = 1 : MCruns
    
    rng('shuffle');
    
    noise_nat = 'noNoise';
    SNR = 40;
    to = cell2mat(dte.t);

    if strcmp(noise_nat,'white')
        e = idinput(length(cell2mat(y_c)),'rgs')*...
                                          std(cell2mat(y_c))*10^(-SNR/20);
        yo = cell2mat(y_c) + e;
    elseif strcmp(noise_nat,'colored')
        e = filter(1,[1 -.9],randn(length(cell2mat(y_c)),1));
        v = e*std(cell2mat(y_c))*10^(-SNR/20)/std(e);
        e = v;
        yo = cell2mat(y_c) + e;
    elseif strcmp(noise_nat,'noNoise')
        yo = cell2mat(y_c);
    end

    for w = 1 : length(y_c)
       ind = find((to>=dte.t{w}(1)).*(to<=dte.t{w}(end)));
       dte.y{w} = yo(ind);
    end
    
    cruze{mc} = prop.*rPar.*(rand(1,length(rPar))-.5);
    
    init = rPar + cruze{mc};
    
    [par{mc},J{mc},e_flag{mc}] = fminsearch(@ghCost,init,options,dte);
    
    disp(['Monte Carlo run ',num2str(mc),...
                             ' Quadratic Error: ',num2str(J{mc})]);
end


%% MOnte Carlo figures

figure(2); hold on;
subplot(1,2,1);
histogram(cell2mat(J),20);
subplot(1,2,2);
%labels = {'Function Terminated','Excided Iterations','Converted'};
flags = cell2mat(e_flag)'+2;
pie(count(num2str(flags),{'1','2','3'}));

figure(3); hold on;
aux = cell2mat(cruze);
for i = 1 : 7
   subplot(4,2,i); 
   histogram(aux(:,i)+rPar(i),20); hold on;
   title(['$$\overline{x} =',num2str(round(mean(aux(:,i)+rPar(i)),3)),...
       '\pm',num2str(round(std(aux(:,i)+rPar(i)),3)),'$$'],...
       'Interpreter','latex');
   legend([{'Cruzing'}]);
end

figure(4); hold on;
aux_i = cell2mat(par);
for i = 1 : 7
   subplot(4,2,i);
   e_prop = 100.*(rPar(i) - aux_i(:,i))/rPar(i); 
   histogram(e_prop,20);
   title(['$$\overline{x} =',num2str(round(mean(e_prop),3)),...
       '\pm',num2str(round(std(e_prop),3)),'$$'],'Interpreter','latex');
   legend([{'Error'}]); 
end