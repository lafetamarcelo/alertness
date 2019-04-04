function [dts,dta] = alertness_sim(W,noise_nat,SNR,resolution,ind)
    
    days = length(W(:,1));
    
    %
    dts.y = cell(days,1); dts.t = cell(days,1); 
    dts.dhom = cell(days,1); dts.s = cell(days,1);
    
    dta.dhom = cell(days,1); dta.so = cell(days,1);
    dta.yo = cell(days,1); dta.yd = cell(days,1);
    dta.td = cell(days,1); dta.yn = cell(days,1);
    dta.tn = cell(days,1); dta.nhom = cell(days,1);
    dta.SNR = zeros(length(W(:,1)),4);
    
    %% Initialize the real parameters
    omega = pi/12; %omega
    tau = 1/0.0353; %tau
    M = 2.52; %M
    cphase = -16.835*pi/12; %cPhase
    offset = 2.4; %offset
    y_oo = 14.3; %y_oo
    tau_e = 2.6247; %tau_e

    %% Create the initial and final time stamp
    acum = 0;
    for i = 1 : length(W(:,1))
       time.init(i) = acum;
       time.final(i) = acum + W(i,2);
       acum = time.final(i) + W(i,1);
    end
    
    %% Create the alertness model
    
    A = [0, 0, 0, 0;...
          1, 0, 0, -omega^2/tau;...
          0, 1, 0, -omega^2;...
          0, 0, 1, -1/tau];
    
    sys = ss(A,[],eye(4),0);
    initial = 14.3;
    
    %% Simulate the alertness model
    
    for i = 1 : length(W(:,1))
        
    timing = linspace(0,W(i,2),resolution) + time.init(i);    
    
    omega_t = omega*timing(1);
    k1 = M*cos(omega_t + cphase);
    k2 = -M*sin(omega_t + cphase)*omega;
    h0 = initial - k1 - offset;
    
    B = [(offset*omega^2)/tau;...
         (k2 + h0*omega^2*tau + offset*omega^2*tau)/tau;...
         (k1 + offset + k2*tau)/tau;...
          h0 + k1 + offset];
      
    % Simulate the state output  
    [state,t_adv] = lsim(sys,zeros(resolution,1),timing - timing(1),B);
    
    % Include noise to signal
    yo = zeros(resolution,4);
    for kn = 1 : 4
        if strcmp(noise_nat,'white') && (kn ~= 1)
            e = idinput(length(state(:,kn)),'rgs')*std(state(:,kn))*10^(-SNR/20);
            yo(:,kn) = state(:,kn) + e;
            %yo(:,kn) = awgn(state(:,kn),SNR,'measured');
            %e = yo(:,kn) - state(:,kn);
            dta.SNR(i,kn) = snr(state(:,kn),e);
        elseif strcmp(noise_nat,'colored') && (kn ~= 1)
            e = filter(1,[1 -.9],randn(length(state(:,kn)),1));
            v = e*std(state(:,kn))*10^(-SNR/20)/std(e);
%             figure(kn);
%             subplot(2,1,1);
%             plot(randn(length(state(:,kn)),1)); hold on;
%             plot(e); hold off;
%             subplot(2,1,2);
%             plot(v); hold off;
            yo(:,kn) = state(:,kn) + v;
        else
            yo(:,kn) = state(:,kn);
        end
    end
    
    if strcmp(noise_nat,'whiteOut')
        aux = [0 0 0 1]*yo';
        e = idinput(length(aux),'rgs')*std(aux)*10^(-SNR/20);
        %y_ = awgn(aux,SNR,'measured');
        y_ = [0 0 0 1]*yo' + e';
    else
        y_ = [0 0 0 1]*yo';
    end
    
    %
    dts.y{i} = y_(ind{i})';
    dts.s{i} = yo(ind{i},:);
    dts.t{i} = t_adv(ind{i}) + timing(1); 
    
    dta.yd{i} = y_';
    dta.yo{i} = ([0, 0, 0, 1]*state')';
    dta.SNRout{i} = snr(dta.yo{i},dta.yd{i} - dta.yo{i});
    dta.so{i} = state;
    dta.td{i} = timing';
    
    %figure(1)
    %plot(timing,y_,'Color',[.7 .7 .7]); hold on;
    %plot(dta.td{i},[0 0 0 1]*state','r--');
    %scatter(dts.t{i},dts.y{i});
    
    dts.dhom{i} = dts.y{i} - M*cos(omega*dts.t{i}+cphase);
    dta.dhom{i} = dta.yo{i} - M*cos(omega*dta.td{i}+cphase);
    
    if i ~= length(W(:,1))
        hom(i,1) = dta.dhom{i}(end);
    
        to = linspace(0,time.init(i+1)-time.final(i),resolution)';
        Sn = y_oo*(1-exp(-to./tau_e)) + hom(i,1)*exp(-to./tau_e);
        hom(i,2) = Sn(end);
        
        dta.tn{i} = to + time.final(i);
        dta.yn{i} = Sn + M*cos(omega*dta.tn{i} + cphase);
        dta.nhom{i} = Sn;
        
        %figure(1); hold on;
        %plot(dta.tn{i},dta.yn{i},'r--');
        
        initial = dta.yn{i}(end);
    end  
    
    end
    
    dts.initial = time.init;
    dts.final = time.final;
    dta.initial = time.init;
    dta.final = time.final;
end