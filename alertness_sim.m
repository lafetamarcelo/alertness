function [dts,dta] = alertness_sim(ppH,days,space,noise_nat,SNR)
    
    dts.y = cell(days,1); dts.t = cell(days,1); 
    dts.dhom = cell(days,1); dta.dhom = cell(days,1);
    dta.yo = cell(days,1); dta.yd = cell(days,1);
    dta.td = cell(days,1); dta.yn = cell(days,1);
    dta.tn = cell(days,1); dta.nhom = cell(days,1);
    
    %
    omega = pi/12; %omega
    tau = 1/0.0353; %tau
    M = 2.52; %M
    cphase = -16.835*pi/12; %cPhase
    offset = 2.4; %offset
    y_oo = 14.3; %y_oo
    tau_e = 2.6247; %tau_e

    % 
    W = autoGen(days);
    
    acum = 0;
    for i = 1 : length(W(:,1))
       time.init(i) = acum;
       time.final(i) = acum + W(i,2);
       acum = time.final(i) + W(i,1);
    end
    
    % Create the alertness model
    
    A = [0, 0, 0, 0;...
          1, 0, 0, -omega^2/tau;...
          0, 1, 0, -omega^2;...
          0, 0, 1, -1/tau];
    
    sys = ss(A,[],eye(4),0);
    initial = 14.3;
    
    for i = 1 : length(W(:,1))
        
    timing = linspace(0,W(i,2),10000) + time.init(i);    
    
    omega_t = omega*timing(1);
    k1 = M*cos(omega_t + cphase);
    k2 = -M*sin(omega_t + cphase)*omega;
    h0 = initial - k1 - offset;
    
    B = [(offset*omega^2)/tau;...
         (k2 + h0*omega^2*tau + offset*omega^2*tau)/tau;...
         (k1 + offset + k2*tau)/tau;...
          h0 + k1 + offset];
      
    [state,t_adv] = lsim(sys,zeros(10000,1),timing - timing(1),B);
    
    %Include noise to signal
    if strcmp(noise_nat,'white')
        e = idinput(length(state),'rgs')*std(state)*10^(-SNR/20);
        yo = state + e;
    elseif strcmp(noise_nat,'colored')
        e = filter(1,[1 -.9],randn(length(state),4));
        v = e*std(state)*10^(-SNR/20)/std(e);
        e = v;
        yo = state + e;
    elseif strcmp(noise_nat,'noNoise')
        yo = state;
    end
    y_ = [0 0 0 1]*yo';
    
    % Sample out
    if strcmp(space,'random')
        ppW = round(W(i,2)*ppH,0);
        cruze = 2*rand(ppW,1)./ppH - 1/ppH;
        i_cruze = round(cruze*length(y_)/W(i,2),0);
        
        ind = sort(abs(round(linspace(1,length(y_),ppW),0) + i_cruze'));
    elseif strcmp(space,'linear')
        ind = round(linspace(1,length(y_),ppW),0); 
    end
    
    if ind(end) > length(y_); ind(end) = length(y_); end
    if ind(1) <= 0; ind(1) = 1; end
    
    %
    dts.y{i} = y_(ind)';
    dts.t{i} = t_adv(ind) + timing(1);  
    dta.yd{i} = y_';
    dta.yo{i} = ([0, 0, 0, 1]*state')';
    dta.td{i} = timing';
    
    %figure(1)
    %plot(timing,y_,'Color',[.7 .7 .7]); hold on;
    %plot(dta.td{i},[0 0 0 1]*state','r--');
    %scatter(dts.td{i},dts.yd{i});
    
    dts.dhom{i} = dts.y{i} - M*cos(omega*dts.t{i}+cphase);
    dta.dhom{i} = dta.yo{i} - M*cos(omega*dta.td{i}+cphase);
    
    if i ~= length(W(:,1))
        hom(i,1) = dta.dhom{i}(end);
    
        to = linspace(0,time.init(i+1)-time.final(i),1000)';
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