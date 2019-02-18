function dts = allert_sim(time,initialization,model)

dts.yd = cell(length(time.initial),1);
dts.td = cell(length(time.initial),1);
dts.yn = cell(length(time.initial),1);
dts.tn = cell(length(time.initial),1);
hom = zeros(length(time.initial),2);

%% Parameters determination
M = model.M;
cphase = model.cphase;
omega = model.omega;
offset = model.offset;
tau = model.tau;
tau_e = model.tau_e;
y_oo = model.y_oo;

%%
A = [0, 0, 0, 0;...
    1, 0, 0, -omega^2/tau;...
    0, 1, 0, -omega^2;...
    0, 0, 1, -1/tau];

sys = ss(A,[],[0, 0, 0, 1],0);

for i = 1 : length(time.initial)
    
    if i == 1
        tm = time.initial(i);
    else
        tm = time.wake(i);
    end
   
    omega_t = omega*tm;
    k1 = M*cos(omega_t + cphase);
    k2 = -M*sin(omega_t + cphase)*omega;
    h0 = initialization - k1 - offset;
    
    B = [(offset*omega^2)/tau;...
       (k2 + h0*omega^2*tau + offset*omega^2*tau)/tau;...
       (k1 + offset + k2*tau)/tau;...
        h0 + k1 + offset];
    
    %[dts.yd{i},dts.td{i}] = impulse(sys,time.sleep(i)-tm);
    %[dts.yd{i},dts.td{i}] = initial(sys,B,time.sleep(i)-tm);
    timing = linspace(0,time.sleep(i)-tm,10000)';
    [dts.yd{i},dts.td{i}] = lsim(sys,zeros(10000,1),timing,B);
    
    %[simout, tout] = varStepSimS(sys,B,0,timing,'15');
    
    %dts.yd{i} = simout(:,1);
    %dts.td{i} = tout + tm;
    dts.td{i} = dts.td{i} + tm;
    dts.dhom{i} = dts.yd{i} - M*cos(omega*dts.td{i}+cphase);
    
    %% DEBUG plot
    figure(1); subplot(3,1,1); hold on;
    plot(dts.td{i},dts.yd{i},'r--');
    figure(1); subplot(3,1,2); hold on;
    plot(dts.td{i},dts.dhom{i},'r--');
    figure(1); subplot(3,1,3); hold on;
    plot(dts.td{i},dts.yd{i}-dts.dhom{i},'r--');
    
    if i ~= length(time.initial)
        hom(i,1) = dts.yd{i}(end) - ...
                          M*cos(omega*dts.td{i}(end) + cphase);
    
        to = linspace(time.sleep(i),time.wake(i+1),100)' - time.sleep(i);
        Sn = y_oo*(1-exp(-to./tau_e)) + hom(i,1)*exp(-to./tau_e);
        hom(i,2) = Sn(end);
    
        dts.tn{i} = to + time.sleep(i);
        dts.yn{i} = Sn + M*cos(omega*dts.tn{i} + cphase);
        dts.nhom{i} = Sn;
        
        %% DEBUG
        figure(1); subplot(3,1,1); hold on;
        plot(dts.tn{i},dts.yn{i},'b--');
        figure(1); subplot(3,1,2); hold on;
        plot(dts.tn{i},dts.nhom{i},'b--');
        figure(1); subplot(3,1,3); hold on;
        plot(dts.tn{i},dts.yn{i}-dts.nhom{i},'b--');
        
        initialization = dts.yn{i}(end);
    end
end

end