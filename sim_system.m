function dts = sim_system(parameters,time,initial)

%% Parameters light

omega = parameters.est.omega;
tau = parameters.est.tau;
M = parameters.est.M;
cphase = parameters.est.cphase;
offset = parameters.est.offset;
y_oo = parameters.est.y_oo;
tau_e = parameters.est.tau_e;

%%
dts.y = cell(length(time.init),1);
dts.circ = cell(length(time.init),1);
dts.hom = cell(length(time.init),1);
dts.t = cell(length(time.init),1);

hom = zeros(length(time.init)-1,2);

sys = ss(parameters.est.struc.A,[],[0, 0, 0, 1],0);

for k = 1 : length(time.init)
    
    omega_t = omega*time.init(k);
    k1 = M*cos(omega_t + cphase);
    k2 = -M*sin(omega_t + cphase)*omega;
    h0 = initial - k1 - offset;
    
    B = [(offset*omega^2)/tau; ...
         (k2 + h0*omega^2*tau + offset*omega^2*tau)/tau; ...
         (k1 + offset + k2*tau)/tau; ...
          h0 + k1 + offset];
      
    timing = linspace(0,time.final(k)-time.init(k),10000)';
    [dts.yd{k},t_adv] = lsim(sys,zeros(10000,1),timing,B);
    dts.td{k} = t_adv + time.init(k);
    
    dts.dhom{k} = dts.yd{k} - M*cos(omega*dts.td{k}+cphase);
    
    if k ~= length(time.init)
        hom(k,1) = dts.dhom{k}(end);
    
        to = linspace(0,time.init(k+1)-time.final(k),100)';
        Sn = y_oo*(1-exp(-to./tau_e)) + hom(k,1)*exp(-to./tau_e);
        hom(k,2) = Sn(end);
    
        dts.tn{k} = to + time.final(k);
        dts.yn{k} = Sn + M*cos(omega*dts.tn{k} + cphase);
        dts.nhom{k} = Sn;
        
        initial = dts.yn{k}(end);
    end
    
end

end