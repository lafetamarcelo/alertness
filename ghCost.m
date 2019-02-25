function J = ghCost(par,dte)

dts.yd = cell(length(dte.y),1);
dts.td = cell(length(dte.y),1);
dts.yn = cell(length(dte.y),1);
dts.tn = cell(length(dte.y),1);
hom = zeros(length(dte.y),2);
y_comp = cell(length(dte.y),1);

omega = par(1); tau = par(2); M = par(3); cphase = par(4);
offset = par(5); y_oo = par(6); tau_e = par(7);

A = [0, 0, 0, 0;...
     1, 0, 0, -omega^2/tau;...
     0, 1, 0, -omega^2;...
     0, 0, 1, -1/tau];

init = dte.y{1}(1);

sys = ss(A,[],[0, 0, 0, 1],0);

for w = 1 : length(dte.y)
    
    if w == 1
       tm = dte.t{w}(1);
    else
       tm = dte.init(w);
    end
    
    omega_t = omega*tm;
    k1 = M*cos(omega_t + cphase);
    k2 = -M*sin(omega_t + cphase)*omega;
    h0 = init - k1 - offset;
    
     B = [(offset*omega^2)/tau;...
          (k2 + h0*omega^2*tau + offset*omega^2*tau)/tau;...
          (k1 + offset + k2*tau)/tau;...
           h0 + k1 + offset];
    
    timing = linspace(0,dte.final(w)-tm,10000)';
    [dts.yd{w},dts.td{w}] = lsim(sys,zeros(10000,1),timing,B);
    
    dts.td{w} = dts.td{w} + tm;
    dts.dhom{w} = dts.yd{w} - M*cos(omega*dts.td{w}+cphase);
    
    ind = zeros(length(dte.t{w}),1);
    for i = 1 : length(dte.t{w})
       error = abs(dts.td{w} - dte.t{w}(i));
       [~,ind(i)] = min(error);
    end
    
    y_comp{w} = dts.yd{w}(ind);
    
    if w ~= length(dte.init)
        
        hom(w,1) = dts.dhom{w}(end);
    
        to = linspace(dte.final(w),dte.init(w+1),100)' - dte.final(w);
        Sn = y_oo*(1-exp(-to./tau_e)) + hom(w,1)*exp(-to./tau_e);
        hom(w,2) = Sn(end);
    
        dts.tn{w} = to + dte.final(w);
        dts.yn{w} = Sn + M*cos(omega*dts.tn{w} + cphase);
        dts.nhom{w} = Sn;
        init = dts.yn{w}(end);
        
    end
    
end
 
E = cell2mat(y_comp) - cell2mat(dte.y);
J = (E'*E)/length(E);
end