function dts = pred_system(parameters,dte)

%% Parameters light
omega = parameters.est.omega;
tau = parameters.est.tau;
M = parameters.est.M;
cphase = parameters.est.cphase;
offset = parameters.est.offset;

sys = ss(parameters.est.struc.A,[],[0, 0, 0, 1],0);
%%
for w = 1 : length(dte.t)
    dts.y{w}(1) = dte.y{w}(1);
    dts.t{w}(1) = dte.t{w}(1);
    for k = 1 : length(dte.t{w})-1
        
        omega_t = omega*dte.t{w}(k);
        k1 = M*cos(omega_t + cphase);
        k2 = -M*sin(omega_t + cphase)*omega;
        h0 = dte.y{w}(k) - k1 - offset;
        
        B = [(offset*omega^2)/tau; ...
             (k2 + h0*omega^2*tau + offset*omega^2*tau)/tau; ...
             (k1 + offset + k2*tau)/tau; ...
              h0 + k1 + offset];
        
        timing = linspace(0,dte.t{w}(k+1)-dte.t{w}(k),100)';
        [y_s,t_s] = lsim(sys,zeros(100,1),timing,B);
        
        dts.t{w}(k+1) = t_s(end) + dte.t{w}(k);
        dts.y{w}(k+1) = y_s(end);
    end
end

end