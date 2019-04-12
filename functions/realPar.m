function par = realPar(dte)

omega = double(pi/12); M = 2.52; tau = 1/0.0353; 
cphase = double(-16.835*pi/12); dc_level = 2.4;
y_oo = 14.3; tau_e = 2.6247;

Par.RealStruc.A = [0, 0, 0, 0;...
                   1, 0, 0, -omega^2/tau; ...
                   0, 1, 0, -omega^2;...
                   0, 0, 1, -1/tau];

for i = 1 : wind_s
    omega_t = omega*(dte.t{i}(1) - 24*(i-1));
    h0 = dte.y{i}(1) - M*cos(omega_t  + cphase) - dc_level;
    k1 = M*cos(cphase+omega_t);
    k2 = -M*sin(cphase+omega_t)*omega;
    Par.RealPar(i,:) = [omega, tau, cphase, M, h0, dc_level, y_oo, tau_e];
    Par.RealStruc.B{i} = [ (dc_level*omega^2)/tau;...
                       (k2 + h0*tau*omega^2 + dc_level*tau*omega^2)/tau;...
                       (k1+dc_level+k2*tau)/tau;...
                        h0 + k1 + dc_level];
end

parameters.real.omega = mean(Par.RealPar(:,1));
parameters.real.tau = mean(Par.RealPar(:,2));
parameters.real.cphase = mean(Par.RealPar(:,3));
parameters.real.M = mean(Par.RealPar(:,4));
parameters.real.h0 = Par.RealPar(:,5);
parameters.real.offset = mean(Par.RealPar(:,6));
parameters.real.y_oo = 14.3;
parameters.ral.tau_e = 1/0.381;

parameters.real.struc = Par.RealStruc;

par = parameters.real;

end