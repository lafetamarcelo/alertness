function parameters = par_retr(A, z0, dte)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

wind_s = length(dte.y);

parameters.Names{1} = 'Omega'; parameters.Names{2} = 'Tau';
parameters.Names{3} = 'Phase'; parameters.Names{4} = 'M';
parameters.Names{5} = 'h(0)';

parameters.TexNames{1} = '$$\omega_n$$'; 
parameters.TexNames{2} = '$$\tau$$';
parameters.TexNames{3} = '$$\phi_{c}$$'; 
parameters.TexNames{4} = '$$M$$';
parameters.TexNames{5} = '$$h(0)$$'; 
parameters.TexNames{6} = '$$A_{DC}$$'; 
parameters.TexNames{7} = '$$y_{\infty}$$';
parameters.TexNames{8} = '$$\tau_{e}$$';

%%
Par.EstStruc.A = A;

% Direct 
%tauhat = -1/A(4,4);
omegahat = (-A(3,4))^.5;

% Combined
% omegahat = (-A(2,4)*tauhat_)^.5;
 tauhat = 1/(-A(2,4)/omegahat^2);

offset = z0{1}(1)*tauhat/omegahat^2;

X = [omegahat^2, 0, 1/tauhat;...
    0, 1/tauhat, 1;...
    1, 1, 0];

comp = [-offset*omegahat^2;...
        -offset/tauhat;...
        -offset];

for k = 1 : wind_s
    b_comp = z0{k}(2:end)' + comp;
    h0_k1_k2hat = X\b_comp;
    h0hat = h0_k1_k2hat(1);

    tan_cphase = h0_k1_k2hat(3)/(-omegahat)/h0_k1_k2hat(2);     
    
    if(tan_cphase < 0)
        cphasehat = atan(tan_cphase)-pi;
        Mhat = h0_k1_k2hat(2)/cos(cphasehat);
    else
        cphasehat = atan(tan_cphase);
        Mhat = h0_k1_k2hat(2)/cos(cphasehat);
    end
    
    cphase_h =  cphasehat - omegahat*(dte.t{k}(1)-24*(k-1));
    
    Par.EstPar(k,:) = [omegahat, tauhat, cphase_h,...
        Mhat, h0hat, offset];
    
    Par.EstStruc.B{k} = z0{k};
    
end


parameters.est.omega = mean(Par.EstPar(:,1));
parameters.est.tau = mean(Par.EstPar(:,2));
parameters.est.cphase = mean(Par.EstPar(:,3));
parameters.est.M = mean(Par.EstPar(:,4));
parameters.est.h0 = Par.EstPar(:,5);
parameters.est.offset = mean(Par.EstPar(:,6));

parameters.est.struc = Par.EstStruc;

parameters.est.omega_ = abs((-A(2,4)*tauhat)^.5);
parameters.est.tau_ = abs(1/(-A(2,4)/omegahat^2));
%% Real parameters

omega = double(pi/12); M = 2.52; tau = 1/0.0353; 
cphase = double(-16.835*pi/12); dc_level = 2.4;
y_oo = 14.3; tau_e = 1/0.381;

Par.RealStruc.A = [0, 0, 0, 0;...
                   1, 0, 0, -omega^2/tau; ...
                   0, 1, 0, -omega^2;...
                   0, 0, 1, -1/tau];

for i = 1 : wind_s
    omega_t = omega*dte.t{i}(1);
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
parameters.real.tau_e = 1/0.381;

parameters.real.struc = Par.RealStruc;

end

