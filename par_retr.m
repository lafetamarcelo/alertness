function Par = par_retr(A, z0, dte)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

wind_s = length(dte.y);

Par.Names{1} = 'Omega'; Par.Names{2} = 'Tau';
Par.Names{3} = 'Phase'; Par.Names{4} = 'M';
Par.Names{5} = 'h(0)';
Par.TexNames{1} = '$$\omega_n$$'; Par.TexNames{2} = '$$\tau$$';
Par.TexNames{3} = '$$\phi_{c}$$'; Par.TexNames{4} = '$$M$$';
Par.TexNames{5} = '$$h(0)$$'; 
Par.TexNames{6} = '$$A_{DC}$$'; Par.TexNames{7} = '$$y_{\infty}$$';
Par.TexNames{8} = '$$\tau_{e}$$';


Par.EstStruc.A = A;

tauhat = -1/A(4,4);
omegahat = (-A(3,4))^.5;

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
%% Real parameters

omega = double(pi/12); M = 2.52; tau = 1/0.0353; 
cphase = double(-16.835*pi/12); dc_level = 0;
y_oo = 14.3; tau_e = 2.6247;

Par.RealStruc.A = [0, 0, 0, 0;...
                   1, 0, 0, -omega^2/tau; ...
                   0, 1, 0, -omega^2;...
                   0, 0, 1, -1/tau];

for i = 1 : wind_s
    omega_t = omega*(dte.t{i}(1) - 24*(k-1));
    h0 = dte.y{i}(1) - M*cos(omega_t  + cphase) - dc_level;
    k1 = M*cos(cphase+omega_t);
    k2 = -M*sin(cphase+omega_t)*omega;
    Par.RealPar(i,:) = [omega, tau, cphase, M, h0, dc_level, y_oo, tau_e];
    Par.RealStruc.B{i} = [ (dc_level*omega^2)/tau;...
                       (k2 + h0*tau*omega^2 + dc_level*tau*omega^2)/tau;...
                       (k1+dc_level+k2*tau)/tau;...
                        h0 + k1 + dc_level];
end

