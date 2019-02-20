function parameters = est_regr(dte,struc,version)
                                                                       
    reg_filt = ss(struc.A',struc.C',eye(4),0);
    wind_s = length(dte.y);
    
    Phi = cell(wind_s,2+wind_s);
    Y = cell(wind_s,1);
    Yr = cell(wind_s,1);
    %% Regression matrix
    
    for w = 1 : wind_s
        
        % Impulse filtering - Input [B] informations
        delt = dte.t{w}(end)-dte.t{w}(1);
        step_t = .0001;%*min(diff(dte.t{w}));
        puls_time = 0:step_t:delt;
        puls_out = impulse(reg_filt,puls_time);
        
        ind = zeros(length(dte.t{w}),1);
        for i = 1 : length(dte.t{w})
           error = abs(puls_time - (dte.t{w}(i) - dte.t{w}(1)));
           [~,ind(i)] = min(error);
           
           % simulated time error plotting
           error_graph(i) = min(error);
        end
        
        figure(2); subplot(2,1,1); hold on;
        plot(error_graph,'Color',[.7 .7 .7]);
        
        Phi{w,2} = puls_out(ind(2:end),1);
        X = puls_out(ind(2:end),2:4);
        
        % y(t) filtering - Output [L] informations
        t_calc = dte.t{w}-dte.t{w}(1);
        [sig_out,t_out] =  sim_vs(reg_filt,0,dte.y{w},t_calc,version);
        %sig_out = lsim(reg_filt,dte.y{w},puls_time(ind));
        
        ind = zeros(length(dte.t{w}),1);
        for i = 1 : length(dte.t{w})
           error = abs(t_out - (dte.t{w}(i) - dte.t{w}(1)));
           [~,ind(i)] = min(error);
           
           % simulated time error plotting
           error_graph(i) = error(ind(i));
        end
        
        figure(2); subplot(2,1,2); hold on;
        plot(error_graph,'Color',[.4 .4 .4]);
        
        Phi{w,1} = [zeros(length(ind(2:end)),1),sig_out(ind(2:end),2:4)];
        Y{w} = dte.y{w}(2:end);
        Yr{w} = sig_out(ind(2:end),5);
        
        % Shifting the window information
        for j = 1 : wind_s
            if w == j
                Phi{w,j+2} = X(:,1:3); 
            else
                Phi{w,j+2} = zeros(length(Phi{w,1}),3);
            end
        end
        
    end

    %% Determine the parameters
    
    Phi_t = cell2mat(Phi);
    
    Theta = [0;Phi_t(:,2:end)\cell2mat(Y)];
 
    J = (cell2mat(Y) - cell2mat(Phi)*Theta)'*(cell2mat(Y) ...
                                - cell2mat(Phi)*Theta)/length(cell2mat(Y));
                            
    Thetar = pinv(cell2mat(Phi))*cell2mat(Yr);
 
    Jr = (cell2mat(Yr) - cell2mat(Phi)*Thetar)'*(cell2mat(Yr) ...
                               - cell2mat(Phi)*Thetar)/length(cell2mat(Yr));
    %% Parameters matrix reconstruction
    
    L = Theta(1:4);
    
    B = cell(wind_s,1);
    for i = 1 : wind_s
        B{i}(1) = Theta(5); 
        B{i}(2) = Theta(6+3*(i-1)); 
        B{i}(3) = Theta(7+3*(i-1));
        B{i}(4) = Theta(8+3*(i-1));
    end
    
    A = struc.A + L*struc.C;
    
    %% Just for DEBUG
    
%     load('rTHETA.mat')
%     rL = (model.A - struc.A)*pinv(model.C);
%     
%     % Parameters determination
%     M = model.p.M;
%     cphase = model.p.cphase;
%     omega = model.p.omega;
%     offset = model.p.offset;
%     tau = model.p.tau;
%     tau_e = model.p.tau_e;
%     y_oo = model.p.y_oo;
%     
%     
%     rTheta = rL;
%     y_r = cell(length(dte.y),1);
%     for i = 1 : length(dte.y)
%         omega_t = omega*dte.t{i}(1);
%         k1 = M*cos(omega_t + cphase);
%         k2 = -M*sin(omega_t + cphase)*omega;
%         h0 = dte.y{i}(1) - k1 - offset;
% 
%         rB = [(offset*omega^2)/tau;...
%            (k2 + h0*omega^2*tau + offset*omega^2*tau)/tau;...
%            (k1 + offset + k2*tau)/tau;...
%             h0 + k1 + offset];
%         
%         if i~=1
%             rTheta = [rTheta;rB(2:end)];
%         else
%             rTheta = [rTheta;rB];
%         end
%         
%         y_r{i} = dte.y{i}(2:end);
%     end 
%     
%     J_comp = (cell2mat(Y) - cell2mat(Phi)*rTheta)'*(cell2mat(Y) ...
%                                - cell2mat(Phi)*rTheta)/length(cell2mat(Y));
%     
%     J_real = (cell2mat(y_r) - cell2mat(Phi)*rTheta)'*(cell2mat(y_r) ...
%                               - cell2mat(Phi)*rTheta)/length(cell2mat(y_r));
%                           
%     rPhi = cell2mat(y_r)*pinv(rTheta);
%     
%     error_Phi = cell2mat(Phi) - rPhi;
    %% Retrieve the model parameters
                  
    par = par_retr(A,B,dte);    
    
    %% Simulate the system for each homeostatic level
    
    h.y = cell(length(dte.y)-1,1);
    h.t = cell(length(dte.y)-1,1);
    
    sys = ss(par.est.struc.A,[],[0, 0, 0, 1],0);
    sys_rev = ss(-par.est.struc.A,[],[0, 0, 0, 1],0);
    
    for k = 1: length(dte.y)
        
        B = par.est.struc.B{k};
        
        timing = linspace(0,dte.final(k)-dte.t{k}(1),10000)';
        [y_adv,t_adv] = lsim(sys,zeros(10000,1),timing,B);
        t_adv = t_adv + dte.t{k}(1);
        
        if (dte.t{k}(1)-dte.init(k)) ~= 0
            timing_rev = linspace(0,dte.t{k}(1)-dte.init(k),100)';
            [y_rev,t_rev] = lsim(sys_rev,zeros(100,1),timing_rev,B);
            t_rev = t_rev + dte.init(k);

            dts.y{k} = [fliplr(y_rev')'; y_adv(2:end)];
            dts.t{k} = [t_rev; t_adv(2:end)];
        else
            dts.y{k} = y_adv;
            dts.t{k} = t_adv;
        end
        
        dts.circ{k} = par.est.M * cos(par.est.omega*dts.t{k} + par.est.cphase);
        dts.hom{k} = dts.y{k} - dts.circ{k};
        
        if k ~= 1
            h.y{k-1}(1,1) = dts.hom{k-1}(end);
            h.y{k-1}(1,2) = dts.hom{k}(1);
            h.t{k-1}(1,1) = dts.t{k}(1) - dts.t{k-1}(end); 
            %%DEBUG figure
            subplot(3,1,3); hold on;
            scatter(dts.t{k-1}(end),h.y{k-1}(1,1),'xr','LineWidth',1.8);
            scatter(dts.t{k}(1),h.y{k-1}(1,2),'or','LineWidth',1.8);
        end
        
        %% Just for DEBUG
        figure(3); subplot(3,1,1); hold on;
        plot(dts.t{k},dts.y{k},'k-','LineWidth',1.2);
        subplot(3,1,2); hold on;
        plot(dts.t{k},dts.circ{k},'k-','LineWidth',1.2);
        subplot(3,1,3); hold on;
        plot(dts.t{k},dts.hom{k},'k-','LineWidth',1.2);
    end
    
    %% Estimate the night parameters
    
    % Estimate options
    maxITER = 200;
    options = optimoptions(@lsqnonlin,'Jacobian','on',...
        'Display','off','TolFun',1e-6,'MaxIter',maxITER,'TolX',1e-6);
    initial = [10,25];
    
    func = @cost_function;
    y_lsq = cell2mat(h.y);
    
    [x,~,~,~,~] = lsqnonlin(func,initial,[],[],options,...
           cell2mat(h.t),y_lsq(:,1),y_lsq(:,2));
    
    par.est.y_oo = x(1); par.est.tau_e = x(2);
    parameters = par;
end