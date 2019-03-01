function parameters = est_regr(dte,struc,version)
                                                                       
    reg_filt = ss(struc.A',struc.C',eye(4),0);
    wind_s = length(dte.y);
    
    Phi = cell(wind_s,2+wind_s);
    Phi_iv = cell(wind_s,2+wind_s);
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
        
        %figure(2); subplot(2,1,1); hold on;
        %plot(error_graph,'Color',[.7 .7 .7]);
        
        Phi{w,2} = puls_out(ind(2:end),1);
        Phi_iv{w,2} = puls_out(ind,1);
        X = puls_out(ind(2:end),2:4);
        X_iv = puls_out(ind,2:4);
        
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
        
        %figure(2); subplot(2,1,2); hold on;
        %plot(error_graph,'Color',[.4 .4 .4]);
        
        Phi{w,1} = [zeros(length(ind(2:end)),1),sig_out(ind(2:end),2:4)];
        Phi_iv{w,1} = [zeros(length(ind),1),sig_out(ind,2:4)];
        
        Y{w} = dte.y{w}(2:end);
        Yr{w} = sig_out(ind(2:end),5);
        
        % Shifting the window information
        for j = 1 : wind_s
            if w == j
                Phi{w,j+2} = X(:,1:3);
                Phi_iv{w,j+2} = X_iv(:,1:3);
            else
                Phi{w,j+2} = zeros(length(Phi{w,1}),3);
                Phi_iv{w,j+2} = zeros(length(Phi_iv{w,1}),3);
            end
        end
        
    end

    %% Determine the parameters
    
    Phi_t = cell2mat(Phi);
    
    Theta = [0;Phi_t(:,2:end)\cell2mat(Y)];
 
    J = (cell2mat(Y) - cell2mat(Phi)*Theta)'*(cell2mat(Y) ...
                                - cell2mat(Phi)*Theta)/length(cell2mat(Y));

    %% Introduction to residual modeling 
    
    Ep = cell2mat(Y) - Phi_t*Theta;
    
    figure(10); hold on;
    subplot(2,1,1); plot(Ep,'LineWidth',1.2);
    subplot(2,1,2); 
    plot(cell2mat(Y),'LineWidth',1.2); hold on; 
    plot(Phi_t*Theta,'r--','LineWidth',1.2);
    hold off;
    
    itt = 50;
    J_it = zeros(itt,1);
    Ph = cell(itt,1);
    Th = cell(itt,1);
    Th_e = cell(itt,1);
    for it = 1 : itt
        
        figure(10);
        subplot(2,1,1); 
        plot(Ep,'LineWidth',1.2);
        
        ind_i = 1;
        Phi_e = cell(length(dte.y),1);
     
        for w = 1 : length(dte.y)
            
            t_calc = dte.t{w}(2:end) - dte.t{w}(2);
            ind_f = ind_i + length(t_calc) - 1;
            [sig_out,~] =  sim_vs(reg_filt,0,-Ep(ind_i:ind_f),...
                                                           t_calc,version);
            ind_i = ind_f + 1;
            
            %for j = 1 : length(dte.y)
            %    if j ~= w
            %        Phi_e{w,j} = zeros(length(t_calc),4);
            %    else
            %        Phi_e{w,j} = sig_out(:,1:4);
            %    end
            %end
            Phi_e{w} = sig_out(:,1:4);
        end
        
        Phi_it = [cell2mat(Phi) , cell2mat(Phi_e)];
        
        Ph{it} = Phi_it;
        
        Theta_it = [0 ; Phi_it(:,2:end)\cell2mat(Y)];
        
        J_it(it) = (cell2mat(Y) - Phi_it*Theta_it)'*...
                     (cell2mat(Y) - Phi_it*Theta_it)/length(cell2mat(Y));
        
        Th{it} = Theta_it(1:length(Theta));
        Th_e{it} = Theta_it(length(Theta)+1:end);
        
        figure(10);
        subplot(2,1,2);
        plot(cell2mat(Y),'LineWidth',1.2); hold on; 
        plot(Phi_it(:,1:length(Theta))*Th{it},'r--','LineWidth',1.2);
        hold off;
        subplot(2,1,1); hold on;
        plot(Phi_it(:,length(Theta)+1:end)*Th_e{it},'r.','LineWidth',1.2);
        hold off;
        

        Ep = cell2mat(Y) - Phi_it*Theta_it;
    end
                     
    figure(19);
    plot(1:50,J_it);
    
    [~,index] = min(J_it);
    
    m = size(cell2mat(Phi),2);
    figure(20);
    plot(cell2mat(Y)); hold on;
    plot(Phi_it*Theta_it);
    plot(cell2mat(Phi)*Th{index}(1:m));
    
    Theta = Th{itt};
    
    %% Introduction to the regularization
%    
%    reg = logspace(0,10,1000);
%    
%    Phi_t = Phi_t(:,2:end);
%    
%    I = eye(size(Phi_t,2));
%    Jreg = zeros(length(reg),1);
%    for k = 1 : length(reg)
%                            
%      Theta = (Phi_t'*Phi_t + reg(k)*I)\Phi_t'*cell2mat(Y);                    
%      Theta = [0; Theta];                      
%      
%      L = Theta(1:4);
%      B = cell(wind_s,1);
%      for i = 1 : wind_s
%         B{i}(1) = Theta(5); 
%         B{i}(2) = Theta(6+3*(i-1)); 
%         B{i}(3) = Theta(7+3*(i-1));
%         B{i}(4) = Theta(8+3*(i-1));
%      end
%      A = struc.A + L*struc.C;
%      
%      par = par_retr(A,B,dte);
%      
%      y_c = cell(length(dte.valid.y),1);
%      for w = 1 : length(dte.valid.y)
%         time.init = dte.valid.t{w}(1);
%         time.final = dte.valid.final(w);
%         initial = dte.valid.y{w}(1);
% 
%         dts = sim_system(par,time,initial);
%         
%         ind = zeros(length(dte.valid.y{w}),1);
%         for i = 1 : length(ind)
%             error = abs(dte.valid.t{w}(i) - dts.td{1});
%             [~, ind(i)] = min(error);
%         end
%         
%         y_c{w} = dts.yd{1}(ind);
%      end
%      E = cell2mat(y_c) - cell2mat(dte.valid.y);
%      
%      Jreg(k) = (E'*E)/length(E);
%    end 
%    
%    figure(15); 
%    plot(reg,Jreg); hold on; set(gca, 'XScale', 'log', 'YScale', 'log')
%    
%    [~, ind_reg] = min(Jreg);
%    
%    reg_candy = reg(ind_reg);
%    
%    Theta = (Phi_t'*Phi_t + reg_candy*I)\Phi_t'*cell2mat(Y);                    
%    Theta = [0; Theta];                      
   
   %% IV intro                        
    
%     IV = cell2mat(Phi_iv)*Theta;
%     
%     figure(3); hold on;
%     plot(cell2mat(dte.y),'LineWidth',1.2,'Color',[.8 .8 .8]);
%     plot(IV,'r--');
    
%     Phi_iv = Phi;
%     in = 1;
%     
%     for w = 1 : length(dte.y)
%         
%         out = in + length(dte.y{w}) - 1;
%         iv_sig = IV(in:out);
%         in = out + 1;
%         
%         t_calc = dte.t{w}-dte.t{w}(1);
%         [sig_out,~] =  sim_vs(reg_filt,0,iv_sig,t_calc,version);
%         
%         Phi_iv{w,1} = [zeros(length(sig_out(2:end,1)),1),sig_out(2:end,2:4)];
%     end
%     
%     Phi_iv = cell2mat(Phi_iv);
%     
%     Theta = (Phi_iv(:,2:end)'*Phi_t(:,2:end))\Phi_iv(:,2:end)'*cell2mat(Y);
%     Theta = [0;Theta];
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
            %subplot(3,1,3); hold on;
            %scatter(dts.t{k-1}(end),h.y{k-1}(1,1),'xr','LineWidth',1.8);
            %scatter(dts.t{k}(1),h.y{k-1}(1,2),'or','LineWidth',1.8);
        end
        
        %% Just for DEBUG
%         figure(3); subplot(3,1,1); hold on;
%         plot(dts.t{k},dts.y{k},'k-','LineWidth',1.2);
%         subplot(3,1,2); hold on;
%         plot(dts.t{k},dts.circ{k},'k-','LineWidth',1.2);
%         subplot(3,1,3); hold on;
%         plot(dts.t{k},dts.hom{k},'k-','LineWidth',1.2);
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