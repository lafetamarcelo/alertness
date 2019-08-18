function parameters = est_regr(dte,struc,version,method)
                                                                       
    reg_filt = ss(struc.A',struc.C',eye(4),0);
    wind_s = length(dte.y);
    
    Phi = cell(wind_s,2+wind_s);
    Y = cell(wind_s,1);
    t_iv = cell(wind_s,1);
    
    %% Regression matrix
    
    for w = 1 : wind_s
        
        % Impulse filtering - Input [B] informations
        delt = dte.t{w}(end)-dte.t{w}(1);
        step_t = .0001;
        puls_time = 0:step_t:delt;
        puls_out = impulse(reg_filt,puls_time);
        
        ind = zeros(length(dte.t{w}),1);
        for i = 1 : length(dte.t{w})
           error = abs(puls_time - (dte.t{w}(i) - dte.t{w}(1)));
           [~,ind(i)] = min(error);
           
           % simulated time error plotting
           %error_graph(i) = min(error);
        end
        
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
           %error_graph(i) = error(ind(i));
        end
        
        %figure(2); subplot(2,1,2); hold on;
        %plot(error_graph,'Color',[.4 .4 .4]);
        
        Phi{w,1} = [zeros(length(ind(2:end)),1),-sig_out(ind(2:end),2:4)];
        
        Y{w} = dte.y{w}(2:end);
        t_iv{w} = dte.t{w}(2:end);
        
        % Shifting the window information
        for j = 1 : wind_s
            if w == j
                Phi{w,j+2} = [X(:,1), -X(:,2), X(:,3)];
            else
                Phi{w,j+2} = zeros(size(Phi{w,1},1),3);
            end
        end
        
    end

    %% Determine the parameters
    
    Phi_t = cell2mat(Phi);
    
    Theta = lsqnonneg(Phi_t,cell2mat(Y));
    
    %Theta = [0;Phi_t(:,2:end)\cell2mat(Y)];
    %Theta = Phi_t\cell2mat(Y);
    
    J = (cell2mat(Y) - cell2mat(Phi)*Theta)'*(cell2mat(Y) ...
                                - cell2mat(Phi)*Theta)/length(cell2mat(Y));
    
    E = (cell2mat(Y) - cell2mat(Phi)*Theta)'*(cell2mat(Y) ...
                                                    - cell2mat(Phi)*Theta);
                            
    L_ = Theta(1:4);
    
    B_ = cell(wind_s,1);
    for i = 1 : wind_s
        B_{i}(1) = Theta(5); 
        B_{i}(2) = Theta(6+3*(i-1)); 
        B_{i}(3) = -Theta(7+3*(i-1));
        B_{i}(4) = Theta(8+3*(i-1));
    end
    
    A_ = struc.A - L_*struc.C;
    
    par_b = par_retr(A_,B_,dte);
    
    %% Introduction to residual modeling 
    
    if strcmp(method,'extended')
    
    Ep = cell2mat(Y) - Phi_t*Theta;
    
    for iter = 1 : 50
    
        ind_i = 1;
        Phi_e = cell(length(dte.y),1);

        for w = 1 : length(dte.y)

            t_calc = dte.t{w}(2:end) - dte.t{w}(2);
            ind_f = ind_i + length(t_calc) - 1;
            %[sig_out,~] =  sim_vs(reg_filt,0,Ep(ind_i:ind_f),...
            %                                               t_calc,version);

            e = Ep(ind_i:ind_f);
            sig_out = [[zeros(1,1); e(1:end-1)], ...
                        [zeros(2,1); e(1:end-2)]];%,...
                         %[zeros(3,1); e(1:end-3)],...
                         %[zeros(4,1); e(1:end-4)]];

            ind_i = ind_f + 1;

            Phi_e{w} = sig_out;%(:,1:4);
        end

        Phi_it = [cell2mat(Phi) , cell2mat(Phi_e)];

        Theta_it = [0 ; Phi_it(:,2:end)\cell2mat(Y)];

        J_it(iter,1) = (cell2mat(Y) - Phi_it*Theta_it)'*...
                       (cell2mat(Y) - Phi_it*Theta_it)/length(cell2mat(Y));
        J_it(iter,2) = (cell2mat(Y) - cell2mat(Phi)*Theta_it(1:length(Theta)))'*...
                       (cell2mat(Y) - cell2mat(Phi)*Theta_it(1:length(Theta)))/length(cell2mat(Y));
                 
        Ep = cell2mat(Y) - Phi_it*Theta_it;
    
    end   
   
    Theta = Theta_it(1:length(Theta));
    
    end

   %% IV intro                        
   if strcmp(method,'instrumental')
       
       y_iv = cell2mat(Phi)*Theta;
       y_reg = cell2mat(Y);
       Phi_iv = Phi;
       
       tstp = cell2mat(t_iv);
       
       for w = 1 : length(dte.y)
           
            t_calc = dte.t{w}(2:end)-dte.t{w}(2);
            index = find((tstp >= t_iv{w}(1)).*(tstp <= t_iv{w}(end)));
            Y{w} = y_reg(index(2:end));
            
            [sig_out,t_out] =  sim_vs(reg_filt,0,y_iv(index),...
                t_calc,version);
            %sig_out = lsim(reg_filt,dte.y{w},puls_time(ind));
            
            ind = zeros(length(dte.t{w})-1,1);
            for i = 2 : length(dte.t{w})
               error = abs(t_out - (dte.t{w}(i) - dte.t{w}(2)));
               [~,ind(i-1)] = min(error);

               % simulated time error plotting
               %error_graph(i) = error(ind(i));
            end

            Phi_iv{w,1} = [zeros(length(ind(2:end)),1),...
                                                sig_out(ind(2:end),2:4)];
            Phi{w,1} = Phi{w,1}(2:end,:);                                
            
            
            
            for j = 2 : size(Phi_iv,2)
               Phi_iv{w,j} = Phi_iv{w,j}(2:end,:);
               Phi{w,j} = Phi{w,j}(2:end,:);
            end
            
       end
       
       Z = cell2mat(Phi_iv);
       Phi_t = cell2mat(Phi);
       
       Theta = (Z(:,2:end)'*Phi_t(:,2:end))\Z(:,2:end)'*cell2mat(Y);
       
       Theta = [0; Theta];
       %Theta = (Z'*Phi_t)\Z'*cell2mat(Y);
   end
    

   %% Parameters matrix reconstruction
    
    L = Theta(1:4);
    
    B = cell(wind_s,1);
    for i = 1 : wind_s
        B{i}(1) = Theta(5); 
        B{i}(2) = Theta(6+3*(i-1)); 
        B{i}(3) = -Theta(7+3*(i-1));
        B{i}(4) = Theta(8+3*(i-1));
    end
    
    A = struc.A - L*struc.C;
    
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
        end
    end
    
    %% Estimate the night parameters
    
    % Estimate options
    maxITER = 200;
    options = optimoptions(@lsqnonlin,'Jacobian','on',...
        'Display','off','TolFun',1e-6,'MaxIter',maxITER,'TolX',1e-6);
    initial = [14,2.5];
    
    func = @cost_function;
    y_lsq = cell2mat(h.y);
    lb = [8 , 1.8]; ub = [28, 3.8];
    
    [x,~,~,~,~] = lsqnonlin(func,initial,lb,ub,options,...
           cell2mat(h.t),y_lsq(:,1),y_lsq(:,2));
    
    par.est.y_oo = x(1); par.est.tau_e = x(2);
    
    parameters = par;
    
end