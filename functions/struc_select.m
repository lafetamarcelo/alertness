function [struc,dte] = struc_select(method,dte)

%load('ploting.mat','dta');

if strcmp(method,'trivial')
    
    zeta = .1;
    omegan = pi/12;
    syscpoly = poly([-1*0.0353,...
        -zeta*omegan+omegan*sqrt(zeta^2-1),...
        -zeta*omegan-omegan*sqrt(zeta^2-1)]);
    A0 = [[zeros(1,2);eye(2)], flipud(-syscpoly(2:end)')];

    alphacand = poly(eig(A0));

    struc.n = 4;
    struc.A = [zeros(1,struc.n);[eye(struc.n-1),...
                                        flipud(-alphacand(2:end)')]]; 
    struc.C = [zeros(1,struc.n-1),1];
    
elseif strcmp(method,'barycenter')
    
    str_att.n = 4;
    str_att.C = [zeros(1,str_att.n-1), 1];
    
    zeta = .0;
    wi = pi/12;
    prop = .4;
    condition = 'false';
    
    % initialize the testing parameters
    wc = linspace(0,prop*wi,15) + wi - prop*.5*wi;
    
    while strcmp(condition,'false')
    
    J = zeros(length(wc),1);
    %aux_ii = cell(length(wc),1);
    
    for k = 1 : length(wc)
        
        omegan = wc(k);
        syscpoly = poly([-1*0.0353,...
            -zeta*omegan+omegan*sqrt(zeta^2-1),...
            -zeta*omegan-omegan*sqrt(zeta^2-1)]);
        A0 = [[zeros(1,2);eye(2)], flipud(-syscpoly(2:end)')];

        alphacand = poly(eig(A0));

        str_att.n = 4;
        str_att.A = [zeros(1,str_att.n);[eye(str_att.n-1),...
                                        flipud(-alphacand(2:end)')]];
        
        %
        parameters = est_regr(dte,str_att,'16','ls');

        % Simulate the alertness level
        time.init = [dte.t{1}(1), dte.init(2:end)];
        time.final = dte.final;
        initial = dte.y{1}(1);

        dts = sim_system(parameters,time,initial);
        
        Y = cell(length(dte.y),1);
        for w = 1 : length(dte.y)
            ind = zeros(length(dte.y{w}),1);
            for j = 1 : length(dte.y{w})
                error = abs(dts.td{w}-dte.t{w}(j)); 
                [~,ind(j)] = min(error);
            end
            Y{w} = dts.yd{w}(ind);
        end
        
        error = cell2mat(dte.y) - cell2mat(Y);
        
        J(k) = error'*error;
        disp(['  Curiosity point - ',num2str(k)]);
    end
    
    mu = min(10/std(J),250);
    for ii = 1:length(wc), aux(ii) = wc(ii)*exp(-mu*J(ii)); end
	
	barywc = sum(aux)./sum(exp(-mu*J));
    
    location = abs(barywc - wc);
    [~,index] = min(location);
    
    if strcmp('on','on')
        figure(6); hold on;
        plot(wc,J','--','LineWidth',1.5); hold on;
    end
    
    if (index <= 5) || (index > 10)
        wc = wc + (barywc - wc(8));
    else
        condition = 'true';
        figure(6); hold on;
        line([barywc barywc],[min(J),max(J)],'LineWidth',1.2,...
                                                      'Color',[1 0 0]);
    end
     
    end
    %% 
    
    %%
    syscpoly = poly([-1*0.0353,...
                    -zeta*barywc+barywc*sqrt(zeta^2-1),...
                            -zeta*barywc-barywc*sqrt(zeta^2-1)]);
    A_bary = [[zeros(1,2);eye(2)], flipud(-syscpoly(2:end)')]; 
    alphabary = poly(eig(A_bary));
    
    struc.n = 4;
    struc.A = [zeros(1,struc.n);[eye(struc.n-1),...
                                        flipud(-alphabary(2:end)')]];
    struc.C = [zeros(1,struc.n-1),1];
    
elseif strcmp(method,'barycenter damped')
    
    str_att.n = 4;
    str_att.C = [zeros(1,str_att.n-1), 1];
    
    wi = pi/12;
    tau_c = linspace(15,30,30);
    
    J = zeros(length(tau_c),1);
    %aux_ii = cell(length(wc),1);
    
    for k = 1 : length(tau_c)
        
        str_att.A = [0 0 0 0;...
                     1 0 0 -wi^2/tau_c(k);...
                     0 1 0 -wi^2;...
                     0 0 1 -1/tau_c(k)];
        
        %
        parameters = est_regr(dte,str_att,'16','ls');

        % Simulate the alertness level
        time.init = [dte.t{1}(1), dte.init(2:end)];
        time.final = dte.final;
        initial = dte.y{1}(1);

        dts = sim_system(parameters,time,initial);
        
        Y = cell(length(dte.y),1);
        for w = 1 : length(dte.y)
            ind = zeros(length(dte.y{w}),1);
            for j = 1 : length(dte.y{w})
                error = abs(dts.td{w}-dte.t{w}(j)); 
                [~,ind(j)] = min(error);
            end
            Y{w} = dts.yd{w}(ind);
        end
        
        %error = cell2mat(dte.y) - cell2mat(Y);
        ye = cell2mat(dte.y); yhat = cell2mat(Y);
        %J(k) = error'*error;
        J(k) = min(norm(ye-yhat,2)/...
			norm(ye-kron(mean(ye),ones(length(ye),1)),2),1);
        disp(['Curiosity point - ',num2str(k)]);
    end
    
    mu = min(10/std(J),250);
    for ii = 1:length(tau_c); aux(ii) = tau_c(ii)*exp(-mu*J(ii)); end
	
	barytau_c = sum(aux)./sum(exp(-mu*J));
    
    if strcmp('on','on')
        figure(6); hold on;
        plot(tau_c,J,'--','LineWidth',1.5); hold on;
        plot([barytau_c barytau_c],[.001, 1],'Color',[.8 .8 .8],...
            'LineWidth',1.5); hold on;
        set(gca, 'YScale', 'log'); hold off;
    end
    
    
    %% 
    struc.n = 4;
    struc.A = [0 0 0 0;...
               1 0 0 -wi^2/barytau_c;...
               0 1 0 -wi^2;...
               0 0 1 -1/barytau_c];
    struc.C = [zeros(1,struc.n-1),1];
    
elseif strcmp(method, 'grid filtering')
    
    omegan = pi/12;
    taun = 28;
    
    str_att.n = 4;
    str_att.C = [zeros(1,str_att.n-1),1];
    
    prop = 0.5;
    
    omega_c = linspace(0,100,10);
    tau_c = linspace(.1,100,10);
    
    J = zeros(length(omega_c),length(tau_c));
    Je = zeros(length(omega_c),length(tau_c));
    for i = 1 : length(omega_c)
        for j = 1 : length(tau_c)
            str_att.A = [0 0 0 0;...
                         1 0 0 -omega_c(i)^2/tau_c(j);...
                         0 1 0 -omega_c(i)^2;...
                         0 0 1 -1/tau_c(j)];
            
            figure(4); hold on;
            sSmodel = ss(str_att.A',str_att.C',eye(4),0);
            pzmap(sSmodel);
            
            parameters = est_regr(dte,str_att,'16','ls');

            % Simulate the alertness level
            time.init = [dte.t{1}(1), dte.init(2:end)];
            time.final = dte.final;
            initial = dte.y{1}(1);

            dts = sim_system(parameters,time,initial);
            
            figure(5);
            Y = cell(length(dte.y),1);
            for w = 1 : length(dte.y)
                ind = zeros(length(dte.y{w}),1);
                for k = 1 : length(dte.y{w})
                    error = abs(dts.td{w}-dte.t{w}(k)); 
                    [~,ind(k)] = min(error);
                end
                Y{w} = dts.yd{w}(ind);
                
                % Ploting results
                plot(dta.td{w},dta.yo{w},'Color',[.6 .6 .6],...
                                                'LineWidth',1.2); hold on;
                if w ~= length(dta.yo)
                    plot(dta.tn{w},dta.yn{w},'--','Color',[.6 .6 .6],...
                                                          'LineWidth',1.2);
                end
                plot(dte.t{w},Y{w},'b-x','LineWidth',1.4)    
            end
            hold off;
            
            error = cell2mat(dte.y) - cell2mat(Y);
            Je(i,j) = error'*error;
            
            ye = cell2mat(dte.y); yhat = cell2mat(Y);
            J(i,j) = min(norm(ye-yhat,2)/...
                norm(ye-kron(mean(ye),ones(length(ye),1)),2),1);
            disp(['Curiosity point - ',num2str(i),' -- ',num2str(j)]); 
        end
    end
    
    figure(5); 
    surf(omega_c,tau_c,J);
    
    figure(6);
    surf(omega_c,tau_c,Je);
    
elseif strcmp(method, 'sGolay')
    
    framelen = 11;
    order = 3;
      
    sig = cell(length(dte.y),1);
    x = cell(length(dte.y),1);
    steady = cell(length(dte.y),1);
    cmplt = cell(length(dte.y),1);
    m = (framelen-1)/2;
    B = sgolay(order,framelen);
    
    for i = 1 : length(dte.y)
        x{i} = detrend(dte.y{i});
        trend = dte.y{i} - x{i};
        sig{i} = sgolayfilt(x{i}-x{i}(1),3,11) + x{1}(1) + trend;
        
        steady{i} = conv(x{i}-x{i}(1),B(m+1,:),'same') + x{i}(1) + trend;
        
        ybeg = B(1:m,:)*x{i}(1:framelen);
        
        lx = length(x{i});
        yend = B(framelen-m+1:framelen,:)*x{i}(lx-framelen+1:lx);

        cmplt{i} = steady{i};
        cmplt{i}(1:m) = ybeg + trend(1:m);
        cmplt{i}(lx-m+1:lx) = yend + trend(lx-m+1:lx);
    end
    
    figure(2); hold on;
    scatter(cell2mat(dte.t),cell2mat(dte.y),25,[.8 .8 .8]);
    for i = 1 : length(dte.y)
        %plot(dte.t{i},sig{i},'k-','LineWidth',1.2); hold on;
        %plot(dte.t{i},steady{i},'r--','LineWidth',1.2);
        scatter(dte.t{i},cmplt{i},'bx','LineWidth',1);
        dte.y{i} = cmplt{i};
    end
    %legend([{'Real'},{'samples'},{'sGolay'},{'Steady'}])
    
    fs1 = 1;  fs2 = 2;
    [p,q] = rat(fs2/fs1);
    normFc = .98 / max(p,q);
    order = 256 * max(p,q);
    beta = 12;

    lpFilt = firls(order, [0 normFc normFc 1],[1 1 0 0]);
    lpFilt = lpFilt .* kaiser(order+1,beta)';
    lpFilt = lpFilt / sum(lpFilt);

    % multiply by p
    lpFilt = p * lpFilt;
    
    for i = 1 : length(dte.y)
        
        %detrending data
        a(1) = (cmplt{i}(end)-cmplt{i}(1)) / (dte.t{i}(end)-dte.t{i}(1));
        a(2) = cmplt{i}(1);

        % detrend the signal
        xdetrend = cmplt{i} - polyval(a,dte.t{i}-dte.t{i}(1));
        [sig{i},time] = resample(xdetrend,dte.t{i}-dte.t{i}(1),...
                                                   fs2,p,q,lpFilt);
        sig{i} = sig{i} + polyval(a,time);
        time = time + dte.t{i}(1);
        
        figure(2); hold on;
        plot(time, sig{i},'r-','LineWidth',1.2);
        dte.y{i} = sig{i}; dte.t{i} = time;
    end
    
    hold off;
    
    zeta = .0;
    omegan = pi/12;
    syscpoly = poly([-1*0.0353,...
        -zeta*omegan+omegan*sqrt(zeta^2-1),...
        -zeta*omegan-omegan*sqrt(zeta^2-1)]);
    A0 = [[zeros(1,2);eye(2)], flipud(-syscpoly(2:end)')];

    alphacand = poly(eig(A0));

    struc.n = 4;
    struc.A = [zeros(1,struc.n);[eye(struc.n-1),...
                                        flipud(-alphacand(2:end)')]]; 
    struc.C = [zeros(1,struc.n-1),1];
    
    
end

end