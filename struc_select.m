function struc = struc_select(method,dte)

if strcmp(method,'trivial')
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
        parameters = est_regr(dte,str_att,'15');

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
        plot(wc,J,'--','LineWidth',1.5); hold on;
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
    
end

end