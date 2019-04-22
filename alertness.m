%% Create an alertness class for estimation and prediction:
% Inputs:   - estimation data dte:: y{cells}, t{cells}
%           - simulation data dts:: initial{array}, t{cells}
% Outputs:  - model :: physical_parameters, dts, structure

classdef alertness < handle
    
   properties % Define the public information
    %Initialize the parameters such as Folkard models
    omega = pi/12; cTau = 1/0.0353; M = 2.52; 
    cPhase = -16.835*pi/12; dc = 2.4; 
    y_ = 14.3; eTau = 1/0.381;
    % Determine the basic data sets used
    dte %The estimation data set
    dtv %The validation data set
    dts %The simulated data set
    % Determine the pre defined structure
    A = [zeros(1,4);[eye(3), zeros(3,1)]]; % observer matrix
    C = [zeros(1,3), 1]; % observer matrix
    Ao % The dynamic constant alertness matrix
    % Set the first estimation information
    structure = 'trivial'; % select the structure selector
    identification = 'least squares'; % select the identification technique
    algorithm = 'MOLI'; % select the algorithm to be used
   end
    
   properties (Access = private)
    degree = 4;
    % Possible structure selectors
    strucPoss = [{'Barycenter'}, {'Let loose'}, {'sGolay filtering'},...
        {'Grid search'}, {'Folkard'}];
    % Possible identification algorithms
    identPoss = [{'least squares'}, {'instrumental variable'},...
        {'extended'}, {'LPV'}, {'N.N.L.S.'}, {''}];
    % Possible algorithms
    algPoss = [{'MOLI'}, {'Drived search'}, {'Folkard'}];
    % Model parameters
    L % The output regressor MOLI parameter matrix
    B % the input regressor MOLI parameter matrix
    % Regression model data
    theta % The regressor parameters
    Phi % The regressor information
    Y % Considered output regressor
   end
   
   methods
       % Initialize the alertness parameters 
       function obj = alertness(struc, ident, alg)
          if ismember(struc, obj.strucPoss)
            obj.structure = struc;
          else
            error('Structure model is not possible.');
          end
          if ismember(ident, obj.identPoss)
            obj.identification = ident;
          else
            error('Identification technique not possible.');
          end
          if ismember(alg, obj.algPoss)
            obj.algorithm = alg;
          else
            error('Algorithm does not exist!')
          end    
       end
       
       %% Flow control functions
       function obj = fit(obj, dte, dtv)
           % Verify the input and validation data-set
           %obj.check(dte); obj.check(dtv);
           obj.dte = dte; obj.dtv = dtv;
           
           % Select the model structure
           obj.strucSelect(obj.dte);
           
           % Estimate the day parameters
           obj.estimate(obj.dte);
           
           % Validate the model
           obj.dts = obj.simulate(obj.dtv);
       end
       
       % Select the model structure
       function obj = strucSelect(obj, dte)
          disp('Initialize the observer filter search...')
          % Determine the observer structure
          if strcmp('Barycenter', obj.structure)
            obj.barycenter(dte);
          elseif strcmp('Folkard', obj.structure)
            obj.trivial_struc();
          elseif strcmp('sGolay filtering', obj.structure)
            obj.sGolayFlow(dte);
          elseif strcmp('Grid search', obj.structure)
            obj.gridSearch(dte);
          elseif strcmp('Let loose', obj.structure)
            % Does nothing!!
          else
            error('Observer search technique does not match the possibilities'); 
          end
          disp('Done!')
       end
       % Estimate the day parameters
       function obj = estimate(obj, dte)
           % Select the algorithm to be used
           if strcmp('MOLI', obj.algorithm)
               obj.MOLI(dte);
               obj.nonlinear_estimate(obj.dte);
           elseif strcmp('Drived search', obj.algorithm)
               obj.drivedSearch(dte);
               obj.nonlinear_estimate(obj.dte);
           elseif strcmp('Folkard', obj.algorithm)
               obj.Ao = obj.A;
           else
              error('No algorithm with such name!');
           end
       end
       %Estimate the night parameters
       function obj = nonlinear_estimate(obj, dte)
           wSize = length(dte.y); h = zeros(wSize-1, 3);
           
           for wii = 1 : wSize - 1
            % Simulate the system backwards
            if dte.t{wii+1}(1) ~= dte.init(wii+1)
                tOmega = obj.omega*dte.t{wii+1}(1);
                k1 = obj.M*cos(tOmega + obj.cPhase);
                k2 = -obj.M*sin(tOmega + obj.cPhase)*obj.omega;
                h0 = dte.y{wii+1}(1) - k1 - obj.dc;
                Beta = [(obj.dc*obj.omega^2)/obj.cTau;...
                    (k2/obj.cTau + h0*obj.omega^2 + obj.dc*obj.omega^2);...
                    (k1 + obj.dc + k2*obj.cTau)/obj.cTau;...
                    h0 + k1 + obj.dc];
                
                alertnessBack = ss(-obj.Ao, Beta, obj.C, 0);
                finalTime = dte.t{wii+1}(1) - dte.init(wii+1);
                [initA, ~] = impulse(alertnessBack, finalTime);
            else
                initA = flip(dte.y{wii+1});
            end
            % Simulate the system forward
            if dte.t{wii}(end) ~= dte.final(wii)
                tOmega = obj.omega*dte.t{wii}(end);
                k1 = obj.M*cos(tOmega + obj.cPhase);
                k2 = -obj.M*sin(tOmega + obj.cPhase)*obj.omega;
                h0 = dte.y{wii}(end) - k1 - obj.dc;
                Beta = [(obj.dc*obj.omega^2)/obj.cTau;...
                    (k2/obj.cTau + h0*obj.omega^2 + obj.dc*obj.omega^2);...
                    (k1 + obj.dc + k2*obj.cTau)/obj.cTau;...
                    h0 + k1 + obj.dc];
                
                alertnessForw = ss(obj.Ao, Beta, obj.C, 0);
                finalTime = dte.final(wii) - dte.t{wii}(end);
                [finalA, ~] = impulse(alertnessForw, finalTime);
            else
                finalA = dte.y{wii};
            end
            % Atribute the vicinity data
            h(wii,1) = initA(end) - obj.M*cos(obj.omega*dte.init(wii+1) ...
                       + obj.cPhase); 
            h(wii,2) = finalA(end) - obj.M*cos(obj.omega*dte.final(wii) ...
                       + obj.cPhase);
            h(wii,3) = dte.init(wii+1) - dte.final(wii);
           end
           
           %Estimate night parameters
           maxITER = 200;
           options = optimoptions(@lsqnonlin, 'Jacobian', 'on',...
               'Display', 'off', 'TolFun', 1e-6, 'MaxIter', maxITER,...
               'TolX', 1e-6);
           initial = [14, 2.5];
           %optimFunc = obj.cost_function;
           % Setting the upper and lower bound for search
           lb = [7 , 1.8]; ub = [17 , 3.8];
           % Estimate with non linear least squares
           [x,~,~,~,~] = lsqnonlin(@obj.costFunc, initial, lb, ub,...
                                options, h);
           % Set the obtained parameters
           obj.y_ = x(1); obj.eTau = x(2);
       end
       % Simulate the alertness model
       function dt = simulate(obj, dts)
           wSize = length(dts.time);
           
           initialHom = dts.initial;
           
           for wii = 1 : wSize
               tOmega = obj.omega*dts.time(wii, 1);
               k1 = obj.M*cos(tOmega + obj.cPhase);
               k2 = -obj.M*sin(tOmega + obj.cPhase)*obj.omega;
               h0 = initialHom - k1 - obj.dc;
               Beta = [(obj.dc*obj.omega^2)/obj.cTau;...
                    (k2/obj.cTau + h0*obj.omega^2 + obj.dc*obj.omega^2);...
                    (k1 + obj.dc + k2*obj.cTau)/obj.cTau;...
                    h0 + k1 + obj.dc];
               alertnessModel = ss(obj.Ao, Beta, obj.C, 0);
               
               simT = linspace(0, dts.time(wii,2)-dts.time(wii,1), 1000)';
               [dt.dayOut{wii}, simTime] = impulse(alertnessModel, simT);
               dt.dayTime{wii} = simTime + dts.time(wii, 1);
               dt.dayHom{wii} = dt.dayOut{wii} - ...
                         obj.M*cos(obj.omega*dt.dayTime{wii} + obj.cPhase);
               
               if wii ~= wSize
                  hom(wii,1) = dt.dayHom{wii}(end); 
                  
                  simT = linspace(0, dts.time(wii+1,1)-dts.time(wii,2), 100)';
                  dt.nightHom{wii} = obj.y_*(1-exp(-simT./obj.eTau)) + ...
                        hom(wii,1)*exp(-simT./obj.eTau);
                  dt.nightTime{wii} = simT + dts.time(wii,2);
                  dt.nightOut{wii} = dt.nightHom{wii} + ...
                      obj.M*cos(obj.omega*dt.nightTime{wii} + obj.cPhase);
                  
                  initialHom = dt.nightOut{wii}(end);
               end 
               
           end
       end
       
       %% Parameter estimation algorithms
       % MOLI algorithms
       function obj = MOLI(obj, dte)
       % Estimate the parameter vectors L and B
           wSize = length(dte.y);
           % Determine the regressor
           obj.regMOLI(dte);
           
           % Solve the regressor problem
           disp('Solving the regression problem...');
           
           if strcmp('N.N.L.S.', obj.identification)
               obj.theta = lsqnonneg(obj.Phi,obj.Y); 
           elseif strcmp('least squares', obj.identification)
               obj.theta = obj.Phi\obj.Y;
           elseif strcmp('intrumental variable', obj.identification)
               obj.theta = intrumVar(dte);
           elseif strcmp('extended L.S.', obj.identification)
               obj.theta = extendedLS(dte);
           else
               error('No possible identification approach!')
           end
           
           % Determine the MOLI parameters matrices L[-] -- B[-]
           obj.L = obj.theta(1:4);
           for k = 1 : wSize
               obj.B{k}(1) = obj.theta(5);
               obj.B{k}(2) = obj.theta(6 + (k-1)*3);
               obj.B{k}(3) = -obj.theta(7 + (k-1)*3);
               obj.B{k}(4) = obj.theta(8 + (k-1)*3);
           end
           obj.Ao = obj.A - obj.L * obj.C;
           
           % Retrieve the physical parameters
           obj.parameters_retrieve(dte);
       end
       
       function obj = regMOLI(obj, dte)
       % Create the MOLI regressor        
           % Create the Morse Observer
           observerMorse = ss(obj.A',obj.C',eye(4),0);
           % Determine the constants and initialize variables
           wSize = length(dte.y); obj.Phi = cell(wSize, wSize+2);
           selectOut = cell(wSize, 1);
           for wii = 1 : wSize
               refTime = dte.t{wii} - dte.t{wii}(1);
               % Observe the input impulse response -- B[.] info
               impTime = 0:0.0001:refTime(end);
               impStates = impulse(observerMorse, impTime);
               % Error analysis and select the simulation data
               impInd = zeros(length(dte.t{wii}),1);
               for k = 1 : length(impInd)
                   e = abs(impTime - (dte.t{wii}(k) - dte.t{wii}(1)));
                   [~, impInd(k)] = min(e);
               end
               
               % Observe the output filtered response -- L[.] info
               [outStates, outTime] = obj.simuLinkSolve(observerMorse,...
                                                   0, dte.y{wii}, refTime);
               outInd = zeros(length(dte.t{wii}),1);
               for k = 1 : length(outInd)
                   e = abs(outTime - (dte.t{wii}(k) - dte.t{wii}(1)));
                   [~, outInd(k)] = min(e);
               end
               
               % Do not provide the first simulated data
               outInd = outInd(2:end); impInd = impInd(2:end);
               
               % Create the regressor matrix
               obj.Phi{wii, 1} = -outStates(outInd, 1:4); % L[.] info
               obj.Phi{wii, 2} = impStates(impInd,1); % DC info
               
               for k = 1 : wSize % B[.] info
                    if wii == k
                        obj.Phi{wii,k+2} = [ impStates(impInd,2), ...
                                            -impStates(impInd,3), ...
                                             impStates(impInd,4)];
                    else
                        obj.Phi{wii,k+2} = zeros(length(outInd), 3);
                    end
               end
               % Estipulate the needed output data
               selectOut{wii} = dte.y{wii}(outInd);
           end
           
           % Build the regressor matrix
           obj.Phi = cell2mat(obj.Phi);
           obj.Y = cell2mat(selectOut);
       end
       
       function obj = parameters_retrieve(obj, dte)
       % Determine the physical / biological parameters    
           wSize = length(dte.y); h0 = zeros(wSize, 1);
           M_ = zeros(wSize,1); cPhase_ = zeros(wSize,1);
           
           eigA = eig(obj.Ao);
           
           % Use the Folkard information to determine \omega and \tau_c
           tauCand = -1./real(eigA); tauCand = tauCand(abs(tauCand) > 1e-02); 
           omegaCand = imag(eigA); omegaCand = omegaCand(abs(omegaCand) > 1e-02);
           tauError = abs(obj.cTau - tauCand);
           omegaError = abs(obj.omega - omegaCand);
           [~, tauInd] = min(tauError);
           [~, omegaInd] = min(omegaError);
           
           obj.cTau = tauCand(tauInd);
           obj.omega = omegaCand(omegaInd);
           
           if isempty(obj.cTau), obj.cTau = 0.01; end
           if isempty(obj.omega), obj.omega = 0.01; end
           
           %obj.cTau = -1/real(eigA(3));
           %obj.omega = imag(eigA(1));
           obj.dc = obj.B{1}(1)*obj.cTau/obj.omega^2;
           
           % Set the regressor for retrieving
           infoX = [obj.omega^2, 0, 1/obj.cTau; ...
                    0, 1/obj.cTau, 1; ...
                    1, 1, 0];
           compDC = [-obj.dc*obj.omega^2; ... 
                     -obj.dc/obj.cTau; ...
                     -obj.dc];
           for k = 1 : wSize
              compB = obj.B{k}(2:end)' + compDC; 
              aux = infoX\compB;
              h0(k) = aux(1);
              
              cPhaseTan = aux(3)/(-obj.omega)/aux(2);
              
              if (cPhaseTan < 0)
                 cPh = atan(cPhaseTan)-pi;
              else
                 cPh = atan(cPhaseTan);
              end
              M_(k) = aux(2)/cos(cPh);
              cPhase_(k) = cPh - obj.omega*(dte.t{k}(1) - 24*(k-1));
           end
           obj.M = mean(M_);
           obj.cPhase = mean(cPhase_);
       end
       
       % Guided search approach
       function obj = drivedSearch(obj, dte)
           obj.Ao = obj.A;
           
           % First shift the phase 
           obj.shiftPhase(dte);
           
       end
       
       function obj = shiftPhase(obj, dte)
       % Shifts the algorithm phase
       
        % Manipulating the time data for simulation
        dtSim.time = [[dte.t{1}(1); dte.init(2:end)'], dte.final'];
        dtSim.initial = dte.y{1}(1);
        dtSim.y = dte.y; dtSim.t = dte.t;
        
        % Determine the possible phase values 
        searchQuant = 50;
        phasePoss = linspace(-2*pi, 2*pi, searchQuant);
        
        J = zeros(searchQuant, 1);
        for pii = 1 : searchQuant
           obj.cPhase = phasePoss(pii);
            
           % Simulate the output
           dt = obj.simulate(dtSim);
           
           % Check the fitness
           J(pii) = obj.fitnessStruc(dtSim, dt, 'squared error');
           
           figure(3); hold on;
           scatter(cell2mat(dte.t), cell2mat(dte.y), 'ko', 'LineWidth', 1.4);
           for i = 1 : length(dte.y)
            plot(dt.dayTime{i}, dt.dayOut{i}, 'r-', 'LineWidth', 1.4);
            if i ~= length(dte.y)
             plot(dt.nightTime{i}, dt.nightOut{i}, 'b-', 'LineWidth', 1.4);
            end
           end
           hold off;
           
           disp(['Curiosity point - ', num2str(pii)]);
        end
        
        figure(5); hold on; 
        plot(phasePoss, J, 'k-', 'LineWidth', 1.4);
        hold off;
        
        [~, indP] = min(J); obj.cPhase = phasePoss(indP);
       end
       %% Filtering/Observer methods
       % Barycenter technique
       function obj = barycenter(obj, dte)
          
          obj.degree = 4; 
          % Manipulating the time data for simulation
          dtSim.time = [[dte.t{1}(1); dte.init(2:end)'], dte.final'];
          dtSim.initial = dte.y{1}(1);
          dtSim.y = dte.y; dtSim.t = dte.t;
          
          % Determine the curiosity points
          f1 = 2/(2*pi); f2 = 0.5/(2*pi);
          fc = linspace(f1, f2, 30); % 2. until 80. Hertz -- step 5. Hz
          wc = 2*pi.*fc; % transforming into radians
          
          % Determine the butterworth poles
          rBut = zeros(obj.degree, 1); iBut = zeros(obj.degree, 1); 
          for k = 1 : obj.degree
            rBut(k) = -sin(.5*pi*(2*k-1) / obj.degree);
            iBut(k) = cos(.5*pi*(2*k-1) / obj.degree);
          end
          polesBut = complex(rBut, iBut);
          polyBut = real(poly(polesBut));
          
          % Test each cut--off frequency
          wcSize = length(wc); J = zeros(wcSize,1);
          for wii = 1 : wcSize
              wCand = wc(wii); wProp = zeros(obj.degree, 1);
              for k = 1 : obj.degree
                  wProp(k) = wCand^k;
              end
              
              alphaCand = -flipud(polyBut(2:end)'.*wProp);
              obj.C = [zeros(1, obj.degree-1), wProp(end)];
              obj.A = [[zeros(1, obj.degree-1); eye(obj.degree-1)],alphaCand];
              
              % Estimate the model
              obj.estimate(dte);
              obj.nonlinear_estimate(dte);
              
              % Simulate the data
              dtb = obj.simulate(dtSim);
              
              % Determine the cost funtion
              J(wii) = obj.fitnessStruc(dtSim, dtb, 'squared error');
              
              disp(['Curiosity point ', num2str(wii)]);
          end
          
          mu = min(10/std(J), 250);
          for ii = 1:wcSize, aux(ii) = wc(ii)*exp(-mu*J(ii)); end
          
          barywc = sum(aux)./sum(exp(-mu*J));
          
          figure(4); hold on; 
          plot(wc, J);
          scatter(barywc, 1, 'rx', 'LineWidth', 1.6);
          set(gca, 'YScale', 'log');
          
          for k = 1 : obj.degree
             wProp(k) = barywc^k;
          end
            
          alpha = -flipud(polyBut(2:end)'.*wProp);
          obj.C = [zeros(1, obj.degree-1), wProp(end)];
          obj.A = [[zeros(1, obj.degree-1); eye(obj.degree-1)],alpha];
          
       end
       
       % sGolay technique
       function obj = sGolayFlow(obj, dte)
       
           framelen = 11; order = 3;

           % Resampling the data 
           dtr = obj.resampleData(dte);

           % Determining the sGolay filtered out 
           dtsG = obj.sGolay(dtr, framelen, order);

           % Estimate the model parameters
           obj.estimate(dtsG);

           % NonLinear estimate
           obj.nonlinear_estimate(dtsG);
           
           obj.dte = dtsG;

           % Manipulating the time data for simulation
           dtSim.time = [[dte.t{1}(1); dte.init(2:end)'], dte.final'];
           dtSim.initial = dte.y{1}(1);
           dtSim.y = dte.y; dtSim.t = dte.t;

           % Simulate the data
           obj.dts = obj.simulate(dtSim);
       
       end
       
       function dtsG = sGolay(obj, dte, framelen, order)
        wSize = length(dte.y); dtsG = dte; 
        
        gSignal = cell(wSize,1);
        steady = cell(wSize, 1);
        cmplt = cell(wSize, 1);
        
        x = dte.y;
        m = (framelen - 1) / 2;
        sG = sgolay(order, framelen);
            
        for i = 1 : wSize
            gSignal{i} = sgolayfilt(x{i}, order, framelen);

            steady{i} = conv(x{i}, sG(m+1,:), 'same');

            ybeg = sG(1:m,:) * x{i}(1:framelen);

            lx = length(x{i});
            yend = sG(framelen-m+1:framelen, :) * x{i}(lx-framelen+1 : lx);

            cmplt{i} = steady{i};
            cmplt{i}(1 : m) = ybeg;
            cmplt{i}(lx-m+1 : lx) = yend;
            
            dtsG.t{i} = dte.t{i}(dte.t{i} <= obj.dte.t{i}(end)); 
            dtsG.y{i} = cmplt{i}(1 : length(dtsG.t{i}));
        end 
       
       figure(1); hold on;
       plot(cell2mat(dte.t), cell2mat(cmplt), 'k-', 'LineWidth', 1.2);
       plot(cell2mat(dte.t), cell2mat(gSignal), 'r--', 'LineWidth', 1.2);
       hold off;
       
       end
       
       function dtr = resampleData(~, dte)
        wSize = length(dte.y); dtr = dte;
        
        % Prepare the filter structure
        fs1 = 1;  fs2 = 2;
        [p, q] = rat(fs2/fs1);
        normFc = .98 / max(p, q);
        order = 256 * max(p, q);
        beta = 12;

        lpFilt = firls(order, [0 normFc normFc 1],[1 1 0 0]);
        lpFilt = lpFilt .* kaiser(order + 1, beta)';
        lpFilt = lpFilt / sum(lpFilt);

        % multiply by p
        lpFilt = p * lpFilt;
        
        sigRes = cell(wSize, 1);
        for i = 1 : wSize
            %detrending data
            a(1) = (dte.y{i}(end)-dte.y{i}(1)) / (dte.t{i}(end)-dte.t{i}(1));
            a(2) = dte.y{i}(1);

            % detrend the signal
            xdetrend = dte.y{i} - polyval(a,dte.t{i}-dte.t{i}(1));
            [sigRes{i},time] = resample(xdetrend, dte.t{i}-dte.t{i}(1),...
                                                       fs2, p, q, lpFilt);
            sigRes{i} = sigRes{i} + polyval(a,time);
            time = time + dte.t{i}(1);

            figure(1); hold on;
            plot(time, sigRes{i},'b--','LineWidth',1.2);
            dtr.y{i} = sigRes{i}; dtr.t{i} = time;
        end
        hold off;
       end
       
       % Folkard structure
       function obj = trivial_struc(obj) 
          alpha = [0;...
                   -obj.omega^2/obj.cTau;...
                   -obj.omega^2;...
                   -1/obj.cTau];
          obj.A = [[zeros(1,3);eye(3)], alpha]; 
       end
       
       % Grid search
       function obj = gridSearch(obj, dte)
           
           % Manipulating the time data for simulation
           dtSim.time = [[dte.t{1}(1); dte.init(2:end)'], dte.final'];
           dtSim.initial = dte.y{1}(1);
           dtSim.y = dte.y; dtSim.t = dte.t;
           
           % Determine possible parameters
           prop = .4; wQuant = 10; tQuant = 10; zeta = .0;
           wComp = prop * obj.omega; tComp = prop * obj.cTau;
           wc = linspace(obj.omega - wComp, obj.omega + wComp, wQuant);
           tc = linspace(obj.cTau - tComp, obj.cTau + tComp, tQuant);
           
           J = zeros(wQuant, tQuant);
           for i = 1 : wQuant
              for j = 1 : tQuant
                 
                 omegan = wc(i); taun = tc(j); 
                 % Create the filtering structure
                 syscpoly = poly([-1*1/taun,...
                                  -zeta*omegan+omegan*sqrt(zeta^2-1),...
                                  -zeta*omegan-omegan*sqrt(zeta^2-1)]);
                 A0 = [[zeros(1,2);eye(2)], flipud(-syscpoly(2:end)')];

                 alphacand = poly(eig(A0));

                 obj.A = [zeros(1, obj.degree);[eye(obj.degree - 1),...
                                              flipud(-alphacand(2:end)')]];
                                                  
                 % Estimate the model
                 obj.estimate(dte);
                 obj.nonlinear_estimate(dte);
                 
                 % Simulate the data
                 dtb = obj.simulate(dtSim);

                 % Determine the cost funtion
                 J(i, j) = obj.fitnessStruc(dtSim, dtb, 'squared error');
              end
           end
           
           figure(4);
           surf(wc, tc, J);
           hold off;
           
           minimun = min(min(J));
           [wInd, tInd] = find(J == minimun);
           
           wBarry = wc(wInd); tBarry = tc(tInd);
           
           syscpoly = poly([-1*1/tBarry,...
                                  -zeta*wBarry+wBarry*sqrt(zeta^2-1),...
                                  -zeta*wBarry-wBarry*sqrt(zeta^2-1)]);
           A0 = [[zeros(1,2);eye(2)], flipud(-syscpoly(2:end)')];

           alphabarry = poly(eig(A0));

           obj.A = [zeros(1, obj.degree);[eye(obj.degree - 1),...
                                          flipud(-alphabarry(2:end)')]];
           
       end
       %% Tools functions
       function obj = check(obj,dt)
           % Check Size
           tSize = length(dt.t); ySize = length(dt.y);
           iSize = length(dt.init); fSize = length(dt.final);
           if (tSize < 1) || (ySize < 1) || (iSize < 1) || (fSize < 1)
              error('Input data is not valid!');
           end
           
           % Check type
       end
       
       function [xOut, tOut] = simuLinkSolve(~, model, x0, y, t)
           m = model;
           % reform the data 
           t = reshape(t,length(t),1);
           y = reshape(y,max(size(y)),min(size(y)));
           % prepare the data for simulation
           simin_y = [t, y];
           % determine the minimum step and the time simulated
           Tmax = .1*min(diff(t));
           if Tmax < 1e-05; Tmax = .001; end 
           tsim = t(end);
           
           % setup the problem options
           simout = [];
           options = simset('SrcWorkspace','current',...
                        'Solver','ode23','MaxStep',Tmax,'MinStep','auto');
           % simulate the system
           sim('vsSimModS',[],options);
           
           xOut = simout; tOut = tout;
       end
       
       function [c, J] = costFunc(~, x, h)
          time = h(:,3); y_d = h(:,2); y_a = h(:,1);
          % lsqnonlin error function
          c = x(1)-x(1).*exp(-time./x(2)) + y_d.*exp(-time./x(2)) - y_a;
    
          % Jacobian calculation
          J(:,1) = 1 - exp(-time./x(2));
          J(:,2) = -x(1).*time.*exp(-time./x(2))./(x(2).*x(2)) ...
                               + y_d.*time.*exp(-time./x(2))./(x(2).*x(2)); 
       end
       
       function J = fitnessStruc(~, dte, dts, method)
           % Initialize variables
           wSize = length(dte.y); yS = cell(wSize,1);
           % Check if the predicted data matches comparison in time
           for wii = 1 : wSize
            simSize = length(dts.dayOut{wii}); 
            realSize = length(dte.y{wii});
            % Check the data size for each one
            if simSize ~= realSize
                ind = zeros(realSize,1);
                for k = 1 : realSize
                   error = abs(dts.dayTime{wii} - dte.t{wii}(k)); 
                   [~, ind(k)] = min(error);
                end
                yS{wii} = dts.dayOut{wii}(ind);
            else
                yS{wii} = dts.dayOut{wii};
            end
           end
           % Determine the cost value
           yS = cell2mat(yS); yE = cell2mat(dte.y);
           if strcmp(method, 'squared error') 
            J = (yS - yE)'*(yS - yE);
           elseif strcmp(method, 'squared error normalized')
            J = ((yS - yE)'*(yS - yE))/length(yS);
           elseif strcmp(method, 'BFR')
            
           elseif strcmp(method, 'variance')
               
           else
            error('Fitness method not allowed!')
           end
           
       end
   end
    
end