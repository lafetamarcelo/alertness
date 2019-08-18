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
    dtp %The prediction data set
    % Determine the pre defined structure
    A = [zeros(1,4);[eye(3), zeros(3,1)]]; % observer matrix
    C = [zeros(1,3), 1]; % observer matrix
    Ao % The dynamic constant alertness matrix
    % Set the first estimation information
    preprocess = 'None'; % select the data preprocessing
    structure = 'Folkard'; % select the structure selector
    identification = 'N.N.L.S.'; % select the identification technique
    algorithm = 'MOLI'; % select the algorithm to be used
    showResults = 'None'; % determine the results to be shown
   end
    
   properties (Access = private)
    degree = 4;
    % Possible data preprocessing
    preprPoss = [{'None'}, {'Resample'}, {'Outliers'}];
    % Possible filtering procedures
    filterPoss = [{'None'},{'sGolay filtering'}];
    % Possible structure algorithms 
    strucPoss = [{'Folkard'},{'Let loose'},{'Barycenter'},...
        {'sGolay filtering'}, {'Grid search'}];
    % Possible algorithms
    algPoss = [{'MOLI'}, {'Drived search'}, {'Force Compute'}, {'Folkard'}];
    % Possible identification algorithms
    identPoss = [{'None'} ,{'Least squares'}, {'Extended Least Squares'},...
        {'Minimal squares'}, {'LPV'}, {'N.N.L.S.'}, {'Support Vector'},...
        {'Simulated Annealing'}, {'Particle Swarm'},{'Genetic Algorithm'},...
        {'Particle Swarm'}];
    % Possible results show
    showResPoss = [{'None'}, {'Sleepness'}, {'Validation'}];
    % Boundary info
    upperBound % Upper bound to detect the outliers 
    lowerBound % Lower bound for outliers
    % Model parameters
    L % The output regressor MOLI parameter matrix
    B % the input regressor MOLI parameter matrix
    % Regression model data
    theta % The regressor parameters
    Phi % The regressor information
    compensator % The regressor compensator
    Y % Considered output regressor
    % Plotting info
    figInd = 10 % initial figure plot
    % Extended Least Squares plot
    extThetaError
   end
   
   methods
       % Initializcde the alertness parameters 
       function obj = alertness(prepr, struc, alg, ident, show)
          if ismember(prepr, obj.preprPoss)
            obj.preprocess = prepr;
          else
            error('Pre processing technique is not possible.')
          end
          if ismember(struc, obj.strucPoss)
            obj.structure = struc;
          else
            error('Structure model is not possible.');
          end
          if ismember(alg, obj.algPoss)
            obj.algorithm = alg;
          else
            error('Algorithm does not exist!')
          end
          if ismember(ident, obj.identPoss)
            obj.identification = ident;
          else
            error('Identification technique not possible.');
          end  
          if ismember(show, obj.showResPoss)
            obj.showResults = show;
          else
            error('Show results is not possible.') 
          end
       end
       
       %% Flow control functions

       function obj = fit(obj, dte, dtv)
           % Verify the input and validation data-set
           %obj.check(dte); obj.check(dtv);
           obj.dtv = dtv;
           
           % Pre process the data set
           obj.dte = obj.process(dte);
           %obj.dtv = obj.process(dtv);
           
           % Filter the data set
           %obj.dte = obj.filter(dte);
           
           % Select the model structure
           obj.strucSelect(obj.dte);
           
           % Estimate the day parameters
           obj.estimate(obj.dte);
           
           % Validate the model
           obj.dts = obj.simulate(obj.dtv);
           
           % One step ahead prediction
           %obj.dtp = obj.predict(obj.dte);
           
           % Show the results
           obj.plotResults();
       end
       
       % Select the data pre processing
       function dts = process(obj, dte)
           if strcmp(obj.preprocess, 'Resample')
              % Resampling the data 
              dts = obj.resampleData(dte);
           elseif strcmp(obj.preprocess, 'Outliers')
              % Remove possible outliers
              dts = obj.outlierData(dte);
           elseif strcmp(obj.preprocess, 'None')
              dts = dte;
           else
               error('Pre process technique does not match!');
           end
           
        end
       % Select the model structure
       function obj = strucSelect(obj, dte)
          disp('Initialize the observer filter search...')
          % Determine the observer structure
          if strcmp('Barycenter', obj.structure)
            obj.barycenter(dte);
          elseif strcmp('Folkard', obj.structure)
            obj.trivialStruc();
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
               obj.nonlinear_estimate(dte);
           elseif strcmp('Drived search', obj.algorithm)
               obj.drivedSearch(dte);
               obj.nonlinear_estimate(dte);
           elseif strcmp('Folkard', obj.algorithm)
               obj.Ao = obj.A;
           elseif strcmp('Force Compute', obj.algorithm)
               obj.forceCompute(dte);
               obj.nonlinear_estimate(dte);
           else
              error('No algorithm with such name!');
           end
       end
       
       % Estimate the night parameters
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
       % Predict the alertness model
       function dt = predict(obj, dts)
          wSize = length(dts.y);
           
          for wii = 1 : wSize
           ySize = length(dts.y{wii});
           for k = 1 : ySize
            tOmega = obj.omega*dts.t{wii}(k);
            k1 = obj.M*cos(tOmega + obj.cPhase);
            k2 = -obj.M*sin(tOmega + obj.cPhase)*obj.omega;
            h0 = dts.y{wii}(k) - k1 - obj.dc;
            Beta = [(obj.dc*obj.omega^2)/obj.cTau;...
                 (k2/obj.cTau + h0*obj.omega^2 + obj.dc*obj.omega^2);...
                 (k1 + obj.dc + k2*obj.cTau)/obj.cTau;...
                 h0 + k1 + obj.dc];
            alertnessModel = ss(obj.Ao, Beta, obj.C, 0);
            
            ref = dts.t{wii}(k);
            if k ~= ySize
             Ts = dts.t{wii}(k+1) - ref;
            else
             Ts = dts.final(wii) - ref; 
            end
            
            if Ts ~= 0
             tSim = linspace(0, Ts, 100);
             [yPred, tPred] = impulse(alertnessModel, tSim);
             tPred = tPred + ref;
            
             if k == 1
              dt.dayOut{wii} = yPred;
              dt.dayTime{wii} = tPred;
             else
              dt.dayOut{wii} = [dt.dayOut{wii}; yPred];
              dt.dayTime{wii} = [dt.dayTime{wii}; tPred];
             end
            end
           end
          end
       end
       % Plot results
       function obj = plotResults(obj)
          if strcmp(obj.showResults, 'Sleepness')
            obj.sleepnessScatter(obj.dte, obj.dts);
          elseif strcmp(obj.showResults, 'Validation')
            obj.sleepnessScatter(obj.dte, obj.dts);
            if strcmp(obj.identification, 'Extended Least Squares')
                obj.extendedPlot();
            end
            obj.plotPrediction(obj.dte, obj.dtp);
          elseif strcmp(obj.showResults, 'None')
               
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
           elseif strcmp('Least squares', obj.identification)
               obj.theta = obj.Phi\obj.Y;
           elseif strcmp('Extended Least Squares', obj.identification)
               obj.extendedLS(dte);
           elseif strcmp('Minimal squares', obj.identification)
               obj.boundedLS(dte);
           elseif strcmp('Support Vector', obj.identification)
               obj.theta = obj.svmTry();
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
       
       function obj = extendedLS(obj, dte)
           
          wSize = length(dte.y);
          phiExt = obj.Phi; % Initialize the extended regressor
          
          iterations = 50;
          nx = 1;
          
          % Set up variables
          thetaExt = cell(iterations, 1);
          resid = cell(iterations, 1);
          % Initialize the first theta argument
          thetaExt{1} = phiExt\obj.Y;
          for it = 1 : iterations 
              % Determine the residuls
              resid{it} = obj.Y - phiExt*thetaExt{it};
              
              % Create the extended regressor
              upInd = 0;
              epsilon = zeros(length(obj.Y) ,nx);
              for wii = 1 : wSize
               lowInd = upInd + 1;
               upInd = upInd + length(dte.y{wii}) - 1;
               error = resid{it}(lowInd:upInd);
               for k = 1 : nx
                epsilon(lowInd:upInd, k) = [zeros(k,1);error(1:end-k)];
               end
              end
              
              phiExt = [obj.Phi, epsilon];
              
              thetaExt{it+1} = phiExt\obj.Y;
              
              if it == 1 
                 thetaExt{it} = [thetaExt{it}; zeros(nx, 1)]; 
              end
              
              obj.extThetaError(it) = (thetaExt{it+1}-thetaExt{it})'...
                                            *(thetaExt{it+1}-thetaExt{it});
          end
          
          obj.theta = thetaExt{iterations}(1:end-nx);
       end
       
       function obj = boundedLS(obj, dte)
        % Use the bounded biological parameters to estimate the bounds for 
        % MOLI paramertes
           wSize = length(dte.y); pSize = size(obj.Phi);
           
           % Initialize the vectors
           lbTheta = zeros(pSize(2), 1); 
           ubTheta = zeros(pSize(2), 1);
           
           % Determine the bounded values
           lbOmega = pi/24; ubOmega = pi/6;
           lbTau = 18; ubTau = 38;
           lbM = 1.5; ubM = 3.5;
           %lbPhase = -2*pi; ubPhase = 0.0;
           lbDC = 1.5; ubDC = 3.5;
           
           % Determine the lower and upper bound for L and DC info from B
           lbTheta(1) = .0;
           lbTheta(2) = -ubOmega^2/lbTau;
           lbTheta(3) = -ubOmega^2;
           lbTheta(4) = -1/lbTau;
           lbTheta(5) = lbDC*lbOmega^2/ubTau;
           
           ubTheta(1) = .0;
           ubTheta(2) = -lbOmega^2/ubTau;
           ubTheta(3) = -lbOmega^2;
           ubTheta(4) = -1/ubTau;
           ubTheta(5) = ubDC*ubOmega^2/lbTau;
           
           % Determine each B parameter lower and upper bounds
           for wii = 1 : wSize
              ind = 3*(wii - 1) + 5;
              % Determining the lower bounds
              q0 = -ubM*ubOmega;
              k0 = -ubM;
              Sw = dte.y{wii}(1) - ubM - ubDC;
              lbTheta(ind + 1) = q0/lbTau + Sw*lbOmega^2 + ubDC*lbOmega^2;
              lbTheta(ind + 2) = (k0 + lbDC)/ubTau + q0;
              lbTheta(ind + 3) = Sw + k0 + lbDC;
              
              % Determining the upper bounds
              q0 = ubM*ubOmega;
              k0 = ubM;
              Sw = dte.y{wii}(1) -lbM - lbDC;
              ubTheta(ind + 1) = q0/lbTau + Sw*ubOmega^2 + ubDC*ubOmega^2;
              ubTheta(ind + 2) = (k0 + ubDC)/lbTau + q0;
              ubTheta(ind + 3) = Sw + k0 + ubDC;
           end
           
           % Estimate the parameter vector
           obj.theta = lsqlin(obj.Phi, obj.Y, ...
                                [], [], [], [], lbTheta, ubTheta);
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
           
           obj.Ao(1,4) = .0;
           obj.Ao(2,4) = -obj.omega^2/obj.cTau;
           obj.Ao(3,4) = -obj.omega^2;
           obj.Ao(4,4) = -1/obj.cTau;
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
       
       % Computed solution
       function obj = forceCompute(obj, dte)
        if strcmp(obj.identification, 'Simulated Annealing')
         % Optimal simulation for Annealing algorithm
         obj.simAnnealing(dte);
        elseif strcmp(obj.identification, 'Particle Swarm')
         obj.partSwarm(dte);
        elseif strcmp(obj.identification, 'Genetic Algorithm')
         obj.genAlg(dte);
        elseif strcmp(obj.identification, 'None')
         % Default algorithm
         obj.nonlinLS(dte);
        else
         error('Estimation procedure not possible') 
        end
       end
        
       function obj = simAnnealing(obj, dte)
        % Simulate the annealing algorithm
           
           % Determine the simulation data
           dtSim.time = [[dte.t{1}(1); dte.init(2:end)'], dte.final'];
           dtSim.initial = dte.y{1}(1);
           dtSim.y = dte.y; dtSim.t = dte.t;
           
           % create some variables for the plot
           global iter__  bVal__ dcurves__ ncurves__ ind__;
           iter__ = 0; bVal__ = 0; ind__ = 1;
           
           % Determine the initial states 
           x0(1) = obj.omega + 0.8;  x0(2) = obj.cTau - 8;
           x0(3) = obj.M - 2;        x0(4) = obj.cPhase + 2;
           x0(5) = obj.dc + 2;
           
           method = 'squared error';
           optimFunc = @(x) obj.simCost(x, dtSim, dte, method);
           options = optimoptions(@simulannealbnd,...
               'TemperatureFcn', 'temperatureexp','Display', 'iter');
           lb = [pi/24, 18, 1.5, -2*pi, 1.5];
           ub = [pi/6, 38, 3.5, .0, 3.5];
           x = simulannealbnd(optimFunc, x0, lb, ub, options);
           
           % Determine the parameters
           obj.omega = x(1); obj.cTau = x(2);
           obj.M = x(3);     obj.cPhase = x(4);
           obj.dc = x(5);
           
           obj.Ao = [.0  .0 .0 .0;...
                      1. .0 .0 -obj.omega^2/obj.cTau;...
                     .0  1. .0 -obj.omega^2;...
                     .0  .0 1. -1/obj.cTau];
       end
       
       function obj = partSwarm(obj, dte)
        % Use particle swarm to estimate the best simulation model
           
           % Determine the simulation data
           dtSim.time = [[dte.t{1}(1); dte.init(2:end)'], dte.final'];
           dtSim.initial = dte.y{1}(1);
           dtSim.y = dte.y; dtSim.t = dte.t;
           
           method = 'squared error';
           optimFunc = @(x) obj.simCost(x, dtSim, dte, method);
           options = optimoptions('particleswarm', ...
               'HybridFcn', @patternsearch, ...
               'UseParallel', true, ...
               'Display', 'iter');
           lb = [pi/24, 18, 1.5, -2*pi, 1.5];
           ub = [pi/6, 38, 3.5, .0, 3.5];
           x = particleswarm(optimFunc, 5, lb, ub, options);
           
           % Determine the parameters
           obj.omega = x(1); obj.cTau = x(2);
           obj.M = x(3);     obj.cPhase = x(4);
           obj.dc = x(5);
           
           obj.Ao = [.0 .0 .0 .0;...
                      1. .0 .0 -obj.omega^2/obj.cTau;...
                      .0 1. .0 -obj.omega^2;...
                      .0 .0 1. -1/obj.cTau];
           
       end
       
       function obj = genAlg(obj, dte)
           % Determine the simulation data
           dtSim.time = [[dte.t{1}(1); dte.init(2:end)'], dte.final'];
           dtSim.initial = dte.y{1}(1);
           dtSim.y = dte.y; dtSim.t = dte.t;
           
           method = 'squared error';
           optimFunc = @(x) obj.simCost(x, dtSim, dte, method);
           options = optimoptions(@ga, 'Display', 'iter');
           lb = [pi/24, 18, 1.5, -2*pi, 1.5];
           ub = [pi/6, 38, 3.5, .0, 3.5];
           x = ga(optimFunc, 5, [], [], [], [], lb, ub, [], [], options);
           
           % Determine the parameters
           obj.omega = x(1); obj.cTau = x(2);
           obj.M = x(3);     obj.cPhase = x(4);
           obj.dc = x(5);
           
           obj.Ao = [.0 .0 .0 .0;...
                      1. .0 .0 -obj.omega^2/obj.cTau;...
                      .0 1. .0 -obj.omega^2;...
                      .0 .0 1. -1/obj.cTau];
       end
       
       function obj = nonlinLS(obj, dte)
        % Simulate the response with non linear least squares
           % Determine the simulation data
           dtSim.time = [[dte.t{1}(1); dte.init(2:end)'], dte.final'];
           dtSim.initial = dte.y{1}(1);
           dtSim.y = dte.y; dtSim.t = dte.t;
           
           % Determine the initial states 
           x0(1) = obj.omega; x0(2) = obj.cTau;
           x0(3) = obj.M;     x0(4) = obj.cPhase;
           
           method = 'error vector';
           optimFunc = @(x) obj.simCost(x, dtSim, dte, method);
           options = optimoptions(@lsqnonlin, 'Display', 'iter');
           lb = [pi/24, 18, 1.5, -2*pi];
           ub = [pi/6, 38, 3.5, .0];
           x = lsqnonlin(optimFunc, x0, lb, ub, options);
           
           % Determine the parameters
           obj.omega = x(1); obj.cTau = x(2);
           obj.M = x(3);     obj.cPhase = x(4);
       end
       
       function c = simCost(obj, x, dts, dte, method)
       % Determine the cost function for the force computed
            
            obj.omega = x(1); obj.cTau = x(2);
            obj.M = x(3); obj.cPhase = x(4);
            obj.dc = x(5);
            
            obj.Ao = [.0 .0 .0 .0;...
                      1. .0 .0 -obj.omega^2/obj.cTau;...
                      .0 1. .0 -obj.omega^2;...
                      .0 .0 1. -1/obj.cTau];
                
            % Simulate the system output
            dt = simulate(obj, dts);
            
            % Determine the cost value
            c = fitnessStruc(obj, dte, dt, method);
            
            % Plot the results
            plotAnnealingResults(dt, c);
            
       end
       %% Data preprocessing
       % Data resampler
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
       
       % Outlier detector
       function dtr = outlierData(obj, dte)
        wSize = length(dte.t);
        % Manipulating the time data for simulation
        dtSim.time = [[dte.t{1}(1); dte.init(2:end)'], dte.final'];
        dtSim.initial = dte.y{1}(1);
        dtSim.y = dte.y; dtSim.t = dte.t;
        
        % Select the Folkard initial model to search for outliers
        obj.trivialStruc(); obj.Ao = obj.A;
        
        % Determine the possible phase values 
        searchQuant = 50;
        phasePoss = linspace(-2*pi, 2*pi, searchQuant);
        
        Ts = cell(wSize, 1);
        Ys = cell(wSize, searchQuant);
        for pii = 1 : searchQuant
           obj.cPhase = phasePoss(pii);
            
           % Simulate the outputm
           dt = obj.simulate(dtSim);
           
           for wii = 1 : wSize 
            % Select the nearest to the time
            ind = zeros(length(dte.t{wii}),1);
            for k = 1 : length(ind)
                error = abs(dte.t{wii}(k) - dt.dayTime{wii});
                [~, ind(k)] = min(error);
            end
            % Determine the output data
            Ys{wii, pii} = dt.dayOut{wii}(ind);
            % Determine the time respective
            Ts{wii, 1} = dt.dayTime{wii}(ind);
           end
        end
        
        dtr = dte; % Initialize return data-set
        for wii = 1 : wSize
           auxOut = cell2mat(Ys(wii, :))';
           obj.upperBound{wii,1} = max(auxOut)';
           obj.lowerBound{wii,1} = min(auxOut)';
           
           % Eliminate the data whose is out of bounds
           boundedInd = boolean((dte.y{wii} <= obj.upperBound{wii})...
               .* (dte.y{wii} >= obj.lowerBound{wii})); 
           dtr.y{wii} = dte.y{wii}(boundedInd);
           dtr.t{wii} = dte.t{wii}(boundedInd);
        end
       end
       
       %% Filtering/Observer methods
        % Folkard structure
       function obj = trivialStruc(obj) 
          alpha = [0;...
                   -obj.omega^2/obj.cTau;...
                   -obj.omega^2;...
                   -1/obj.cTau];
          obj.A = [[zeros(1,3);eye(3)], alpha]; 
       end
       
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
       
           framelen = 7; order = 2;

           % Determining the sGolay filtered out 
           dtsG = obj.sGolay(dte, framelen, order);

           % Estimate the model parameters
           obj.estimate(dtsG);
           
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
            
            dtsG.t{i} = dte.t{i}(dte.t{i} <= obj.dte.final(i)); 
            dtsG.y{i} = cmplt{i}(1 : length(dtsG.t{i}));
        end 
       
       figure(1); hold on;
       for i  = 1 : wSize
        plot(dte.t{i}, cmplt{i}, 'k-', 'LineWidth', 1.2);
        plot(dte.t{i}, gSignal{i}, 'r--', 'LineWidth', 1.2);
       end
       hold off;
       
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
       
       %% Tool functions
       function obj = check(obj, dt)
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
            
           elseif strcmp(method, 'error vector')
            J = yS - yE;
           else
            error('Fitness method not allowed!')
           end
           
       end
       
       function beta = svmTry(obj)
          X = [obj.Phi(:,1:4), obj.Phi(:, 6:end)];
          %X = obj.Phi; 
          Ysv = obj.Y;
          svmModel = fitrsvm(X, Ysv, 'Standardize', true);%, 'OptimizeHyperparameters', 'auto');
              
          [~, betaSize] = size(obj.Phi);
          beta = zeros(betaSize, 1);
          phiRec = obj.Phi(svmModel.IsSupportVector, :);
          for aii = 1 : length(phiRec)
             beta = beta + svmModel.Alpha(aii).*phiRec(aii, :)';
          end
          beta = [beta(1:4); svmModel.Bias; beta(5:end)];
       end
       
       function J = distFolk(obj)
          omegaFolk = pi/12; cTauFolk = 1/0.0353; MFolk = 2.52; 
          cPhaseFolk = -16.835*pi/12; dcFolk = 2.4; 
          J = zeros(1,5);
          
          J(1) = abs(100*(omegaFolk - obj.omega)/omegaFolk);
          J(2) = abs(100*(cTauFolk - obj.cTau)/cTauFolk);
          J(3) = abs(100*(MFolk - obj.M)/MFolk);
          J(4) = abs(100*(cPhaseFolk - obj.cPhase)/cPhaseFolk);
          J(5) = abs(100*(dcFolk - obj.dc)/dcFolk);
       end
       
       %% Plotting functions
       function obj = sleepnessScatter(obj, dte, dts)
           
           wSize = length(dte.y);
           
           figure(obj.figInd + 1); hold on;
           set(gcf,'color','w');
           title('Without night estimate', 'Interpreter', 'latex');
           ylim([0 10]);
           for wii = 1 : wSize
              if wii ~= wSize
                  % Include rectangles
                  colors = [.8 .8 .8];
                  dot_1 = dts.nightTime{wii}(1);
                  dot_2 = dts.nightTime{wii}(end) - dts.nightTime{wii}(1);
                  rectangle('Position', [dot_1,0,dot_2,16], 'FaceColor', ...
                      colors, 'EdgeColor', colors);
              end
              
              % Determine the sleepness
              sleep = 10 - 9.*(dte.y{wii} - 1)./13;
              sleepSim = 10 - 9.*(dts.dayOut{wii} - 1)./13;
              
              % Scatter the output
              scatter(dte.t{wii}, sleep, 'ko', 'LineWidth', 1.4);
              plot(dts.dayTime{wii}, sleepSim, '-', 'LineWidth', 1.6, ...
                  'Color', [.8 .0 .0]);
           end
           legend([{'Measure'}, {'Model simulated'}], ...
               'Interpreter', 'latex');
           
           figure(obj.figInd + 2); hold on;
           set(gcf,'color','w');
           title('With night estimate', 'Interpreter', 'latex');
           ylim([0 10]);
           for wii = 1 : wSize
              if wii ~= wSize
                  % Include rectangles
                  colors = [.8 .8 .8];
                  dot_1 = dts.nightTime{wii}(1);
                  dot_2 = dts.nightTime{wii}(end) - dts.nightTime{wii}(1);
                  rectangle('Position', [dot_1,0,dot_2,16], 'FaceColor',...
                      colors, 'EdgeColor', colors);
              end
               
              % Determine the sleepness
              sleep = 10 - 9.*(dte.y{wii} - 1)./13;
              sleepSim = 10 - 9.*(dts.dayOut{wii} - 1)./13;
              
              % Scatter the output
              scatter(dte.t{wii}, sleep, 'ko', 'LineWidth', 1.4);
              plot(dts.dayTime{wii}, sleepSim, '-', 'LineWidth', 1.4, ...
                  'Color', [.8 .0 .0]);
              
              % Solve the night estimate
              if wii ~= wSize
                 sleepSimN = 10 - 9.*(dts.nightOut{wii} - 1)./13;
                 
                 plot(dts.nightTime{wii}, sleepSimN, '--', ...
                     'LineWidth', 1.4, 'Color', [.0 .0 .8]);
              end 
           end
           %legend([{'Measure'}, {'Model simulated'}], ...
           %    'Interpreter', 'latex'); hold off;
          
           obj.figInd = obj.figInd + 2;
       end
       
       function obj = extendedPlot(obj)
       % Plot the extended least squares parameters convergence
         xSize = length(obj.extThetaError);
         
         figure(obj.figInd + 1); hold on;
         set(gca, 'YScale', 'log');
         scatter(1:xSize, obj.extThetaError, 'ro', 'LineWidth', 1.5);
         plot(1:xSize, obj.extThetaError, 'r--', 'LineWidth', 1.5);
         
         obj.figInd = obj.figInd + 1;
       end
       
       function obj = plotPrediction(obj, dte, dts)
        wSize = length(dts.dayOut);
           
        figure(obj.figInd + 1); hold on;
        scatter(cell2mat(dte.t), cell2mat(dte.y), 'ko', 'LineWidth', 1.4);
        
        tReg = cell(wSize, 1);
        for wii = 1 : wSize
         plot(dts.dayTime{wii}, dts.dayOut{wii}, 'r-', 'LineWidth', 1.4);
         tReg{wii} = dte.t{wii}(2:end); 
        end
        
        scatter(cell2mat(tReg), obj.Phi*obj.theta, 'bx', 'LineWidth', 1.4);
        
        obj.figInd = obj.figInd + 1;
       end
   end
    
end