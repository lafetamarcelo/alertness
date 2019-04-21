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
    dte %The estimation data set
    dtv %The validation data set
    dts %The simulated data set
    A = [zeros(1,4);[eye(3), zeros(3,1)]]; % observer matrix
    C = [zeros(1,3), 1]; % observer matrix
    Ao % The dynamic constant alertness matrix
    structure = 'trivial'; % select the structure selector
    identification = 'least squares'; % select the identification technique
   end
    
   properties (GetAccess = private) % Define the private information
    % Possible structure search setting
    strucPoss = ['Folkard','Let loose','Barycenter','sGolay filtering'];
    degree = 4;
    % Possible identification algorithms
    identPoss = ['least squares', 'instrumental variable',...
        'extended', 'LPV', 'N.N.L.S'];
    % Model parameters
    L % The output regressor MOLI matrix parameters
    B
    % Regression model data
    theta % The regressor parameters
    Phi % The regressor information
    Y % Considered output regressor
   end
   
   methods
       % Initialize the alertness parameters 
       function obj = alertness(struc, ident) % initialize the model
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
       end
       
       %% Flow control functions
       function obj = fit(obj, dte, dtv)
           % Verify the input and validation data-set
           %obj.check(dte); obj.check(dtv);
           obj.dte = dte; obj.dtv = dtv;
           
           % Select the model structure
           obj.struc_select(obj.dte);
           
           % Estimate the model parameters
           obj.estimate(dte);
           
           % NonLinear estimate
           obj.nonlinear_estimate(dte);
           
           % Validate the model
           obj.dts = obj.simulate(dtv);
       end
       
       % Select the model structure
       function obj = struc_select(obj, dte)
          disp('Initialize the observer filter search...')
          % Determine the observer structure
          if strcmp('Barycenter', obj.structure)
            obj.barycenter(dte);
          elseif strcmp('Folkard', obj.structure)
            obj.trivial_struc();
          elseif strcmp('sGolay filtering', obj.structure)
            obj.sGolay_resample(obj, dte);
          elseif strcmp('Let loose', obj.structure)
            % Does nothing!!
          else
            error('Observer search technique does not match the possibilities'); 
          end
          disp('Done!')
       end
       % Estimate the regressor model parameters
       function obj = estimate(obj, dte)
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
       %Estimate the night parameters
       function obj = nonlinear_estimate(obj, dte)
           wSize = length(dte.y); h = zeros(wSize-1, 3);
           
           for wii = 1 : wSize - 1
            % Simulate the system backwards
            if dte.t{wii+1}(1) ~= dte.init(wii+1)
                alertnessBack = ss(-obj.Ao, obj.B{wii+1}', obj.C, 0);
                finalTime = dte.t{wii+1}(1) - dte.init(wii+1);
                [initA, ~] = impulse(alertnessBack, finalTime);
            else
                initA = flip(dte.y{wii+1});
            end
            % Simulate the system forward
            if dte.t{wii}(end) ~= dte.final(wii)
                alertnessForw = ss(obj.Ao, obj.B{wii}', obj.C, 0);
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
       
       % Determine/Retrieve the physical parameters
       function obj = parameters_retrieve(obj, dte)
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
       
       %% Filtering/Observer methods
       function obj = barycenter(obj, dte)
          
          obj.degree = 4; 
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
              
              % Manipulating the time data for simulation
              dtSim.time = [[dte.t{1}(1); dte.init(2:end)'], dte.final'];
              dtSim.initial = dte.y{1}(1);
              dtSim.y = dte.y; dtSim.t = dte.t;
              
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
       
       function obj = trivial_struc(obj) 
          alpha = [0;...
                   -obj.omega^2/obj.cTau;...
                   -obj.omega^2;...
                   -1/obj.cTau];
          obj.A = [[zeros(1,3);eye(3)], alpha]; 
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
       
       function [xOut, tOut] = simuLinkSolve(obj, model, x0, y, t)
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
       
       function [c, J] = costFunc(obj, x, h)
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