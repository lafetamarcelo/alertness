%% Create an alertness class for estimation and prediction:
% Inputs:   - estimation data dte:: y{cells}, t{cells}
%           - simulation data dts:: initial{array}, t{cells}
% Outputs:  - model :: physical_parameters, dts, structure

classdef alertness < handle
    
   properties % Define the public information
    physical_parameters %The alertness parameters information
    dte %The estimation data set
    dtv %The simulated data set
    A = [zeros(1,4);[eye(3), zeros(3,1)]]; % observer matrix
    C = [zeros(1,3), 1]; % observer matrix
    structure = 'trivial'; % select the structure selector
    identification = 'least squares'; % select the identification technique
    
   end
    
   properties (GetAccess = private) % Define the private information
    strucPoss = ['Folkard','Let loose','Barycenter','sGolay filtering'];
    identPoss = ['least squares', 'instrumental variable', 'extended', 'LPV'];
   end
   
   methods
       % Initialize the alertness parameters 
       function obj = alertness(struc, ident) % initialize the model
          if ismember(obj.strucPoss, struc)
            obj.structure = struc;
          else
            error('Structure model is not possible.');
          end
          if ismember(obj.identPoss, ident)
            obj.identification = ident;
          else
            error('Identification technique not possible.');
          end
       end
       
       %% Flow control functions
       function obj = fit(obj, dte, dtv)
           % Verify the input and validation data-set
           check(dte); check(dtv);
           obj.dte = dte; obj.dtv = dtv;
           
           % Select the model structure
           struc_select(obj.dte);
           
           % Estimate the model parameters
           estimate(dte);
           
           % Validate the model
           obj.dts = simulate(dtv);
       end
       
       % Select the model structure
       function obj = struc_select(obj, dte)
          disp('Initialize the observer filter search...')
          % Determine the observer structure
          if strcmp('Barycenter', obj.structure)
            barycenter(obj, dte);
          elseif strcmp('Folkard', obj.structure)
            trivial_struc(obj);
          elseif strcmp('sGolay filtering', obj.structure)
            sGolay_resample(obj, dte);
          elseif strcmp('Let loose', obj.structure)
            % Does nothing!!
          else
            error('Observer search technique does not match the possibilities'); 
          end
          disp('Done!')
       end
       
       function obj = estimate(obj, dte)
           % Create the Morse Observer
           observerMorse = ss(obj.A',obj.C',eye(4),0);
           % Determine the constants and initialize variables
           wSize = length(dte.y); Phi = cell(wSize, wSize+2);
           Y = cell(wSize, 1);
           for wii = 1 : wSize
               refTime = dte.t{wii} - dte.t{wii}(1);
               % Observe the input impulse response -- B[.] info
               impTime = 0:.0001:refTime(end);
               impStates = impulse(observerMorse, impTime);
               % Error analysis and select the simulation data
               impInd = zeros(length(dte.t{wii}),1);
               for k = 1 : lenght(impInd)
                   e = impTime - (dte.t{wii}(k) - dte.t{wii}(1));
                   [~, impInd(k)] = min(e);
               end
               
               % Observe the output filtered response -- L[.] info
               [outStates, outTime] = simuLinkSolve(observerMorse, 0, ...
                                                    dte.y{wii}, refTime);
               outInd = zeros(length(dte.t{wii}),1);
               for k = 1 : lenght(outInd)
                   e = outTime - (dte.t{wii}(k) - dte.t{wii}(1));
                   [~, outInd(k)] = min(e);
               end
               
               % Do not provide the first simulated data
               outInd = outInd(2:end); impInd = impInd(2:end);
               
               % Create the regressor matrix
               Phi{wii, 1} = outStates(outInd, 1:end); % L[.] info
               Phi{wii, 2} = ones(length(outInd), 1); % DC info
               
               for k = 1 : wSize % B[.] info
                    if wii == k
                        Phi{wii,k+2} = [ impStates(impInd,1), ...
                                        -impStates(impInd,2), ...
                                         impStates(impInd,3)];
                    else
                        Phi{wii,k+2} = zeros(length(outInd), 3);
                    end
               end
               % Estipulate the needed output data
               Y{wii} = dte.y(outInd);
           end
           
           % Build the regressor matrix
           obj.Phi = cell2mat(Phi);
           obj.Y = cell2mat(Y);
           
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
           
           % Retrieve the physical parameters
           parameters_retrieve(dte);
           
       end
       
       
       function obj = parameters_retrieve(obj, dte)
            
       end
       
       
       function dts = sim(obj, dts)
           
           
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
       
       function [xOut, tOut] = simuLinkSolve(model, x0, ys, ts)
           m = model;
           % reform the data 
           ts = reshape(ts,length(ts),1);
           ys = reshape(ys,max(size(ys)),min(size(ys)));
           % prepare the data for simulation
           simin_y = [ts, ys];
           % determine the minimum step and the time simulated
           Tmax = .1*min(diff(ts));
           if Tmax < 1e-05; Tmax = .001; end 
           tsim = ts(end);
           
           % setup the problem options
           simout = [];
           options = simset('SrcWorkspace','current',...
                        'Solver','ode23','MaxStep',Tmax,'MinStep','auto');
           % simulate the system
           sim('vsSimModS',[],options);
           
       end
   end
   
   
    
end