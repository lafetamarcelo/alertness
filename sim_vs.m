function [ simout, tout] = sim_vs( MODEL, x0, y, t, version)
% function [ simout ] = varStepSim( MODEL, x0, u, t)
%
% Simulate the response of state-space linear model MODEL to the input 
% sequence u (Nxnu) at N irregular sampling instants defined in vector 
% t (Nx1). Variable x0 is the model initial condition.
%
% EXAMPLE:
% t = cumsum(1./logspace(-2,1,8))';
% u = rand(length(t),1);
% MODEL = ss(tf(1,[1 1]));
% simout = varStepSim(MODEL, 0, u, t); 

m = MODEL;

t = reshape(t,length(t),1);
y = reshape(y,max(size(y)),min(size(y)));

if length(y) == 1
    simin_y = [t,zeros(length(t),1)];
    %simin_y = [t,zeros(length(t),1),ones(length(t),1)];
else
    simin_y = [t,y];
end

Tmax = .1*min(diff(t));
tsim = t(end);
%tsim = 2*t(end);

simout = [];
options = simset('SrcWorkspace','current',...
    'Solver','ode23','MaxStep',Tmax,'MinStep','auto');

if strcmp(version, '15')
    sim('vsSimModS_R2015',[],options);
elseif strcmp(version, '16')
    sim('vsSimModS',[],options);
end
%% Just for debug
%plot(tout,simout(:,1),tout,simout(:,2));
%hold on;
%plot(t,u,'.b')
%legend('input','output')

end
