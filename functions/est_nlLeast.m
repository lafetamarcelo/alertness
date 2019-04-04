function par = est_nlLeast(dte, method)
%% Estimate the non linear least squares regressor

         %M    %cphase  %omega %dc  %tau
init_ = [2.52, 2.4071, pi/12, 2.4, 1/0.0353];
lb = [.1, -2*pi, pi/24, .5, 24];
ub = [10, 2*pi, pi/6, 5.5, 38];

if ismember('nonLinear LS', method)

    options = optimoptions(@lsqnonlin, ...
                            'Algorithm','trust-region-reflective',...
                              'MaxIter', 200, ... 
                                'SpecifyObjectiveGradient', true,...
                                   'Display','iter','Jacobian','on');
                          %trust-region-reflective
    x = lsqnonlin(@cost_func, init_, [], [], options, dte);
    %x = lsqcurvefit(@cost_func, init_, [], [], lb, ub, dte);

    M = x(1); cphase = x(2); omega = x(3);
    dc = x(4); tau = x(5);
    
elseif ismember('genetic algorithm', method)

    Fitness = @(x) cost_ga(x, dte);
    options = optimoptions(@ga, 'UseVectorized', false, 'Display','iter');
    [x, ~] = ga(Fitness, 5, [], [], [], [], lb, ub, [], options);

    M = x(1); cphase = x(2); omega = x(3);
    dc = x(4); tau = x(5);
    
elseif ismember('', method)
    
    
    
    
end


y_sim = cell(length(dte.t),1);
t_sim = cell(length(dte.t),1);
figure(2); hold on;
for k = 1 : length(dte.t)
    
    s_t = linspace(dte.t{k}(1),dte.t{k}(end),100);
    e_t = s_t - s_t(1);
    
    c_i = M*cos(omega*s_t + cphase);
    
    ho = dte.y{k}(1) - c_i(1);
    
    s_i = (ho-dc).*exp(-e_t/tau)+dc; 
    
    y_sim{k} = s_i + c_i;
    t_sim{k} = s_t;
    subplot(3,1,1); hold on;
    plot(t_sim{k}, y_sim{k}, 'k', 'LineWidth', 1.2);
    scatter(dte.t{k}, dte.y{k}, 'r', 'LineWidth', 1.4);
    subplot(3,1,2); hold on;
    plot(t_sim{k}, c_i, 'k', 'LineWidth', 1.2);
    subplot(3,1,3); hold on;
    plot(t_sim{k}, s_i, 'k', 'LineWidth', 1.2);
end

par.M = M; par.cphase = cphase; par.omega = omega;
par.dc = dc; par.tau = tau;


end

function [c, J] = cost_func(x, dte)
    M = x(1); cphase = x(2); omega = x(3);
    dc = x(4); tau = x(5);
    
    size_ = length(dte.t);
    
    J_ = cell(size_,5);
    y_hat = cell(size_,1);
    for i = 1 : size_
       c_i = M*cos(omega*dte.t{i} + cphase);
       
       Sa = dte.y{i}(1) - c_i(1);
       t_ = dte.t{i} - dte.t{i}(1);
       s_i = (Sa-dc).*exp(-t_/tau)+dc; 
       
       y_hat{i} = s_i + c_i;
       
       figure(1); hold on;
       plot(dte.t{i}, y_hat{i})
       
       J_{i,1} = cos(omega*dte.t{i} + cphase);
       J_{i,2} = -M*sin(omega*dte.t{i} + cphase);
       J_{i,3} = -M.*dte.t{i}.*sin(omega*dte.t{i} + cphase);
       J_{i,4} = 1 - exp(-t_/tau);
       J_{i,5} = (Sa - dc).*t_.*exp(-t_/tau)./tau^2;
    end
    
    J = cell2mat(J_);
    c = cell2mat(dte.y) - cell2mat(y_hat);
    %c = c'*c;
end

function c = cost_ga(x, dte)
    M = x(1); cphase = x(2); omega = x(3);
    dc = x(4); tau = x(5);
    
    size_ = length(dte.t);
    
    y_hat = cell(size_,1);
    for i = 1 : size_
       c_i = M*cos(omega*dte.t{i} + cphase);
       
       Sa = dte.y{i}(1) - c_i(1);
       t_ = dte.t{i} - dte.t{i}(1);
       s_i = (Sa-dc).*exp(-t_/tau)+dc; 
       
       y_hat{i} = s_i + c_i;
       
       figure(1); hold on;
       plot(dte.t{i}, y_hat{i})
    end
    
    c = cell2mat(dte.y) - cell2mat(y_hat);
    c = c'*c;
end