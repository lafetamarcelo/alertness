function [c,J] = cost_function(x,time,y_d,y_a)
    
    %lsqnonlin
    c = x(1)-x(1).*exp(-time./x(2)) + y_d.*exp(-time./x(2)) - y_a;
    
    %Jacobian
    J(:,1) = 1 - exp(-time./x(2));
    J(:,2) = -x(1).*time.*exp(-time./x(2))./(x(2).*x(2)) ...
        + y_d.*time.*exp(-time./x(2))./(x(2).*x(2));
    
    %c = sum(c.*c);
    
    %figure(5)
    %subplot(2,1,1); hold on;
    %scatter(x(1),c,'r','LineWidth',1.5);
    %subplot(2,1,2); hold on;
    %scatter(x(2),c,'b','LineWidth',1.5);
end