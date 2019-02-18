function struc = struc_select(method)

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
end

end