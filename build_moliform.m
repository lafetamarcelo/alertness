
% Symbolic nominal model
clc;
clear;
% syms('omega','k1','k2','tau','h0','real');
% As = [-1/tau,0,0; 0, 0,-omega^2; 0, 1, 0];
% Bs = [h0; 1; 0];
% Cs = [1, k1, k2];

syms('omega','k1','k2','tau','h0','offset','real');
As = [-1/tau,0,0, 0;...
    0, 0,-omega^2, 0;...
    0, 1, 0, 0;
    0 0 0 0];
Bs = [h0; 1; 0; offset];
Cs = [1, k1, k2, 1];

% syms('omega','k1','k2','ks','dtau','tauw','h0','p','real');
% As = [-1/(dtau*p + tauw),0,0; 0, 0,-omega^2; 0, 1, 0];
% Bs = [h0, ks; k2, 0; k1, 0];
% Cs = [1, 0, 1];

ny = size(Cs,1);
nx = size(As,1);
l = size(As,1);		% list of observability indices --- l=3: ?nica forma.

invP = vpa(zeros(nx));
ct = 1;
for i = 1:ny
	ni = l(i);
	ci = Cs(i,:);
	
	for ini = 1:ni
		invP(ct,:) = ci*As^(ini-1);
		ct = ct + 1;
	end
end


P = simplify(inv(invP));

Q = vpa(zeros(nx));
aux = 0; ct = 1;
for i = 1:ny
	ni = l(i);
	aux = aux + ni;
	for ini = 1:ni
		Q(:,ct) = (As^(ini-1))*P(:,aux);
		ct = ct + 1;
	end
end

Ac = simplify(Q\As*Q)   %simplify(Q+diff(inv(Q),'p')*Q)
Bc = simplify(Q\Bs)
Cc = simplify(Cs*Q)


% syms('omega','k1','k2','real');
% As = [0,1;omega^2,0];
% Bs = [0;1];
% Cs = [k1, k2];
% 
% O = [Cs; Cs*As];
% P = inv(O);
% Q = [P(:,2), As*P(:,2)];
% 
% Ac = simplify(Q\As*Q)
% Bc = simplify(Q\Bs)
% Cc = simplify(Cs*Q)