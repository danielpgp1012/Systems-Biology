function [dxdt] = ODEq1(x,t,k)
%ODEq1 returns time derivative of all components in system
%   
k1=k(1);
k_1=k(2);
k2=k(3);
k3=k(4);
k_3=k(5);

S=x(1);
E=x(2);
ES=x(3);
P=x(4);
ESI=x(5);
EI=x(6);

dEdt=-k1*E*S+k_1*ES+k2*ES-k3*E*I+k_3*EI;
dSdt=-k1*E*S+k_1*ES-k1*EI*S+k_1
end

