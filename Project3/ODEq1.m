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
ESI=x(4);
EI=x(5);
P=x(6);


dSdt=-k1*E*S+k_1*ES-k1*EI*S+k_1*ESI;
dEdt=-k1*E*S+k_1*ES+k2*ES-k3*E*I+k_3*EI;
dESdt=k1*E*S-(k_1+k2)*ES-k3*I*ES+k_3*ESI;
dESIdt=k3*I*E-(k_3+k_1)*ESI+S*EI*k1;
dEIdt=ESI*k_1-(k1*S+k_3)*EI+k3*E*I;
dPdt=k2*ES;

dxdt=[dSdt;dEdt;dESdt;dESIdt;dEIdt;dPdt];
% Mass Balance
E_tot=E+ES+EI+ESI;
Kd=E*S/ES;
Ki=ES/ESI;
Ki=E*I/EI;
end

