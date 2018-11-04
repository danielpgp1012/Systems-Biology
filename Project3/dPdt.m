function [dxdt] = dPdt(t,S,I,Kd,Ki,vmax)
%ODEq1 returns time derivative of all components in system
%   
dxdt=vmax*S./((S+Kd)*(1+I/Ki));
end

