function f = acetate_constrain(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
flux = x(25);
f = (11.5*0.158*exp(-0.0859*flux))/(11.5+0.158*(exp(-0.0859*flux)-1))-flux;
end

