function obj = objfun(x,A)
%Objective Function Returns biomass over sum of fluxes
%   Detailed explanation goes here
obj=-x(17)/sum(x.^2);
end

