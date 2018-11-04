function [dCdt] = ODEq2(t,C,vmax,kcat,Km,ka,ki,v1)
%ODEq1 returns time derivative of all components in system
%   
%initialization of rates
v=zeros(10,1);

%initialization of concentrations
MAPKKK=C(1);
MAPKKK_P=C(2);
MAPKK=C(3);
MAPKK_P=C(4);
MAPKK_PP=C(5);
MAPK=C(6);
MAPK_P=C(7);
MAPK_PP=C(8);

if (nargin<8)
    v(1)=vmax(1)*MAPKKK*(1+ka*(MAPK_PP))/((Km(1)+MAPKKK)*(1+ki*MAPK_PP));
else
    v(1)=v1*(1+ka*(MAPK_PP))/(1+ki*MAPK_PP);
end
v(2) = vmax(2)*MAPKKK_P/(Km(2)+MAPKK_P);
v(3) = kcat(3)*MAPKKK_P*MAPKK/(Km(3)+MAPKK);
v(4) = kcat(4)*MAPKKK_P*MAPKK_P/(Km(4)+MAPKK_P);
v(5) = vmax(5)*MAPKK_PP/(Km(5)+MAPKK_PP);
v(6) = vmax(6)*MAPKK_P/(Km(6)+MAPKK_P);
v(7) = kcat(7)*MAPKK_PP*MAPK/(Km(7)+MAPK);
v(8) = kcat(8)*MAPKK_PP*MAPK_P/(Km(8)+MAPK_P);
v(9) = vmax(9)*MAPK_PP/(Km(9)+MAPK_PP);
v(10)= vmax(10)*MAPK_P/(Km(10)+MAPK_P);

% Reaction Coefficients
%  V1  v2  v3  v4  v5  v6  v7  v8  v9  v10
A=[-1,  1,  0,  0,  0,  0,  0,  0,  0,  0;...
    1, -1,  0,  0,  0,  0,  0,  0,  0,  0;...
    0,  0, -1,  0,  0,  1,  0,  0,  0,  0;...
    0,  0,  1, -1,  1, -1,  0,  0,  0,  0;...
    0,  0,  0,  1, -1,  0,  0,  0,  0,  0;...
    0,  0,  0,  0,  0,  0, -1,  0,  0,  1;...
    0,  0,  0,  0,  0,  0,  1, -1,  1, -1;...
    0,  0,  0,  0,  0,  0,  0,  1, -1,  0];


dCdt=A*v;

end

