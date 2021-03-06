function [x_new,t] = a_to_b(x,p)
%Function that returns new state from old state
%   x: State as the number of molecules 
%
%generate two random numbers between 0 and 1
r1=rand(1);
r2=rand(1);
%calculate propensities
cf=p.kf;
cb=p.kb;
na=6.02e23;
% Na and Nb
Na=x(1)*p.V*na;
Nb=x(2)*p.V*na;

a_1=Na*cf; %a_1=cf*Na
a_2=Nb*cb; %a_2=cb*Nb
a_0=a_1+a_2; %a_0=sum(a_i)

if r1==0
    r1=eps;
end

t = -log(r1)/a_0;

%Step 3: Check if reaction goes f or b
if a_2<r2*a_0
    %forward reaction occurs
    Ext=t*p.kf*x(1);
    x_new=[x(1)-Ext,x(2)+Ext];
else
    %backward rxn occurs
    Ext=t*p.kb*x(2);
    x_new=[x(1)+Ext,x(2)-Ext];
    
end
