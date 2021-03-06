%% Question 1
%k1=k(1);
%k_1=k(2);
%k2=k(3);
%k3=k(4);
%k_3=k(5);
Kd=0.25;
Ki=1;
vmax=50;
legendstr=cell(4,1);
s=linspace(0,2);
t=0;
ODE=@(s,I) dPdt(t,s,I,Kd,Ki,vmax);
for I=0:3
v=ODE(s,I);
plot(s,v)
hold on
legendstr{I+1}=(strcat("[ I ]= ",num2str(I), " mM"));
end
legend(legendstr{:})
xlabel('S(mM)')
ylabel('v(mmol/min)')
title('Production rate as a Function of Substrate')
hold off

%% Question 2 b
%Constant Initialization
vmax = [2.5; 0.25; 0; 0; 0.75; 0.75; 0; 0; 0.5; 0.5]; %nMs-1

kcat = [0; 0; 0.025; 0.025; 0; 0; 0.025; 0.025]; % s-1

Km = [10; 8; 15; 15; 15; 15; 15; 15; 15; 15]; %in nM

ka=0;
ki=0;

C0 = struct('Names',{{'MAPKKK','MAPKKK-P','MAPKK','MAPKK-P','MAPKK-PP','MAPK',...
    'MAPK-P','MAPK-PP'}},'ccs',[100;0;300;0;0;300;0;0]);

C  = struct('Names',{{'MAPKKK','MAPKKK-P','MAPKK','MAPKK-P','MAPKK-PP','MAPK',...
    'MAPK-P','MAPK-PP'}},'ccs',[]);
%set time span
t_span=[0, 500];

%Solve System
[t,C.ccs] = ode15s(@(t,C) ODEq2(t,C,vmax,kcat,Km,ka,ki),t_span,C0.ccs);
%% Figure
figure
plot(t,C.ccs)
legend(C.Names)
xlabel('time (s)')
ylabel('Concentration (nM)')
ylim([-15 315])
%% c
v1 = 0:0.01:1;
MAPK_PPss = zeros(length(v1),1);
t_span=[0, 900]; %Increase time span in case SS is extended for lower v1s
for i=1:length(v1)
vmax(1) = v1(i);
[t,y] = ode15s(@(t,C) ODEq2(t,C,vmax,kcat,Km,ka,ki),t_span,C0.ccs);  

MAPK_PPss(i)=y(end,8); %take steady state value
end
%% Figure Q2d
figure
plot (v1,MAPK_PPss)
title('SwitchLike Dose Response')
xlabel('V_{max,1}')
ylabel('[MAPK-PP]')
ylim([-15 315])

%% Q2d

vmax(1)=2.5;
ka=0;
ki=0.1;

t_span=([0, 10000]);

[t,y] = ode15s(@(t,C) ODEq2(t,C,vmax,kcat,Km,ka,ki),t_span,C0.ccs);  
figure 
plot(t,y(:,end))
title ('Oscillatory Response')
xlabel('Time (s)')
ylabel('[MAPK-PP]')
%legend(C.Names);

%%
vs=[1,0.5,0.2,0.1];
legstr=cell(length(vs),1);
figure
for i=1:length(vs)
   vmax(1)=vs(i);
   [t,y] = ode15s(@(t,C) ODEq2(t,C,vmax,kcat,Km,ka,ki),t_span,C0.ccs);
   plot(t,y(:,end))
   hold on
   legstr{i}=strcat('V_{max,1}=',num2str(vs(i)));
end
legend(legstr{:})
xlabel('time(s)')
ylabel('[MAPK-PP]')
title('Oscillatory Response at Different V_{1,max}')
ylim([-15 300])
