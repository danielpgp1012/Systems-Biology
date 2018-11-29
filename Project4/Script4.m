%% Question 1

%% b
p=struct('V',1.7e-12,'kf',2,'kb',1); %Volume in L, 
x_0=[1e-9,0]; %cc of A and B in M
No_timesteps=10000;
No_repetitions=10;
time=cell(No_repetitions,1);
X=cell(No_repetitions,1);

for i=1:No_repetitions
    time{i}=zeros(No_timesteps,1);
    X{i}=zeros(No_timesteps,2);
    X{i}(1,:)=x_0; %initial conditions
    for j=2:No_timesteps
        [X{i}(j,:),dt]=a_to_b(X{i}(j-1,:),p);      
        time{i}(j)=time{i}(j-1)+dt;
    end
    
end

%% c
for i=1:No_repetitions
stairs(time{i},X{i}*1e9)
hold on 
end
ylabel('Concentration [nM]')
xlabel('Time (s)')

%% d
time_plot=linspace(0,2); %plot until steady state is reached
X_A=zeros(length(time_plot),No_repetitions);
X_B=zeros(length(time_plot),No_repetitions);
for i=1:No_repetitions
   X_A(:,i)=transpose(spline(time{i},transpose(X{i}(:,1)),time_plot)); 
   X_B(:,i)=transpose(spline(time{i},transpose(X{i}(:,2)),time_plot));
   
end
X_A_avg=mean(X_A,2);
X_B_avg=mean(X_B,2);

X_A_std=transpose(std(transpose(X_A)));
X_B_std=transpose(std(transpose(X_B)));

X_A_plot=[X_A_avg-X_A_std,X_A_avg,X_A_avg+X_A_std];
X_B_plot=[X_B_avg-X_B_std,X_B_avg,X_B_avg+X_B_std];

figure 
plot1=plot (time_plot,X_A_plot(:,1)*1e9,'--',time_plot,X_A_plot(:,2)*1e9,'-',time_plot,X_A_plot(:,3)*1e9,'--');
hold on
plot2=plot (time_plot,X_B_plot(:,1)*1e9,'--',time_plot,X_B_plot(:,2)*1e9,'-',time_plot,X_B_plot(:,3)*1e9,'--');



ylabel('Concentration (nM)')
xlabel('time (s)')


%% e Mass Action SOlution
dXdt=@(t,x) [-p.kf*x(1)+p.kb*x(2);p.kf*x(1)-p.kb*x(2)];
tspan=[0,2];
[t,x]=ode45(dXdt,tspan,x_0);
gca;
plot3=plot(t,x*1e9);

title('Concentration Evolution with V=1.7pL')
legend([plot1(2) plot2(2) plot3(1) plot3(2)],'A','B','A_{MA}','B_{MA}')
hold off

%% f
%% b
p.V=1.7e-13; %Volume in L, 
x_0=[1e-9,0]; %cc of A and B in M
No_timesteps=10000;
No_repetitions=10;
time=cell(No_repetitions,1);
X=cell(No_repetitions,1);

for i=1:No_repetitions
    time{i}=zeros(No_timesteps,1);
    X{i}=zeros(No_timesteps,2);
    X{i}(1,:)=x_0; %initial conditions
    for j=2:No_timesteps
        [X{i}(j,:),dt]=a_to_b(X{i}(j-1,:),p);      
        time{i}(j)=time{i}(j-1)+dt;
    end
    
end

%% c
figure
for i=1:No_repetitions
stairs(time{i},X{i}*1e9)
hold on 
end
ylabel('Concentration [nM]')
xlabel('Time (s)')

%% c zoomed in
figure
for i=1:No_repetitions
stairs(time{i},X{i}*1e9)
hold on 
end
ylabel('Concentration [nM]')
xlabel('Time (s)')
xlim([0 2])

%% d
time_plot=linspace(0,2); %plot until steady state is reached
X_A=zeros(length(time_plot),No_repetitions);
X_B=zeros(length(time_plot),No_repetitions);
for i=1:No_repetitions
   X_A(:,i)=transpose(spline(time{i},transpose(X{i}(:,1)),time_plot)); 
   X_B(:,i)=transpose(spline(time{i},transpose(X{i}(:,2)),time_plot));
   
end
X_A_avg=mean(X_A,2);
X_B_avg=mean(X_B,2);

X_A_std=transpose(std(transpose(X_A)));
X_B_std=transpose(std(transpose(X_B)));

X_A_plot=[X_A_avg-X_A_std,X_A_avg,X_A_avg+X_A_std];
X_B_plot=[X_B_avg-X_B_std,X_B_avg,X_B_avg+X_B_std];

figure 
plot1=plot (time_plot,X_A_plot(:,1)*1e9,'--',time_plot,X_A_plot(:,2)*1e9,'-',time_plot,X_A_plot(:,3)*1e9,'--');
hold on
plot2=plot (time_plot,X_B_plot(:,1)*1e9,'--',time_plot,X_B_plot(:,2)*1e9,'-',time_plot,X_B_plot(:,3)*1e9,'--');



ylabel('Concentration (nM)')
xlabel('time (s)')

%% e Mass Action SOlution
dXdt=@(t,x) [-p.kf*x(1)+p.kb*x(2);p.kf*x(1)-p.kb*x(2)];
tspan=[0,2];
[t,x]=ode45(dXdt,tspan,x_0);
gca;
plot3=plot(t,x*1e9);
ylim([0 1])
title('Concentration Evolution with V=0.17pL')
legend([plot1(2) plot2(2) plot3(1) plot3(2)],'A','B','A_{MA}','B_{MA}')