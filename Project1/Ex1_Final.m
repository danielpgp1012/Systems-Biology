close all
%% Load model
load('ecoli_core_model.mat')
%% Load the COBRA Toolbox
addpath(genpath(fullfile('..','ext')))

%%  Load dependencies
addpath(genpath(fullfile('..','matTFA')))
addpath(genpath(fullfile('..','thermoDatabases')))
addpath(genpath(fullfile('..','models')))
addpath(genpath(fullfile('..','plotting')))
 
%% Part 1
%1.1


model = changeObjective(model,model.rxns{13}); %Selct biomass

%1.2
model = changeRxnBounds(model,'EX_glc(e)',-10.0,'l');

%1.3
FBAsolution = optimizeCbModel(model,'max');

%1.4 & 1.5
%We will try to find the exchange reactions where the reactant has at least
%1 C

[selExc,selUpt]=findExcRxns(model,0); %outputs exchange reactions

uptakes=model.rxns(selExc); %selects substrates from exchange reactions

substratesModel=extractSubNetwork(model,uptakes); %create submodel only
% with given substrates of substrates

cReactions=findCarbonRxns(substratesModel,1); %find substrates that contain
% AT LEAST 1 carbon

biomass.aerobic=zeros(length(cReactions),1); %initialize aerobic biomass vector

biomass.anaerobic=biomass.aerobic; %"" anaerobic ""

model=changeRxnBounds(model,cReactions,0,'l'); %set all carbon substrate uptakes to 0

modelaerobic=model; %initialize aerobic model

modelanaerobic=model; %initialize anaerobic model

modelaerobic=changeRxnBounds(modelaerobic,'EX_o2(e)',-20,'b'); %set oxygen lower boundary uptake to -20mmol/DW/h
modelanaerobic=changeRxnBounds(modelanaerobic,'EX_o2(e)',0,'b'); %Set anaerobic conditions

for i=1:length(cReactions) %Change substrate each iteration and set its lower boundary to 10. 
modelaerobic=changeRxnBounds(modelaerobic,cReactions{i},-10.0,'l');
modelanaerobic=changeRxnBounds(modelanaerobic,cReactions{i},-10.0,'l');
%Optimize models 
FBAsolutionaerobic=optimizeCbModel(modelaerobic,'max');
FBAsolutionanaerobic=optimizeCbModel(modelanaerobic,'max');

%Store biomass in vectors
biomass.aerobic(i)=FBAsolutionaerobic.obj;
biomass.anaerobic(i)=FBAsolutionanaerobic.obj;

%Store glucose model
if (strcmp(cReactions{i},'EX_glc(e)'))
    modelo2save=modelaerobic;
    modelasave=modelanaerobic;
    FBAsolutionaerobicsave=FBAsolutionaerobic;
    FBAsolutionanaerobicsave=FBAsolutionanaerobic;
end

%Set all carbon reactions to 0
modelaerobic=changeRxnBounds(modelaerobic,cReactions,0.0,'l');
modelanaerobic=changeRxnBounds(modelanaerobic,cReactions,0.0,'l');
end

%% Display data
%bar graph: top->aerobic biomass, bot->anaerobic biomass
fig1=figure;
c=categorical(cReactions);
h1=subplot(2,1,1);
bar(c,biomass.aerobic)
title('Aerobic Max. Biomass')
ylabel('Max Biomass')
hold on
h2=subplot(2,1,2);
bar(c,biomass.anaerobic)
title('Anaerobic Max. Biomass')
ylabel('Max Biomass')
ylim([0 max(biomass.anaerobic)+0.05])
hold off

% Tables
T1=table(cReactions,biomass.aerobic,biomass.anaerobic); %yield is normalized to consumption
T1.Properties.VariableNames = {'Substrates' 'Aerobic_Yield' 'Anaerobic_Yield'}

%csv table to observe in ESCHER
T_ae=table(modelo2save.rxns,FBAsolutionaerobicsave.full);
writetable(T_ae, 'aerobic_flux_data.csv');
T_an=table(modelasave.rxns,FBAsolutionanaerobicsave.full);
writetable(T_an, 'anerobic_flux_data.csv');
%% Part 2
%2.1
model=changeRxnBounds(model,cReactions,1000,'u');
model=changeRxnBounds(model,'EX_o2(e)',1000,'u');
model=changeRxnBounds(model,cReactions,0.0,'l'); %Set all carbon reactions to 0
model=changeRxnBounds(model,'EX_glc(e)',-10,'l'); %Set glucose uptake to -10mmol/gDW/h
model=changeRxnBounds(model,'EX_o2(e)',-20,'l'); %Set oxygen uptake to -20mmol/gDW/h
cutoff=0.1; %Specify given cutoff of 10% of max biomass with given conditions of glucose
Tol=FBAsolutionaerobicsave.obj*cutoff; %Tol is set to be the minimum threshold
[grRatio, grRateKO, grRateWT, hasEffect] = singleGeneDeletion(model); %Compute effect of genes in biomass production

essentialGenes=model.genes(find(grRateKO<Tol)); %Find genes for which the biomass 'grRateKO' is lower than threshold 
%being deleted

%2.2

% Vary Substrates 
biomass.genesaerobic=zeros(length(cReactions),1); %vector containing number of essential genes for each substrate
biomass.genesanaerobic=biomass.genesaerobic; %vector containing number of essential genes for each substrate
model=changeRxnBounds(model,cReactions,0.0,'l'); %set all substrate uptake to 0

threshold.aerobic=cutoff*biomass.aerobic; %vector with cuttoff for each substrate
threshold.anaerobic=cutoff*biomass.anaerobic;

for i=1:length(cReactions)
model=changeRxnBounds(model,cReactions{i},-10,'l'); %Update substrate with every iteration. Update its max uptake to -10 mmol/gDW/h
[grRatio, grRateKO, grRateWT, hasEffect] = singleGeneDeletion(model); %Do gene deletion calculation
essentialGenes=model.genes(find(grRateKO<threshold.aerobic(i))); %Find those genes that after deletion decrease the biomass production to less than 10% of the max
biomass.genesaerobic(i)=length(essentialGenes); %Store total number of genes in vector

%for anaerobic condition we only need to check glucose and fructose
if ((strcmp(cReactions{i},'EX_glc(e)')) || (strcmp(cReactions{i},'EX_fru(e)'))) 
model=changeRxnBounds(model,'EX_o2(e)',0,'l');
[grRatio, grRateKO, grRateWT, hasEffect] = singleGeneDeletion(model); %hasEffect: Boolean
essentialGenes=model.genes(find(grRateKO<threshold.anaerobic(i)));
biomass.genesanaerobic(i)=length(essentialGenes);
model=changeRxnBounds(model,'EX_o2(e)',-20,'l'); %set oxygen consumption rate to -20 again
end

model=changeRxnBounds(model,cReactions,0,'l'); %set all substrate consumption rates to 0

end
%2.3 Display data
%Bar Graph 
fig2=figure;
f1=subplot(2,1,1);
bar(c,biomass.genesaerobic)
title('Aerobic and Anaerobic Essential Genes for Each Substrate')
hold on
f2=subplot(2,1,2);
bar(c,biomass.genesanaerobic)
title('Anaerobic Max. Biomass')
hold off
%table
t2=table(cReactions,biomass.genesaerobic,biomass.genesanaerobic);
t2.Properties.VariableNames={'Substrate' 'Aerobic_Genes' 'Anaerobic_Genes'}

%% Part 3: Flux Variabilty
%3.1
%Initialize cells
minFluxx.aerobic=[];
minFluxx.anaerobic=[];
maxFluxx.aerobic=[];
maxFluxx.anaerobic=[];
%Change solver
solverOK=changeCobraSolver('pdco','QP'); %QP type of problem
%Set glucose and oxygen uptakes
model=changeRxnBounds(model,'EX_glc(e)',-10,'l');
model=changeRxnBounds(model,'EX_o2(e)',-20,'l');
%List relevant Reactions
Relevant_rxns={'EX_ac(e)';'EX_co2(e)';'EX_lac-D(e)';'ACONTb';'MDH';'EX_o2(e)';'EX_glc(e)'};
%Calculate flux variability using Cobra Function
[minFlux,maxFlux,Vmin,Vmax] = fluxVariability(model,80,'max',Relevant_rxns);
%Store results in vectors
minFluxx.aerobic=minFlux;
maxFluxx.aerobic=maxFlux;

%3.2
%Set anaerobic conditions
model=changeRxnBounds(model,'EX_o2(e)',0,'l');

[minFlux,maxFlux,Vmin,Vmax] = fluxVariability(model,80,'max',Relevant_rxns);
minFluxx.anaerobic=minFlux;
maxFluxx.anaerobic=maxFlux;

%3 Table Display
T_3=table(Relevant_rxns,minFluxx.aerobic,maxFluxx.aerobic,minFluxx.anaerobic,maxFluxx.anaerobic);
T_3.Properties.VariableNames = {'Substrates' 'MinFlux_Aerobic' 'MaxFlux_Aerobic' 'MinFlux_Anaerobic' 'MaxFlux_Anaerobic'}


%% PART 4 - Glucose
%4.1
UpGlucose=linspace(0,30,40); %Create vector from 0 to 30 of glucose
model=changeRxnBounds(model,'EX_o2(e)',-20,'b'); %change oxygen max uptake to -20 mmol/gDW/h

Biomass_vector=zeros(length(UpGlucose),1); %Set initialize biomass vector   

%Cycle through all glucose max uptake rates and observe effect on biomass
%production
for i=1:length(UpGlucose)
model = changeRxnBounds(model,'EX_glc(e)',-UpGlucose(i),'b');
sol=optimizeCbModel(model,'max');
Biomass_vector(i)=sol.obj;
end

fig3=figure;
plot(-UpGlucose,Biomass_vector)
xlabel('Glucose Uptake [mmol/gDW/h]')
ylabel('Biomass [mmol/gDW/h]')
title('Biomass Production as function of Glucose Uptake Rate')
%%
%4.2
UpOxygen=linspace(0,30,40); %oxygen vector 
Objective=zeros(length(UpGlucose),length(UpOxygen)); %initialize biomass matrix 

for i=1:length(UpOxygen) %nested for loop looking at biomass output at each glucose and oxygen combination
    model = changeRxnBounds(model,'EX_o2(e)',-UpOxygen(i),'b');
    for j=1:length(UpGlucose)
        model = changeRxnBounds(model,'EX_glc(e)',-UpGlucose(j),'b');
        sol=optimizeCbModel(model,'max');
        Objective(i,j)=sol.obj;
    end
end

fig4=figure;
surfl(UpGlucose,UpOxygen,Objective); %3D plot
xlabel('Glucose Uptake [mmol/gDW/h]')
zlabel('Biomass [mmol/gDW/h]')
ylabel('Oxygen Uptake [mmol/gDW/h]');
%%
%4.3
%Compute shadowprices with built-in function: For figures to work fine,
%need to edit 'phenotypePhasePlane' and remove the numbering after figure
%()
[growthRates, shadowPrices1, shadowPrices2] = phenotypePhasePlane(model, 'EX_glc(e)', 'EX_o2(e)',50,30,30);
%first graph corresponds to glucose

%%CSV Table Generation to observe in ESCHER

Glucose.values=[14.29, 8.163, 5.306, 0.8163]; %Chose datapoints within specific shadowprices
Oxygen.values=[9.796,15.1,16.33,14.29]; %Corresponding data points to glucose vector
Glucose.ShadowPrices=[2, 5, 13, 0]; %Shadow Prices of glucose Values *1e-2

%for loop
for i=1:length(Oxygen.values)
    model=changeRxnBounds(model,{'EX_o2(e)','EX_glc(e)'},-[Oxygen.values(i),Glucose.values(i)],'b');
    Solution_sp=optimizeCbModel(model,'max');
    names=strcat("ShadowPrice_Glucose",num2str(Glucose.ShadowPrices(i)),".csv"); 
    Table_sp=table(model.rxns,Solution_sp.full);
    writetable(Table_sp,names)
end

%% PART 4 - Succinate
%4.1
model=changeRxnBounds(model,cReactions,0,'b'); %set all carbon substrate uptakes to 0
model=changeRxnBounds(model,cReactions,1000,'u'); %set upper boundary to 1000
UpSuccinate=linspace(0,30,40); %Create vector from 0 to 30 of glucose
model=changeRxnBounds(model,'EX_o2(e)',-20,'b'); %change oxygen max uptake to -20 mmol/gDW/h

Biomass_vector_s=zeros(length(UpSuccinate),1); %Set initialize biomass vector   
sol_s={};
%Cycle through all glucose max uptake rates and observe effect on biomass
%production
for i=1:length(UpSuccinate)
    model = changeRxnBounds(model,'EX_succ(e)',-UpSuccinate(i),'b');
    sol_s=optimizeCbModel(model,'max');
    Biomass_vector_s(i)=sol_s.obj;
end

fig5=figure
plot(-UpSuccinate,Biomass_vector_s)
xlabel('Succinate Uptake [mmol/gDW/h]')
ylabel('Biomass [mmol/gDW/h]')
title('Biomass Production as function of Succinate Uptake Rate')
%%
%4.2
model=changeRxnBounds(model,cReactions,0,'b'); %set all carbon substrate uptakes to 0
model=changeRxnBounds(model,cReactions,1000,'u'); %set upper boundary to 1000
UpOxygen=linspace(0,30,40); %oxygen vector 
Objective_s=zeros(length(UpSuccinate),length(UpOxygen)); %initialize biomass matrix 

for i=1:length(UpOxygen) %nested for loop looking at biomass output at each glucose and oxygen combination
    model = changeRxnBounds(model,'EX_o2(e)',-UpOxygen(i),'b');
    for j=1:length(UpSuccinate)
        model = changeRxnBounds(model,'EX_succ(e)',-UpSuccinate(j),'b');
        sol_s=optimizeCbModel(model,'max');
        Objective_s(i,j)=sol_s.obj;
    end
end

fig6=figure;
surfl(UpSuccinate,UpOxygen,Objective_s); %3D plot
xlabel('Succinate Uptake [mmol/gDW/h]')
zlabel('Biomass [mmol/gDW/h]')
ylabel('Oxygen Uptake [mmol/gDW/h]');
%%
%4.3
%Compute shadowprices with built-in function: For figures to work fine,
%need to edit 'phenotypePhasePlane' and remove the numbering after figure
%()
[growthRates, shadowPrices1, shadowPrices2] = phenotypePhasePlane(model, 'EX_succ(e)', 'EX_o2(e)',50,30,30);
%first graph corresponds to glucose

%%CSV Table Generation to observe in ESCHER

Succinate.values=[14.29, 8.163, 5.306, 0.8163]; %Chose datapoints within specific shadowprices
Oxygen.values=[9.796,15.1,16.33,14.29]; %Corresponding data points to glucose vector
Succinate.ShadowPrices=[2, 5, 13, 0]; %Shadow Prices of glucose Values *1e-2

%for loop
for i=1:length(Oxygen.values)
    model=changeRxnBounds(model,{'EX_o2(e)','EX_succ(e)'},-[Oxygen.values(i),Succinate.values(i)],'b');
    Solution_sp_s=optimizeCbModel(model,'max');
    names=strcat("ShadowPrice_Succinate",num2str(Succinate.ShadowPrices(i)),".csv"); 
    Table_sp_s=table(model.rxns,Solution_sp_s.full);
    writetable(Table_sp_s,names)
end
model=changeRxnBounds(model,'EX_succ(e)',1000,'u');
%% PART 4 - Pyruvate
%4.1
model=changeRxnBounds(model,cReactions,0,'b'); %set all carbon substrate uptakes to 0
model=changeRxnBounds(model,cReactions,1000,'u'); %set upper bounds to 1000
UpPyruvate=linspace(0,30,40); %Create vector from 0 to 30 of glucose
model=changeRxnBounds(model,'EX_o2(e)',-20,'b'); %change oxygen max uptake to -20 mmol/gDW/h

Biomass_vector_p=zeros(length(UpPyruvate),1); %Set initialize biomass vector   

%Cycle through all glucose max uptake rates and observe effect on biomass
%production
for i=1:length(UpPyruvate)
model = changeRxnBounds(model,'EX_pyr(e)',-UpPyruvate(i),'b');
sol_p=optimizeCbModel(model,'max');
Biomass_vector_p(i)=sol_p.obj;
end

fig7=figure;
plot(-UpPyruvate,Biomass_vector_p)
xlabel('Pyruvate Uptake [mmol/gDW/h]')
ylabel('Biomass [mmol/gDW/h]')
title('Biomass Production as function of Pyruvate Uptake Rate')
%%
%4.2
model=changeRxnBounds(model,cReactions,0,'b'); %set all carbon substrate uptakes to 0
model=changeRxnBounds(model,cReactions,1000,'u'); %set upper bounds to 1000
UpOxygen=linspace(0,30,40); %oxygen vector 
Objective_p=zeros(length(UpOxygen),length(UpPyruvate)); %initialize biomass matrix 

for i=1:length(UpOxygen) %nested for loop looking at biomass output at each glucose and oxygen combination
    model = changeRxnBounds(model,'EX_o2(e)',-UpOxygen(i),'b');
    for j=1:length(UpPyruvate)
        model = changeRxnBounds(model,'EX_pyr(e)',-UpPyruvate(j),'b');
        sol_p=optimizeCbModel(model,'max');
        Objective_p(i,j)=sol_p.obj;
    end
end

fig8=figure;
surfl(UpPyruvate,UpOxygen,Objective_p); %3D plot
xlabel('Pyruvate Uptake [mmol/gDW/h]')
zlabel('Biomass [mmol/gDW/h]')
ylabel('Oxygen Uptake [mmol/gDW/h]');
%%
%4.3
%Compute shadowprices with built-in function: For figures to work fine,
%need to edit 'phenotypePhasePlane' and remove the numbering after figure
%()
[growthRates, shadowPrices1, shadowPrices2] = phenotypePhasePlane(model, 'EX_pyr(e)', 'EX_o2(e)',50,30,30);
%first graph corresponds to glucose
model=changeRxnBounds(model,cReactions,0,'b'); %set all carbon substrate uptakes to 0
model=changeRxnBounds(model,cReactions,1000,'u'); %set upper bounds to 1000
%%CSV Table Generation to observe in ESCHER

Pyruvate.values=[14.29, 8.163, 5.306, 0.8163]; %Chose datapoints within specific shadowprices
Oxygen.values_p=[9.796,15.1,16.33,14.29]; %Corresponding data points to glucose vector
Pyruvate.ShadowPrices=[2, 5, 13, 0]; %Shadow Prices of glucose Values *1e-2

%for loop
for i=1:length(Oxygen.values)
    model=changeRxnBounds(model,{'EX_o2(e)','EX_succ(e)'},-[Oxygen.values_p(i),Pyruvate.values(i)],'b');
    Solution_sp_p=optimizeCbModel(model,'max');
    names=strcat("ShadowPrice_Pyruvate",num2str(Pyruvate.ShadowPrices(i)),".csv"); 
    Table_sp_p=table(model.rxns,Solution_sp_p.full);
    writetable(Table_sp_p,names)
end

