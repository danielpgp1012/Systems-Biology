%% Initialise Cobrat Toolbox
initCobraToolbox;

%% Load the model
load ('ecoli_core_model.mat');

%% PART 1
%Q1
model = changeObjective(model,model.rxns{13},1);

%Q2
model = changeRxnBounds(model,'EX_glc(e)',-10.0,'l');

%Q3
FBAsolution = optimizeCbModel(model,'max');

%Q4&5
[selExc, selUpt] = findExcRxns(model, 0); %all exchange reactions -> sel Exc = boolean vector 1 if exc, 0 if not.
uptakes=model.rxns(selExc); 
substratesModel=extractSubNetwork(model,uptakes); %output = Cobra model of subnetwork
cReactions=findCarbonRxns(substratesModel,1); %all reactions of the substrates with at least 1 carbon

biomass.aerobic=zeros(length(cReactions),1); %array of 0
biomass.anaerobic=biomass.aerobic;
modelaerobic=model;
modelanaerobic=model;
modelaerobic=changeRxnBounds(modelaerobic,'EX_o2(e)',-20,'l');
modelanaerobic=changeRxnBounds(modelanaerobic,'EX_o2(e)',0,'l');

for i=1:length(cReactions)
modelaerobic=changeRxnBounds(modelaerobic,cReactions{i},-10.0,'l');
modelanaerobic=changeRxnBounds(modelanaerobic,cReactions{i},-10.0,'l');

FBAsolutionaerobic=optimizeCbModel(modelaerobic,'max');
FBAsolutionanaerobic=optimizeCbModel(modelanaerobic,'max');

biomass.aerobic(i)=FBAsolutionaerobic.obj;
biomass.anaerobic(i)=FBAsolutionanaerobic.obj;

%saving the model with glucose as substrate for after the one-by-one evalution
if (strcmp(cReactions{i},'EX_glc(e)'))
    modelo2save=modelaerobic;
    modelasave=modelanaerobic;
    FBAsolutionaerobicsave=FBAsolutionaerobic;
    FBAsolutionanaerobicsave=FBAsolutionanaerobic;
end

%setting the reactions back to 0 for the next loop
modelaerobic=changeRxnBounds(modelaerobic,cReactions,0.0,'l');
modelanaerobic=changeRxnBounds(modelanaerobic,cReactions,0.0,'l');

end

%Q6
%For the ESCHER maps
To2=table(modelo2save.rxns,FBAsolutionaerobicsave.x); 
writetable(To2, 'aerobic_flux_data.csv')
Tan=table(modelasave.rxns,FBAsolutionanaerobicsave.x) ;
writetable(Tan, 'anaerobic_flux_data.csv')

%Matlab plots
c=categorical(cReactions);
h1=subplot(2,1,1);
bar(c,biomass.aerobic)
title('Aerobic Max. Biomass')
h2=subplot(2,1,2);
bar(c,biomass.anaerobic)
title('Anaerobic Max. Biomass')

%% PART 4

UpGlucose=linspace(0,30,30);
model=changeRxnBounds(model,'EX_o2(e)',-10,'l');
Biomass_vector=zeros(length(UpGlucose),1);

for i=1:length(UpGlucose)
model = changeRxnBounds(model,'EX_glc(e)',-UpGlucose(i),'l');
sol=optimizeCbModel(model,'max');
Biomass_vector(i)=sol.obj;
end

figure
plot(-UpGlucose,Biomass_vector)
xlabel('Glucose Uptake [mmol/gDW/h]')
ylabel('Biomass [mmol/gDW/h]')

UpOxygen=linspace(0,30,30);
Objective=zeros(length(UpGlucose),length(UpOxygen));
for i=1:length(UpGlucose)
    model = changeRxnBounds(model,'EX_glc(e)',-UpGlucose(i),'l');
    for j=1:length(UpOxygen)
        model = changeRxnBounds(model,'EX_o2(e)',-UpOxygen(j),'l');
        sol=optimizeCbModel(model,'max');
        Objective(i,j)=sol.obj;
    end
end

figure
%C = UpGlucose.*UpOxygen;
surfl(-UpGlucose,-UpOxygen,Objective);
xlabel('Glucose Uptake [mmol/gDW/h]')
zlabel('Biomass [mmol/gDW/h]')

ylabel('Oxygen Uptake [mmol/gDW/h]');
%Q3
%TODO

    


