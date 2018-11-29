%% Load the COBRA Toolbox
addpath("C:\Users\danie\Documents\GitHub\cobratoolbox")
addpath('C:\tomlab\tomlab20181122_1')

initCobraToolbox
%% convert model
%model=sbmlTestModelToMat(pwd, pwd);

%% load model
load ('iCbu641.mat');

%% Change objective
model=changeObjective(model,'Biomass',1);

%% Test biomass at given glycerol uptake
%[selExc,selUpt]=findExcRxns(model,0); %outputs exchange reactions

%uptakes=model.rxns(selUpt); %selects substrates from exchange reactions

%substratesModel=extractSubNetwork(model,uptakes); %create submodel only
% with given substrates of substrates

%cReactions=findCarbonRxns(substratesModel,1); %find substrates that contain

load ('cNames') ; %Filter carbon-containing substrates
indices=find(ismember(model.rxnNames,cNames));
cReactions = model.rxns(indices);
model=changeRxnBounds(model,cReactions,0,'l'); %change boundaries
glycerol_rxn=cReactions(find(ismember(cNames,'Glycerol exchange')));
model=changeRxnBounds(model,glycerol_rxn,-10,'l');

%% find solution
solFBA=optimizeCbModel(model);

%% for loop to plot curv
uptake_flux=linspace(0,60);
biomass=zeros(length(uptake_flux),1);
x0=zeros(length(model.rxns),length(uptake_flux));
for i=1:length(uptake_flux)
    model=changeRxnBounds(model,glycerol_rxn,-uptake_flux(i),'l');
    solFBA=optimizeCbModel(model);
    biomass(i)=solFBA.f;
    x0(:,i)=solFBA.full;
end

figure 
plot (uptake_flux,biomass)
ylabel('Growth rate (h^{-1})')
xlabel('Glycerol Uptake flux mmol/gDWh')

%% block formic and hydrogen production

BlockRxnNames={'Hydrogen exchange ','Formate exchange '};
BlockRxnIndices=[find(ismember(model.rxnNames,'Hydrogen exchange ')),find(ismember(model.rxnNames,'Formate exchange '))];
BlockRxnCodes=model.rxns(BlockRxnIndices);
 model=changeRxnBounds(model,BlockRxnCodes,0,'u');
 

biomass2=zeros(length(uptake_flux),1);
x0=zeros(length(model.rxns),length(uptake_flux));
for i=1:length(uptake_flux)
    model=changeRxnBounds(model,glycerol_rxn,-uptake_flux(i),'l');
    solFBA=optimizeCbModel(model);
    biomass2(i)=solFBA.f;
    x0(:,i)=solFBA.full;
end

figure 
plot (uptake_flux,biomass2)
ylabel('Growth rate (h^{-1})')
xlabel('Glycerol Uptake flux mmol/gDWh')
 
%% fix biomass and find 
load ('model_NL.mat')
%% Find Enzyme Reactions
T=strfind(model_NL.rxnNames,'ase');
indeces_enzymes=zeros(length(T),1);
k=1;
for i=1:length(T)
    if ~isempty(T{i})
        indeces_enzymes(k)=i;
        k=k+1;
    end
end
indeces_enzymes(k:end)=[];
model_NL.enzyme_rxn=indeces_enzymes;
model_NL.objFunction='x(17)/sum(x(NLPproblem.enzyme_rxn).^2)';
model_NL.b=zeros(701,1);
%%

uptake_flux=linspace(0,60,30);
model_NL=changeRxnBounds(model_NL,cReactions,0,'l'); %change boundaries
biomass_NL=zeros(length(uptake_flux),1);
%model_NL.x0=ones(891,1)*0.5;
model_NL.objFunction='x(17)';
 model_NL.x0=x0(:,50);
 model_NL=changeObjective(model_NL,'Biomass',1);
for i=1:length(uptake_flux)
   
    %model_NL=changeRxnBounds(model_NL,glycerol_rxn,-uptake_flux(i),'b');
    model_NL=changeRxnBounds(model_NL,glycerol_rxn,-20,'b');
    solFBA=solveCobraNLP(model_NL,'MaxFunctionEvaluations',10e4,'OptimalityTolerance',1e-4);
    biomass_NL(i)=solFBA.full(17);
    model_NL.x0=solFBA.full;
end


%% 
figure 
plot (uptake_flux_NL,biomass_NL)
ylabel('Growth rate (h^{-1})')
xlabel('Glycerol Uptake flux mmol/gDWh')


%% Minimize ATP  production as well
model_NL.w=0.04; %weight factor that correlates with experimental results
model_NL.objFunction='x(17)/((1-NLPproblem.w)*x(1)^2+NLPproblem.w*sum(x(NLPproblem.enzyme_rxn).^2))';

uptake_flux=linspace(0,60,30);
model_NL=changeRxnBounds(model_NL,cReactions,0,'l'); %change boundaries
biomass_NL2=zeros(length(uptake_flux),1);
%model_NL.x0=ones(891,1)*0.5;
 model_NL.x0=x0(:,30);
 model_NL=changeObjective(model_NL,'Biomass',1);
for i=1:length(uptake_flux)
   
    model_NL=changeRxnBounds(model_NL,glycerol_rxn,-uptake_flux(i),'b');
    solFBA=solveCobraNLP(model_NL,'MaxFunctionEvaluations',3e4,'OptimalityTolerance',1e-3);
    biomass_NL2(i)=solFBA.full(17);
    model_NL.x0=solFBA.full;
end



%% Conopt Solver Setup
%x_0 = x0(:,40);
%HessPattern = zeros(size(x,1));
model_NL=changeRxnBounds(model_NL,glycerol_rxn,0,'b');
x_0 = x0(:,30);
NLP = conAssign('objfun', 'objGradient', [], [], model_NL.lb, model_NL.ub, 'NLP', x_0,0, [], ...
                           model_NL.A, model_NL.b, model_NL.b);

glycerol_index = find(ismember(model_NL.rxns,glycerol_rxn));
uptake_flux=linspace(5,60);
biomass_NL3=zeros(length(uptake_flux),1);


for i=1:length(uptake_flux)
   
    NLP.x_L(glycerol_index) = -uptake_flux(i);
    NLP.x_U(glycerol_index) = -uptake_flux(i);
    solConopt = tomRun('conopt',NLP,1);
    biomass_NL3(i)=solConopt.x_k(17);
    NLP.x_0 = solConopt.x_k;
end
%% Comparison plot
uptake_flux_NL=linspace(0,60,30);
uptake_flux=linspace(0,60);
uptake_flux_NL2 = linspace(5,60);

figure
plot(uptake_flux,biomass,uptake_flux,biomass2,...
    uptake_flux_NL2,biomass_NL3);
legend('Max Biomass','Max Biomass no H_2','Max Biomass/sum(v_i^2)}')
ylabel('Growth rate (h^{-1})')
xlabel('Glycerol Uptake flux mmol/gDWh')
