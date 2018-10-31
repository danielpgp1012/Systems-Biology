%%clear
restoredefaultpath

%% Load solver
%Ensure to change path in load_cplex
cplex_loaded = load_cplex;

if ~cplex_loaded
    error('CPLEX not found')
end

%% Load the COBRA Toolbox
addpath("C:\Users\danie\Documents\GitHub\cobratoolbox")
initCobraToolbox
%%  Load dependencies
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\thermoDatabases");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\models");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\plotting");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA\thermo");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA\utilities");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA\sampling");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA\io");
%addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\ext\Cobra205_vDec2014");

%addpath(genpath(fullfile('..','matTFA')))
%addpath(genpath(fullfile('..','thermoDatabases')))
%addpath(genpath(fullfile('..','models')))
%addpath(genpath(fullfile('..','plotting')))

%% Change Cobra Solver
changeCobraSolver('ibm_cplex');
%changeCobraSolver('glpk');

%% Load Model
load("small_ecoli.mat");

%% Load Thermo Database
tmp = load('thermo_data.mat');
ReactionDB = tmp.DB_AlbertyUpdate;
clear tmp

%% FBA solution q1
model_red=changeObjective(model_red,'Ec_biomass_iJO1366_WT_53p95M');
substrates.FBA={'DM_glc_e','DM_lac-D_e','DM_ac_e','DM_etoh_e'};
oxygen={'DM_o2_e'};
model_red=changeRxnBounds(model_red,oxygen,-1000,'l'); %set aerobic conditions
model_red=changeRxnBounds(model_red,oxygen,1000,'u');

model_red=changeRxnBounds(model_red,substrates.FBA,0,'l');
model_red=changeRxnBounds(model_red,substrates.FBA,1000,'u');

biomass.FBA_aerobic=zeros(length(substrates.FBA),1);

biomass.TFA_aerobic=zeros(length(substrates.FBA),1);
biomass.TFA_anaerobic=zeros(length(substrates.FBA),1); 


for i=1:length(substrates.FBA)
    %aerobic first
    model_red=changeRxnBounds(model_red,substrates.FBA{i},-10,'l');
    solFBA = optimizeCbModel(model_red);
    biomass.FBA_aerobic(i)=solFBA.f; %assign objective
    model_red=changeRxnBounds(model_red,substrates.FBA{i},0,'l');
end

%anaerobic
biomass.FBA_anaerobic=zeros(length(substrates.FBA),1);
model_red=changeRxnBounds(model_red,oxygen,0,'l'); 
for i=1:length(substrates.FBA)
    %aerobic first
    model_red=changeRxnBounds(model_red,substrates.FBA{i},-10,'l');
    solFBA = optimizeCbModel(model_red);
    biomass.FBA_anaerobic(i)=solFBA.f; %assign objective
    model_red=changeRxnBounds(model_red,substrates.FBA{i},0,'l');
end
model_red=changeRxnBounds(model_red,oxygen,-1000,'l');
model_red=changeRxnBounds(model_red,substrates.FBA,-1000,'l');
clear solFBA
%% TFA

%prepped_model = prepModelforTFA(model_red, ReactionDB, model_red.CompartmentData); %Adds reaction data and compartment data
                                                                                   %creates constraints

%% Convert to TFA
%min_obj = 0;
%tmp = convToTFA(prepped_model, ReactionDB, [], 'DGo', [], min_obj); %

%% Add net flux variables, which are equal to forwards flux - backwards flux

%this_tmodel = addNetFluxVariables(tmp); %We split every reaction into forward and backward. It is inconvenient.
%clear tmp
%clear prepped_model
%clear ReactionDB
%% Remove constraints for reverse and forward reactions
%metnames_r=strcat('R_',this_tmodel.rxns);
%this_tmodel.var_lb(find(ismember(this_tmodel.varNames,metnames_r))) = 0;
%this_tmodel.var_ub(find(ismember(this_tmodel.varNames,metnames_r))) = 1000;

%% save model
%save this_tmodel
%% Load model
%Start here
load this_tmodel


%% TFA q1
% Get the variables representing the net fluxes
substrates_TFA={'R_DM_glc_e';'R_DM_lac-D_e';'R_DM_ac_e';'R_DM_etoh_e'};
substrate_indices=find_cell(substrates_TFA,this_tmodel.varNames);
this_tmodel.var_ub(substrate_indices) = 0;

%Change oxygen boundary
oxygen_name={'NF_DM_o2_e'};
oxygen_index=find_cell(oxygen_name,this_tmodel.varNames);
this_tmodel.var_lb(oxygen_index) = -1000;
%aerobic conditions
for i=1:length(substrates_TFA)
    this_tmodel.var_ub(substrate_indices(i)) = 10;
    soltFA = solveTFAmodelCplex(this_tmodel);
    biomass.TFA_aerobic(i)=soltFA.val;
    this_tmodel.var_ub(substrate_indices(i)) = 0;
end

%anaerobic conditions
this_tmodel.var_lb(oxygen_index) = 0;

for i=1:length(substrates_TFA)
    this_tmodel.var_ub(substrate_indices(i)) = 10;
    soltFA = solveTFAmodelCplex(this_tmodel);
    biomass.TFA_anaerobic(i)=soltFA.val;
    this_tmodel.var_ub(substrate_indices(i)) = 0;
end
this_tmodel.var_lb(oxygen_index)=-1000;
save biomass
%% Generate Comparison Table q1
T1=table(substrates_TFA,round(biomass.FBA_aerobic,3),round(biomass.FBA_anaerobic,3),round(biomass.TFA_aerobic,3),round(biomass.TFA_anaerobic,3)); %yield is normalized to consumption
T1.Properties.VariableNames = {'Substrates' 'FBA_Aerobic_Yield' 'FBA_Anaerobic_Yield' 'TFA_Aerobic_Yield' 'TFA_Anaerobic_Yield'}
%We see that the solutions are close enough

%% Load table q2
load('metabolomics.mat','metabolomics') %cell containing: 1: LC_IDs, 2:Names 3: Concentrations 4: Standard Deviations (50% for those not defined)
%% Set variables from table
IDs=metabolomics{1};
metabolite_names=metabolomics{2};
concentrations=metabolomics{3}*1e-3; 
std_dev_concentrations=metabolomics{4}*1e-3;
%% Find indices
IDs_in_model=IDs(ismember(IDs,this_tmodel.varNames)); %list of names in the model... size of IDs in model

indices_inmetabolomics=ismember(IDs,this_tmodel.varNames); %1 if in model, 0 if it is not... size of IDS

metabolomics_indices=find_cell(IDs_in_model,this_tmodel.varNames); %indices in model of metabolomics ... 

IDs_in_model_lb=IDs_in_model(log(concentrations(indices_inmetabolomics)-std_dev_concentrations(indices_inmetabolomics))>this_tmodel.var_lb(metabolomics_indices)); %Select indices where conc-stddev>lb in model

indices_inmetabolomics_lb=ismember(IDs,IDs_in_model_lb); %1 if in model and concentration-stdev>lb, otherwise 0

metabolomics_indices_lb=find_cell(IDs_in_model_lb,this_tmodel.varNames); %indices to change lower bounds
%% Change model
%Changes bounds of concentrations of metabolomics that are in the model
this_model_metabolomics=this_tmodel;

%% change 
this_model_metabolomics.var_ub(metabolomics_indices) = log(concentrations(indices_inmetabolomics)+std_dev_concentrations(indices_inmetabolomics));

this_model_metabolomics.var_lb(metabolomics_indices_lb) = log(concentrations(indices_inmetabolomics_lb)-std_dev_concentrations(indices_inmetabolomics_lb));

this_model_metabolomics.var_ub(substrate_indices) = 0;


%Change oxygen boundary
oxygen_name={'NF_DM_o2_e'};
oxygen_index=find_cell(oxygen_name,this_tmodel.varNames);
this_model_metabolomics.var_lb(oxygen_index) = -1000;


%aerobic conditions
biomass.TFA_aerobic_2=zeros(length(substrates_TFA),1);
biomass.TFA_anaerobic_2=zeros(length(substrates_TFA),1);

for i=1:length(substrates_TFA)
    this_model_metabolomics.var_ub(substrate_indices(i)) = 10;
    soltFA = solveTFAmodelCplex(this_model_metabolomics);
    biomass.TFA_aerobic_2(i)=soltFA.val;
    this_model_metabolomics.var_ub(substrate_indices(i)) = 0;
end

%anaerobic conditions
this_model_metabolomics.var_lb(oxygen_index) = 0;

for i=1:length(substrates_TFA)
    this_model_metabolomics.var_ub(substrate_indices(i)) = 10;
    soltFA = solveTFAmodelCplex(this_model_metabolomics);
    biomass.TFA_anaerobic_2(i)=soltFA.val;
    this_model_metabolomics.var_ub(substrate_indices(i)) = 0;
end
clear metabolomics
clear concentrations
clear std_dev_concentrations
clear solTFA
save this_model_metabolomics;

%% Table 2
T2=table(substrates_TFA,round(biomass.FBA_aerobic,3),round(biomass.FBA_anaerobic,3),(round(biomass.TFA_aerobic_2,3)),(round(biomass.TFA_anaerobic_2,3))); %yield is normalized to consumption
T2.Properties.VariableNames = {'Substrates' 'FBA_Aerobic_Yield' 'FBA_Anaerobic_Yield' 'TFA_Aerobic_Yield' 'TFA_Anaerobic_Yield'}

%% Table 3
%metabolomics used in model
fprintf("The metabolites used in the model are the following:\n");
metabolite_names(indices_inmetabolomics)
fprintf("The metabolites not used in the model are:\n");
metabolite_names(~(indices_inmetabolomics))

%% Question 3
load biomass
%% Aerobic Variable initialization
%Initialization of storage variables
fluxes_FBA_aerobic=cell(length(substrates.FBA),1);
fluxes_TFA_aerobic=cell(length(substrates_TFA),1);
for cid=1:length(fluxes_FBA_aerobic)
fluxes_FBA_aerobic{cid}.min=zeros(length(model_red.rxns),1);
fluxes_FBA_aerobic{cid}.max=zeros(length(model_red.rxns),1);
fluxes_TFA_aerobic{cid}.min=zeros(length(this_tmodel.rxns),1);
fluxes_TFA_aerobic{cid}.max=zeros(length(this_tmodel.rxns),1);
end

%% Aerobic FBA

model_red=changeRxnBounds(model_red,substrates.FBA,0,'l'); %set substrate uptake to 0
for subs=1:length(substrates.FBA)
    model_red=changeRxnBounds(model_red,'Ec_biomass_iJO1366_WT_53p95M',biomass.FBA_aerobic(subs),'b'); %find biomass index
    model_red=changeRxnBounds(model_red,substrates.FBA{subs},-10,'l'); %set specific substrate uptake to -10
    for i=1:length(model_red.rxns) %iterate over all reactions
        model_red=changeObjective(model_red,model_red.rxns{i}); %change objective to ith reaction
        solFBA = optimizeCbModel(model_red); %maximize ith reaction
        fluxes_FBA_aerobic{subs}.max(i)=solFBA.f; %store it in fluxes_FBA cell
        solFBA = optimizeCbModel(model_red,'min'); %minimize ith reaction
        fluxes_FBA_aerobic{subs}.min(i)=solFBA.f; %store it
        %fprintf('Substrate: %s,Metabolite: %s, min: %0.3f , max: %0.3f\n', substrates.FBA{subs},model_red.rxns{i},fluxes_FBA{subs}.min(i),fluxes_FBA{subs}.max(i))
    end
    model_red=changeRxnBounds(model_red,substrates.FBA{subs},0,'l'); %set substrate uptake back to 0
    
end


%% Anaerobic variable initialization
fluxes_FBA_anaerobic=cell(length(substrates.FBA),1);
fluxes_TFA_anaerobic=cell(length(substrates_TFA),1);
for cid=1:length(fluxes_FBA_anaerobic)
fluxes_FBA_anaerobic{cid}.min=zeros(length(model_red.rxns),1);
fluxes_FBA_anaerobic{cid}.max=zeros(length(model_red.rxns),1);
fluxes_TFA_anaerobic{cid}.min=zeros(length(this_tmodel.rxns),1);
fluxes_TFA_anaerobic{cid}.max=zeros(length(this_tmodel.rxns),1);
end
%% Anaerobic FBA

model_red=changeRxnBounds(model_red,oxygen,0,'l'); %set anaerobic conditions
model_red=changeRxnBounds(model_red,substrates.FBA,0,'l');
for subs=1:length(substrates.FBA)-2
    model_red=changeRxnBounds(model_red,'Ec_biomass_iJO1366_WT_53p95M',biomass.FBA_anaerobic(subs),'b');
    model_red=changeRxnBounds(model_red,substrates.FBA{subs},-10,'l');
    for i=1:length(model_red.rxns)
        
        model_red=changeObjective(model_red,model_red.rxns{i});
        solFBA = optimizeCbModel(model_red);
        fluxes_FBA_anaerobic{subs}.max(i)=solFBA.f;
        solFBA = optimizeCbModel(model_red,'min');
        fluxes_FBA_anaerobic{subs}.min(i)=solFBA.f;
        %fprintf('Substrate: %s,Metabolite: %s, min: %0.3f , max: %0.3f\n', substrates.FBA{subs},model_red.rxns{i},fluxes_FBA{subs}.min(i),fluxes_FBA{subs}.max(i))
    end
    model_red=changeRxnBounds(model_red,substrates.FBA{subs},0,'l');
    
end
clear solFBA

%% TFA

%Find names of compounds in model
metnames=strcat('NF_',this_tmodel.rxns); 


biomass_index=find_cell('F_Ec_biomass_iJO1366_WT_53p95M',this_tmodel.varNames); %find biomass index

%% Aerobic TFA
this_tmodel.var_ub(substrate_indices)=0; %set substrates to 0

this_tmodel.var_lfb(oxygen_index)=-1000; %set aerobic conditions

for subs=1:length(substrates.FBA)
    this_tmodel.var_lb(biomass_index)=biomass.TFA_aerobic(subs); %change lower and upper bounds of biomass to max
    
    this_tmodel.var_ub(substrate_indices(subs))=10; %change substrate uptake to -10 mmol/gh
    
    [fluxes_TFA_aerobic{subs}.min,fluxes_TFA_aerobic{subs}.max]=min_max(this_tmodel,metnames);
    
    this_tmodel.var_ub(substrate_indices(subs))=0;
    
end
%% save fluxes
save ('fluxes_TFA_aerobic')

%% Anaerobic TFA

this_tmodel.var_lb(oxygen_index)=0;
for subs=1:length(substrates.FBA)-2
    this_tmodel.var_lb(biomass_index)=biomass.TFA_anaerobic(subs); %change lower and upper bounds of biomass to max
    
    this_tmodel.var_ub(substrate_indices(subs))=10; %change substrate uptake to -10 mmol/gh
    
    [fluxes_TFA_aerobic{subs}.min,fluxes_TFA_aerobic{subs}.max]=min_max(this_tmodel,metnames);
    
    this_tmodel.var_ub(substrate_indices(subs))=0;
    
end
%% Save anaerobic Flux

save ('fluxes_TFA_anaerobic')

%% Directionalities
%initialize cells
sign_FBA_aerobic=cell(length(substrates_TFA),1);
sign_FBA_anaerobic=cell(length(substrates_TFA),1);
sign_TFA_aerobic=cell(length(substrates_TFA),1);
sign_TFA_anaerobic=cell(length(substrates_TFA),1);
sign_table_aerobic=cell(length(substrates_TFA),1);
sign_table_anaerobic=cell(length(substrates_TFA),1);
%iterate over substrates
for i=1:length(sign_FBA_aerobic)
    %find signs of FBA and TFA separately
    sign_FBA_aerobic{i}=sign(sign(fluxes_FBA_aerobic{i}.max)+sign(fluxes_FBA_aerobic{i}.min)); %get overall directionality: 0 is bidirectional
    sign_FBA_anaerobic{i}=sign(sign(fluxes_FBA_anaerobic{i}.max)+sign(fluxes_FBA_anaerobic{i}.min));
        
    sign_TFA_aerobic{i}=sign(sign(fluxes_FBA_aerobic{i}.max)+sign(fluxes_FBA_aerobic{i}.min));
 
    sign_TFA_anaerobic{i}=sign(sign(fluxes_FBA_aerobic{i}.max)+sign(fluxes_FBA_aerobic{i}.min));
    
    %compare FBA and TFA signs: if product of two signs is smaller or equal
    %to 0, then TFA has different direction from FBA. Set sign of TFA. If
    %FBA is bidirectional, the product will be 0 and we'll take TFA sign.
    %If TFA is bidirectional, the product will be 0 and we'll take TFA sign
    %(0). If product is larger than 1, make sign be 0 as both TFA and FBA
    %will have same direction
  
    sign_table_aerobic{i}(sign_FBA_aerobic{i}.*sign_TFA_aerobic{i}<=0)=sign_TFA_aerobic{i}(sign_FBA_aerobic{i}.*sign_TFA_aerobic{i}<=0);
    sign_table_aerobic{i}(sign_FBA_aerobic{i}.*sign_TFA_aerobic{i}>0)=0;
    
    sign_table_anaerobic{i}(sign_FBA_anaerobic{i}.*sign_TFA_anaerobic{i}<=0)=sign_TFA_anaerobic{i}(sign_FBA_anaerobic{i}.*sign_TFA_anaerobic{i}<=0);
    sign_table_anaerobic{i}(sign_FBA_anaerobic{i}.*sign_TFA_anaerobic{i}>0)=0;
    
    %Table 
    T_aerobic=table(this_tmodel.rxns,transpose(sign_table_aerobic{i}));
    writetable(T_aerobic,strcat(substrates_TFA{i},'directionalityaerobic.csv'));
    if (i<3) %skip 0 biomass for anaerobic conditions
        T_anaerobic=table(this_tmodel.rxns,transpose(sign_table_anaerobic{i}));
        
        writetable(T_anaerobic,strcat(substrates_TFA{i},'directionalityanaerobic.csv'));
    end
end

%% 4
% Load Variables
load biomass
%% Initialize Storage Variables
metnames=strcat('LC_',this_tmodel.mets);

mets_aerobic=cell(length(substrate_indices),1);
mets_anaerobic=cell(length(substrate_indices),1);
for i=1:length(substrate_indices)
   mets_aerobic{i}.max=zeros(length(this_tmodel.mets),1);
   mets_aerobic{i}.min=zeros(length(this_tmodel.mets),1);
   mets_anaerobic{i}.max=zeros(length(this_tmodel.mets),1);
   mets_anaerobic{i}.min=zeros(length(this_tmodel.mets),1);
end

%% TFA Aerobic no metabolomics
this_tmodel.var_ub(substrate_indices)=0; %set substrates to 0

this_tmodel.var_lb(oxygen_index)=-1000; %set aerobic conditions




for subs=1:length(substrates.FBA)
    this_tmodel.var_lb(biomass_index)=biomass.TFA_aerobic(subs); %change lower and upper bounds of biomass to max
    
    this_tmodel.var_ub(substrate_indices(subs))=10; %change substrate uptake to -10 mmol/gh
    

    [mets_aerobic{subs}.min,mets_aerobic{subs}.max]=min_max(this_tmodel,metnames);
         
    
    this_tmodel.var_ub(substrate_indices(subs))=0; 
end

%% TFA Anaerobic no metabolomics
this_tmodel.var_lb(oxygen_index)=0;

this_tmodel.var_ub(substrate_indices)=0;

for subs=1:length(substrates.FBA)-2
   this_tmodel.var_lb(biomass_index)=biomass.TFA_anaerobic(subs); %change lower and upper bounds of biomass to max
    
   this_tmodel.var_ub(substrate_indices(subs))=10; %change substrate uptake to -10 mmol/gh
    
   [mets_anaerobic{subs}.min,mets_anaerobic{subs}.max]=min_max(this_tmodel,metnames);

   this_tmodel.var_ub(substrate_indices(subs))=0; 
end

%% metabolomics
%% load model
load this_model_metabolomics
%% Initialize storage variables
mets_aerobic_met=cell(length(substrate_indices),1);
mets_anaerobic_met=cell(length(substrate_indices),1);
for i=1:length(substrate_indices)
   mets_aerobic_met{i}.max=zeros(length(this_tmodel.mets),1);
   mets_aerobic_met{i}.min=zeros(length(this_tmodel.mets),1);
   mets_anaerobic_met{i}.max=zeros(length(this_tmodel.mets),1);
   mets_anaerobic_met{i}.min=zeros(length(this_tmodel.mets),1);
end

%% TFA aerobic with metabolomics
this_model_metabolomics.var_ub(substrate_indices)=0; %set substrates to 0

this_model_metabolomics.var_lb(oxygen_index)=-1000; %set aerobic conditions

for subs=1:length(substrates.FBA)-3 %only evaluate glucose because no biomass was seen for the rest
    this_model_metabolomics.var_lb(biomass_index)=biomass.TFA_aerobic_2(subs); %change lower and upper bounds of biomass to max
    
    this_model_metabolomics.var_ub(substrate_indices(subs))=10; %change substrate uptake to -10 mmol/gh
    

    [mets_aerobic_met{subs}.min,mets_aerobic_met{subs}.max]=min_max(this_model_metabolomics,metnames);
         
    
    this_model_metabolomics.var_ub(substrate_indices(subs))=0; 
end

%% TFA anaerobic with metabolomics

this_model_metabolomics.var_ub(substrate_indices)=0; %set substrates to 0

this_model_metabolomics.var_lb(oxygen_index)=0; %set aerobic conditions

for subs=1:length(substrates.FBA)-3 %only evaluate glucose because no biomass was seen for the rest
    this_model_metabolomics.var_lb(biomass_index)=biomass.TFA_anaerobic_2(subs); %change lower and upper bounds of biomass to max
    
    this_model_metabolomics.var_ub(substrate_indices(subs))=10; %change substrate uptake to -10 mmol/gh
    

    [mets_anaerobic_met{subs}.min,mets_anaerobic_met{subs}.max]=min_max(this_model_metabolomics,metnames);
         
    
    this_model_metabolomics.var_ub(substrate_indices(subs))=0; 
end

%% Table formation
T4_aerobic=table(metnames(ismember(metnames,this_tmodel.varNames)),mets_aerobic{1}.min,mets_aerobic{1}.max,mets_aerobic_met{1}.min,mets_aerobic_met{1}.max,mets_aerobic{2}.min,mets_aerobic{2}.max,mets_aerobic{3}.min,mets_aerobic{3}.max,mets_aerobic{4}.min,mets_aerobic{4}.max);
T4_aerobic.Properties.VariableNames = {'Metabolites' 'min_glc' 'max_glc' 'min_mets_glc' 'max_mets_glc' 'min_lacD' 'max_lacD' 'min_ac' 'max_ac' 'min_etoh' 'max_etoh'}
writetable(T4_aerobic,'Aerobic_Comparison_Q4.csv');
%% anaerobic table
T4_anaerobic=table(metnames(ismember(metnames,this_tmodel.varNames)),mets_anaerobic{1}.min,mets_anaerobic{1}.max,mets_anaerobic_met{1}.min,mets_anaerobic_met{1}.max,mets_anaerobic{2}.min,mets_anaerobic{2}.max,mets_anaerobic{3}.min,mets_anaerobic{3}.max,mets_anaerobic{4}.min,mets_anaerobic{4}.max);
T4_anaerobic.Properties.VariableNames = {'Metabolites' 'min_glc' 'max_glc' 'min_mets_glc' 'max_mets_glc' 'min_lacD' 'max_lacD' 'min_ac' 'max_ac' 'min_etoh' 'max_etoh'}
writetable(T4_aerobic,'Anaerobic_Comparison_Q4.csv');