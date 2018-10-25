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

%% FBA
model_red=changeObjective(model_red,'Ec_biomass_iJO1366_WT_53p95M');
substrates.FBA={'DM_glc_e','DM_lac-D_e','DM_ac_e','DM_etoh_e'};
oxygen={'DM_o2_e'};
model_red=changeRxnBounds(model_red,oxygen,-1000,'l');
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

%% TFA
%Check for blocked reactions

prepped_model = prepModelforTFA(model_red, ReactionDB, model_red.CompartmentData); %Adds reaction data and compartment data
                                                                                   %creates constraints

%% Convert to TFA
min_obj = 0;
tmp = convToTFA(prepped_model, ReactionDB, [], 'DGo', [], min_obj); %

%% Add net flux variables, which are equal to forwards flux - backwards flux
% NF_rxn = F_rxn - B_rxn
this_tmodel = addNetFluxVariables(tmp); %We split every reaction into forward and backward. It is inconvenient.
%addNetFluxVariables will merge fwd and backward reactions
%A matrix contains fluxes, energies, deltaGs, and integers. It has much
%more constraints
%Variables: 
%F - Forward
%B - Backward
%NF - Net flux
%varub - Contains boundaries for ALL variables
%to constraint flux we need to constrain NF
%dG delta G
%dGo deltaG standard
%objtype: 
%f: which variable we are optimizing
%objtype: -1 : maximize. 1 : minimize (need to test)
%Vartypes
%C - Continuous: can be any real value
%B - Can only be binary
%
%% save model
save this_tmodel
%% Load model
%Start here
load this_tmodel

%% Solve tFA
soltFA = solveTFAmodelCplex(this_tmodel);

%% Perform TVA
% Get the variables representing the net fluxes
substrates_TFA={'NF_DM_glc_e';'NF_DM_lac-D_e';'NF_DM_ac_e';'NF_DM_etoh_e'};
substrate_indices=find_cell(substrates_TFA,this_tmodel.varNames);
this_tmodel.var_lb(substrate_indices) = 0;

%Change oxygen boundary
oxygen_name={'NF_DM_o2_e'};
oxygen_index=find_cell(oxygen_name,this_tmodel.varNames);
this_tmodel.var_lb(oxygen_index) = -1000;
%aerobic conditions
for i=1:length(substrates_TFA)
    this_tmodel.var_lb(substrate_indices(i)) = -10;
    soltFA = solveTFAmodelCplex(this_tmodel);
    biomass.TFA_aerobic(i)=soltFA.val;
    this_tmodel.var_lb(substrate_indices(i)) = 0;
end

%anaerobic conditions
this_tmodel.var_lb(oxygen_index) = 0;

for i=1:length(substrates_TFA)
    this_tmodel.var_lb(substrate_indices(i)) = -10;
    soltFA = solveTFAmodelCplex(this_tmodel);
    biomass.TFA_anaerobic(i)=soltFA.val;
    this_tmodel.var_lb(substrate_indices(i)) = 0;
end

%% Generate Comparison Table
T1=table(substrates_TFA,round(biomass.FBA_aerobic,3),round(biomass.FBA_anaerobic,3),round(biomass.TFA_aerobic,3),round(biomass.TFA_anaerobic,3)); %yield is normalized to consumption
T1.Properties.VariableNames = {'Substrates' 'FBA_Aerobic_Yield' 'FBA_Anaerobic_Yield' 'TFA_Aerobic_Yield' 'TFA_Anaerobic_Yield'}
%We see that the solutions are close enough

%% 2 
load('metabolomics.mat') %cell containing: 1: LC_IDs, 2:Names 3: Concentrations 4: Standard Deviations (50% for those not defined)
IDs=metabolomics{1};
metabolite_names=metabolomics{2};
concentrations=metabolomics{3};
std_dev_concentrations=metabolomics{4};
metabolomics_indices=find_cell(IDs,this_tmodel.varNames); %Find indices in model that correspond to metabolomics
indices_inmetabolomics=ismember(IDs,this_tmodel.varNames); %1 if metabolomic is in model and 0 if it is not
indices_inmetabolomics_lb=ismember(indices_inmetabolomics,indices_inmetabolomics(concentrations<=std_dev_concentrations));
metabolomics_indices_lb=find_cell(IDs(indices_inmetabolomics_lb),this_tmodel.varNames);
%Changes bounds of concentrations of metabolomics that are in the model
this_tmodel.var_ub(metabolomics_indices) = log(concentrations(indices_inmetabolomics)+std_dev_concentrations(indices_inmetabolomics));

this_tmodel.var_lb(metabolomics_indices_lb) = log(concentrations(indices_inmetabolomics_lb)-std_dev_concentrations(indices_inmetabolomics_lb));

%set substrate uptake to 0
this_tmodel.var_lb(substrate_indices) = 0;

%Change oxygen boundary
oxygen_name={'NF_DM_o2_e'};
oxygen_index=find_cell(oxygen_name,this_tmodel.varNames);
this_tmodel.var_lb(oxygen_index) = -1000;
%aerobic conditions
biomass.TFA_aerobic_2=zeros(length(substrates_TFA),1);
biomass.TFA_anaerobic_2=zeros(length(substrates_TFA),1);
for i=1:length(substrates_TFA)
    this_tmodel.var_lb(substrate_indices(i)) = -10;
    soltFA = solveTFAmodelCplex(this_tmodel);
    biomass.TFA_aerobic2(i)=soltFA.val;
    this_tmodel.var_lb(substrate_indices(i)) = 0;
end

%anaerobic conditions
this_tmodel.var_lb(oxygen_index) = 0;

for i=1:length(substrates_TFA)
    this_tmodel.var_lb(substrate_indices(i)) = -10;
    soltFA = solveTFAmodelCplex(this_tmodel);
    biomass.TFA_anaerobic2(i)=soltFA.val;
    this_tmodel.var_lb(substrate_indices(i)) = 0;
end
%% Table 2
T2=table(substrates_TFA,round(biomass.FBA_aerobic,3),round(biomass.FBA_anaerobic,3),transpose(round(biomass.TFA_aerobic2,3)),transpose(round(biomass.TFA_anaerobic2,3))); %yield is normalized to consumption
T2.Properties.VariableNames = {'Substrates' 'FBA_Aerobic_Yield' 'FBA_Anaerobic_Yield' 'TFA_Aerobic_Yield' 'TFA_Anaerobic_Yield'}

%% Table 3
%metabolomics used in model
fprintf("The metabolites used in the model are the following:\n");
metabolite_names(indices_inmetabolomics)
fprintf("The metabolites not used in the model are:\n");
metabolite_names(~(indices_inmetabolomics))

%% Question 3
indices_metabolites=find(ismember(this_tmodel,'NF'))