%%clear
restoredefaultpath

%% Load solver
%Ensure to change path in load_cplex
cplex_loaded = load_cplex;

if ~cplex_loaded
    error('CPLEX not found')
end

%% Load the COBRA Toolbox
%addpath("C:\Users\danie\Documents\GitHub\cobratoolbox")
%initCobraToolbox
%%  Load dependencies
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\thermoDatabases");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\models");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\plotting");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA\thermo");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA\utilities");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA\sampling");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA\io");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\ext\Cobra205_vDec2014");

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
substrates={'DM_glc_e','DM_lac-D_e','DM_ac_e','DM_etoh_e'};
oxygen={'DM_o2_e'};
model_red=changeRxnBounds(model_red,oxygen,-1000,'l');
model_red=changeRxnBounds(model_red,oxygen,1000,'u');
%Find exchange reactions
[selExc,selUpt]=findExcRxns(model_red,0); %outputs exchange reactions

uptakes=model_red.rxns(selExc); %selects substrates from exchange reactions

substratesModel=extractSubNetwork(model_red,uptakes); %create submodel only
% with given substrates of substrates

cReactions=findCarbonRxns(substratesModel,1); %find substrates that contain
% AT LEAST 1 carbon
cReactions=cReactions(find(contains(cReactions,'_e')));

%Set all uptakes to 0
%model_red=changeRxnBounds(model_red,cReactions,0,'l');
%model_red=changeRxnBounds(model_red,cReactions,1000,'u');


biomass.FBA_aerobic=zeros(length(substrates),1);

biomass.TFA_aerobic=zeros(length(substrates),1);
biomass.TFA_aerobic=zeros(length(substrates),1);
for i=1:length(substrates)
    %aerobic first
    model_red=changeRxnBounds(model_red,substrates{i},-10,'l');
    solFBA = optimizeCbModel(model_red);
    biomass.FBA_aerobic(i)=solFBA.obj; %assign objective
    model_red=changeRxnBounds(model_red,substrates{i},0,'l');
end

%anaerobic
biomass.FBA_anaerobic=zeros(length(substrates),1);
model_red=changeRxnBounds(model_red,oxygen,0,'l'); 
for i=1:length(substrates)
    %aerobic first
    model_red=changeRxnBounds(model_red,substrates{i},-10,'l');
    solFBA = optimizeCbModel(model_red);
    biomass.FBA_anaerobic(i)=solFBA.obj; %assign objective
end

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
%% Solve tFA
soltFA = solveTFAmodelCplex(this_tmodel);

%% Perform TVA
% Get the variables representing the net fluxes
ssubstrate_names={'NF_DM_glc_e','NF_DM_ac_e','NF_DM_lac-D_e','NF_DM_etoh_e'};
substrate_indices=getAllVar(this_tmodel,Substrate_names);
this_tmodel.var_lb(Substrate_indices) = 0;
for i=1:length(substrate_names)
    this_tmodel.var_lb(Substrate_indices{i}) = -10;
    soltFA = solveTFAmodelCplex(this_tmodel);
    biomass.TFA_aerobic(i)=soltFA.val;
end
{'NF_DM_o2_e'}