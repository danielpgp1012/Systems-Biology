%%clear
restoredefaultpath

%% Load solver

cplex_loaded = load_cplex;

if ~cplex_loaded
    error('CPLEX not found')
end

%% Load the COBRA Toolbox
addpath(genpath(fullfile('..','ext')))

%%  Load dependencies
addpath(genpath(fullfile('..','matTFA')))
addpath(genpath(fullfile('..','thermoDatabases')))
addpath(genpath(fullfile('..','models')))
addpath(genpath(fullfile('..','plotting')))

%% Change Cobra Solver
changeCobraSolver('cplex_direct')

%% Load Model
model = load("GSmodel_Ecoli.mat");
