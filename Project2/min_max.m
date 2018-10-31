function [min_vals,max_vals] = min_max(model,vars)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%  Load dependencies
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\thermoDatabases");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\models");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\plotting");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA\thermo");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA\utilities");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA\sampling");
addpath("C:\Users\danie\Documents\EPFL MA1\Fall 1\Principles and Applications of Systems Biology\MATFA\matTFA\io");
cplex_loaded = load_cplex;

if ~cplex_loaded
    error('CPLEX not found')
end

%%
indices_mets=find(ismember(model.varNames,vars)); %find indices in model
max_vals=zeros(length(indices_mets),1); %initialize variables containing max and mins
min_vals=zeros(length(indices_mets),1);

if (length(indices_mets)<length(vars))
   warning('Some variables were not found in model');  %notify the user if some variables are not in model
end 
model.f(model.f==1)=0; %set all objectives to 0
for i=1:length(indices_mets)
        model.f(indices_mets(i))=1; %change objective to each variable
        soltFA = solveTFAmodelCplex(model); %maximize ith variable
        if isempty(soltFA.val)
            soltFA.val=0;
             warning('no solution');
        end
        max_vals(i)=soltFA.val; %store in max vector
        model.objtype=1; %change from maximize to minimize
        soltFA = solveTFAmodelCplex(model); %minimize ith variable
        if isempty(soltFA.val)
            soltFA.val=0;
            warning('no solution'); %notify user if no solution was found
        end
        min_vals(i)=soltFA.val; %store min value
        model.f(indices_mets(i))=0; %set objective to 0
        model.objtype=-1; %set maximze
        
        %fprintf('Substrate: %s,Metabolite: %s, min: %0.3f , max: %0.3f\n', substrates.FBA{subs},model_red.rxns{i},fluxes_FBA{subs}.min(i),fluxes_FBA{subs}.max(i))
end

