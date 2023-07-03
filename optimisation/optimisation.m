cd ../cobratoolbox % path to cobratoolbox directory 
initCobraToolbox;

cd ../optimisation % path to the optimisation script 

%% read cb model
model = readCbModel('VMH_Bacteroides_caccae_ATCC_43185.mat');               
model = changeObjective(model, model.rxns(find(model.c)));

% make necessary modification to accomodate B caccae heterotrophia 
model_modification;  

%% protelysis reaction
model = addReaction(model, 'protease_rxn','reactionName','Proteolysis pathway modelling', ...
    'reactionFormula','1 wheyprotein[e] + 23 atp[c] -> 21 adp[c] + 8.3 ala_L[e] + 2.3 arg_L[e] + 13.1 asp_L[e] + 3.6 cys_L[e] + 18.4 glu_L[e] + 3.8 gly[e] + 2 his_L[e] + 6.5 ile_L[e] + 14.2 leu_L[e] + 10.3 lys_L[e] + 2.4 met_L[e] + 3.4 phe_L[e] + 5.9 pro_L[e] + 6.5 ser_L[e] + 6.5 thr_L[e] + 1.7 trp_L[e] + 3.1 tyr_L[e] + 7  val_L[e]', ...
    'lowerBound', 0, 'upperBound', 10,'subsystem','Proteolytic pathway');
    
%% add a second protelysis reaction specific to model the rest of the protein in the model
model = addReaction(model, 'EX_MedProtein_e', 'reactionName', 'WheyProtein Exchange Reaction', ...
    'ReactionFormula', 'wheyprotein[e] <=>','lowerBound', -5, 'upperBound', 0,'subsystem', 'Exchange/demand reaction');
    
model = addReaction(model, 'med_protease_rxn','reactionName','Proteolysis pathway modelling for medium', ...
    'reactionFormula','1 protein[e] + 23 atp[c] -> 21 adp[c] + 5.7 ala_L[e] + 0.8 arg_L[e] + 3.6 asp_L[e] + 6.6 glu_L[e] + 11.3 gly[e] + 0.7 his_L[e] + 2 ile_L[e] + 2.1 leu_L[e] + 0.9 lys_L[e] + 0.4 met_L[e] + 1.1 phe_L[e] + 6 pro_L[e] + 1.1 ser_L[e] + 0.8 thr_L[e] + 0.2 tyr_L[e] + 2.6  val_L[e]', ...
    'lowerBound', 0, 'upperBound', 10,'subsystem','Proteolytic pathway');

%% bounds the model's reaction with -15 and 15 max rates
mll=model.lb;
muu=model.ub;
iml=find(mll<0);
imu=find(muu>0);
mll(iml)=-15;
muu(imu)=15;
model.lb=mll;
model.ub=muu;

%% shut all EX reaction
for i = 1:length(model.rxns)
    if strcmp(model.rxns{i}(1:2),'EX') && model.lb(i)==-15 %%% ATTENTION A BIEN MODIFIER ICI AUSSI
        model.lb(i) = 0;
    end
end

%% get medium information
T = readtable('EX_flux_bounds_initial.csv');
ExFlux = table2array(T(:,2));
ExName = table2cell(T(:,1));

%% switch back on only metabolites present in medium
for j=1:length(ExName)
    if ExFlux(j) < 0 || ExFlux(j) == 0
        model = changeRxnBounds(model,ExName(j),ExFlux(j),'l');
    else
        model = changeRxnBounds(model,ExName(j),ExFlux(j),'u');
    end
end

%% added to make sure those metabolites couldn't be import
model = changeRxnBounds(model, 'EX_succ(e)',0,'l');
model = changeRxnBounds(model, 'EX_ac(e)',0,'l');
model = changeRxnBounds(model, 'EX_succ(e)',0,'l');
model = changeRxnBounds(model, 'EX_ppa(e)',0,'l');
model = changeRxnBounds(model, 'EX_lac_L(e)',0,'l');
model = changeRxnBounds(model, 'EX_lac_D(e)',0,'l');
model = changeRxnBounds(model, 'EX_isoval(e)',0,'u');

%% Make two different models depending on the substrates 
modelG20 = model; 
modelG20 = changeRxnBounds(modelG20, 'EX_biomass(e)',0.8,'u');
modelG20 = changeRxnBounds(modelG20, 'EX_biomass(e)',0.58,'l');

modelG0 = model; 
modelG0 = changeRxnBounds(modelG0, 'EX_biomass(e)',0.6,'u');
modelG0 = changeRxnBounds(modelG0, 'EX_biomass(e)',0.47,'l');
modelG0 = changeRxnBounds(modelG0,'protease_rxn',0,'u');
modelG0 = changeRxnBounds(modelG0,'med_protease_rxn',0,'u');
modelG0 = changeRxnBounds(modelG0,'EX_WheyProtein_e',0,'l');
modelG0 = changeRxnBounds(modelG0,'EX_MedProtein_e',0,'l'); 

% Computation of upper bounds from data
T_growth = readtable('growth.csv');
data = struct();
data.G0.Biomass.Mean = T_growth.MeanG0;
data.G20.Biomass.Mean = T_growth.MeanG20; 
data.G0.Biomass.Stdev = T_growth.SDG0;
data.G0.Biomass.Stdev(data.G0.Biomass.Stdev == 0) = 0.001; 
data.G20.Biomass.Stdev = T_growth.SDG20; 
data.G20.Biomass.Stdev(data.G20.Biomass.Stdev == 0) = 0.001;
data.G0.Biomass.timePoints = T_growth.TimePoints; 
data.G20.Biomass.timePoints = T_growth.TimePoints; 

growth = table2struct(T_growth,"ToScalar",true);

%% Integration biomass
int_0_T_biomassG20 = 0;
int_0_T_biomassG0 = 0; 
for i=1:length(growth.TimePoints)
    curr_time = growth.TimePoints(i);
    if curr_time < 12 % time integration of biomass from time 0 to time 12
        int_0_T_biomassG20 = int_0_T_biomassG20 + (growth.MeanG20(i+1) + growth.MeanG20(i))/2 * (growth.TimePoints(i+1) - growth.TimePoints(i)); % trapezoidal integration scheme
        int_0_T_biomassG0 = int_0_T_biomassG0 + (growth.MeanG0(i+1) + growth.MeanG0(i))/2 * (growth.TimePoints(i+1) - growth.TimePoints(i)); % trapezoidal integration scheme
    end
end

%% get metabolite data from G20 medium
T_metabolites = readtable('metabolite_dosage.csv');
uniq_met = unique(T_metabolites.Metabolites); 
data.TimePoints = T_metabolites.Temps;

for i = 1:length(uniq_met)
    data.G0.(char(uniq_met{i})).Mean= []; 
    data.G0.(char(uniq_met{i})).Stdev = [];
    data.G20.(char(uniq_met{i})).Mean = []; 
    data.G20.(char(uniq_met{i})).Stdev = []; 
end 

m = 0; 
n = 0; 
for u = 1:length(uniq_met)
    for i = 1:size(T_metabolites,1)
        if isequal(T_metabolites.Substrate(i), {'G0' }) && isequal(T_metabolites.Metabolites(i), uniq_met(u))
            data.G0.(char(uniq_met{u})).Mean(end+1,1) = T_metabolites.Mean(i);
            data.G0.(char(uniq_met{u})).Stdev(end+1,1) = T_metabolites.Stdev(i);
            if data.G0.(char(uniq_met{u})).Stdev(end,1) == 0
               data.G0.(char(uniq_met{u})).Stdev(end,1) = 0.001; 
            end 
            data.G0.(char(uniq_met{u})).ID = T_metabolites.Model_idx(i);
        end
        if isequal(T_metabolites.Substrate(i), {'G20' }) && isequal(T_metabolites.Metabolites(i), uniq_met(u))
           data.G20.(char(uniq_met{u})).Mean(end+1,1) = T_metabolites.Mean(i);
           data.G20.(char(uniq_met{u})).Stdev(end+1,1) = T_metabolites.Stdev(i);
           if data.G20.(char(uniq_met{u})).Stdev(end,1) == 0
               data.G20.(char(uniq_met{u})).Stdev(end,1) = 0.001; 
            end 
           data.G0.(char(uniq_met{u})).ID = T_metabolites.Model_idx(i);
           n = n+1;  
        end
    end
end

for u = 1:length(uniq_met)
    curr_met = data.G0.(char(uniq_met{u})).ID;
    max_fluxG0 = (data.G0.(char(uniq_met{u})).Mean(3) - data.G0.(char(uniq_met{u})).Mean(1))/int_0_T_biomassG0;
    max_fluxG20 = (data.G20.(char(uniq_met{u})).Mean(3) - data.G20.(char(uniq_met{u})).Mean(1))/int_0_T_biomassG20;
    if max_fluxG0<0 
        modelG0 = changeRxnBounds(modelG0,curr_met,max_fluxG0,'l');
    end
    if max_fluxG0>0
        modelG0 = changeRxnBounds(modelG0,curr_met,max_fluxG0,'u');
    end 
    if  max_fluxG20<0 
        modelG20 = changeRxnBounds(modelG20,curr_met,max_fluxG20,'l');
    end
    if max_fluxG20>0
        modelG20 = changeRxnBounds(modelG20,curr_met,max_fluxG20,'u');
    end
end

[~,loc] = ismember('EX_WheyProtein_e', modelG20.rxns); 
modelG20 = changeRxnBounds(modelG20,'EX_MedProtein_e',modelG20.lb(loc)/2,'l'); 
modelG0 = changeRxnBounds(modelG0,'EX_MedProtein_e',modelG20.lb(loc)/2,'l'); 
modelG0 = changeRxnBounds(modelG0,'EX_WheyProtein_e',0,'l'); 

modelG20 = changeRxnBounds(modelG20,'EX_lcts(e)',0,'u');
modelG0 = changeRxnBounds(modelG0,'EX_succ(e)',10,'u');
modelG20 = changeRxnBounds(modelG20,'EX_lcts(e)',0,'u');
modelG0 = changeRxnBounds(modelG0,'EX_succ(e)',10,'u');

%% Get metabolites production and consumption for the observation function 
Mets_for_obs = {'Glucose','Lactose','Acetate','Propionate'}; 

model_list=[modelG0 modelG20]; 

%% Objective function:
restrected_objective_func = @(bounds) objective_func(bounds,data,model_list);
%% Optimization
bornes_ini = [5,5,-5,0.6,0.8]; % [ general bounds for amino acids, glycine, aspartate, G0 initial growth rate, G20 initial growth rate]

lb = [0,0,-15,0,0];
ub = [15,15,0,5,5];
A = [];
b = [];
Aeq = [];
beq = [];

fid = fopen('opti_parameters_G0_G20_co_opti.txt', 'w');
fprintf(fid, 'GenUB\tGly\tAsp\tmuG0\tmuG20\n');
fclose(fid);


options = optimset('Display', 'iter', 'MaxFunEvals',100000, 'MaxIter',100000);

[res,fval,exitflag,output]= fmincon(restrected_objective_func,bornes_ini,A,b, Aeq,beq, lb,ub,[],options);

%Growth rate prediction using the optimal bounds
[result, ac_fba, ppa_fba,glc_fba,lcts_fba,b1,flux, model1,model2,timeVec] = model_func(res,model_list);

x = 'B_caccae_GLM';
writeCbModel(model1,'format','mat', 'fileName', x)
x2 = 'B_caccae_P20';
writeCbModel(model2,'format','mat', 'fileName', x2)


res_biom = b1;
res_ac = ac_fba;
res_ppa = ppa_fba;
res = res;

writematrix(res_biom);
writematrix(res_ac);
writematrix(res_ppa);
writematrix(res);
