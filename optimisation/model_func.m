function [result, ac_fba, ppa_fba,glc_fba,lcts_fba,biomass,flux, model1,model2,timeVec] = model_func(bounds,model_list)
result = [];
flux = bounds;
    ac_fba= [];
    ppa_fba= [];
    glc_fba= [];
    lcts_fba= [];
    
    b1 = [];
TG0 = readtable('Concentration_table_G0.csv'); %% read initial concentration table
Xini = table2array(TG0(:,2)); %% initial concentration of amino acids
Xsub = table2cell(TG0(:,1)); %% uptake reaction list of amino acids

TG20 = readtable('Concentration_table_G20.csv'); %% read initial concentration table
Xini(:,2) = table2array(TG20(:,2)); %% initial concentration of amino acids
Xsub(:,2) = table2cell(TG20(:,1)); %% uptake reaction list of amino acids

% FBA and dFBA 
initBiomass = .041; %% initial biomass concentration
timeStep = 0.05; nSteps = 240;
exclUptakeRxns = {'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)'};


% Estimated rate and regulatation parameters     
x_glc = [-18.2468,95.1863]; %% glucose
x_prot = [ -0.2301,15.7123,-0.0011,15.0534]; %% proteins
x_lcts = [-0.3840,9.5025]; %% lactose 
k = [x_lcts(2), x_prot(2),x_prot(4)]; 
k_glc = x_glc(2);

model_list(1) = changeRxnBounds(model_list(1),'EX_biomass(e)',bounds(4),'u'); %% G0
model_list(2) = changeRxnBounds(model_list(2),'EX_biomass(e)',bounds(5),'u'); %% G20

biomass =[];
for i = 1:length(model_list)
    % amino acids with unsignificant change in concentration over time
    model_list(i) = changeRxnBounds(model_list(i),'EX_ala_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_ala_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_arg_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_arg_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_ile_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_ile_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_leu_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_leu_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_lys_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_lys_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_met_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_met_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_gln_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_gln_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_his_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_his_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_ser_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_ser_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_phe_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_phe_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_val_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_val_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_thr_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_thr_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_tyr_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_tyr_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_asn_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_asn_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_glu_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_glu_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_pro_L(e)',bounds(1),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_pro_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_gly_L(e)',-bounds(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_asp_L(e)',bounds(1),'u');
    
    % amino acids with significant change in concentration over time
    model_list(i) = changeRxnBounds(model_list(i),'EX_gly_L(e)',bounds(2),'u');
    model_list(i) = changeRxnBounds(model_list(i),'EX_asp_L(e)',bounds(3),'l');
    
    % other metabolites 
    model_list(i) = changeRxnBounds(model_list(i),'EX_glc_D(e)',x_glc(1),'l');  
    model_list(i) = changeRxnBounds(model_list(i),'EX_lcts(e)',x_lcts(1),'l');
    model_list(2) = changeRxnBounds(model_list(2),'EX_WheyProtein_e',x_prot(1),'l');
    model_list(i) = changeRxnBounds(model_list(i),'EX_MedProtein_e',x_prot(1)/2,'l');
    
    FBAsolutionmodel = optimizeCbModel(model_list(i), 'max');
    result(i) = FBAsolutionmodel.f;
    [concentrationMatrix,excRxnNames,timeVec,biomassVec,flux] = mydynamicFBA(model_list(i),Xsub(:,i),Xini(:,i), initBiomass, timeStep, nSteps, Xsub,exclUptakeRxns,k,k_glc);

    [~,acetate] = ismember('EX_ac(e)', excRxnNames);
    [~,ppa] = ismember('EX_ppa(e)', excRxnNames);
    [~,glucose] = ismember('EX_glc_D(e)', excRxnNames);
    [~,lactose] = ismember('EX_lcts(e)', excRxnNames); 
    
    try 
    ac_fba(:,end+1) = full(concentrationMatrix(acetate,:));
    ppa_fba(:,end+1) = full(concentrationMatrix(ppa,:));
    glc_fba(:,end+1) = full(concentrationMatrix(glucose,:));
    lcts_fba(:,end+1) = full(concentrationMatrix(lactose,:));
    biomass(:,end+1) = biomassVec;
    end 
end

model1 = model_list(1);
model2 = model_list(2); 

     writematrix(bounds, 'temp_res.txt');
     writematrix(ac_fba, 'temp_ac.txt');
     writematrix(ppa_fba, 'temp_ppa.txt');
     writematrix(glc_fba, 'temp_glc.txt');
     writematrix(lcts_fba, 'temp_lcts.txt');
     writematrix(b1, 'temp_biomass.txt');  
    fid = fopen('opti_parameters.txt', 'a+');
    fprintf(fid, '%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', bounds(1),bounds(2),bounds(3),bounds(4),bounds(5));
    fclose(fid);
end
