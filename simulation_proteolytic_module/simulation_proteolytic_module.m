%% 1) GET Models and optimized parameters  : 

model = readCbModel('B_caccae_P20.mat');  
model2 = readCbModel('B_caccae_GLM.mat');  
dFBA_param = csvread('res.txt');

model = changeRxnBounds(model, 'EX_succ(e)', 10,'u');
model2 = changeRxnBounds(model2, 'EX_succ(e)', 10,'u');

%% initialising dFBA 
initBiomass = .041; timeStep = 0.05; nSteps =  480; % dFBA parameters
gr_res = readtable('growth.csv');

% Initial concentration for the 3 media 
T1 = readtable('Concentration_table_G20.csv');
Xini = table2array(T1(:,2));
Xsub = table2cell(T1(:,1));
T2 = readtable('Concentration_table_G0.csv'); 
Xini2 = table2array(T2(:,2));
Xsub2 = table2cell(T2(:,1));
T3 = readtable('Concentration_table_G2.csv');
Xini3 = table2array(T3(:,2));
Xsub3 = table2cell(T3(:,1));

%% file for graph 
met_res = readtable('met_res.csv');

% General parameters for dfba 
exclUptakeRxns = {'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)'}; % rxn to exclude 

%% regulation parameters 
x_glc = [-18.2468,95.1863];
x_prot = [ -0.2301,15.7123,-0.0011,15.0534]; 
x_lcts = [-0.3840,9.5025];
k = [x_lcts(2), x_prot(2),x_prot(4)]; 
k_glc = x_glc(2); 


% G20 
[concentrationMatrix,excRxnNames,timeVec,biomassVec,flux,concentrations] = mydynamicFBA(model,Xsub,Xini, initBiomass, timeStep, nSteps, Xsub,exclUptakeRxns,k,k_glc);
conc = full(concentrationMatrix);

% G0
[concentrationMatrix2,excRxnNames2,timeVec2,biomassVec2,flux2,concentrations2] = mydynamicFBA(model2,Xsub2,Xini2, initBiomass, timeStep, nSteps, Xsub,exclUptakeRxns,k,k_glc);
conc2 = full(concentrationMatrix2);
% G2 
[concentrationMatrix3,excRxnNames3,timeVec3,biomassVec3,flux3,concentrations3] = mydynamicFBA(model,Xsub3,Xini3, initBiomass, timeStep, nSteps, Xsub,exclUptakeRxns,k,k_glc);
conc3 = full(concentrationMatrix3);


currDate = strrep(datestr(datetime), ':', '_'); 
mkdir(currDate)
cd(currDate)
%% Graph generation - setting and bio data 
SG0 = 'b*';
SG20 = 'r*';
SG2 = 'k*';
SG0_sim = 'b--';
SG20_sim = 'r--';
SG2_sim = 'k--';
x = [0,6,12,24];

%% Generate graph 

x2 = table2array(gr_res(:,1)); %% time point
b1 = table2array(gr_res(:,2)); % G0
b2 = table2array(gr_res(:,4)); % G20
b3 = table2array(gr_res(:,3)); % G2
b1err = table2array(gr_res(:,5)); % G0
b2err = table2array(gr_res(:,7)); % G20
b3err = table2array(gr_res(:,6)); % G2

errorbar(x2,b2,b2err,SG20)
hold on
errorbar(x2,b1,b1err,SG0)
hold on 
errorbar(x2,b3,b3err,SG2)
plot(timeVec,biomassVec,SG20_sim,timeVec3,biomassVec3,SG2_sim,timeVec2,biomassVec2,SG0_sim )

title('Biomass','FontSize', 20, 'FontWeight', 'normal');
laby = sprintf(' Biomass (g.L⁻¹)');
ylabel(laby,'FontSize', 18);
yl = ylim;
ylim([0 yl(2)]);
xlabel('Time (h)','FontSize', 18);
legend('P20','P2', 'GLM','Location','northwest')
set(gca, 'TickLabelInterpreter', 'latex')
set(gca, 'XTickLabel', get(gca, 'XTick'), 'YTickLabel', get(gca, 'YTick'), 'FontSize', 12)
print('g20_g0_g2_biomass','-dpng')
close;


%% Acetate
[~, loc] = ismember('EX_ac(e)',excRxnNames);
[~, loc2] = ismember('EX_ac(e)',excRxnNames2);
[~, loc3] = ismember('EX_ac(e)',excRxnNames3);
acetateG20 = conc(loc, :);
acetateG0 = conc2(loc2, :);
acetateG2 = conc3(loc3, :);
y = [table2array(met_res(9+48:12+48,6))];
err = [table2array(met_res(9+48:12+48,7))];
y2 = [table2array(met_res(49:3+49,6))];
err2 = [table2array(met_res(49:3+49,7))];
y3 = [table2array(met_res(53:56,6))];
err3 = [table2array(met_res(53:56,7))];
errorbar(x,y,err,SG20)
hold on
errorbar(x,y2,err2,SG0)
hold on
errorbar(x,y3,err3,SG2)
hold on 
plot(timeVec, acetateG20,SG20_sim,timeVec2, acetateG0, SG0_sim,timeVec3, acetateG2, SG2_sim)
title(char(table2array(met_res(1+48,1))));
laby = sprintf('%s (mM)', char(table2array(met_res(1+48,1))));
ylabel(laby);
yl = ylim;
ylim([0 yl(2)]);
xlabel('Time (h)');
legend('G20','SG0','G2')
print('g20_g0_g2_acetate','-dpng')
close; 


%% Propionate
[~, loc] = ismember('EX_ppa(e)',excRxnNames);
[~, loc2] = ismember('EX_ppa(e)',excRxnNames2);
[~, loc3] = ismember('EX_ppa(e)',excRxnNames3);
propionateG20 = conc(loc, :);
propionateG0 = conc2(loc2, :);
propionateG2 = conc3(loc3, :);
y = [table2array(met_res(9+60:12+60,6))];
err = [table2array(met_res(9+60:12+60,7))];
y2 = [table2array(met_res(61:3+61,6))];
err2 = [table2array(met_res(61:3+61,7))];
y3 = [table2array(met_res(65:68,6))];
err3 = [table2array(met_res(65:68,7))];
errorbar(x,y,err,SG20)
hold on 
errorbar(x,y2,err2,SG0)
hold on
errorbar(x,y3,err3,SG2)
hold on
plot(timeVec, propionateG20, SG20_sim, timeVec2, propionateG0, SG0_sim,timeVec3, propionateG2, SG2_sim)
title(char(table2array(met_res(1+60,1))));
laby = sprintf('%s (mM)', char(table2array(met_res(1+60,1))));
ylabel(laby);
yl = ylim;
ylim([0 yl(2)]);
xlabel('Time (h)');
legend('G20','G0','G2')
print('g20_g0_g2_propionate','-dpng')
close; 

%% LACTOSE
[~, loc] = ismember('EX_lcts(e)',excRxnNames);
[~, loc2] = ismember('EX_lcts(e)',excRxnNames2);
[~, loc3] = ismember('EX_lcts(e)',excRxnNames3);
lcts = conc(loc, :);
lcts2 = conc2(loc2, :);
lcts3 = conc3(loc3, :);
y = [table2array(met_res(9:12,6))];
err = [table2array(met_res(9:12,7))];
y2 = [table2array(met_res(1:4,6))];
err2 = [table2array(met_res(1:4,7))];
y3 = [table2array(met_res(5:8,6))];
err3 = [table2array(met_res(5:8,7))];
errorbar(x,y,err,SG20)
hold on
errorbar(x,y2,err2,SG0)
hold on 
errorbar(x,y3,err3,SG2)
hold on 
plot(timeVec, lcts, SG20_sim,timeVec2, lcts2, SG0_sim,timeVec3, lcts3, SG2_sim )
title(char(table2array(met_res(1,1))));
laby = sprintf('%s (mM)', char(table2array(met_res(1,1))));
ylabel(laby);
yl = ylim;
ylim([0 yl(2)]);
xlabel('Time (h)');
legend('G20','G0','G2')
print('g20_g0_g2_lactose','-dpng')
close; 

%% glucose
[~, loc] = ismember('EX_glc_D(e)',excRxnNames);
[~, loc2] = ismember('EX_glc_D(e)',excRxnNames2);
[~, loc3] = ismember('EX_glc_D(e)',excRxnNames3);
glc = conc(loc, :);
glc2 = conc2(loc2, :);
glc3 = conc3(loc3, :);
y = [table2array(met_res(9+12:12+12,6))];
err = [table2array(met_res(9+12:12+12,7))];
y2 = [table2array(met_res(13:16,6))];
err2 = [table2array(met_res(13:16,7))];
y3 = [table2array(met_res(17:20,6))];
err3 = [table2array(met_res(17:20,7))];
errorbar(x,y,err,SG20)
hold on
errorbar(x,y2,err2,SG0)
hold on
errorbar(x,y3,err3,SG2)
hold on
plot(timeVec, glc, SG20_sim,timeVec2, glc2, SG0_sim,timeVec3, glc3, SG2_sim )
title(char(table2array(met_res(1+12,1))));
laby = sprintf('%s (mM)', char(table2array(met_res(1+12,1))));
ylabel(laby);
yl = ylim;
ylim([0 yl(2)]);
xlabel('Time (h)');
legend('G20','G0','G2')
print('g20_g0_g2_glucose','-dpng')
close; 

%% protéines
[~, loc] = ismember('EX_WheyProtein_e',excRxnNames);
[~, loc3] = ismember('EX_WheyProtein_e',excRxnNames3);
WP = conc(loc, :);
WP3 = conc3(loc3, :);
[~, loc] = ismember('EX_MedProtein_e',excRxnNames);
[~, loc2] = ismember('EX_MedProtein_e',excRxnNames2);
[~, loc3] = ismember('EX_MedProtein_e',excRxnNames3);
MedP = conc(loc, :);
MedP2 = conc2(loc2, :);
MedP3 = conc3(loc3, :);

%% Protein graph 
whole_protein = WP + MedP;
whole_protein3 = WP3 + MedP3;
y = [table2array(met_res(81:84,6))];
err = [table2array(met_res(81:84,7))];
y2 = [table2array(met_res(73:76,6))];
err2 = [table2array(met_res(73:76,7))];
y3 = [table2array(met_res(77:80,6))];
err3 = [table2array(met_res(77:80,7))];
errorbar(x,y,err,SG20)
hold on
errorbar(x,y2,err2,SG0)
hold on
errorbar(x,y3,err3,SG2)
hold on
plot(timeVec, whole_protein,SG20_sim,timeVec2, MedP2, SG0_sim,timeVec3,whole_protein3,SG2_sim)
title(char(table2array(met_res(73,1))));
laby = sprintf('%s (mM)', char(table2array(met_res(73,1))));
ylabel(laby);
yl = ylim;
ylim([0 yl(2)]);
xlabel('Time (h)');
legend('G20','GO','G2')
print('g20_g0_g2_protein','-dpng')
close; 

%% plot all produced metabolites
met_produced_name = {}; 
met_produced_conc = []; 
met_consumed_name = {}; 
met_consumed_conc = [];
for i = 1:length(excRxnNames) 
    if conc(i,1) < 1000 && conc(i,end) > conc(i,1) && strcmp(excRxnNames{i}(1:2),'EX')
        met_produced_name(end+1) = excRxnNames(i);
        met_produced_conc(end+1,:) = conc(i,:);
    end 
    if conc(i,1) < 1000 && conc(i,end) < conc(i,1) && strcmp(excRxnNames{i}(1:2),'EX')
        met_consumed_name(end+1) = excRxnNames(i);
        met_consumed_conc(end+1,:) = conc(i,:);
    end 
end

subplot(1,2,1)
plot(timeVec,met_produced_conc)
xlabel('Time (h)')
legend(replace(met_consumed_name, {'EX_','_'},{'',' '}),'Location','eastoutside')
title('P20 produced metabolites')
text(-4,40,'A','FontSize',12,'FontWeight','bold') 
subplot(1,2,2)
plot(timeVec,met_consumed_conc)
xlabel('Time (h)')
legend(replace(met_consumed_name, {'EX_','_'},{'',' '}),'Location','eastoutside')
title('P20 consumed metabolites')
text(-4,30,'B','FontSize',12,'FontWeight','bold') 
print('P20_consumed_and_produced_metabolites','-dpng')
close;


met_produced_nameG0 = {}; 
met_produced_concG0 = []; 
met_consumed_nameG0 = {}; 
met_consumed_concG0 = [];

for i = 1:length(excRxnNames2) 
    if conc2(i,1) < 1000 && conc2(i,end) > conc2(i,1) && strcmp(excRxnNames2{i}(1:2),'EX')
        met_produced_nameG0(end+1) = excRxnNames2(i);
        met_produced_concG0(end+1,:) = conc2(i,:);
    end 
    if conc2(i,1) < 1000 && conc2(i,end) < conc2(i,1) && strcmp(excRxnNames2{i}(1:2),'EX')
        met_consumed_nameG0(end+1) = excRxnNames2(i);
        met_consumed_concG0(end+1,:) = conc2(i,:);
    end 
end
plot(timeVec2,met_produced_concG0)
xlabel('Time (h)')
legend(strrep(met_produced_nameG0, 'EX_',''),'Location','eastoutside')
title('G0 produced metabolites')
print('G0_produced_metabolites','-dpng')
close;
plot(timeVec2,met_consumed_concG0)
xlabel('Time (h)')
legend(strrep(met_consumed_nameG0, 'EX_',''),'Location','eastoutside')
title('G0 consumed metabolites')
print('G0_consumed_metabolites','-dpng')
close;


met_produced_nameG2 = {}; 
met_produced_concG2 = []; 
met_consumed_nameG2 = {}; 
met_consumed_concG2 = [];
for i = 1:length(excRxnNames3) 
    if conc3(i,1) < 1000 && conc3(i,end) > conc3(i,1) && strcmp(excRxnNames3{i}(1:2),'EX')
        met_produced_nameG2(end+1) = excRxnNames3(i);
        met_produced_concG2(end+1,:) = conc3(i,:);
    end 
    if conc3(i,1) < 1000 && conc3(i,end) < conc3(i,1) && strcmp(excRxnNames3{i}(1:2),'EX')
        met_consumed_nameG2(end+1) = excRxnNames3(i);
        met_consumed_concG2(end+1,:) = conc3(i,:);
    end     
end

plot(timeVec3,met_produced_concG2)
xlabel('Time (h)')
legend(strrep(met_produced_nameG2, 'EX_',''),'Location','eastoutside')
title('G2 produced metabolites')
print('G2_produced_metabolites','-dpng')
close;
plot(timeVec3,met_consumed_concG2)
xlabel('Time (h)')
legend(strrep(met_consumed_nameG2, 'EX_',''),'Location','eastoutside')
title('G2 consumed metabolites')
print('G2_consumed_metabolites','-dpng')
close;

%% flux that max out
for i = 1:length(flux)
    if model.lb(i) == flux(i,1) && model.lb(i) ~= 0 && strcmp(model.rxns{i}(1:2),'EX')
        sprintf('flux %s is limitating\n', char(model.rxns(i)))
    end
end 

for i = 1:length(flux3)
    if model.lb(i) == flux3(i,1) && model.lb(i) ~= 0 && strcmp(model.rxns{i}(1:2),'EX')
        sprintf('flux %s is limitating\n', char(model.rxns(i)))
    end
end 

for i = 1:length(flux2)
    if model2.lb(i) == flux2(i,1) && model2.lb(i) ~= 0 && strcmp(model2.rxns{i}(1:2),'EX')
        sprintf('flux %s is limitating\n', char(model2.rxns(i)))
    end
end 


%% R SQUARE FOR METABOLITE
met = {'Acetate','Propionate','Glucose','Lactose','Succinate'};
rxn_model = {'EX_ac(e)','EX_ppa(e)','EX_glc_D(e)','EX_lcts(e)','EX_succ(e)'}; 
Rsq = struct(); 
S = struct(); 
timepoint = [0,6,12,24]; 
YcalcG0 = [0;0;0;0];
YcalcG2 = [0;0;0;0];
YcalcG20 = [0;0;0;0];
YmeanTotG0 = []; 
YmeanTotG2 = [];
YmeanTotG20 = [];
YcalcTotG0 = []; 
YcalcTotG2 = []; 
YcalcTotG20 = []; 
YmeanG0 = [];
YmeanG2 = [];
YmeanG20 = [];
for i = 1:length(rxn_model) %% FOR EACH METABOLITE
    [pres,~] = ismember(met_res.Reaction,rxn_model(i));
    
    curr_met = met_res(pres,:);
    
    [a,~] = ismember(curr_met.Media,'G0');
    [b,~] = ismember(curr_met.Media,'G2');
    [c,~] = ismember(curr_met.Media,'G20');
    
    YmeanG0 = curr_met.Mean(a);
    YmeanG2 = curr_met.Mean(b);
    YmeanG20 = curr_met.Mean(c);
    
    YmeanTotG0(:,end+1) = YmeanG0;
    YmeanTotG2(:,end+1) = YmeanG2;
    YmeanTotG20(:,end+1) = YmeanG20;
    
    S.G0.(char(met(i))).Data = YmeanG0;
    S.G2.(char(met(i))).Data = YmeanG2;
    S.G20.(char(met(i))).Data = YmeanG20;
    
    YG0 = curr_met(a,3:5);
    YG2 = curr_met(b,3:5);
    YG20 = curr_met(c,3:5);
    
    reaction = rxn_model(i);
    
    [~, loc] = ismember(reaction,excRxnNames2);  %% G0
    [~, loc2] = ismember(reaction,excRxnNames3); %% G2
    [~, loc3] = ismember(reaction,excRxnNames);  %% G20
    
    for k = 1:length(timepoint)
        v = timeVec == timepoint(k);
        YcalcG0(k) = conc2(loc,v);  %% G0
        YcalcG2(k) = conc3(loc2,v); %% G2
        YcalcG20(k) = conc(loc3,v); %% G20
    end
        YcalcTotG0(:,end+1) = YcalcG0; 
        YcalcTotG2(:,end+1) = YcalcG2; 
        YcalcTotG20(:,end+1) = YcalcG20;
        
        S.G0.(char(met(i))).pred = YcalcG0;
        S.G2.(char(met(i))).pred = YcalcG2;
        S.G20.(char(met(i))).pred = YcalcG20;
        
    for j = 1:3 %% FOR EACH REPLICAT
        YiG0 = table2array(YG0(:,j));
        YiG2 = table2array(YG2(:,j));
        YiG20 = table2array(YG20(:,j));
    end
end
 
%% RELATIVE ERROR FOR GROWTH
timepoint = table2array(gr_res(:,1));

biomG0 = table2array(gr_res(:,2)); 
biomG2 = table2array(gr_res(:,3)); 
biomG20 =  table2array(gr_res(:,4)); 

predG0 = zeros(length(timepoint),1); 
predG2 = zeros(length(timepoint),1); 
predG20 = zeros(length(timepoint),1); 

for i = 1:length(timepoint)
    v = timeVec == timepoint(i);
    predG0(i,1) = biomassVec2(1,v);
    predG2(i,1) = biomassVec3(1,v);
    predG20(i,1) =  biomassVec(1,v);
end 


%% G0/G2/G20 
% loglog scale data vs predicted 
y = [biomG0; biomG2; biomG20;YmeanTotG0(:); YmeanTotG2(:); YmeanTotG20(:)]; % Biological Data (mean)
x = [predG0;predG2; predG20;YcalcTotG0(:); YcalcTotG2(:); YcalcTotG20(:)];  % Predicted Values
Y = [ones(length(y),1) y];
b = x\Y; 
Max = max(x);
Min = min(x); 
Rsq1 = 1 - sum((x - y).^2)/sum((y - mean(y)).^2); % le rsquare 
subplot(2,1,1);
loglog([Min,Max],[Min,Max])
hold on 
loglog(y,x,'*b')
hold on 
loglog(sort(y),b(1)+b(2)*sort(y),':r') % b1 et b2 : indique le biais 
xlabel('Predicted')
ylabel('Actual')
str = sprintf('R-squared = %f',Rsq1); 
ylim([min(y),max(y)]);
xlim([min(x),max(x)]);
ylimits = ylim;
xlimits = xlim; 
str = sprintf('R-squared = %f',Rsq1);
text(0.99*xlimits(1),0.99*ylimits(2),str,'HorizontalAlignment','left','VerticalAlignment','top')
res = (x - y) ./ y; 
subplot(2,1,2); 
semilogx(x,res,'*b')
hold on
semilogx(y,0.*x)
hold on 
semilogx(y,0.1.*ones(length(x)),':r')
hold on
semilogx(y,-0.1.*ones(length(x)),':r')
xlabel('Predicted')
ylabel('Residual')
ylim([-1 2.5])
xlim([10^-1.4,10^1.5])
fig1 = figure(1);
print('G0_G2_G20_loglog_data_vs_predicted','-dpng')
close;

%% G2 G0/G2/G20 BIOMASS
y = [biomG0; biomG2; biomG20]; % Biological Data 
x = [predG0;predG2; predG20];  % Predicted Values
Y = [ones(length(y),1) y];
b = Y\x; % 
Max = max(x);
Min = min(x); 
Rsq1 = 1 - sum((x - y).^2)/sum((y - mean(y)).^2);
subplot(2,1,1);
loglog([Min,Max],[Min,Max])
hold on 
loglog(y,x,'*b')
hold on 
loglog(sort(y),b(1)+b(2)*sort(y),':r')
xlabel('Predicted (log)')
ylabel('Actual (log)')
str = sprintf('R-squared = %f',Rsq1); 
ylim([min(y),max(y)]);
xlim([min(x),max(x)]);
ylimits = ylim;
xlimits = xlim; 
str = sprintf('R-squared = %f',Rsq1);
text(0.99*xlimits(1),0.99*ylimits(2),str,'HorizontalAlignment','left','VerticalAlignment','top')
res = (x - y) ./ y; 
subplot(2,1,2); 
semilogx(x,res,'*b')
hold on
semilogx(y,0.*x)
hold on 
semilogx(y,0.1.*ones(length(x)),':r')
hold on
semilogx(y,-0.1.*ones(length(x)),':r')
xlabel('Predicted')
ylabel('Residual')
ylim([-1 1])
fig1 = figure(1);
fig1.WindowState = 'maximized';
print('G0_G2_G20_loglog_data_vs_predicted_biomass','-dpng')
close;

%% G0/G20
% loglog scale data vs predicted 
y = [biomG0; biomG20;YmeanTotG0(:); YmeanTotG20(:)]; 
x = [predG0; predG20;YcalcTotG0(:);YcalcTotG20(:)];
Y = [ones(length(y),1) y];
b = x\Y; % 
Max = max(x);
Min = min(x); 
Rsq1 = 1 - sum((x - y).^2)/sum((y - mean(y)).^2); % le rsquare 
subplot(2,1,1);
loglog([Min,Max],[Min,Max])
hold on 
loglog(y,x,'*b')
hold on 
xlabel('Predicted (log)')
ylabel('Actual (log)')
str = sprintf('R-squared = %f',Rsq1); 
ylim([min(y),max(y)]);
xlim([min(x),max(x)]);
ylimits = ylim;
xlimits = xlim; 
str = sprintf('R-squared = %f',Rsq1);
text(0.99*xlimits(1),0.99*ylimits(2),str,'HorizontalAlignment','left','VerticalAlignment','top')
res = (x - y) ./ y; 
subplot(2,1,2); 
semilogx(x,res,'*b')
hold on
semilogx(y,0.*x)
hold on 
semilogx(y,0.1.*ones(length(x)),':r')
hold on
semilogx(y,-0.1.*ones(length(x)),':r')
xlabel('Predicted')
ylabel('Residual')
ylim([-1 2.5])
xlim([10^-1.4,10^1.5])
print('G0_G20_loglog_data_vs_predicted','-dpng')
close;

%% G2 
% loglog scale data vs predicted 
y = [ biomG2;YmeanTotG2(:)]; 
x = [predG2;YcalcTotG2(:)];
Y = [ones(length(y),1) y];
b = x\Y; % 
Max = max(x);
Min = min(x); 
Rsq1 = 1 - sum((x - y).^2)/sum((y - mean(y)).^2); % le rsquare 
subplot(2,1,1);
loglog([Min,Max],[Min,Max])
hold on 
loglog(y,x,'*b')
hold on 
xlabel('Predicted (log)')
ylabel('Actual (log)')
str = sprintf('R-squared = %f',Rsq1); 
ylim([min(y),max(y)]);
xlim([min(x),max(x)]);
ylimits = ylim;
xlimits = xlim; 
str = sprintf('R-squared = %f',Rsq1);
text(0.99*xlimits(1),0.99*ylimits(2),str,'HorizontalAlignment','left','VerticalAlignment','top')
res = (x - y) ./ y; 
subplot(2,1,2); 
semilogx(x,res,'*b')
hold on
semilogx(y,0.*x)
hold on 
semilogx(y,0.1.*ones(length(x)),':r')
hold on
semilogx(y,-0.1.*ones(length(x)),':r')
xlabel('Predicted')
ylabel('Residual')
ylim([-1 1])
xlim([10^-1.4,10^1.5])
print('G2_loglog_data_vs_predicted','-dpng')
close;

%% ONLY FOR BIOMASS G0 and G20
y = [biomG0; biomG20]; 
x = [predG0; predG20];
Y = [ones(length(y),1) y];
b = Y\x; % 
Max = max(x);
Min = min(x); 
Rsq1 = 1 - sum((x - y).^2)/sum((y - mean(y)).^2); % le rsquare 
subplot(2,1,1);
loglog([Min,Max],[Min,Max])
hold on 
loglog(y,x,'*b')
hold on 
xlabel('Predicted (log)')
ylabel('Actual (log)')
str = sprintf('R-squared = %f',Rsq1); 
ylim([min(y),max(y)]);
xlim([min(x),max(x)]);
ylimits = ylim;
xlimits = xlim; 
str = sprintf('R-squared = %f',Rsq1);
text(0.99*xlimits(1),0.99*ylimits(2),str,'HorizontalAlignment','left','VerticalAlignment','top')
res = (x - y) ./ y; 
subplot(2,1,2); 
semilogx(x,res,'*b')
hold on
semilogx(y,0.*x)
hold on 
semilogx(y,0.1.*ones(length(x)),':r')
hold on
semilogx(y,-0.1.*ones(length(x)),':r')
xlabel('Predicted')
ylabel('Residual')
ylim([-1 1])
print('G0_G20_loglog_data_vs_predicted_biomass','-dpng')
close;

%% ONLY FOR BIOMASS G2
y = [biomG2]; 
x = [predG2];
Y = [ones(length(y),1) y];
b = Y\x; % 
Max = max(x);
Min = min(x); 
Rsq1 = 1 - sum((x - y).^2)/sum((y - mean(y)).^2); % le rsquare 
subplot(2,1,1);
loglog([Min,Max],[Min,Max])
hold on 
loglog(y,x,'*b')
hold on 
xlabel('Predicted (log)')
ylabel('Actual (log)')
str = sprintf('R-squared = %f',Rsq1); 
ylim([min(y),max(y)]);
xlim([min(x),max(x)]);
ylimits = ylim;
xlimits = xlim; 
str = sprintf('R-squared = %f',Rsq1);
text(0.99*xlimits(1),0.99*ylimits(2),str,'HorizontalAlignment','left','VerticalAlignment','top')
res = (x - y) ./ y; 
subplot(2,1,2); 
semilogx(x,res,'*b')
hold on
semilogx(y,0.*x)
hold on 
semilogx(y,0.1.*ones(length(x)),':r')
hold on
semilogx(y,-0.1.*ones(length(x)),':r')
xlabel('Predicted')
ylabel('Residual')
ylim([-1 1])
xlim([10^-1.39, 10^-0.1])
print('G2_loglog_data_vs_predicted_biomass','-dpng')
close;

%% G0 G2 G20 EACH METABOLITE + BIOMASS 
y = [biomG0; biomG2; biomG20]; % Biological Data (mean)
x = [predG0;predG2; predG20];  % Predicted Values
Y = [ones(length(y),1) y];
b = Y\x; % 
Max = max(x);
Min = min(x); 
Rsq1 = 1 - sum((x - y).^2)/sum((y - mean(y)).^2);  
subplot(4,3,1);
loglog([Min,Max],[Min,Max])
hold on 
loglog(y,x,'*b')
hold on 
loglog(sort(y),b(1)+b(2)*sort(y),':r')
xlabel('Predicted (log)')
ylabel('Actual (log)')
title('Biomass')
str = sprintf('R-squared = %f',Rsq1); 
ylim([min(y),max(y)]);
xlim([min(x),max(x)]);
ylimits = ylim;
xlimits = xlim; 
str = sprintf('R-squared = %f',Rsq1);
text(0.99*xlimits(1),0.99*ylimits(2),str,'HorizontalAlignment','left','VerticalAlignment','top')
res = (x - y) ./ y; 
subplot(4,3,4); 
semilogx(x,res,'*b')
hold on
semilogx(y,0.*x)
hold on 
semilogx(y,0.1.*ones(length(x)),':r')
hold on
semilogx(y,-0.1.*ones(length(x)),':r')
xlabel('Predicted')
ylabel('Residual')
ylim([-1 1])
pos_vs = [2 3 7 8 9];
pos_res = [5 6 10 11 12]; 
met = {'Acetate','Propionate','Glucose','Lactose','Succinate'};
for i = 1:length(met)
    y = [S.G0.(char(met(i))).Data;S.G2.(char(met(i))).Data;S.G20.(char(met(i))).Data];
    x = [S.G0.(char(met(i))).pred;S.G2.(char(met(i))).pred;S.G20.(char(met(i))).pred];
    Y = [ones(length(y),1) y];
    if i == 4 
        b = Y\x;
    else 
    b = x\Y;
    end 
    Max = max(x);
    Min = min(x);
    Rsq1 = 1 - sum((x - y).^2)/sum((y - mean(y)).^2); 
    subplot(4,3,pos_vs(i))
    loglog([Min,Max],[Min,Max])
    hold on
    loglog(y,x,'*b')
    hold on
    loglog(sort(y),b(1)+b(2)*sort(y),':r') 
    xlabel('Predicted (log)')
    ylabel('Actual (log)')
    title(char(met(i)))
    ylimits = ylim;
    xlimits = xlim; 
    str = sprintf('R-squared = %f',Rsq1);
    text(0.99*xlimits(1),0.99*ylimits(2),str,'HorizontalAlignment','left','VerticalAlignment','top')
    res = (x - y) ./ y;
    subplot(4,3,pos_res(i))
    semilogx(x,res,'*b')
    hold on
    semilogx(y,0.*x)
    hold on
    semilogx(y,0.1.*ones(length(x)),':r')
    hold on
    semilogx(y,-0.1.*ones(length(x)),':r')
    xlabel('Predicted')
    ylabel('Residual')
end
fig1 = figure(1);
fig1.WindowState = 'maximized';
print('G0_G2_G20_loglog_data_vs_predicted_all_param','-dpng')
close;
%% G0 G20 EACH METABOLITE + BIOMASS 
y = [biomG0; biomG20]; % Biological Data (mean)
x = [predG0; predG20];  % Predicted Values
Y = [ones(length(y),1) y];
b = x\Y; 
Max = max(x);
Min = min(x); 
Rsq1 = 1 - sum((x - y).^2)/sum((y - mean(y)).^2); % le rsquare 
subplot(4,3,1);
loglog([Min,Max],[Min,Max])
hold on 
loglog(y,x,'*b')
hold on 
xlabel('Predicted (log)')
ylabel('Actual (log)')
str = sprintf('R-squared = %f',Rsq1); 
ylim([min(y),max(y)]);
xlim([min(x),max(x)]);
ylimits = ylim;
xlimits = xlim; 
str = sprintf('R-squared = %f',Rsq1);
text(0.99*xlimits(1),0.99*ylimits(2),str,'HorizontalAlignment','left','VerticalAlignment','top')
title('Biomass')
res = (x - y) ./ y; 
subplot(4,3,4); 
semilogx(x,res,'*b')
hold on
semilogx(y,0.*x)
hold on 
semilogx(y,0.1.*ones(length(x)),':r')
hold on
semilogx(y,-0.1.*ones(length(x)),':r')
xlabel('Predicted')
ylabel('Residual')
ylim([-1 1])
pos_vs = [2 3 7 8 9];
pos_res = [5 6 10 11 12]; 
met = {'Acetate','Propionate','Glucose','Lactose','Succinate'};
for i = 1:length(met)
    y = [S.G0.(char(met(i))).Data;S.G20.(char(met(i))).Data];
    x = [S.G0.(char(met(i))).pred;S.G20.(char(met(i))).pred];
    Y = [ones(length(y),1) y];
    b = x\Y;
    Max = max(x);
    Min = min(x);
    Rsq1 = 1 - sum((x - y).^2)/sum((y - mean(y)).^2); % le rsquare
    subplot(4,3,pos_vs(i))
    loglog([Min,Max],[Min,Max])
    hold on
    loglog(y,x,'*b')
    xlabel('Predicted (log)')
    ylabel('Actual (log)')
    title(char(met(i)))
    str = sprintf('R-squared = %f',Rsq1);
    ylim([min(y),max(y)]);
    xlim([min(x),max(x)]);
    ylimits = ylim;
    xlimits = xlim; 
    str = sprintf('R-squared = %f',Rsq1);
    text(0.99*xlimits(1),0.99*ylimits(2),str,'HorizontalAlignment','left','VerticalAlignment','top')
    res = (x - y) ./ y;
    subplot(4,3,pos_res(i))
    semilogx(x,res,'*b')
    hold on
    semilogx(x,0.*x)
    hold on
    semilogx(x,0.1.*ones(length(x)),':r')
    hold on
    semilogx(x,-0.1.*ones(length(x)),':r')
    xlabel('Predicted')
    ylabel('Residual')
    ylim([-1 1])
end
fig1 = figure(1);
fig1.WindowState = 'maximized';
print('G0_G20_loglog_data_vs_predicted_all_param','-dpng')
close;

%% G2 EACH METABOLITE + BIOMASS 
y = [biomG2]; % Biological Data (mean)
x = [predG2];  % Predicted Values
Y = [ones(length(y),1) y];
b = x\Y; % 
Max = max(x);
Min = min(x); 
Rsq1 = 1 - sum((x - y).^2)/sum((y - mean(y)).^2);
subplot(4,3,1);
loglog([Min,Max],[Min,Max])
hold on 
loglog(y,x,'*b')
hold on 
xlabel('Predicted (log)')
ylabel('Actual (log)')
title('Biomass')
str = sprintf('R-squared = %f',Rsq1); 
ylimits = ylim;
xlimits = xlim; 
str = sprintf('R-squared = %f',Rsq1);
text(0.99*xlimits(1),0.99*ylimits(2),str,'HorizontalAlignment','left','VerticalAlignment','top')
res = (x - y) ./ y; 
subplot(4,3,4); 
semilogx(x,res,'*b')
hold on
semilogx(y,0.*x)
hold on 
semilogx(y,0.1.*ones(length(x)),':r')
hold on
semilogx(y,-0.1.*ones(length(x)),':r')
xlabel('Predicted')
ylabel('Residual')
ylim([-1 1])
pos_vs = [2 3 7 8 9];
pos_res = [5 6 10 11 12]; 
met = {'Acetate','Propionate','Glucose','Lactose','Succinate'};
for i = 1:length(met)
    y = [S.G2.(char(met(i))).Data];
    x = [S.G2.(char(met(i))).pred];
    Y = [ones(length(y),1) y];
    b = x\Y;
    Max = max(x);
    Min = min(x);
    Rsq1 = 1 - sum((x - y).^2)/sum((y - mean(y)).^2); % le rsquare
    subplot(4,3,pos_vs(i))
    loglog([Min,Max],[Min,Max])
    hold on
    loglog(y,x,'*b')
    hold on
    xlabel('Predicted (log)')
    ylabel('Actual (log)')
    title(char(met(i)))
    str = sprintf('R-squared = %f',Rsq1);
    ylimits = ylim;
    xlimits = xlim; 
    str = sprintf('R-squared = %f',Rsq1);
    text(0.99*xlimits(1),0.99*ylimits(2),str,'HorizontalAlignment','left','VerticalAlignment','top')
    res = (x - y) ./ y;
    subplot(4,3,pos_res(i))
    semilogx(x,res,'*b')
    hold on
    semilogx(x,0.*x)
    hold on
    semilogx(x,0.1.*ones(length(x)),':r')
    hold on
    semilogx(x,-0.1.*ones(length(x)),':r')
    xlabel('Predicted')
    ylabel('Residual')
    ylim([-1 1])
end
fig1 = figure(1);
fig1.WindowState = 'maximized';
print('G2_loglog_data_vs_predicted_all_param','-dpng')
close;

%% TEST ON P20 and REACTION PROTEOLYSIS 
% Comparaison between WP and MD => G0 
% Comparaison with no reaction 
modelP20_without_prot = model; 
modelP20_without_prot = changeRxnBounds(modelP20_without_prot, 'protease_rxn', 0,'u'); 
modelP20_without_prot = changeRxnBounds(modelP20_without_prot, 'med_protease_rxn', 0,'u'); 
modelP20_without_prot = changeRxnBounds(modelP20_without_prot, 'EX_biomass(e)', 0,'l');
[concentrationMatrix_P20_without_prot,excRxnNames_P20_without_prot,timeVec_P20_without_prot,biomassVec_P20_without_prot,flux_P20_without_prot,concentrations_P20_without_prot] = mydynamicFBA(modelP20_without_prot,Xsub,Xini, initBiomass, timeStep, nSteps, Xsub,exclUptakeRxns,k,k_glc);
conc_without_prot = full(concentrationMatrix_P20_without_prot);

modelP20_without_wp = model; 
modelP20_without_wp = changeRxnBounds(modelP20_without_wp,'protease_rxn', 0,'u'); 
modelP20_without_wp = changeRxnBounds(modelP20_without_wp, 'EX_biomass(e)', 0,'l');
[concentrationMatrix_P20_without_wp,excRxnNames_P20_without_wp,timeVec_P20_without_wp,biomassVec_P20_without_wp,flux_P20_without_wp,concentrations_P20_without_wp] = mydynamicFBA(modelP20_without_wp,Xsub,Xini, initBiomass, timeStep, nSteps, Xsub,exclUptakeRxns,k,k_glc);
conc_without_wp = full(concentrationMatrix_P20_without_wp); 

modelP20_without_med = model; 
modelP20_without_med = changeRxnBounds(modelP20_without_med,'med_protease_rxn', 0,'u'); 
modelP20_without_med = changeRxnBounds(modelP20_without_med, 'EX_biomass(e)', 0,'l');
[concentrationMatrix_P20_without_med,excRxnNames_P20_without_med,timeVec_P20_without_med,biomassVec_P20_without_med,flux_P20_without_med,concentrations_P20_without_med] = mydynamicFBA(modelP20_without_med,Xsub,Xini, initBiomass, timeStep, nSteps, Xsub,exclUptakeRxns,k,k_glc);
conc_without_med = full(concentrationMatrix_P20_without_med); 

subplot(1,2,1)
errorbar(x2,b2,b2err,SG20)
hold on
plot(timeVec, biomassVec, 'r--', timeVec_P20_without_prot, biomassVec_P20_without_prot, 'k*-', timeVec_P20_without_wp, biomassVec_P20_without_wp, 'b--', timeVec_P20_without_med, biomassVec_P20_without_med, 'm--', 'MarkerSize', 1)
legend('Experimental data','Sim with both proteolysis modules', 'Sim with no proteolysis module', 'Sim without the whey protein module','Sim without the medium protein module','Location','northwest','FontSize',12)
xlabel('Time (h)', 'FontSize', 15);
ylabel('Biomass (g.L-1)', 'FontSize', 15);
title('Biomass', 'FontSize', 18,'FontWeight', 'normal' ) 
%print('test_reaction_prot_biomass','-dpng')

[~, loc] = ismember('EX_ac(e)',excRxnNames_P20_without_prot);
acetateG20_without_prot = conc_without_prot(loc, :);
[~, loc] = ismember('EX_ac(e)',excRxnNames_P20_without_wp);
acetateG20_without_wp = conc_without_wp(loc, :);
[~, loc] = ismember('EX_ac(e)',excRxnNames_P20_without_med);
acetateG20_without_med = conc_without_med(loc, :);
[~, loc] = ismember('EX_ac(e)',excRxnNames); 
acetateG20 = conc(loc, :);

yac = [table2array(met_res(9+48:12+48,6))];
errac = [table2array(met_res(9+48:12+48,7))];

[~, loc] = ismember('EX_ppa(e)',excRxnNames_P20_without_prot);
propionateG20_without_prot = conc_without_prot(loc, :);
[~, loc] = ismember('EX_ppa(e)',excRxnNames_P20_without_wp); 
propionateG20_without_wp = conc_without_wp(loc, :);
[~, loc] = ismember('EX_ppa(e)',excRxnNames_P20_without_med);
propionateG20_without_med = conc_without_med(loc, :);
[~, loc] = ismember('EX_ppa(e)',excRxnNames); 
propionateG20 = conc(loc, :);

ypp = [table2array(met_res(9+60:12+60,6))];
errpp = [table2array(met_res(9+60:12+60,7))];
subplot(1,2,2)
errorbar(x,yac,errac,'r*') %% ERROR ACEATE
hold on 
errorbar(x,ypp,errpp,'ro') %% ERROR PROPIONATE 
hold on 
plot(timeVec, acetateG20, 'r--',timeVec_P20_without_prot, acetateG20_without_prot,'k*-',timeVec_P20_without_wp,acetateG20_without_wp,'b--',timeVec_P20_without_med,acetateG20_without_med,'m--', 'MarkerSize', 1)
hold on 
plot(timeVec,propionateG20,'r-.',timeVec_P20_without_prot,propionateG20_without_prot, 'k*-',timeVec_P20_without_wp, propionateG20_without_wp, 'b-.',timeVec_P20_without_med,propionateG20_without_med,'m--', 'MarkerSize', 1)
legend('Experimental Ac','Experimental Ppa','Ac with both proteolysis modules', 'Ac with no proteolysis module', 'Ac without the whey protein module','Ac without the medium protein module', 'Ppa with both proteoysis modules','Ppa with no proteolysis module','Ppa without the whey protein module','Ppa without the medium protein module' ,'Location','northwest','FontSize',12)
xlabel('Time (h)', 'FontSize', 15);
ylabel('Metabolites (mM)', 'FontSize', 15);
title('Acetate and Propionate', 'FontSize', 18,'FontWeight', 'normal') 
print('test_reaction_prot','-dpng')

% comparaison with both reactions => G20 

[~, loc] = ismember('EX_ac(e)',excRxnNames);
[~, loc2] = ismember('EX_ac(e)',excRxnNames2);
[~, loc3] = ismember('EX_ac(e)',excRxnNames3);
acetateG20 = conc(loc, :);
acetateG0 = conc2(loc2, :);
acetateG2 = conc3(loc3, :);
yac = [table2array(met_res(9+48:12+48,6))];
errac = [table2array(met_res(9+48:12+48,7))];
y2ac = [table2array(met_res(49:3+49,6))];
err2ac = [table2array(met_res(49:3+49,7))];
y3ac = [table2array(met_res(53:56,6))];
err3ac = [table2array(met_res(53:56,7))];
[~, loc] = ismember('EX_ppa(e)',excRxnNames);
[~, loc2] = ismember('EX_ppa(e)',excRxnNames2);
[~, loc3] = ismember('EX_ppa(e)',excRxnNames3);
propionateG20 = conc(loc, :);
propionateG0 = conc2(loc2, :);
propionateG2 = conc3(loc3, :);
ypp = [table2array(met_res(9+60:12+60,6))];
errpp = [table2array(met_res(9+60:12+60,7))];
y2pp = [table2array(met_res(61:3+61,6))];
err2pp = [table2array(met_res(61:3+61,7))];
y3pp = [table2array(met_res(65:68,6))];
err3pp = [table2array(met_res(65:68,7))];


errorbar(x,yac,errac,'r*')
hold on
errorbar(x,y2ac,err2ac,'b*')
hold on
errorbar(x,y3ac,err3ac,'k*')
hold on 
errorbar(x,ypp,errpp,'ro')
hold on 
errorbar(x,y2pp,err2pp,'bo')
hold on
errorbar(x,y3pp,err3pp,'ko')
hold on
plot(timeVec, acetateG20,'r--',timeVec3, acetateG2, 'k--',timeVec2, acetateG0, 'b--')
hold on 
plot(timeVec, propionateG20, 'r-.',timeVec3, propionateG2, 'k-.',timeVec2, propionateG0, 'b-.')
legend('Acetate P20','Acetate P2','Acetate GLM','Propionate P20','Propionate P2','Propionate GLM','Location','northwest')
xlabel('Time (h)','FontSize', 18);
ylabel('Metabolites (mM)','FontSize', 18);
title('Acetate and Propionate','FontSize', 20, 'FontWeight', 'normal') 
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'XTickLabel', get(gca, 'XTick'), 'YTickLabel', get(gca, 'YTick'), 'FontSize', 12);
print('g20_g0_g2_propionate_acetate','-dpng')
close;
