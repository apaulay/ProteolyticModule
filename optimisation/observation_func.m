function [calc_distance] = observation_func(bounds,data,model_list)
calc_distance = [];
[~, ac_fba , ppa_fba ,glc_fba,lcts_fba,b1_time_vect,flux, model1,model2,timeVec] = model_func(bounds,model_list);

%% biological data 
data_exp = [data.G0.Biomass.Mean(1:12); data.G20.Biomass.Mean(1:12)];
ac = [data.G0.Acetate.Mean(2:3); data.G20.Acetate.Mean(2:3)]; 
ppa = [data.G0.Propionate.Mean(2:3); data.G20.Propionate.Mean(2:3)]; 
glc = [data.G0.Glucose.Mean(2:3); data.G20.Glucose.Mean(2:3)]; 
lcts = [data.G0.Lactose.Mean(2:3); data.G20.Lactose.Mean(2:3)]; 

%% stdev of the biological data 
data_exp_stdev = [data.G0.Biomass.Stdev(1:12); data.G20.Biomass.Stdev(1:12)];
ac_stdev = [data.G0.Acetate.Stdev(2:3); data.G20.Acetate.Stdev(2:3)]; 
ppa_stdev = [data.G0.Propionate.Stdev(2:3); data.G20.Propionate.Stdev(2:3)]; 
glc_stdev = [data.G0.Glucose.Stdev(2:3); data.G20.Glucose.Stdev(2:3)]; 
lcts_stdev = [data.G0.Lactose.Stdev(2:3); data.G20.Glucose.Stdev(2:3)]; 

%% Get simulation data 
ac_sim = [0, 0;0, 0];
ppa_sim = [0, 0;0, 0];
glc_sim = [0, 0;0, 0];
lcts_sim = [0, 0;0, 0];

if ~isempty(ac_fba)
for i = 1:length(timeVec)
    if timeVec(i) == 6 %% equivalent T6
        ac_sim(:,1) = ac_fba(i,:);
        ppa_sim(:,1) = ppa_fba(i,:);
        glc_sim(:,1) = glc_fba(i,:);
        lcts_sim(:,1) = lcts_fba(i,:);
    end
    if timeVec(i) == 12
        ac_sim(:,2) = ac_fba(i,:);
        ppa_sim(:,2) = ppa_fba(i,:);
        glc_sim(:,2) = glc_fba(i,:);
        lcts_sim(:,2) =  lcts_fba(i,:);
    end
end
end

ac_sim = reshape(ac_sim.',[],1);
ppa_sim = reshape(ppa_sim.',[],1);
glc_sim = reshape(glc_sim.',[],1);
lcts_sim =reshape(lcts_sim.',[],1);

time = data.G0.Biomass.timePoints(1:13);
biomass_sim = [];
[~,B] = size(b1_time_vect); %% check that both dfba have been performed. 
for i = 1:length(time(1:end-1))
    if B == 2
        loc = timeVec == time(i);
        idx = find(loc, 1, 'first');
        [b,~] = size(b1_time_vect);
        if idx <= b
            biomass_sim(i,:) = b1_time_vect(idx,:);
        else
            biomass_sim(i,:) = [0,0];
        end
    else
        biomass_sim(i,:) = [0,0];
    end
end 

biomass_sim = reshape(biomass_sim,[],1); 

calc_distance = norm((biomass_sim - data_exp).^2./(data_exp_stdev.^2)) + norm((ac_sim - ac).^2./(ac_stdev.^2)) + norm((ppa_sim - ppa).^2./(ppa_stdev.^2)) + norm((glc_sim - glc).^2./(glc_stdev.^2)) + norm((lcts_sim - lcts).^2./(lcts_stdev.^2));
    
    fid = fopen('opti_observation.txt', 'a+');
    fprintf(fid, '%.4f\n',calc_distance);
    fclose(fid);
end
