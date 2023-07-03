function x = objective_func(bounds,data,model_list)
% bounds contains the parameters to optimize, 
% data contains the growth and metabolites data 
% model_list contains the list of the models
    x = observation_func(bounds,data,model_list);
end
