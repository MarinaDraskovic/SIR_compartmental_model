% SIR model China.
% Cost function - Sum of the squared errors on active and recovered cases.
% Minimized with PSO.
% Added equation y=k*I
% where k=population size*fraction of reported/observed cases.
% Equations of SIR model with fraction of population (not number of people)
% Estimated beta (transmission parameter), gamma (recovery rate), k.
% https://github.com/epimath/param-estimation-SIR/blob/master/EisenbergIdentifiabilityLab.pdf
% Returns [beta gamma k], test, yest
%______________________________________________
function [paramests, test, yest] = SIR_model_paramest(active_cases_data, recovered_cases_data, times)
    % A function that returns S(0), I(0), R(0)
    % S(0) = 1- I(0), I(0) = data(1)/k, R(0) = 0. 
    % Parameters: beta, gamma, k.
    x0_fun = @(params) [1-active_cases_data(1)/params(3); active_cases_data(1)/params(3); 0];

    % SIR model function (No N because we are working with fractions.)
    fun_SIR = @(t, y, params)[-params(1)*y(1)*y(2); (params(1)*y(1)*y(2))-(params(2)*y(2)); params(2)*y(2)];

    % Parameter estimation
    lb = [0,0];
    ub = [0.5,0.5];
    [paramests, ~] = particleswarm(@(p) costFun(fun_SIR, times, abs(p), active_cases_data, recovered_cases_data, x0_fun), 3, lb, ub);
    paramests = abs(paramests);

    % Solve SIR model
    [test,yest] = ode45(fun_SIR, times, x0_fun(paramests), [], paramests);
    
end

%% Cost function
function ret = costFun(fun_SIR, times, params, adata, rdata, x0_fun)
    [~,y] = ode45(fun_SIR, times, x0_fun(params), [], params);
    data=[adata; rdata];
    y=[(y(:,2)*params(3)); (y(:,3)*params(3))];
    ret = sum((data - y).^2);
end

