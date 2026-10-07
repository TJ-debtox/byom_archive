%% BYOM function startvals.m
%
%  Syntax: [Y0,ECx,logc] = startvals(DATA,W)
%
% This function calculates starting values for fitting a dose-response
% curve, and ranges for plotting model curves, from the data set entered.
% 
% As input, it gets:
%
% * _DATA_  is the data matrix
% * _W_     the weight matrix (used for survival data only)
% 
% Output is
% * _Y0_    initial value for control response
% * _ECx_   initial value for ECx (uses the EC50)
% * _logc_  range for plotting the model curve on log-scale
%
% * Author: Tjalling Jager 
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_doseresp.html>
% 
%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start
function [Y0,ECx,logc] = startvals(DATA,W)

if DATA(1,1) < 0 % then we have survival data
    mn_resp = mean(DATA(2:end,2:end)./W,2,'omitnan'); % calculate mean response over all replicates (ignore NaNs)
    % Note: recalculate the mean response to a survival probability
    Y0      = mn_resp(1); % initial value for control response is mean of control (as probability)
else
    mn_resp = mean(DATA(2:end,2:end),2,'omitnan'); % calculate mean response over all replicates (ignore NaNs)
    Y0      = mn_resp(1); % initial value for control response (before scaling)
    mn_resp = mn_resp / mn_resp(1); % scale the mean response to that in the control treatment
end

ind_1   = find(mn_resp > 0.5,1,'last'); % find highest concentration that has less than 50% effect 
ind_2   = find(mn_resp < 0.5,1,'first'); % find lowest concentration that has more than 50% effect 
c       = DATA(2:end,1); % take concentrations from data set (first column)

if c(ind_1) == 0 % catch cases where there is more than 50% effect in the first treatment
    % and the first concentration is zero
    ECx = c(ind_2); % just use the second one (first treatment with >50% effect)
else
    ECx = 10.^(mean(log10(c([ind_1 ind_2])))); % initial value for ECx around 50% effect (geometric mean)
end

logc    = [log10(c(2)) - 2.4*mean(diff(log10(c(c>0)))) log10(1.2*c(end))]; % useful range for plotting the model curve on log-scale
