%% BYOM function simplefun.m (the model as explicit equations)
%
%  Syntax: Xout = simplefun(t,X0,par,c,glo)
%
% This function calculates the output of the model system. It is linked to
% the script files named byom_doseresp_*.m. Therefore, _t_ is used for
% concentrations and _c_ for time! (Note: BYOM normally works with 'time'
% on the x-axis). As input, it gets:
%
% * _t_   is the vector with exposure concentrations
% * _X0_  is a vector with the initial values for states (not used)
% * _par_ is the parameter structure
% * _c_   is the exposure time at which the dose-response is made
% * _glo_ is the structure with information (normally global)
%
% Variable _t_ is handed over as a vector, and scenario name _c_ as single
% number, by <call_deri.html call_deri.m> (you do not have to use them in
% this function). Output _Xout_ (as matrix) provides the output for each
% state at each _t_.
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
function Xout = simplefun(t,X0,par,c,glo)

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

ECx  = par.ECx(1);  % concentration for x% effect (x in glo.x_EC)
Y0   = par.Y0(1);   % response in control
beta = par.beta(1); % slope factor of the dose response

%% Calculate the model output
% This is the actual model, the log-logistic dose-response curve, specified
% as explicit function.

x = glo.x_EC; % the effect level (as percentage)
y = Y0 ./ (1+(x/(100-x))*(t/ECx).^beta); % log-logistic dose-response function

Xout = [y]; % combine all outputs into a matrix