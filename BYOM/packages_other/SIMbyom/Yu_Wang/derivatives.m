%% SIMbyom function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the model system. It is
% linked to the script file byom_sim_yuwang. As input, it gets:
%
% * _t_   is the time point, provided by the ODE solver
% * _X_   is a vector with the previous value of the states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
%
% Time _t_ and scenario name _c_ are handed over as single numbers by
% call_deri (you do not have to use them in this function). Output _dX_ (as
% vector) provides the differentials for each state at _t_.
%
% * Author: Tjalling Jager
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_simbyom.html>
% 
%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM.

%% Start

function dX = derivatives(t,X,par,c,glo)

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

x = X(1); % state 1 
y = X(2); % state 2 
z = X(3); % state 3

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values (probably not needed for the SIMbyom package).

alp   = par.alp(1);   % alpha
bet   = par.bet(1);   % beta
gam   = par.gam(1);   % gamma
del   = par.del(1);   % delta

%% Calculate the derivatives
% This is the actual model, specified as ODE

dx = alp*(y-x);
dy = bet*x-gam.*x*z;
dz = exp(x*y)-del*z;

dX = [dx;dy;dz]; % collect all derivatives in one vector dX
