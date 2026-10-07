%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the model system. As input,
% it gets:
%
% * _t_   is the time point, provided by the ODE solver
% * _X_   is a vector with the previous value of the states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
%
% Time _t_ and scenario name _c_ are handed over as single numbers by
% call_deri.m (you do not have to use them in this function). Output _dX_
% (as vector) provides the differentials for each state at _t_.
%
% * Author: Tjalling Jager (email: tjalling_at_debtox.info) 
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives(t,X,par,c,glo)

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

Ci = X(1); % state 1 is the internal concentration at previous time point

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

ke   = par.ke(1);  % elimination rate constant, d-1
ku   = par.ku(1);  % uptake rate constant factor, L/kg/d

%% Calculate the derivatives
% This is the actual model, specified as ODE:
%
% $$ \frac{d}{dt}C_i=k_u C_w- k_e C_i $$
%
% Note: forcing a sudden stop of exposure, as done below, is not the best
% solution. Such a discontinuity is difficult for the ODE solver. However,
% here it works well enough. In general, breaking up the time vector, and
% solving them piecewise is better (this is done in various BYOM packages
% for TKTD analysis).

if t>7 % for the depuration phase ...
    c = 0; % make the external concentration zero
end
dCi = ku * c - ke * Ci; % first order bioconcentration

dX = [dCi]; % collect derivatives in one vector dX
