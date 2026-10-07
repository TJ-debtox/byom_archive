%% BYOM function simplefun.m (the model as explicit equations)
%
%  Syntax: Xout = simplefun(t,X0,par,c,glo)
%
% This function calculates the output of the model system. This is simple
% one-compartment TK. As input, it gets:
%
% * _t_   is the time vector
% * _X0_  is a vector with the initial values for states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
%
% Time _t_ is handed over as a vector, and scenario name _c_ as single
% number, by <call_deri.html call_deri.m> (you do not have to use them in
% this function). Output _Xout_ (as matrix) provides the output for each
% state at each _t_.
%
% * Author: Tjalling Jager 
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function Xout = simplefun(t,X0,par,c,glo)

%% Unpack initial states
% The state variables enter this function in the vector _X_0.

% Ci0 = X0(1); % state 1 is the internal concentration at t=0

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

ke   = par.ke(1);     % elimination rate constant, d-1
Piw  = par.ku(1)/ke;  % bioconcentration factor, L/kg

%% Calculate the model output
% This is the actual model, specified as explicit function(s). Here, we
% assume that depuration phase starts at t=4.

Ci = c*Piw*(1 - exp(-ke*t)); % internal concentration, accumulation

Ci(t>4) = c*Piw*(1 - exp(-ke*4)) * exp(-ke*(t(t>4)-4)); % modify for depuration phase

Xout = [Ci]; % combine them into a matrix