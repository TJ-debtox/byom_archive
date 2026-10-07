%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the model system. It is
% linked to the script files for Calanus2016, fit body residues with
% two-compartment TK. As input, it gets:
%
% * _t_   is the time point, provided by the ODE solver
% * _X_   is a vector with the previous value of the states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
%
% Time _t_ and scenario name _c_ are handed over as single numbers by
% <call_deri.html call_deri.m> (you do not have to use them in this
% function). Output _dX_ (as vector) provides the differentials for each
% state at _t_.
%
% * Author: Tjalling Jager
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives(t,X,par,c,glo)

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

CS = X(1); % state 1 is the scaled internal concentration in structure at previous time point
CL = X(2); % state 2 is the scaled internal concentration in lipids at previous time point

% these concentrations are both expressed on volume basis.

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

kS   = par.keS(1);    % elimination rate constant structure
PSW  = par.kuS(1)/kS; % bioconcentration factor structure-water
kL   = par.keL(1);    % elimination rate constant lipids
PLS  = par.kuL(1)/kL; % bioconcentration factor lipid-structure

%% Calculate the derivatives
% This is the actual model, specified as a system of two ODEs:

FLS = glo.FLS; % lipid content (VL/VS)

if t>4 && c==1.2 % depuration phase starts at t=4
    c=0;
end
% This is not such a great solution, but it works well enough. Better to
% break up the time vector and solve in two parts (this strategy is used in
% the TKTD packages for BYOM).

dCS = kS*(PSW*c-CS) - kL*FLS*(PLS*CS-CL); % first order bioconcentration
dCL = kL*(PLS*CS-CL); % first order bioconcentration

dX = [dCS;dCL;0]; % collect all derivatives in one vector dX
