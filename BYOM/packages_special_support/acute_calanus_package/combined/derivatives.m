%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the model system. It is
% linked to the script files for Calanus2016, GUTS combined with
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

%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives(t,X,par,c,glo)

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

S  = X(1); % state 1 is the survival probability at previous time point
CS = X(2); % state 2 is the scaled internal concentration in structure at previous time point
CL = X(3); % state 3 is the scaled internal concentration in lipids at previous time point

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

keS  = par.keS(1);     % elimination rate constant structure
PSW  = par.kuS(1)/keS; % bioconcentration factor structure-water
keL  = par.keL(1);     % elimination rate constant lipids
PLS  = par.kuL(1)/keL; % bioconcentration factor lipid-structure

% mS   = par.mS(1);   % no-effect concentration, ref. to structure (used in call_deri)
% bS   = par.bS(1);   % killing rate, ref. to structure (used in call_deri)
% Fs   = par.Fs(1);   % fraction spread of NEC distribution, (-) (used in call_deri)
hb   = par.hb(1);   % background hazard rate

FLS  = glo.FLS; % VL/VS
Tdep = glo.Tdep; % start of depuration

%% Calculate the derivatives
% This is the actual model, specified as a system of ODEs.

if t>Tdep && c==1.2 % start depuration at t=4 (for fitting)
    c=0;
end
% This is not such a great solution, but it works well enough. Better to
% break up the time vector and solve in two parts (this strategy is used in
% the TKTD packages for BYOM).

dCS = keS*(PSW*c-CS) - keL*FLS*(PLS*CS-CL); % first order bioconcentration, structure
dCL = keL*(PLS*CS-CL); % first order bioconcentration, lipids

dS = -hb * S; % only background hazard rate
% mortality due to the chemical is included in call_deri!

dX = [dS;dCS;dCL;0]; % collect derivatives in one vector (last one is a dummy)
