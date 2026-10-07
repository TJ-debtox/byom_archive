%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the full GUTS model system.
% Note that the survival probability dues to chemical stress is all
% calculated in <call_deri.html call_deri.m>. As input, it gets:
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
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_guts.html>

%  Copyright (c) 2012-2023, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives(t,X,par,c,glo)

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

S  = X(1); % state is the survival probability at previous time point
Ci = X(2); % state is the internal concentration at previous time point
Di = X(3); % state is the scaled damage at previous time point

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

ke   = par.ke(1);   % elimination rate constant
kr   = par.kr(1);   % damage repair rate constant
Kiw  = par.Kiw(1);  % bioconcentration factor
% mi   = par.mi(1);   % median of threshold distribution (used in call_deri)
% bi   = par.bi(1);   % killing rate (used in call_deri)
% Fs   = par.Fs(1);   % fraction spread of threshold distribution, (-) (used in call_deri)
hb   = par.hb(1);   % background hazard rate

%% Extract correct exposure for THIS time point
% Allow for external concentrations to change over time, either
% continuously, or in steps, or as a static renewal with first-order
% disappearance. For constant exposure, the code in this section is skipped
% (and could also be removed).

if glo.timevar(1) == 1 % if we are warned that we have a time-varying concentration ...
    c = read_scen(-1,c,t,glo); % use read_scen to derive actual exposure concentration
    % Note: this has glo as input to the function to save time!
    % the -1 lets read_scen know we are calling from derivatives (so need one c)
end

%% Calculate the derivatives
% This is the actual model, specified as a system of ODEs.

dCi = ke * (Kiw * c - Ci); % first order bioconcentration
if glo.fastrep == 1 % is we assume that damage repair is infinitely fast
    dDi = dCi; % same change in damage as in internal conc.
else
    dDi = kr * (Ci - Di); % first order damage build-up from Ci (scaled)
end

dS = -hb * S; % only background hazard rate
% mortality due to the chemical is included in call_deri!

dX = [dS;dCi;dDi]; % collect derivatives in one vector


