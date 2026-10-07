%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the model system. It is
% linked to the script file <byom_debtox_capitella.html
% byom_debtox_capitella.m>. As input, it gets:
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
% * Date: September 2020
% * Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives(t,X,par,c,glo)

% NOTE: glo is now no longer a global here, but passed on in the function
% call from call_deri. That saves 20% calculation time!

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

cV = X(1); % state 1 is the scaled internal concentration
e  = X(2); % state 2 is the scaled reserve density
L  = X(3); % state 3 is body length
% Rc = X(4); % state 4 is cumulative reproduction

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

ke = par.ke(1);  % elimination rate constant, d-1
c0 = par.c0(1);  % no-effect concentration
cT = par.cT(1);  % tolerance concentration

% L0 = par.L0(1);  % initial body length (mm) (only used in call_deri.m)
Lp = par.Lp(1);  % length at puberty (mm)
Lm = par.Lm(1);  % maximum length (mm)
rB = par.rB(1);  % von Bertalanffy growth rate
Rm = par.Rm(1);  % maximum reproduction rate
g  = par.g(1);   % energy-investment ratio
f  = par.f(1);   % scaled functional response

% modifications for the Capitella data set
Tlag = par.Tlag(1); % lag time for growth to start
f0   = par.f0(1);   % scaled functional response control
if c == 0 % if we are dealing with the control ...
    f = f0; % for control, use lower value for f!
end

%% Extract correct exposure for THIS time point
% Allow for external concentrations to change over time, either
% continuously, or in steps, or as a static renewal with first-order
% disappearance. For constant exposure, the code in this section is skipped
% (and could also be removed). Note that glo.timevar is used to signal that
% there is a time-varying concentration. This option is set in call_deri.

if glo.timevar(1) == 1 % if we are warned that we have a time-varying concentration ...
    c = read_scen(-1,c,t,glo); % use read_scen to derive actual exposure concentration
    % Note: this has glo as input to the function to save time!
    % the -1 lets read_scen know we are calling from derivatives (so need one c)
end

%% Calculate the derivatives
% This is the actual model, specified as a system of ODEs:

kM = rB * 3 * (1+g)/g;    % maintenance rate coefficient
v  = Lm * kM * g;         % energy conductance
s  = (1/cT)*max(0,cV-c0); % stress factor

switch glo.moa % stress effect according to selected mechanism
    case 1 % effect on assimilation/feeding
        f = f*max(0,1-s);
    case 2 % effect on maintenance costs
        kM = kM * (1+s);
        Rm = Rm * (1+s);
    case 3 % effect on growth costs
        g  = g * (1+s);
        kM = kM /(1+s);
    case 4 % effect on repro costs
        Rm = Rm / (1+s);
    case 5 % hazard to early embryo
        Rm = Rm * exp(-s);
    case 6 % costs for growth and repro
        g  = g * (1+s);
        kM = kM /(1+s);
        Rm = Rm / (1+s); 
end

de = (f-e)*v/L; % change in scaled energy density
dL = (kM*g/(3*(e+g))) * (e*v/(kM*g) - L); % change in length
if L<Lp % before length at puberty is reached ...
    R = 0; % no reproduction
else
    R = (Rm/(Lm^3-Lp^3))*((v*(L^2)/kM + L^3)*e/(e+g) - Lp^3); % repro rate
    R = max(0,R); % make sure it is always positive
end
dRc = R; % change in cumulative repro is the repro rate

% OPTION 1: DEBtox default, growth dilution and surface:volume ratio
dcV = ke*(Lm/L)*(c-cV) - cV*(3/L)*dL; % change in scaled internal conc.
% % OPTION 2: growth dilution only
% dcV = ke*(c-cV) - cV*(3/L)*dL; % change in scaled internal conc.
% % OPTION 3: no growth dilution and no surface:volume ratio
% dcV = ke*(c-cV); % change in scaled internal conc.

if t < Tlag % before lag time, no change in states ...
    dX = [0;0;0;0]; % collect all derivatives in one vector dX
else
    dX = [dcV;de;dL;dRc]; % collect all derivatives in one vector dX
end
