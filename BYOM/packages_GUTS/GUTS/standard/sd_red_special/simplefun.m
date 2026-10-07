%% BYOM function simplefun.m (the model as explicit equations)
%
%  Syntax: Xout = simplefun(t,X0,par,c,glo)
%
% This function calculates the output of the reduced GUTS-SD model system
% (including survival). As input, it gets:
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
% * Date: June 2022
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_guts.html>

%  Copyright (c) 2012-2023, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function [Xout,haz] = simplefun(t,X0,par,c,glo)

% The calculation of survival in this file only works when the exposure
% concentration is constant over time! For time-varying concentrations, use
% the GUTS version in the 'reduced' folder or the ODE version in the
% 'experimental' directory.

%% Unpack initial states
% The state variables enter this function in the vector _X_0. However, the
% analytical solution is (for now) limited to the situation where the
% initial scaled damage level is zero.

% S0  = X0(glo.locS); % survival probability at t=0
% Dw0 = X0(glo.locD); % scaled damage (referenced to water) at t=0

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

kd   = par.kd(1);   % dominant rate constant
mw   = par.mw(1);   % median threshold
bw   = par.bw(1);   % killing rate
hb   = par.hb(1);   % background hazard rate

%% Calculate the model output
% This is the actual model, specified as explicit function(s):

S  = exp(-hb*t);         % initialise S with background mortality for survival
Dw = (1-exp(-kd*t)) * c; % scaled damage over time
% Note: the initial damage at t=0 is ALWAYS assumed to be zero.

% This is basically the same calculation as in DEBtool fomort
if c > mw % if concentration to test is above the threshold ...
    t0 = -log(1-mw/c)/kd;   % no-effect-time: point where effects start
    f  = (bw/kd)*max(0,exp(-kd*t0) - exp(-kd*t))*c - bw*(max(0,c-mw))*max(0,t-t0);
    S  = min(1,exp(f)) .* S; % multiply with background mortality
end
% This should work without errors: if c<=mw then there is no effect anyway,
% apart from background mortality. This avoids the awkward 1e-8 of the
% original below. Since we only work with one concentration at a time,
% there is also no need for t1.

% t0 = -log(max(1e-8,1-mw/max(1e-8,c)))/kd; % no-effect-time
% t1 = ones(length(t),1);                   % column-matrix of ones
% f  = (bw/kd)*max(0,t1*exp(-kd*t0) - exp(-kd*t)).*(t1*c) - ...
%     bw*(t1*(max(0,c-mw)')).*max(0,t - t1*t0);
% 
% S = min(1,exp(f)) .* S; % multiply with background mortality
   
Xout(:,[glo.locS glo.locD]) = [S Dw]; % combine them into a matrix

haz = bw * max(0,Dw-mw); % calculate hazard for each time point