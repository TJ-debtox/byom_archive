%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the model system. It is
% linked to the script files <byom_debkiss_compound.html
% byom_debkiss_compound.m>. This is the DEBkiss model in compound
% parameters. As input, it gets:
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
% * Back to index <walkthrough_debkiss.html>

%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives(t,X,par,c,glo)

% NOTE: glo is now no longer a global here, but passed on in the function
% call from call_deri. That saves calculation time!

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

L  = X(1); % state is body length
% Rc = X(2); % state is cumulative reproduction (not used)

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

% Unpack globals
kap  = glo.kap;     % approximation for kappa (-)
yP   = glo.yP;      % product of yVA and yAV (-)

% Unpack model parameters for the basic life history
L0   = par.L0(1);   % body length at start (mm)
Lp   = par.Lp(1);   % body length at puberty (mm)
Lm   = par.Lm(1);   % maximum body length (mm)
rB   = par.rB(1);   % von Bertalanffy growth rate constant (1/d)
Rm   = par.Rm(1);   % maximum reproduction rate (#/d)
f    = par.f(1);    % scaled functional response (-)

% Unpack extra parameters for specific cases
Lf   = par.Lf(1);   % body length at half-saturation feeding (mm)
Lj   = par.Lj(1);   % body length at which acceleration stops (mm)
Tlag = par.Tlag(1); % lag time for start development (d)

%% Calculate the derivatives
% This is the actual model, specified as a system of ODEs. This is the
% DEBkiss model, with toxicant effects and starvation module, expressed in
% terms of compound parameters, as presented in the publication (Jager,
% 2020).

L = max(1e-3*L0,L); % make sure that body length is not negative or almost zero (extreme shrinking may do that)
% This should not be needed as shrinking is limited at the bottom of this
% function.

if Lf > 0 % to include feeding limitation for juveniles ...
    f  = f / (1+(Lf^3)/(L^3)); % hyperbolic relationship for f with body volume
end
if Lj > 0 % to include acceleration until metamorphosis ...
    f = f * min(1,L/Lj);
end

% Calcululate the actual derivatives
dL = rB * (f*Lm - L); % ODE for body length

fR = f; % if there is no starvation, f for reproduction is the standard f
% starvation rules can modify the outputs here
if dL < 0 % then we are looking at starvation and need to correct things
    fR = (f - kap * (L/Lm))/(1-kap); % new f for reproduction alone
    if fR >= 0  % then we are in the first stage of starvation: 1-kappa branch can help pay maintenance
        dL = 0; % stop growth, but don't shrink
    else        % we are in stage 2 of starvation and need to shrink to pay maintenance
        fR = 0; % nothing left for reproduction
        dL = (rB/yP) * ((f*Lm/kap) - L); % shrinking rate
    end
end
        
R  = 0; % reproduction rate is zero, unless ... 
if L >= Lp % if we are above the length at puberty, reproduce
    R = max(0,Rm * (fR*Lm*(L^2) - Lp^3)/(Lm^3 - Lp^3));
end

dRc = R; % cumulative reproduction rate

if L <= 0.5 * L0 % if an animal has size less than half the start size ...
    dL = 0; % don't let it grow or shrink any further (to avoid numerical issues)
end

dX = zeros(size(X)); % initialise with zeros in right format
if t >= Tlag % when we are past the lag time ...
    dX = [dL;dRc]; % collect all derivatives in one vector dX
end
