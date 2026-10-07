%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the GUTS model system with
% extra states for immobility. The implementation of the full model is used
% as the basis here, so there is TK and damage dynamics, which are
% calculated separately. For IT, only internal concentration and damage are
% calculated by this function. For SD, also the hazard rates and state
% probabilities (incl. background hazard). As input, it gets:
%
% * _t_   is the time point, provided by the ODE solver
% * _X_   is a vector with the previous value of the states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information
%
% Time _t_ and scenario name _c_ are handed over as single numbers by
% <call_deri.html call_deri.m> (you do not have to use them in this
% function). Output _dX_ (as vector) provides the differentials for each
% state at _t_.
%
% * Author: Tjalling Jager
% * Date: January 2024
% * Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2024, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives(t,X,par,c,glo)

% NOTE: glo is now no longer a global here, but passed on in the function
% call from call_deri. That saves calculation time!

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

Ci = X(glo.locC);  % internal concentration at previous time point 
Di = X(glo.locD);  % scaled damage at previous time point
pA = X(glo.loc_a); % active probability at previous time point
pI = X(glo.loc_i); % immobile probability at previous time point
% State glo.loc_d is for the mortality probability at previous time point (not used in the ODEs)
% State glo.loc_id is for the sum of dead and immobile, calculated in call_deri (not used in the ODEs)

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

ke   = par.ke(1);   % elimination rate constant
kr   = par.kr(1);   % damage repair rate constant
Kiw  = par.Kiw(1);  % bioconcentration factor (set to 1 for 'scaled')
hb   = par.hb(1);   % background hazard rate

if glo.sel == 1 % for SD, need more parameters as health states are calculated below
    bii  = par.bii(1); % killing rate immobility
    bid  = par.bid(1); % killing rate for death
    bir  = par.bir(1); % killing rate for recovery    
    mii  = par.mii(1); % threshold for immobility
    mid  = par.mid(1); % threshold for death
end
% Note: when there is only one threshold used in the analysis, it is copied
% to the two parameters mii and mid in call_deri.

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

% % include the option for a Michaelis-Menten type saturation of exposure
% if isfield(par,'cK')
%     cK   = par.cK(1); % half-saturation concentration
%     if cK > 0 % if we want saturation ...
%         c = c * (cK/(c+cK)); % modify exposure concentrations
%     end
% end

%% Calculate the derivatives
% This is the actual model, specified as a system of ODEs.

if glo.fastslow(1) == 1 % for fast toxicokinetics
    dCi = 0; % no change in internal concentration
    Ci  = c; % internal concentration equals external
else
    dCi = ke * (Kiw * c - Ci); % first order bioconcentration
    % Note: when Kiw = 1, we can view this as scaled internal concentration
end

if glo.fastslow(2) == 1 % if we assume that damage repair is infinitely fast
    dDi = 0;  % no change in scaled damage
    if glo.damconfig == 3 % if switch is set to 3 ...
        Di  = c; % scaled damage equals EXternal concentration
    else
        Di  = Ci; % scaled damage equals INternal concentration
    end
else
    if glo.damconfig == 3 % if switch is set to 3 ...
        dDi = kr * (c - Di); % first order damage build-up from c (EXternal!)
    else
        dDi = kr * (Ci - Di); % first order damage build-up from Ci (scaled INternal)
    end
end

% We need rules to select which state variable (Ci or Di) to use for which
% effect (immobility or death).
Did = Di; % always use damage for death mechanism
if glo.damconfig ~= 2 % if switch is set to 1 or 3 ...
    Dii = Ci;  % use internal concentration as damage for immobilisation
else           % if it is set to 2 ...
    Dii = Did; % use damage Di for immobilisation as well
end

if glo.sel == 1 % for SD, need health states calculated here
    
    % hazard rates for all transitions
    hi = bii * max(0,Dii-mii); % active to immobile
    hd = bid * max(0,Did-mid); % immobile (and active) to dead
    hr = bir * max(0,mii-Dii); % immobile to active
    
    hi  = min(111,hi); % maximise the hazard rate to 99% immobility in 1 hour
    hd  = min(111,hd); % maximise the hazard rate to 99% mortality in 1 hour
    hr  = min(111,hr); % maximise the hazard rate to 99% recovery in 1 hour

    if glo.onehit == 0
        % Calculate changes to probability to be in each state. This
        % assumes a two-step mechanism: only immobilised animals will die
        % (2-hit).
        dpH = -(hi+hb)*pA + hr*pI;
        dpI = hi*pA - (hd+hb+hr)*pI;
        dpD = hb*pA + (hb+hd)*pI;
    else    
        % Otherwise, use two independent mechanisms (active
        % individuals can also die from the chemical: 1-hit).
        dpH = -(hi+hb+hd)*pA + hr*pI;
        dpI = hi*pA - (hd+hb+hr)*pI;
        dpD = (hb+hd)*pA + (hb+hd)*pI;
    end
    
else % for IT, the death-immobile-active states are calculated in call_deri
    
    dpH = 0;
    dpI = 0;
    dpD = 0;
    
end

dX = zeros(size(X)); % initialise with zeros in right format
dX([glo.locC glo.locD glo.loc_a glo.loc_i glo.loc_d glo.loc_id]) = [dCi;dDi;dpH;dpI;dpD;0]; % collect derivatives in one vector


