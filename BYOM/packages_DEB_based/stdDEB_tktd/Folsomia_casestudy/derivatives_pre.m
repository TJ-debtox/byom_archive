%% BYOM function derivatives_pre.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the standard DEB model
% system, but only for the pre-calculations in call_deri: to find egg costs
% and initial values. This saves approx. 25% calculation time in a basic
% fit (controls only), relative to using the full derivatives function. The
% system is specified in terms of powers (energy fluxes in J/d) and body
% length. In this package, no toxicant stress is included. As input, it
% gets:
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
% * Date: April 2022
% * Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2023, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives_pre(t,X,par,c,glo)

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

E  = X(1); % state 1 is the reserve
EH = X(2); % state 2 is the maturity
L  = X(3); % state 3 is the volumetric length

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

% dV   = glo.dV;     % dry weight density of structure (g/cm3)
% delM = glo.delM;   % shape correction coefficient (-)
Lmat   = glo.Lmat;   % turn structure field into a regular matrix

TA   = par.TA(1);   % Arrhenius temperature
spAm = par.spAm(1); % max. surface-specific assimilation rate (J/(cm2 d))
spM  = par.spM(1);  % volume-specific somatic maintenance costs (J/(cm3 d))
spT  = par.spT(1);  % surface-specific maintenance costs (J/(cm2 d))
kJ   = par.kJ(1);   % maturity maintenance rate constant (1/d)
EG   = par.EG(1);   % volume-specific costs for growth (J/cm3)
EHb  = par.EHb(1);  % maturity level at birth (J)
EHj  = par.EHj(1);  % maturity level at metamorphosis (J)
EHp  = par.EHp(1);  % maturity level at puberty (J)
v    = par.v(1);    % energy conductance (cm/d)
kap  = par.kap(1);  % allocation fraction to soma (-)
f    = par.f(1);    % scaled food density (-)
kapH = 1; % maturation efficiency (-), this is only there to have a handle for toxic effects

%% Modify parameters where needed

% Temperature corrections with Arrhenius relationship
FT   = exp(TA/glo.Tref - TA/glo.T); % factor by which to modify rates/times
spAm = spAm * FT;   % modify specific assimilation
spM  = spM  * FT;   % modify volume-specific maintenance
spT  = spT  * FT;   % modify surface-specific maintenance
kJ   = kJ   * FT;   % modify maturity maintenance
v    = v    * FT;   % modify energy conductance

% NOTE: glo.Lmat will collect the length at the various life stage events
% (birth, metamorphosis, puberty and starting length for the simulation).
% It is initialised as vector of NaN's. The code is set up to catch these
% switches with the events function in call_deri. Code below will run
% through the model as if nothing happens, as long as the corresponding
% element of glo.Lmat is not defined. This makes sure that the events
% function can find the event without worrying about the possible
% discontinuities at the event. Note that call_deri stops and restarts at
% events.

% Acceleration with the abj model
if EHj > 0 && ~isnan(Lmat(1)) % then we have acceleration between birth and metamorphosis
    % Checking for Lb in glo.Lmat is to make sure that this piece of code
    % is skipped for the initial embryo simulations.
    del = 1; % start with no acceleration
    if EH > EHb % after birth (this should be superfluous)
        Lb = Lmat(1); % length at birth; this field is now set by call_deri
        if EH < EHj || isnan(Lmat(2))
            % Checking for Lj is to make sure that this piece of code is
            % not run when *looking* for the metamorphosis event. This
            % saves the ODE solver from dealing with a switch when it is
            % not needed.
            del = L/Lb; % accelerate!
        else % we have already reached metamorphosis
            Lj  = Lmat(2); % length at metamorphosis; this field is now set by call_deri
            del = Lj/Lb; % acceleration factor at metamorphosis and onwards
        end
    end
    v    = v    * del; % change reserve mobilisation
    spAm = spAm * del; % change surface-specific assimilation rate
    % Note: also change surface-specific maintance costs?
end

% Make sure embryos have no 'food'
if EH < EHb || isnan(Lmat(1))
    f = 0; % for embryos, there is never any food
    % Checking for Lb is to make sure that this piece of code is not run
    % when *looking* for the birth event.
end

%% Calculate the derivatives
% This is the actual model, specified as a system of ODEs. Note: flux pJ is
% calculated under starvation, but that value is then not used. It may be
% needed for future modules (e.g., respiration and ageing).
    
% Calculate the powers
pA = f * spAm * L^2;        % assimilation power
pS = spM * L^3 + spT * L^2; % maintenance power (somatic and surface-specific)
pC = E * (EG * v * L^2 + pS) / (kap * E + EG * L^3); % mobilisation power
pG = kap * pC - pS;         % growth power
pJ = kJ * EH;               % maturity maintenance power
pR = (1-kap) * pC - pJ;     % maturation/reproduction power

% No starvation module for embryo/control simulations! An events function
% is in place that stops the simulation as soon as l=e (which means that
% starvation is imminent).

% State variables
dE = pA - pC;               % change in reserve
dL = (1/(3*L^2)) * (pG/EG); % change in volumetric length
dL = max(0,dL);             % do not allow shrinking
if EH < EHp || isnan(Lmat(3)) % for embryos and juveniles ...
    % Checking for Lp is to make sure that this piece of code is not run
    % when *looking* for the puberty event. This saves the ODE solver from
    % dealing with a switch when it is not needed. Note that call_deri
    % stops at puberty.
    dEH = max(0,kapH * pR); % maturation (don't allow to become negative)
else                        % for adults ...
    dEH = 0;                % no maturation
end
% NOTE: I do not allow dL and dEH to decrease. When dL starts to decrease,
% the embryo will never reach birth, since maturity will then not increase
% either. This is caught by the events function as scaled reserve e will
% decrease below scaled length l. For getting an initial length L=L0, there
% will be no starvation, as that animal will (be assumed to) experience
% constant f.

dX = [dE;dEH;dL]; % collect all derivatives in one vector
