%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the standard DEB model
% system with TKTD module from DEBtox 2019. The system is specified in
% terms of powers (energy fluxes in J/d) and body length. In this package,
% no toxicant stress is included. As input, it gets:
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
% * Date: February 2023
% * Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2023, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives(t,X,par,c,glo)

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

E  = X(1); % state 1 is the reserve
EH = X(2); % state 2 is the maturity
L  = X(3); % state 3 is the volumetric length
% Rc = X(4); % state 4 is the cumulative reproduction (not used in this function)
Dw = X(5); % state 5 is the scaled damage (referenced to water)
S  = X(6); % state 6 is survival probability

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

% dV    = glo.dV;     % dry weight density of structure (g/cm3)
delM    = glo.delM;   % shape correction coefficient (-)
yP      = glo.yP;     % product of yVE and yEV (-)
Lmat    = glo.Lmat;   % turn structure field into a regular matrix
E0mat   = glo.E0mat;  % turn structure field into a regular matrix
FBV     = glo.FBV;    % dry weight of egg as fraction of dry body weight (-)
KRV     = glo.KRV;    % part. coeff. repro buffer and structure (kg/kg)

% glo.Lwm_ref is the reference maximum physical length for scaling rate constants
Lm_ref  = glo.Lwm_ref * delM; % value in glo is *physical length* (need *volumetric* here)

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
kapR = par.kapR(1); % reproduction efficiency (-)
f    = par.f(1);    % scaled food density (-)
hb   = par.hb(1);   % background hazard rate (d-1)
a    = par.a(1);    % Weibull background hazard coefficient (-)
kapH = 1; % maturation efficiency (-), this is only there to have a handle for toxic effects

% unpack model parameters for the response to toxicants
kd   = par.kd(1);   % dominant rate constant (d-1)
zb   = par.zb(1);   % effect threshold energy budget ([C])
bb   = par.bb(1);   % effect strength energy-budget effects (1/[C])
zs   = par.zs(1);   % effect threshold survival ([C])
bs   = par.bs(1);   % effect strength survival (1/([C] d))

f_rem = f; % remember the input f for egg costs, before we apply starvation

hb    = a * (hb^a) * t^(a-1); % option for Weibull mortalty when a is not 1

%% Extract correct exposure for THIS time point
% Allow for external concentrations to change over time, either
% continuously, or in steps, or as a static renewal with first-order
% disappearance. For constant exposure, the code in this section is skipped
% (and could also be removed). Note that glo.timevar is used to signal that
% there is a time-varying concentration. This option is set in call_deri.

E0_calc = glo.E0_calc(1,:); % use first row by default
if glo.timevar(1) == 1 % if we are warned that we have a time-varying concentration ...
    if size(glo.E0_calc,1)>1            % if multiple rows have been defined
        i_d = floor(c/100);             % extract the data set number from the treatment identifier
        E0_calc = glo.E0_calc(i_d+1,:); % use row for this study
    end
    c = read_scen(-1,c,t,glo); % use read_scen to derive actual exposure concentration
    % Note: this has glo as input to the function to save time!
    % the -1 lets read_scen know we are calling from derivatives (so need one c)
end

%% Modify parameters where needed

% Temperature corrections with Arrhenius relationship
FT   = exp(TA/glo.Tref - TA/glo.T); % factor by which to modify rates/times
spAm = spAm * FT;   % modify specific assimilation
spM  = spM  * FT;   % modify volume-specific maintenance
spT  = spT  * FT;   % modify surface-specific maintenance
kJ   = kJ   * FT;   % modify maturity maintenance
v    = v    * FT;   % modify energy conductance

% Calculate stress factor and hazard rate
s = bb * max(0,Dw-zb); % stress level for metabolic effects
h = bs * max(0,Dw-zs); % hazard rate for effects on survival

h = min(111,h); % maximise the hazard rate to 99% mortality in 1 hour
% Note: this helps in extreme conditions, as the system becomes stiff for
% very high hazard rates. This is especially needed for EPx calculations,
% where the MF is increased until there is effect on all endpoints!

% Define specific stress factors s*, depending on the mode of action as
% specified in the vector with switches glo.moa.
Si = glo.moa * s;  % vector with specific stress factors from switches for mode of action
sA = min(1,Si(1)); % assimilation/feeding (maximise to 1 to avoid negative values for 1-sA)
sM = Si(2);        % maintenance (somatic and maturity)
sG = Si(3);        % growth and maturation costs
sR = Si(4);        % reproduction costs  
sH = Si(5);        % also include hazard to reproduction

% Apply the specific stress factors s* to the primary parameters
spAm = spAm * (1-sA); % effect on assimilation, alternative: /(1+s) increase in handling time
spM  = spM * (1+sM);  % effect on somatic maintenance
kJ   = kJ * (1+sM);   % effect on maturity maintenance
EG   = EG * (1+sG);   % effect on growth costs
kapH = kapH / (1+sG); % effect on maturation costs
kapR = exp(-sH) * kapR / (1+sR); % add effect on repro hazards and costs
% NOTE: the repro costs can also be applied to E0mat (first column).
% However, perhaps it is better to think of this pMoA as an increase in the
% overhead costs of egg production, rather than an increase in initial
% reserve per egg.

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

% Starvation rules may override these fluxes
if pG < 0        % then we have starvation
    if pC > pS + pJ          % mobilisation is still enough to pay somatic and maturity maintenance
        pG = 0;              % stop growth
        pR = pC - pS - pJ;   % maturation/reproduction gets what's left
    elseif pC >= pS          % mobilisation can only pay for somatic maintenance costs
        pG = 0;              % stop growth
        pR = 0;              % stop maturation/reproduction
        % pJ = pC - pS;        % maturity maintenance gets what's left
    else                     % mobilisation cannot even pay for somatic maintenance
        % pJ = 0;              % stop maturity maintenance
        pR = 0;              % stop maturation/reproduction
        pG = (pC - pS) / yP; % shrinking of structure to pay somatic maintenance
    end
end
% NOTE: in the latter two cases, maturity maintenance is not paid. We may
% at some point want to include 'rejuvenation' under these conditions.

% State variables
dE = pA - pC;               % change in reserve
dL = (1/(3*L^2)) * (pG/EG); % change in volumetric length
if EH < EHp || isnan(Lmat(3)) % for embryos and juveniles ...
    % Checking for Lp is to make sure that this piece of code is not run
    % when *looking* for the puberty event. This saves the ODE solver from
    % dealing with a switch when it is not needed. Note that call_deri
    % stops at puberty.
    dEH = max(0,kapH * pR); % maturation (don't allow to become negative)
    R   = 0;                % no reproduction
else                        % for adults ...
    dEH = 0;                % no maturation
    switch E0_calc(2) 
        case 1
            e = f_rem; % base egg costs on the food density only (then it remains constant under stress, but differs between data sets with different f)
        case 2
            e = 1; % base egg costs on optimal food density only (then it remains constant, always)
        otherwise
            e = E * v / (spAm * L^3); % base egg costs on actual reserve density so it changes with stress of mother
    end
    E0  = (interp1(E0mat(:,2),E0mat(:,1),e,'linear','extrap')); % derive egg costs for the current e
    R   = max(0,(kapR / E0) * pR); % reproduction rate in numbers of eggs (don't allow to become negative)
end
dRc = R; % change in cumulated reproduction

% NOTE: in principle, we could also work with a life-stage flag to tell the
% model which stage we're in. That would save the check on isnan(Lmat(i)).
% However, the accelerating phase is not necessarily a life stage, and we
% may want to end it at puberty! Further, Lb and Lj are needed here as
% well for the abj module.

if ~isnan(Lmat(1)) && dL < 0  % don't shrink indefinitely
    if L < 0.25 * glo.Lmat(1) % at a quarter of length at birth
        dL = 0;               % simply stop shrinking ...
    end                       % to avoid numerical problems
end

% For the damage dynamics, there are four feedback factors x* that obtain a
% value based on the settings in the configuration vector glo.feedb: a
% vector with switches for various feedbacks: [surface:volume on uptake,
% surface:volume on elimination, growth dilution, losses with
% reproduction].
Xi = glo.feedb .* [Lm_ref/L,Lm_ref/L,(3/L)*dL,R*FBV*KRV]; % multiply switch factor with feedbacks
xu = Xi(1); % factor for surf:vol scaling uptake rate 
xe = Xi(2); % factor for surf:vol scaling elimination rate 
xG = Xi(3); % factor for growth dilution
xR = Xi(4); % factor for losses with repro

% If switch for surf:vol scaling is zero, the factor must be 1 and not 0!
% Note: this was previously done with a max-to-1 command. However, that is
% not a good idea in combination with Lm_ref (which is not necessarily
% equal to or larger than the actual Lm).
if glo.feedb(1) == 0
    xu = 1;
end
if glo.feedb(2) == 0
    xe = 1;
end

xG = max(0,xG); % stop reverse growth dilution
% NOTE NOTE: reverse growth dilution (concentration by shrinking) is now
% turned OFF as it leads to runaway situations that lead to failure of the
% ODE solvers. However, this needs some further thought!

dDw = kd * (xu * c - xe * Dw) - (xG + xR) * Dw; % ODE for scaled damage
dS  = -(h + hb) * S; % change in survival probability (incl. background mort.)

dX = [dE;dEH;dL;dRc;dDw;dS]; % collect all derivatives in one vector
