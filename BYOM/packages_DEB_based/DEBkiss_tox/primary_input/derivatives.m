%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the model system. It is
% linked to the script files in the directory 'primary_input'. As input,
% it gets:
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
% * Date: January 2022
% * Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2022, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives(t,X,par,c,glo)

% NOTE: glo is now no longer a global here, but passed on in the function
% call from call_deri. That saves 20% calculation time!

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

X = max(X,0); % ensure the states cannot become negative due to numerical issues

Dw = X(glo.locD); % state is the scaled damage (referenced to water)
WV = X(glo.locL); % state is dry structural body mass
% Rc = X(glo.locR); % state is cumulative reproduction (not used)
S  = X(glo.locS); % state is survival probability

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

% unpack globals
dV   = glo.dV;      % dry weight density of structure
delM = glo.delM;    % shape correction coefficient
yAV  = glo.yAV;     % yield of assimilates on structure (starvation) (-)
yBA  = glo.yBA;     % yield of egg buffer on assimilates (-)
yVA  = glo.yVA;     % yield of structure on assimilates (growth) (-)
KRV  = glo.KRV;     % part. coeff. repro buffer and structure (kg/kg)

% unpack model parameters for the basic life history
sJAm = par.sJAm(1); % specific assimilation rate (mg/mm^2/d)
sJM  = par.sJM(1);  % specific maintenance costs (mg/mm^3/d)
kap  = par.kap(1);  % allocation fraction to soma (-)
WB0  = par.WB0(1);  % initial dry weight of egg (mg)
LpM  = par.LpM(1);  % body length at puberty (mm)
f    = par.f(1);    % scaled functional response (-)
hb   = par.hb(1);   % background hazard rate (d-1)
a    = par.a(1);    % Weibull background hazard coefficient (-)

% unpack extra parameters for specific cases
LfM  = par.LfM(1);  % body length at half-saturation feeding (mm)
LjM  = par.LjM(1);  % actual body length at which acceleration stops (mm)
Tlag = par.Tlag(1); % lag time before everything starts (d)

% unpack model parameters for the response to toxicants
kd   = par.kd(1);   % dominant rate constant (d-1)
zb   = par.zb(1);   % effect threshold energy budget ([C])
bb   = par.bb(1);   % effect strength energy-budget effects (1/[C])
zs   = par.zs(1);   % effect threshold survival ([C])
bs   = par.bs(1);   % effect strength survival (1/([C] d))

hb = a * (hb^a) * t^(a-1); % option for Weibull mortalty when a is not 1

% Note: par.fh is used in call_deri, rather than here. That is done for
% speed! Checking whether fh is in par is rather slow, and derivatives is
% called many times.

% translate parameters on length basis to weight basis
WVp = dV * (LpM * delM)^3; % translate puberty length to dry weight
WVf = dV * (LfM * delM)^3; % translate feeding-saturation length to dry weight
Lj  = LjM * delM;          % translate actual acceleration length to volumetric length

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
% This is the actual model, specified as a system of ODEs. This is the full
% DEBkiss model, expressed in primary parameters, extended with damage
% dynamics and effects. See DEBkiss e-book Chapter 5. Note: some unused
% fluxes are calculated because they may be needed for future modules
% (e.g., respiration and ageing).

L  = (WV/dV)^(1/3); % calculate current volumetric length from dry weight
Lm = kap*sJAm/sJM;  % maximum structural body length in control (for scaling kd)

if WVf > 0 % to include feeding limitation for juveniles ...
    sf   = 1 / (1+WVf/WV); % hyperbolic relationship for f with body weight
    sJAm = sJAm * sf; % reduce the max assimilation rate
    % kd = kd*sf; % also reduce dominant rate by same factor? (open for discussion!)
end
if Lj > 0 % to include acceleration until metamorphosis ...
    f = f * min(1,L/Lj); % this implies lower f for L<Lj
end

% select what to do with maturity maintenance
if glo.mat == 1
    sJJ = sJM * (1-kap)/kap; % add specific maturity maintenance with the suggested value
else
    sJJ = 0; % or ignore it completely (simplest DEBkiss model)
end
WVb = WB0 * yVA * kap; % approximate body mass at birth (to put a stop on shrinking)

% calculate stress factor and hazard rate
s = bb*max(0,Dw-zb); % stress level for metabolic effects
h = bs*max(0,Dw-zs); % hazard rate for effects on survival

h = min(111,h); % maximise the hazard rate to 99% mortality in 1 hour
% Note: this helps in extreme conditions, as the system becomes stiff for
% very high hazard rates. This is especially needed for EPx calculations,
% where the MF is increased until there is effect on all endpoints!

% Define specific stress factors s*, depending on the mode of action as
% specified in the vector with switches glo.moa.
Si = glo.moa * s;  % vector with specific stress factors from switches for mode of action
sA = min(1,Si(1)); % assimilation/feeding (maximise to 1 to avoid negative values for 1-sA)
sM = Si(2);        % maintenance (somatic and maturity)
sG = Si(3);        % growth costs
sR = Si(4);        % reproduction costs               
sH = Si(5);        % also include hazard to reproduction

% apply the specific stress factors s* to the primary parameters
sJAm = sJAm *(1-sA); % effect on assimilation, alternative: /(1+s) increase in handling time
sJM  = sJM * (1+sM); % effect on somatic maintenance
sJJ  = sJJ * (1+sM); % effect on maturity maintenance
yVA  = yVA / (1+sG); % effect on growth costs
yBA  = exp(-sH) * yBA / (1+sR); % add effect on repro costs and hazards

% Note: when including maturation in the model explicitly, it makes sense
% to start by combining 'costs for growth' with 'costs for maturation'.
% This ensures that body length is a good proxy for maturity level.

% calculate the main fluxes
JA = f * sJAm * L^2;          % assimilation flux
JM = sJM * L^3;               % somatic maintenance flux
JV = yVA * (kap*JA-JM);       % growth flux

if WV < WVp                   % below size at puberty ...
    JR = 0;                   % no reproduction flux
    JJ = sJJ * L^3;           % maturity maintenance flux
    % JH = (1-kap) * JA - JJ; % maturation flux (not used!)
else                          % above size at puberty ...
    JJ = sJJ * (WVp/dV);      % maturity maintenance flux
    JR = (1-kap) * JA - JJ;   % reproduction flux
end

% starvation rules may override these fluxes
if kap * JA < JM      % allocated flux to soma cannot pay maintenance
    if JA >= JM + JJ  % but still enough total assimilates to pay both maintenances
        JV = 0;       % stop growth
        if WV >= WVp  % for adults ...
            JR = JA - JM - JJ; % repro buffer gets what's left
        else
            % JH = JA - JM - JJ; % maturation gets what's left (not used!)
        end
    elseif JA >= JM   % only enough to pay somatic maintenance
        JV = 0;       % stop growth
        JR = 0;       % stop reproduction
        % JH = 0;     % stop maturation for juveniles (not used)
        % JJ = JA - JM; % maturity maintenance flux gets what's left (not used)
    else              % we need to shrink
        JR = 0;       % stop reproduction
        % JJ = 0;       % stop paying maturity maintenance (not used)
        % JH = 0;     % stop maturation for juveniles (not used)
        JV = (JA - JM) / yAV; % shrink; pay somatic maintenance from structure
    end
end

% calculate the derivatives
dWV = JV;             % change in body mass
dRc = yBA * JR / WB0; % continuous reproduction flux
dS  = -(h + hb)* S;   % change in survival probability (incl. background mort.)

% For the damage dynamics, there are four feedback factors x* that obtain a
% value based on the settings in the configuration vector glo.feedb: a
% vector with switches for various feedbacks: [surface:volume on uptake,
% surface:volume on elimination, growth dilution, losses with
% reproduction].
Xi = glo.feedb .* [Lm/L,Lm/L,dWV/WV,dRc*(WB0/WV)*KRV]; % multiply switch factor with feedbacks
xu = max(1,Xi(1)); % factor for surf:vol scaling uptake rate 
xe = max(1,Xi(2)); % factor for surf:vol scaling elimination rate 
xG = Xi(3);        % factor for growth dilution
xR = Xi(4);        % factor for losses with repro

% if switch for surf:vol scaling is zero, the factor must be 1 and not 0!
% Note: this was previously done with a max-to-1 command. However, that is
% a bit dangerous when modifying the model.
if Xi(1) == 0
    xu = 1;
end
if Xi(2) == 0
    xe = 1;
end

xG = max(0,xG); 
% NOTE NOTE: reverse growth dilution (concentration by shrinking) is now
% turned OFF as it leads to runaway situations that lead to failure of the
% ODE solvers. However, this needs some further thought!

dDw = kd * (xu * c - xe * Dw) - (xG + xR) * Dw; % ODE for scaled damage

% to avoid shrinking to lead to negative body size ...
if WV < WVb/2 % dont shrink below half of the size at birth ...
    dWV = 0;  % simply stop shrinking ...
end

dX = zeros(size(X)); % initialise with zeros in right format
if t >= Tlag % when we are past the lag time ...
    dX([glo.locD glo.locL glo.locR glo.locS]) = [dDw;dWV;dRc;dS]; % collect all derivatives in one vector dX
end
