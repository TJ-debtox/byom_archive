%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the model system. It is
% linked to the script files for the snails in the DEBkiss1 package. As
% input, it gets:
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
% * Back to index <walkthrough_debkiss_paper.html>

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

WV = X(1); % state 1 is the structural body mass
% cR = X(2); % state 2 is the cumulative reproduction (not used here)
WB = X(3); % state 3 is the egg buffer of assimilates

WV = max(0,WV); % for poor choice of paramaters, WV can become nagtive, which smothers the ODE solver

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

Tlag  = glo.Tlag;   % lag time for start embryo development
dV    = glo.dV;     % dry weight density of structure

if c==1
    Tlag = 0; % the c=1 scenario for embryos does not include a lag time
    % for the juvenile/adult runs, this does not hurt as Tlag=0 anyway
end

sJAm  = par.sJAm(1);  % specific assimilation rate 
sJM   = par.sJM(1);   % specific maintenance costs 
WB0   = par.WB0(1);   % initial weight of egg
WVp   = par.WVp(1);   % body mass at puberty
yAV   = par.yAV(1);   % yield of assimilates on volume (starvation)
yBA   = par.yBA(1);   % yield of egg buffer on assimilates
yVA   = par.yVA(1);   % yield of structure on assimilates (growth)
kap   = par.kap(1);   % allocation fraction to soma
f1    = par.f1(1);    % scaled food level, ad libitum
f2    = par.f2(1);    % scaled food level, regime 2
f3    = par.f3(1);    % scaled food level, regime 3
fB    = par.fB(1);    % scaled food level for the embryo (adapted model)
WVf   = par.WVf(1);   % half-saturation size for initial food limitation

%% Calculate the derivatives
% This is the actual model, specified as a system of ODEs.
% Note: some unused fluxes are calculated because they may be needed for
% future modules (e.g., respiration and ageing).

L = (WV/dV)^(1/3); % volumetric length from dry weight

% select what to do with maturity maintenance
if glo.mat == 1
    sJJ = sJM * (1-kap)/kap; % add specific maturity maintenance with the suggested value
else
    sJJ = 0; % or ignore it
end

WVb = WB0 *yVA * kap; % predicted body mass at birth

if WB > 0 % if we have an embryo (run from byom_snail_embryo.m)
    if c==2 % for the 'adapted' scenario ...
        f = fB; % assimilation at limited rate
    else % for the 'standard' scenario
        f = 1;  % assimilation at maximum rate
    end
else % we have juveniles/adults (run from byom_snail_fit.m or byom_snail_starv.m)
    switch c % run through the three feeding regimes to select proper f
        case 1
            f = f1;
        case 2
            f = f2;
            yAV = 7*0.8; % this is used for the adapted starvation simulation
            % no problem if it is also on for the fits, as there is no
            % starvation occurring there anyway!
        case 3
            f = f3;
    end
    if WVf > 0
        f = f / (1+WVf/WV); % hyperbolic relationship for f with body weight
    end
end

% if WV < WVb && WB <= 0
%     error('Your animal is initiated or has shrunk to less than the size at birth! Increase the start weight or check the food levels.')
% end

JA = f * sJAm * L^2;          % assimilation
JM = sJM * L^3;               % somatic maintenance
JV = yVA * (kap*JA-JM);       % growth

if WV < WVp                   % below size at puberty
    JR = 0;                   % no reproduction
    JJ = sJJ * L^3;           % maturity maintenance flux
    % JH = (1-kap) * JA - JJ; % maturation flux (not used!)
else
    JJ = sJJ * (WVp/dV);      % maturity maintenance
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
        JJ = JA - JM; % maturity maintenance flux gets what's left
    else              % we need to shrink
        JR = 0;       % stop reproduction
        JJ = 0;       % stop paying maturity maintenance
        % JH = 0;     % stop maturation for juveniles (not used)
        JV = (JA - JM) / yAV; % shrink; pay somatic maintenance from structure
    end
end

% we do not work with a repro buffer here, so a negative JR does not make
% sense; therefore, minimise to zero
JR = max(0,JR);

% Calcululate the derivatives
dWV = JV;             % change in body mass
dcR = yBA * JR / WB0; % continuous reproduction flux
if WB > 0             % for embryos ...
    dWB = -JA;        % decrease in egg buffer
else
    dWB = 0;          % for juveniles/adults, no change
end
if WV < WVb/4 && dWB == 0 % dont shrink below quarter size at birth
    dWV = 0; % simply stop shrinking ...
end

if t < Tlag 
    dX = [0;0;0];
else
    dX = [dWV;dcR;dWB];       % collect both derivatives in one vector
end
