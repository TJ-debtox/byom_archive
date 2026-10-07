%% BYOM function simplefun.m (the model as explicit equations)
%
%  Syntax: Xout = simplefun(t,X0,par,c,glo)
%
% This function calculates the output of the reduced GUTS-SD/IT model
% system for binary mixtures (damage addition).
%
% As input, it gets:
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
% Author: Tjalling Jager 
% Date: December 2021
% Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function Xout = simplefun(t,X0,par,c,glo)

%% Unpack initial states
% The state variables enter this function in the vector _X0_. The initial
% scaled damage level does not have to be zero. However, a non-zero value
% would generally be meaningless (look at the full model for examples with
% non-zero initial body residue). Furthermore, non-zero _Dw0_ will likely
% lead to erroneous results for the calculation of LCx and LPx.

% S0  = X0(1); % survival probability at t=0
DwA0 = X0(2); % scaled damage (referenced to external concentrations) at t=0
DwB0 = X0(3); % scaled damage (referenced to external concentrations) at t=0

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

hb   = par.hb(1);   % background hazard rate
kdA  = par.kdA(1);  % dominant rate constant chemical A
kdB  = par.kdB(1);  % dominant rate constant chemical B
mw   = par.mw(1);   % median threshold
bw   = par.bw(1);   % killing rate
Fs   = max(1+1e-6,par.Fs(1)); % fraction spread of the threshold distribution (should not be 1)
WB   = par.WB(1);   % weight factor for chemical B (relative to A)
IAB  = par.IAB(1);  % optional interaction factor

%% Calculate the model output
% This is the actual model, specified as explicit function(s). For the
% GUTS-SD cases, we need a long time vector as we need to numerically
% integrate the hazard rate over time. For IT, this is not strictly
% necessary. However, as soon as one starts with time-varying exposure, it
% may be needed, and it will not hurt calculation speed that much.
% 
% NOTE: this file is set up to analyse and predict effects in laboratory
% toxicity tests with a limited duration. If used to analyse longer time
% scales (or FOCUS profiles) a much longer time vector would be needed!
% 
% NOTE: in general, it is a good idea to have the time points from the
% exposure profile into the model time vector. However, this will not make
% a relevant difference for laboratory toxicity tests.

min_t = 500; % minimum length of time vector
t_rem = t;   % remember the original time vector (as we will add to it)
if length(t) < min_t % make sure there are at least min_t points
    t = unique([t;(linspace(t(1),t(end),min_t))']);
end

% Extract identifiers for each compound (assume a factor of mix_fact was used in coding the identifiers)
mix_fact = glo.mix_fact;
cA = mix_fact*floor(c/mix_fact); % identifier for compound A
cB = c - cA;         % identifier for compound B

% assume regular kinetics ... (so no fast or slow kinetics)
TevA = [0 0]; % events setting: without anything else, assume it is zero! 
TevB = [0 0]; % events setting: without anything else, assume it is zero! 
kcA  = 0;     % assume no disappearance of the test chemical
kcB  = 0;     % assume no disappearance of the test chemical
if isfield(glo,'int_scen') % if it exists: use it to derive current external conc.
    if ismember(cA,glo.int_scen) % is cA in the scenario range global?
        [TevA,kcA] = read_scen(-2,cA,t,glo); % use read_scen to derive actual exposure concentration
        % the -2 lets make_scen know we are calling from simplefun (and need events and kc)
    end
    if ismember(cB,glo.int_scen) % is cB in the scenario range global?
        [TevB,kcB] = read_scen(-2,cB,t,glo); % use read_scen to derive actual exposure concentration
        % the -2 lets make_scen know we are calling from simplefun (and need events and kc)
    end    
end

% Calculate two scaled damage levels over time, using sub-function taken
% from the standard GUTS package! The sub-function is at the end of this
% function.
DwA = calc_Dw(t,TevA,kdA,kcA,DwA0);
DwB = calc_Dw(t,TevB,kdB,kcB,DwB0);

% Add scaled damage vectors with a weight factor; we are doing addition here
DwT = DwA + WB * DwB + IAB .* DwA .* DwB;

if glo.sel == 1 % SD translate damage into survival
    
    hz     = bw * max(0,DwT-mw);    % calculate hazard for each time point
    cumhaz = cumtrapz(t,hz+hb);     % integrate the hazard rate numerically
    S      = min(1,exp(-1*cumhaz)); % calculate survival probability, incl. background

elseif glo.sel == 2 % IT translate damage into survival
    
    beta = log(39)/log(Fs); % shape parameter for logistic from Fs
    mw   = max(mw,1e-100);  % make sure that the threshold is not exactly zero ...
    % Make sure that damage does not decrease over time (dead animals don't
    % become alive). This is only needed for time-varying exposure.
    maxDw = DwT; % copy the vector DwT to maxDw
    ind   = find([0;diff(maxDw)]<0,1,'first'); % first index to places where Dw has decreased
    while ~isempty(ind) % as long as there is a decrease somewhere ...
        maxDw(ind:end) = max(maxDw(ind:end),maxDw(ind-1)); % replace every later time with max of that and previous point
        ind = find([0;diff(maxDw)]<0,1,'first'); % any decrease left?
    end
    S = exp(-hb*t) .* (1 ./ (1+(maxDw/mw).^beta)); % survival probability
    % the survival due to the chemical is multiplied with the background survival

end

Xout(:,[glo.locS glo.locD]) = [S DwA DwB]; % combine all state variables into a matrix

[~,loct] = ismember(t_rem,t); % find where the requested time points are in the long Xout
Xout     = Xout(loct,:);      % only keep the ones we asked for


function Dw = calc_Dw(t,Tev,kd,kc,Dw0)

% Sub-function with the standard analytical solutions from the GUTS
% package.

% initialise internal concentrations and damage with NaNs
Dw    = nan(length(t),1);
Dw(1) = Dw0; % the first element is the starting concentration

diff_rel = 1e-5; % minimum relative difference between the rate constants

if size(Tev,2) == 2 % than we have a pulsed or static renewal scenario

    if abs(1-kd/kc) < diff_rel % the analytical solution does not allow the two rate constants
        % to be exactly the same (or too close) ... when kc=0, the abs
        % gives inf, so the part of the code below is not run.
        kd = kd * (1+diff_rel); % increase kd tiny bit
    end

    for i = 1:size(Tev,1) % run through all event periods
        
        a = kd * 1 * Tev(i,2); % modify the uptake flux to the new start concentration in this period
        % compound parameter a is used in the solution below
        
        if i<size(Tev,1)
            ind_t = (t>Tev(i,1) & t<=Tev(i+1,1)); % find the logical indices for the period
        else
            ind_t = (t>Tev(i,1)); % find the logical indices for the last period
            % if the scenario is longer than the data, this is all zeros,
            % which leads to an empty te, but this does not produce an
            % error so it is fine.
        end
        te = t(ind_t) - Tev(i,1); % find the new part of the time vector, and make it start at zero again
        
        % analytical solution
        Dw(ind_t) = (a/(kd-kc))*exp(-kc * te) + (Dw0 - (a/(kd-kc))) * exp(-kd * te);
        
        % find new starting values at exact moment of new period
        if i < size(Tev,1) % don't do this for last period
            te = Tev(i+1,1) - Tev(i,1); % start time for new period
            Dw0 = (a/(kd-kc))*exp(-kc * te) + (Dw0 - (a/(kd-kc))) * exp(-kd * te);
        end
    end
    
else % we are using linear interpolation in a forcing series
    
    if t(end) > Tev(end,1) % do we ask for more points than in Tev?
        Tev = [Tev; t(end) 0 0]; % add a dummy time point in Tev (from t)
    end

    for i = 1:size(Tev,1)-1 % run through all events (not last one that we added)
        
        ind_t = (t>Tev(i,1) & t<=Tev(i+1,1)); % find the logical indices for the period
        te = t(ind_t) - Tev(i,1); % find the new part of the time vector, and make it start at zero again
        
        % Thanks to Bob Kooi for using Maple to obtain the solution
        a = Tev(i,2); % take as initial concentration the new one in this period
        b = Tev(i,3); % take as slope the new one in this period
        Dw(ind_t) = b * te + a - b / kd + exp(-kd * te) * (Dw0 - a + b / kd);
        
        % find new starting values at exact moment of new period
        if i < size(Tev,1)-1 % don't do this for last period
            te   = Tev(i+1,1) - Tev(i,1); % start time for new period
            Dw0 = b * te + a - b / kd + exp(-kd * te) * (Dw0 - a + b / kd);
        end
    end
    
end


