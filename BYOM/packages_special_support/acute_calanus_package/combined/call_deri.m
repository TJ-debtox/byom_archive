%% BYOM function call_deri.m (calculates the model output)
%
%  Syntax: [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)
%
% This function calls the ODE solver to solve the system of differential
% equations specified in <derivatives.html derivatives.m>, or the explicit
% function(s) in <simplefun.html simplefun.m>. As input, it gets:
%
% * _t_   the time vector
% * _par_ the parameter structure
% * _X0v_   a vector with initial states and one concentration (scenario number)
% * _glo_ is the structure with information (normally global)
%
% The output _Xout_ provides a matrix with time in rows, and states in
% columns. This function calls <derivatives.html derivatives.m>. The
% optional output _TE_ is the time at which an event takes place (specified
% using the events function). The events function is set up to catch
% discontinuities. It should be specified according to the problem you are
% simulating. If you want to use parameters that are (or influence) initial
% states, they have to be included in this function.
%
% * Author: Tjalling Jager
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)

% These outputs need to be defined, even if they are not used
Xout2    = []; % additional uni-variate output, not used in this example
zvd      = []; % additional zero-variate output, not used in this example

%% Initial settings
% This part extracts optional settings for the ODE solver that can be set
% in the main script (defaults are set in prelim_checks). Note that the
% files in this folder always use the ODE solver, and that the events
% function is also always used (to catch when the threshold is exceeded).
% Further in this section, initial values can be determined by a parameter
% (overwrite parts of _X0_), and zero-variate data can be calculated. See
% the example BYOM files for more information.

stiff      = glo.stiff;    % set to 1 or 2 to use a stiff solver instead of the standard one
min_t      = 500;          % minimum length of time vector (affects ODE stepsize as well)
break_time = glo.break_time; % break time vector up for ODE solver (1) or don't (0)

if length(stiff) == 1
    stiff(2) = 2; % by default: moderately tightened tolerances
end

% Unpack the vector X0v, which is X0mat for one scenario
X0 = X0v(2:end); % these are the intitial states for a scenario

%% Calculations
% This part calls the ODE solver to calculate the output (the value of the
% state variables over time). There is generally no need to modify this
% part. The solver ode45 generally works well. For stiff problems, the
% solver might become very slow; you can try ode113 or ode15s instead.

c     = X0v(1); % the concentration (or scenario number)
t     = t(:);   % force t to be a row vector (needed when useode=0)
t_rem = t;      % remember the original time vector (as we will add to it)

glo.timevar = [0 0]; % flag for time-varying exposure (second element is to tell read_scen which interval we need)
TE  = 0; % dummy for time of events in the events function
Tev = [0 c]; % exposure profile events setting: without anything else, assume it is constant

InitialStep = max(t)/100; % specify initial stepsize
MaxStep     = max(t)/10;  % specify maximum stepsize
% For constant concentrations, and when breaking the time vector, we can
% use a default; Matlab uses as default the length of the time vector
% divided by 10.

if isfield(glo,'int_scen') && ismember(c,glo.int_scen) % is c in the scenario range global?
    % then we have a time-varying concentration
    
    glo.timevar = [1 0]; % tell derivatives that we have a time-varying treatment
    % This is a time-saver! The isfield and ismember calls take some
    % time that rapidly multiplies as derivatives is called many times.
    
    Tev = read_scen(-2,c,-1,glo); % use make_scen again to derive actual exposure concentration
    % the -2 lets make_scen know we need events, the -1 that this is also needed for splines!
    min_t = max(min_t,length(Tev(Tev(:,1)<t(end),1))*2); 
    % For very long exposure profiles (e.g., FOCUS profiles), we now
    % automatically generate a larger time vector (twice the number of the
    % relevant points in the scenario).
    if break_time == 0
        InitialStep = t_rem(end)/(10*min_t); % initial step size
        MaxStep     = t_rem(end)/min_t; % maximum step size
        % For the ODE solver, when we have a time-varying exposure set here, we
        % base minimum step size on min_t. Small stepsize is a good idea for
        % pulsed exposures; otherwise, stepsize may become so large that a
        % concentration change is missed completely. When we break the time
        % vector, limiting step size is not needed.
    end
    
end
% For all GUTS cases, we need a long time vector. For the hazard and proper
% model, we need to numerically integrate the hazard rate over time. For
% IT, we need to find the maximum damage. If the time span of the data set
% or simulation is very long (and the number of points in Tev is not),
% consider increasing the size of min_t.
if length(t) < min_t % make sure there are at least min_t points
    t = unique([t;(linspace(t(1),t(end),min_t))']);
end

% Always use the ODE solver to calculate the solution, and always use the
% events function.
% 
% Note: when using an ODE solver and time-varying exposure, it is a good
% idea to break up the time vector and analyse them piecewise. 

% This needs some further study ... The tolerances work out differently for
% different solvers.
switch stiff(2)
    case 1 % for ODE15s, slightly tighter tolerances seem to suffice (for ODE113: not tested yet!)
        RelTol  = 1e-4; % relative tolerance (tightened)
        AbsTol  = 1e-7; % absolute tolerance (tightened)
    case 2 % somewhat tighter tolerances ...
        RelTol  = 1e-5; % relative tolerance (tightened)
        AbsTol  = 1e-8; % absolute tolerance (tightened)
    case 3 % for ODE45, very tight tolerances seem to be necessary
        RelTol  = 1e-9; % relative tolerance (tightened)
        AbsTol  = 1e-9; % absolute tolerance (tightened)
end

% The code below is meant for discontinous time-varying exposure. It will
% also work for constant exposure, but breaking up will probably not be
% very effective for piecewise polynomials (type 1). It will run each
% interval between two exposure 'events' (in Tev) separately, stop the
% solver, and restart. This ensures that the discontinuities are no problem
% anymore.

T = Tev(:,1); % time vector with events
if T(end) > t(end) % scenario may be longer than t(end)
    T(T>t(end)) = []; % remove all entries that are beyond the last time point
    % this may remove one point too many, but that will be added next
end
if T(end) < t(end) % scenario may (now) be shorter than we need
    T = cat(1,T,t(end)); % then add last point from t
end

options = odeset; % start with default options for the ODE solver

% With an events functions ... we have additional output arguments for
% events: TE catches the time of an event, YE the states at the event,
% and IE the number of the event
options = odeset(options,'RelTol',RelTol,'AbsTol',AbsTol,'Events',@eventsfun,'InitialStep',InitialStep,'MaxStep',MaxStep); % add an events function

t = unique([T;t;(T(1:end-1)+T(2:end))/2]); % combine T, t, and halfway-T into new time vector

if break_time == 0

    % Simply use the ODE solver once
    switch stiff(1) % note that YE and IE can be added as additional outputs
        case 0
            [tout,Xout,TE,~,~] = ode45(@derivatives,t,X0,options,par,c,glo);
        case 1
            [tout,Xout,TE,~,~] = ode113(@derivatives,t,X0,options,par,c,glo);
        case 2
            [tout,Xout,TE,~,~] = ode15s(@derivatives,t,X0,options,par,c,glo);
    end
    
else
    % NOTE: predefining *should* increase speed, but the speed gain is not
    % so clear in this case. I expect some gain when using many exposure
    % intervals.
    tout    = nan(length(t)+5,1); % this vector will collect the time output from the solver (5 more than needed)
    Xout    = nan(length(tout),length(X0)); % this matrix will collect the states output from the solver (5 more than needed)
    ind_i   = 1; % index for where we are in tout and Xout
    % Note: I initialise tout and Xout with more than the elements of t
    % because the events may be added as well (only stopping events, which
    % are not used yet). The number is rather arbitrary but should be equal
    % to, or more than, the number of events.

    for i = 1:length(T)-1 % run through all intervals between events
        t_tmp = [T(i);T(i+1)]; % start with start and end time for this period
        t_tmp = unique([t_tmp;t(t<t_tmp(2) & t>t_tmp(1))]); % and add the time points from t that fit in there
        glo.timevar = [glo.timevar(1) i]; % to tell read_scen which interval we are in
        % Transferring the interval is a huge time saver!
        
        % NOTE: adding a time point in between when t_tmp is length 2 is
        % not needed. We've added T and half-way-into T into time vector t.
        % Since we do NOT stop at events, there is no way to have a time
        % vector of two elements!

        % use ODE solver to find solution in this interval
        % Note: the switch-case may be compromising calculation speed; however,
        % it needs to be investigated which solver performs best
        switch stiff(1)
            case 0
                [tout_tmp,Xout_tmp,TE,~,~] = ode45(@derivatives,t_tmp,X0,options,par,c,glo);
            case 1
                [tout_tmp,Xout_tmp,TE,~,~] = ode113(@derivatives,t_tmp,X0,options,par,c,glo);
            case 2
                [tout_tmp,Xout_tmp,TE,~,~] = ode15s(@derivatives,t_tmp,X0,options,par,c,glo);
        end

        % collect output in correct location
        nt = length(tout_tmp); % length of the output time vector
        tout(ind_i:ind_i-1+nt)   = tout_tmp; % add tout_tmp to correct position in tout
        Xout(ind_i:ind_i-1+nt,:) = Xout_tmp; % add Xout_tmp to correct position in Xout

        ind_i = ind_i-1+nt;       % update ind_i for next round
        X0    = Xout_tmp(end,:)'; % update X0 for next round
        % Note: the way ind_i is used, there will be no double time points
        % in tout and Xout.
    end
        
end
    
if isempty(TE) || all(TE == 0) % if there is no event caught
    TE = +inf; % return infinity
end

%% Output mapping
% _Xout_ contains a row for each state variable. It can be mapped to the
% data. If you need to transform the model values to match the data, do it
% here. 
% 
% This set of files is geared towards the use of the ODE-solver solution
% for the damage level in derivatives. Therefore, only the damage level and
% background survival is available at this point, and the survival due to
% the chemical is added for each death mechanism below.

mi   = par.mS(1);             % median of threshold distribution
bi   = par.bS(1);             % killing rate
Fs   = max(1+1e-6,par.Fs(1)); % fraction spread of the threshold distribution
beta = log(39)/log(Fs);       % shape parameter for logistic from Fs
S    = Xout(:,glo.locS);      % take background survival from the model
Di   = Xout(:,glo.locD);      % take the correct state variable for scaled damage

switch glo.sel
    case 1 % stochastic death
        haz    = bi * max(0,Di-mi);          % calculate hazard for each NEC
        cumhaz = cumtrapz(tout,haz);            % integrate the hazard rate numerically
        S      = S .* min(1,exp(-1*cumhaz)); % calculate survival probability, incl. background
    
    case 2 % individual tolerance, log-logistic distribution is used
        mi = max(mi,1e-20); % make sure that the threshold is not zero ...
        % Make sure that Dw does not decrease over time (dead animals dont
        % become alive). Uses built-in Matlab function cummax.
        maxDi = cummax(Di); % copy the vector Di to maxDi and find cumulative maximum
        S     = S .* (1 ./ (1+(maxDi/mi).^beta)); % survival probability
        % the survival due to the chemical is multiplied with the background
        % survival (calculated in derivatives), assuming logistic distrib.

    case 3 % mixed model
        % This is a fast way for GUTS proper. Damage is calculated only
        % once, as it is the same for all individuals. Survival for
        % different NECs is calculated below.
        n          = 200; % number of slices from the threshold distribution
        Fs2        = 999^(1/beta); % fraction spread for 99.9% of the distribution
        z_range    = linspace(mi/(1.5*Fs2),mi*Fs2,n); % range of NECs to cover 99.9%
        prob_range = ((beta/mi)*(z_range/mi).^(beta-1)) ./ ((1+(z_range/mi).^beta).^2); % pdf for the log-logistic (Wikipedia)
        prob_range = prob_range / sum(prob_range); % normalise the densities to exactly one
        S1         = zeros(length(tout),1); % initialise the survival probability over time with zeros
        Di = Xout(:,glo.locD);   % take the correct state variable for scaled damage
        for i = 1:n % run through the different thresholds
            haz    = bi * max(0,Di-z_range(i)); % calculate hazard for each NEC
            cumhaz = cumtrapz(tout,haz);           % integrate the hazard rate numerically
            surv   = min(1,exp(-1*cumhaz));     % calculate survival probability
            S1     = S1 + surv * prob_range(i); % add to S1, weighted for this NECs prob density
        end
        S = S .* S1; % make sure to add background hazard
end

FLS = glo.FLS; % VL/VS
Xout(:,4) = (Xout(:,2)+FLS*Xout(:,3))/(1+FLS); % make the fourth state the total concentration

Xout(:,glo.locS) = S;           % replace correct state by newly calculated survival prob.
[~,loct] = ismember(t_rem,tout);   % find where the requested time points are in the long Xout
Xout = Xout(loct,:);            % only keep the ones we asked for

%% Events function
% Modify this part of the code if _eventson_=1. 
% This subfunction catches the 'events': in this case, it looks for the
% internal concentration where the threshold is exceeded, and the time at
% which depuration starts (for body residues only).
%
% Note that the eventsfun has the same inputs, in the same sequence, as
% <derivatives.html derivatives.m>.

function [value,isterminal,direction] = eventsfun(t,X,par,c,glo)

mS      = par.mS(1); % median for thresholds
nevents = 2;         % number of events that we try to catch

value      = zeros(nevents,1); % initialise with zeros
value(1)   = X(glo.locD) - mS; % thing to follow is damage minus threshold
if c == 1.2
    value(2)   = t - 4; % catch depuration phase
end
isterminal = zeros(nevents,1); % do NOT stop the solver at an event
direction  = zeros(nevents,1); % catch ALL zero crossing when function is increasing or decreasing
