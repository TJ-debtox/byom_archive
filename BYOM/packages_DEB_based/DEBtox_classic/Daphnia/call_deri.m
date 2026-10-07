%% BYOM function call_deri.m (calculates the model output)
%
%  Syntax: [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)
%
% This function calls the ODE solver to solve the system of differential
% equations specified in <derivatives.html derivatives.m>. It is specific
% for the model in the directory 'DEBtox_classic', and modified to deal
% with time-varying exposure by splitting up the ODE solving to avoid hard
% switches.
%
% As input, it gets:
% * _t_   the time vector
% * _par_ the parameter structure
% * _X0v_   a vector with initial states and one concentration (scenario number)
% * _glo_ the structure with various types of information (used to be global)
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
% * Date: September 2023
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_debtox.html>

%  Copyright (c) 2012-2023, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)

% These outputs need to be defined, even if they are not used
Xout2    = []; % additional uni-variate output, not used in this example
zvd      = []; % additional zero-variate output, not used in this example

%% Initial settings
% This part extracts optional settings for the ODE solver that can be set
% in the main script (defaults are set in prelim_checks). Further in this
% section, initial values can be determined by a parameter (overwrite parts
% of X0), and zero-variate data can be calculated. See the example BYOM
% files for more information.
% 
% This package will always use the ODE solver; therefore the general
% settings in the global glo for the calculation (useode and eventson) are
% removed here. This packacge also always uses the events function, so the
% option eventson is not used. Furthermore, simplefun.m is removed.

stiff      = glo.stiff;      % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
break_time = glo.break_time; % break time vector up for ODE solver (1) or don't (0)
min_t      = 500; % minimum length of time vector (affects ODE stepsize only, when needed)

if length(stiff) == 1
    stiff(2) = 1; % by default: normally tightened tolerances
end

% Unpack the vector X0v, which is X0mat for one scenario
X0 = X0v(2:end); % these are the intitial states for a scenario
c  = X0v(1);     % the concentration (or scenario number)

% Start with a check on time-varying exposure scenarios. This is a
% time-saver! The isfield and ismember calls take some time that rapidly
% multiplies as derivatives is called many times.
glo.timevar = [0 0]; % flag for time-varying exposure (second element is to tell read_scen which interval we need)
if isfield(glo,'int_scen') && ismember(c,glo.int_scen) % is c in the scenario range of the global?
    glo.timevar = [1 0]; % tell derivatives that we have a time-varying treatment
end

% if needed, extract parameters from par that influence initial states in X0
L0    = par.L0(1);  % initial body length (mm)
X0(3) = L0;         % put this parameter in the correct location of the initial vector

% % if needed, calculate model values for zero-variate data from parameter set
% if ~isempty(zvd)
%     zvd.ku(3) = par.Piw(1) * par.ke(1); % add model prediction as third value
% end

%% Calculations
% This part calls the ODE solver to calculate the output (the value of the
% state variables over time). There is generally no need to modify this
% part. The solver ode45 generally works well. For stiff problems, the
% solver might become very slow; you can try ode15s instead.

t  = t(:);      % force t to be a row vector (needed when useode=0)
t_rem = t;      % remember the original time vector (as we will add to it)

% The code below is meant for discontinous time-varying exposure. It will
% also work for constant exposure. Breaking up will not be as effective for
% piecewise polynomials (type 1), and slow, so better turn that off
% manually when you try splining. The code will run each interval between
% two exposure 'events' (in Tev) separately, stop the solver, and restart.
% This ensures that the discontinuities are no problem anymore.

TE  = 0; % dummy for time of events in the events function
Tev = [0 c]; % exposure profile events setting: without anything else, assume it is constant
glo.timevar = [0 0]; % flag for time-varying exposure (second element is to tell read_scen which interval we need)

InitialStep = max(t)/100; % specify initial stepsize
MaxStep     = max(t)/10;  % specify maximum stepsize
% For constant concentrations, and when breaking the time vector, we can
% use a default; Matlab uses as default the length of the time vector
% divided by 10.

if glo.timevar(1) == 1 % time-varying concentrations?
    % if it exists: use it to derive current external conc.
    Tev = read_scen(-2,c,t,glo); % use read_scen to derive actual exposure concentration
    % the -2 lets read_scen know we need events
    
    if break_time == 0 && size(Tev,1) > 2
        min_t = max(min_t,length(Tev(Tev(:,1)<t(end),1))*2);
        % The min_t is only used to limit the step size of the ODE
        % solver, since the all calculations are done in derivatives.m.
        InitialStep = t(end)/(10*min_t); % initial step size
        MaxStep     = t(end)/min_t;      % maximum step size
        % For the ODE solver, when we have a time-varying exposure set
        % here, we base minimum step size on min_t. Small stepsize is a
        % good idea for pulsed exposures; otherwise, stepsize may become so
        % large that a concentration change is missed completely. When we
        % break the time vector, limiting step size is not needed.
    end
end

% This is a means to include a delay caused by the brood pounch in species
% like Daphnia. The repro data are for the appearance of neonates, but egg
% production occurs earlier. This global shifts the model output in this
% function below. This way, the time vector in the data does not need to be
% manipulated, and the model plots show the neonate production as expected.
% (NEW)
bp = 0; % by default, no brood-pouch delay
if isfield(glo,'Tbp') && glo.Tbp > 0 % if there is a global specifying a brood-pouch delay ...
    tbp = t(t>glo.Tbp)-glo.Tbp; % extra times needed to calculate brood-pounch delay
    t   = unique([t;tbp]);      % add the shifted time points to account for brood-pouch delay
    bp  = 1;                    % signal rest of code that we need brood-pouch delay
end

% specify options for the ODE solver
options = odeset; % start with default options for the ODE solver
% This needs further study ... events function removed. Events function
% needs to be considered very carefully for this model (and would only be
% useful for SD).
switch stiff(2)
    case 1 % for ODE15s, slightly tighter tolerances seem to suffice (for ODE113: not tested yet!)
        RelTol  = 1e-4; % relative tolerance (tightened)
        AbsTol  = 1e-7; % absolute tolerance (tightened)
    case 2 % somewhat tighter tolerances ...
        RelTol  = 1e-5; % relative tolerance (tightened)
        AbsTol  = 1e-8; % absolute tolerance (tightened)
    case 3 % for ODE45, very tight tolerances seem to be necessary in some cases
        RelTol  = 1e-9; % relative tolerance (tightened)
        AbsTol  = 1e-9; % absolute tolerance (tightened)
end
options = odeset(options,'RelTol',RelTol,'AbsTol',AbsTol,'Events',@eventsfun,'InitialStep',InitialStep,'MaxStep',MaxStep); % set options
% Note: setting tolerances is pretty tricky. For some cases, tighter
% tolerances are needed but not for others. For ODE45, tighter tolerances
% seem to work well, but not for ODE15s.

T = Tev(:,1); % time vector with events
if T(end) > t(end) % scenario may be longer than t(end)
    T(T>t(end)) = []; % remove all entries that are beyond the last time point
    % this may remove one point too many, but that will be added next
end
if T(end) < t(end) % scenario may (now) be shorter than we need
    T = cat(1,T,t(end)); % then add last point from t
end
% If a lag time is used, it is a good idea to add it as an event as ODE45
% does not like such a switch either.
if isfield(par,'Tlag') && par.Tlag(1) > 0
    T = unique([par.Tlag(1);T]);
end

% Breaking up the time vector is faster and more robust when calibrating on
% data from time-varying exposure that contain discontinuities. However,
% when running through very detailed exposure profiles (e.g., those from
% FOCUS), it is probably best to use break_time=0. Testing indicates that
% breaking up will be much slower than simply running the ODE solver across
% the entire profile (and produces very similar output).
% 
% Note: it is possible to initialise tout and Xout with NaNs, and avoid
% these matrices to grow. This may be a bit slower for exposure scenarios
% with few events, but faster for scenarios with many.

t = unique([T;t;(T(1:end-1)+T(2:end))/2]); % combine T, t, and halfway-T into new time vector
% this hopefully prevents the ODE solver from missing exposure pulses, and
% it makes sure that, for break_time=1, we never have a temporary time
% vector of length 2.

if break_time == 0 
    
    % simply use the ODE solver for the entire time vector
    switch stiff(1)
        case 0
            [tout,Xout,TE,~,~] = ode45(@derivatives,t,X0,options,par,c,glo);
        case 1
            [tout,Xout,TE,~,~] = ode113(@derivatives,t,X0,options,par,c,glo);
        case 2
            [tout,Xout,TE,~,~] = ode15s(@derivatives,t,X0,options,par,c,glo);
    end
    
else
    
    % Transferring the interval to derivatives is a huge time saver! It is
    % not safe to use i as index to Tev since we may add Tlag to T (which
    % is not in Tev); so better calculate it explicitly. Using <find> may
    % be a bit slow, that is why this is done outside of the loop calling
    % the ODE solver.
    ind_Tev = zeros(length(T)-1,1);
    for i = 1:length(T)-1 % run through all intervals between events
        ind_Tev(i) = find(Tev(:,1) <= T(i),1,'last');
    end
    
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

    % run the ODE solver piece-wise across all exposure events specified
    for i = 1:length(T)-1 % run through all intervals between events
        t_tmp = [T(i);T(i+1)]; % start with start and end time for this period
        t_tmp = unique([t_tmp;t(t<t_tmp(2) & t>t_tmp(1))]); % and add the time points from t that fit in there
        glo.timevar(2) = ind_Tev(i); % tell derivatives which part of Tev we're in
        
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
% This is not so useful as only the TE found in the last time period is
% returned.

%% Output mapping
% _Xout_ contains a row for each state variable. It can be mapped to the
% data. If you need to transform the model values to match the data, do it
% here. 

if bp == 1 % if we need a brood-pouch delay ... (NEW)
    [~,loct] = ismember(tbp,tout); % find where the extra brood-pouch time points are in the long Xout
    Xbp      = Xout(loct,glo.locR); % only keep the brood-pouch ones we asked for
end
t = t_rem; % return t to the original input vector

% Select the correct time points to return to the calling function
[~,loct] = ismember(t,tout); % find where the requested time points are in the long Xout
Xout     = Xout(loct,:);     % only keep the ones we asked for

if bp == 1 % if we need a brood-pouch delay ... (NEW)
    [~,loct] = ismember(tbp+glo.Tbp,t); % find where the extra brood-pouch time points SHOULD BE in the long Xout
    Xout(:,glo.locR)    = 0;   % clear the reproduction state variable
    Xout(loct,glo.locR) = Xbp; % put in the brood-pouch ones we asked for
end

% Xout(:,1) = Xout(:,1).^3; % e.g., do something on first column, like cube it ...

% NEW: this is to accommodate data sets with time-to-death, which require a
% (pseudo) hazard rate at the time of observation of death (next to the
% survival over time for animals that survive or that are removed or
% go missing). This still needs to be tested in the classic DEBtox context!
if isfield(glo,'hazout') && glo.hazout == 1
    c0s = par.c0s(1); % no-effect concentration survival
    b   = par.b(1);   % killing rate
    hb  = par.hb(1);  % background hazard rate

    cV  = Xout(:,glo.locD); % take the correct state variable for scaled conc./damage

    haz  = b * max(0,cV-c0s); % calculate the hazard rate
    haz  = min(111,haz); % maximise the hazard rate to 99% mortality in 1 hour
    % Note: this must be the same as in derivatives.
    
    Xout2 = haz(:) + hb; % add background hazard (make sure it is a column vector)
    Xout2 = Xout2(loct); % only keep the time points we asked for
end

% % To obtain the output of the derivatives at each time point. The values in
% % dXout might be used to replace values in Xout, if the data to be fitted
% % are the changes (rates) instead of the state variable itself.
% % dXout = zeros(size(Xout)); % initialise with zeros
% for i = 1:length(t) % run through all time points
%     dXout(i,:) = derivatives(t(i),Xout(i,:),par,c,glo); 
%     % derivatives for each stage at each time
% end

%% Events function
% This subfunction catches the 'events': in this case, this one catches
% when the scaled internal concentration/damage exceeds one of the
% thresholds, and catches the switch at puberty.
%
% Note that the eventsfun has the same inputs, in the same sequence, as
% <derivatives.html derivatives.m>.

function [value,isterminal,direction] = eventsfun(t,X,par,c,glo)

% Note: glo needs to be an input to the events function as well, rather
% than a global since the inputs for the events function must match those
% for derivatives.

Lp  = par.Lp(1);  % length at puberty (mm)
c0  = par.c0(1);  % no-effect concentration sub-lethal
c0s = par.c0s(1); % no-effect concentration survival

nevents = 3;      % number of events that we try to catch

value       = zeros(nevents,1); % initialise with zeros
value(1)    = X(3) - Lp;        % follow when body length exceeds length at puberty
value(2)    = X(1) - c0;        % follow when scaled conc./damage exceed the NEC
value(3)    = X(1) - c0s;       % follow when scaled conc./damage exceed the NEC
isterminal  = zeros(nevents,1); % do NOT stop the solver at an event
direction   = zeros(nevents,1); % catch ALL zero crossing when function is increasing or decreasing
