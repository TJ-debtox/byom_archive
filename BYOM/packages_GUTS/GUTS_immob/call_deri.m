%% BYOM function call_deri.m (calculates the model output)
%
%  Syntax: [Xout TE] = call_deri(t,par,X0v,glo)
%
% This function calls the ODE solver to solve the system of differential
% equations specified in <derivatives.html derivatives.m>. This function
% calculates the model output for the GUTS model system with extra states
% for immobility. The full model is used here, so TK and damage dynamics
% are calculated separately. As input, it gets:
%
% * _t_   the time vector
% * _par_ the parameter structure
% * _X0v_   a vector with initial states and one concentration (scenario number)
% * _glo_ the structure with various types of information (used to be global)
%
% The output _Xout_ provides a matrix with time in rows, and states in
% columns. This function calls <derivatives.html derivatives.m>. The
% optional output _TE_ is the time at which an event takes place (specified
% using the events function, which is for now NOT used). 
%
% This function is for the immobility GUTS model, and is modified from the
% standard one in the GUTS package. For SD, <derivatives.html
% derivatives.m>. is used to calculate all state variables, apart from the
% sum of immobiles and dead animals. For IT, <derivatives.html
% derivatives.m>. only calculates internal concentration and damage;
% mortality due to chemical stress is calculated in this function as a form
% of 'output mapping'.
%
% * Author: Tjalling Jager
% * Date: January 2024
% * Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2024, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo) 

Xout2    = []; % additional uni-variate output, not used in this example
zvd      = []; % additional zero-variate output, not used in this example

%% Initial settings
% Note: all the regular global settings in glo (useode, eventson and stiff)
% are ignored for this function. Instead, the code below only uses the ODE
% solver. Events function is NOT used for now, but the settings for stiff
% solver etc. are now inclided. 

stiff      = glo.stiff; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
min_t      = 500;       % minimum length of time vector (affects ODE stepsize as well, when break_time=0)
break_time = glo.break_time; % break time vector up for ODE solver (1) or don't (0)
% Breaking the time vector is generally a good idea when modelling pulsed forcings
if length(stiff) == 1 % backwards compatability
    stiff(2) = 1; % by default: normally tightened tolerances
end

% Unpack the vector X0v, which is X0mat for one scenario
X0 = X0v(2:end); % these are the intitial states for a scenario

%% Calculations
% This part calls the ODE solver to calculate the output (the value of the
% state variables over time). There is generally no need to modify this
% part. The solver ode45 generally works well. For stiff problems, the
% solver might become very slow; you can try ode113 or ode15s instead.

c       = X0v(1); % the concentration (or scenario number)
t       = t(:);   % force t to be a row vector (needed when useode=0)
t_rem   = t;      % remember the original time vector (as we will add to it)

% The code below is meant for discontinous time-varying exposure. It will
% also work for constant exposure. Note that breaking up the time vector
% will generally not be very useful for piecewise polynomials (type 1). It
% will run each interval between two exposure 'events' (in Tev) separately,
% stop the solver, and restart. This ensures that the discontinuities are
% no problem anymore.

TE  = 0;     % dummy for time of events in the events function
Tev = [0 c]; % exposure profile events setting: without anything else, assume it is constant
glo.timevar = [0 0]; % flag for time-varying exposure (second element is to tell read_scen which interval we need)

InitialStep = max(t)/100; % specify initial stepsize
MaxStep     = max(t)/10;  % specify maximum stepsize
% For constant concentrations and SD, we can use a default; Matlab uses as
% default the length of the time vector divided by 10.

if isfield(glo,'int_scen') && ismember(c,glo.int_scen) % is c in the scenario range global?
    glo.timevar = [1 0]; % tell simplefun and derivatives that we have a time-varying treatment
    % This is a time-saver! The isfield and ismember calls take some
    % time that rapidly multiplies as derivatives is called many times.
    Tev = read_scen(-2,c,-1,glo); % use read_scen to derive actual exposure concentration
    % the -2 lets read_scen know we need events, the -1 that this is also needed for splines!
    
    min_t = max(min_t,length(Tev(Tev(:,1)<t(end),1))*2); 
    % For very long exposure profiles (e.g., FOCUS profiles), we now
    % automatically generate a larger time vector for IT (twice the number
    % of the relevant points in the scenario). For SD, it affects the
    % stepsize when not breaking the time vector.
    if break_time == 0 && size(Tev,1) > 2
        InitialStep = t_rem(end)/(10*min_t); % initial step size
        MaxStep     = t_rem(end)/min_t; % maximum step size
        % For the ODE solver, when we have a time-varying exposure set
        % here, we base minimum step size on min_t. Small stepsize is a
        % good idea for pulsed exposures; otherwise, stepsize may become so
        % large that a concentration change is missed completely. When we
        % break the time vector, limiting step size is not needed.
    end      
end
% For SD, we don't need a long time vector as survival is calculated by the
% ODE solver. However, for IT, we do need it to find the maximum of damage!
% In principle, this is not needed for constant exposure, but perhaps still
% a good idea as someone might modify derivatives in some way ...
if glo.sel == 2 && length(t) < min_t % for individual tolerance, 
    t = unique([t;(linspace(t(1),t(end),min_t))']); % make sure there are at least min_t points
end

% For SD, we can use a single threshold, which is then copied to the
% separate thresholds. This is now done here rather than in derivatives,
% which saves a bit of time. For IT, it makes little sense to have a joint
% threshold, but someone might like to try it anyway ...
if isfield(par,'mi') % then use one threshold
    par.mii(1) = par.mi(1);  % copy threshold to immobility threshold
    par.mid(1) = par.mi(1);  % copy threshold to death threshold
end

T = Tev(:,1); % time vector with events
if T(end) > t(end) % scenario may be longer than t(end)
    T(T>t(end)) = []; % remove all entries that are beyond the last time point
    % this may remove one point too many, but that will be added next
end
if T(end) < t(end) % scenario may (now) be shorter than we need
    T = cat(1,T,t(end)); % then add last point from t
end
% always use the ODE solver, and for now do NOT use the events function

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
options = odeset(options,'RelTol',RelTol,'AbsTol',AbsTol,'InitialStep',InitialStep,'MaxStep',MaxStep); % add an events function
% Note: setting tolerances is pretty tricky. For some cases, tighter
% tolerances are needed but not for others. For ODE45, tighter tolerances
% seem to work well, but not for ODE15s.

t = unique([T;t;(T(1:end-1)+T(2:end))/2]); % combine T, t, and halfway-T into new time vector
% this hopefully prevents the ODE solver from missing exposure pulses, and
% it makes sure that, for break_time=1, we never have a temporary time
% vector of length 2.

if break_time == 0 % don't break up the time vector
    
    switch stiff(1)
        case 0
            [tout,Xout] = ode45(@derivatives,t,X0,options,par,c,glo);
        case 1
            [tout,Xout] = ode113(@derivatives,t,X0,options,par,c,glo);
        case 2
            [tout,Xout] = ode15s(@derivatives,t,X0,options,par,c,glo);
    end
    
else % break up the time vector and apply the ODE-solver piecewise
    
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
        % [tout_tmp,Xout_tmp,TE,~,~] = ode45(@derivatives,t_tmp,X0,options,par,c,glo); % when using an events function ...
        switch stiff(1)
            case 0
                [tout_tmp,Xout_tmp] = ode45(@derivatives,t_tmp,X0,options,par,c,glo);
            case 1
                [tout_tmp,Xout_tmp] = ode113(@derivatives,t_tmp,X0,options,par,c,glo);
            case 2
                [tout_tmp,Xout_tmp] = ode15s(@derivatives,t_tmp,X0,options,par,c,glo);
        end
        
        Xout_tmp = max(0,Xout_tmp);  % in extreme cases, states can become ever so slightly negative
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
t = tout; % this is now the new t (to be compared to t_rem) that matches Xout

if isempty(TE) || all(TE == 0) % if there is no event caught
    TE = +inf; % return infinity
end

%% Output mapping
% _Xout_ contains a row for each state variable. It can be mapped to the
% data. If you need to transform the model values to match the data, do it
% here. 
%
% This mapping is here used for getting Ci and Di correct in case of fast
% kinetics or fast damage repair, and for IT to caclulate the output
% probabilities for the different states. For SD, these are already
% calculated in derivatives. For the immobility case, it looks like an ODE
% version is always needed for SD, since we are dealing with a system of 4
% ODEs. Damage may be calculated analytically, but that still leaves 3
% linked ODEs for survival. Therefore, the calculation below is mainly for
% IT.

Di   = Xout(:,glo.locD); % take the correct state variable for scaled damage
Ci   = Xout(:,glo.locC); % take the correct state variable for body residues
 
% Calculate external concentration. This is not needed in all cases (only
% for fast TK or for fast repair AND damconfig=3). Here, I just calculate
% external concentration for fast kinetics or fast damage dynamics.
if any(glo.fastslow==1)
    c_v = c * ones(length(t),1); % if no exposure profile is specified, simply copy c over all time points
    if glo.timevar == 1 % if we are warned that we have a time-varying concentration ...
        c_v = read_scen(-3,c,t,glo); % use read_scen to derive actual exposure concentration vector
        % the -3 lets read_scen know we are calling for fast/slow kinetics (and need a conc. vector)
    end
end

if glo.fastslow(1) == 1 % for fast toxicokinetics
    Ci    = c_v; % fast kinetics: internal concentration equals external concentration
    Ci(1) = X0(glo.locC); % make initial one Ci0 (especially needed for IT)
    Xout(:,glo.locC) = Ci; % replace internal concentrations in output vector
end
if glo.fastslow(2) == 1 % is we assume that damage repair is infinitely fast
    if glo.damconfig == 3 % if switch is set to 3 ...
        Di  = c_v; % scaled damage equals EXternal concentration
    else
        Di  = Ci; % scaled damage equals internal concentration
    end
    Xout(:,glo.locD) = Di; % replace damage in output vector
end

if glo.sel == 2 % individual tolerance, log-logistic distribution is used
    
    mii  = par.mii(1);       % median of threshold distribution for immobility
    mid  = par.mid(1);       % median of threshold distribution for death
    hb   = par.hb(1);        % background hazard rate
    mii  = max(mii,1e-100);  % make sure that the threshold is not exactly zero ...
    mid  = max(mid,1e-100);  % make sure that the threshold is not exactly zero ...

    % We need rules to select which state variable (Ci or Di) to use for which
    % effect (immobility or death).
    Did = Di; % always use damage for death mechanism
    if glo.damconfig ~= 2 % if switch is set to 1 or 3 ...
        Dii = Ci;  % use internal concentration for immobilisation
    else
        Dii = Did; % use damage Di for immobilisation as well as death
    end
    
    if isfield(par,'Fs') % assume that same beta holds for death and immobility response
        Fs     = max(1+1e-6,par.Fs(1)); % fraction spread of BOTH threshold distributions
        beta_d = log(39)/log(Fs);  % shape parameter for logistic from Fs
        beta_i = beta_d;           % take same beta for immobility as for death
    else % use separate ones
        Fsi    = max(1+1e-6,par.Fsi(1)); % fraction spread for immobility
        Fsd    = max(1+1e-6,par.Fsd(1)); % fraction spread for death
        beta_i = log(39)/log(Fsi); % shape parameter for logistic from Fs
        beta_d = log(39)/log(Fsd); % shape parameter for logistic from Fs
    end
    
    % Calculate an additional damage vector maxDid that does not decrease
    % over time (dead animals dont become alive).
    maxDid = Did; % copy the vector Did to maxDid
    ind    = find([0;diff(maxDid)]<0,1,'first'); % first index to places where Dw has decreased
    while ~isempty(ind) % as long as there is a decrease somewhere ...
        maxDid(ind:end) = max(maxDid(ind:end),maxDid(ind-1)); % replace every later time with max of that and previous point
        ind = find([0;diff(maxDid)]<0,1,'first'); % any decrease left?
    end
    
    % probabilities to be in each state
    pD = 1 ./ (1+(maxDid/mid).^-beta_d); % dead probability using maxDid
    pA = min(1 ./ (1+(Dii/mii).^beta_i) , 1-pD); % active probability using Dii
    pI = 1 - pA - pD;                   % immobility probability from the difference
    
    % Note: if we want to exclude the possibility for recovery (e.g., to
    % calculate worst-case EPx), the best way may be to calculate pH here
    % with maxDid. That would need some switch, though.
    
    % correct for background survival
    Sb = exp(-hb*t); % background survival probability    
    pD = pD + (1-Sb) .* (pA + pI);
    pA = pA .* Sb;
    pI = pI .* Sb;
    
    % replace correct states by newly calculated probabilities
    Xout(:,[glo.loc_a glo.loc_i glo.loc_d glo.loc_id]) = [pA pI pD pI+pD]; 
    
else % for stochastic death, derivatives calculates everything apart from the sum I+D (and Dw for fast kinetics)
    
    Xout(:,glo.loc_id) = Xout(:,glo.loc_i) + Xout(:,glo.loc_d); % sum probabilities death and immobility
    
end

[~,loct] = ismember(t_rem,t);   % find where the requested time points are in the long Xout
Xout     = Xout(loct,:);        % only keep the ones we asked for

%% Events function
% This subfunction catches the 'events': in this case, it looks for the
% damage/concentration level where the thresholds are exceeded (this only
% makes sense for SD; for IT, the effects are calculated in call_deri after
% the ODE solver).
%
% Note that the eventsfun has the same inputs, in the same sequence, as
% <derivatives.html derivatives.m>.

function [value,isterminal,direction] = eventsfun(t,X,par,c,glo)

% NOTE: AT THIS MOMENT, THE EVENTS FUNCTION IS NOT USED. THIS NEEDS SOME
% THOUGHT HOW TO STRUCTURE THIS PROPERLY WITH THE VARIOUS OPTIONS FOR FAST
% KINETICS. HOWEVER, CODE ALSO WORKS FINE WITHOUT THIS.

if glo.sel == 2 || any(glo.fastslow == 1) % for individual tolerance or fast kinetics, this is not helpful
    value      = 0;
    isterminal = 0;
    direction  = 0;
    return % for IT, we can return immediately
end

if isfield(par,'mi') % same threshold
    mii  = par.mi(1);  % threshold for immobility
    mid  = par.mi(1);  % threshold for death
else % separate threshold
    mii  = par.mii(1);  % threshold for immobility
    mid  = par.mid(1);  % threshold for death
end

nevents = 2; % number of events that we try to catch
value   = zeros(nevents,1); % initialise with zeros

% We need rules to select which state variable (Ci or Di) to use for which
% effect (immobility or death).
Did = X(glo.locD); % use damage for death
value(1) = Did - mid; % thing to follow is death damage minus threshold

if glo.damconfig ~= 2 % if switch is set to 1 or 3 ...
    Dii = X(glo.locC); % use internal concentration as damage for immobilisation
else
    Dii = Did; % use damage Di for immobilisation as well
end
value(2) = Dii - mii; % thing to follow is immob. damage minus threshold
    
isterminal = zeros(nevents,1); % do NOT stop the solver at an event
direction  = zeros(nevents,1); % catch ALL zero crossing when function is increasing or decreasing
