%% BYOM function call_deri.m (calculates the model output)
%
%  Syntax: [Xout,XE,Xout2,zvd,,stdDEB_start] = call_deri(t,par,X0v,glo)
%
% This function calls the ODE solver to solve the system of differential
% equations specified in <derivatives.html derivatives.m>. This version is
% specific for the standard DEB model with TKTD module. It runs through the
% embryo to tabulate the relationship between egg costs and scaled reserve
% density at birth, and start the 'real' run at birth. As input, it gets:
%
% * _t_   the time vector
% * _par_ the parameter structure
% * _X0v_   a vector with initial states and one concentration (scenario number)
% * _glo_ is the structure with information (normally global)
%
% The output _Xout_ provides a matrix with time in rows, and states in
% columns. This function calls <derivatives.html derivatives.m>. The
% optional output _XE_ is a matrix with the information on the events that
% are caught (matrix with time of event, event nr, and states at events).
% The events function is set up to catch discontinuities. It should be
% specified according to the problem you are simulating. 
%
% * Author: Tjalling Jager
% * Date: June 2022
% * Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2023, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function [Xout,XE,Xout2,zvd,stdDEB_start] = call_deri(t,par,X0v,glo)

% These outputs need to be defined, even if they are not used
Xout2    = []; % additional uni-variate output, not used in this example
zvd      = []; % additional zero-variate output, not used in this example

% These outputs need to be defined, if we return before the end of the function
Xout         = [];
XE           = [];
stdDEB_start = [];

warning('off','backtrace') % at this moment, several warnings can be produced, and we don't need all the extra info

%% Check if we can immediately return!
% This is a quick check whether puberty and/or metamorphosis occur before
% birth. This would lead to problems in the current version of this
% function, as the events function will then always stop before birth, so
% birth is never found. For now, I will simply return with an empty output
% matrix. In transfer, that means that the minloglik will be made +inf.
% Added a check whether the maturity threshold for metamorphosis is higher
% than for puberty. This is not really a problem, since the code further on
% will automatically put metamorphosis at puberty, but it messes up the
% CIs. Also added EHb=0, which might happen when profiling and does not
% lead to meaningful results.

if par.EHp(1) <= par.EHb(1) || (par.EHj(1) <= par.EHb(1) && par.EHj(1) > 0) || par.EHj(1) > par.EHp(1) || par.EHb(1) == 0
    return
end

%% Initial settings
% This part extracts optional settings for the ODE solver that can be set
% in the main script (defaults are set in prelim_checks). Note that the
% entries in X0mat are NOT used for initial values. The initial values are
% calculated from the model parameters, and simulated for birth or a
% specific length L0.
% 
% This package will always use the ODE solver; therefore the general
% settings in the global glo for the calculation (useode and eventson) are
% removed here. This packacge also always uses the events function, so the
% option eventson is not used. Furthermore, simplefun.m is removed.

stiff      = glo.stiff; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
min_t      = 500; % minimum length of time vector (affects ODE stepsize only, when needed)
break_time = glo.break_time; % break time vector up for ODE solver (1) or don't (0)
names_sep  = glo.names_sep;

if length(stiff) == 1
    stiff(2) = 1; % by default: normally tightened tolerances
end

c  = X0v(1); % the concentration (or scenario number)
X0 = X0v(2:end); % these are the intitial states for a scenario
% Start with a check on time-varying exposure scenarios. This is a
% time-saver! The isfield and ismember calls take some time that rapidly
% multiplies as derivatives is called many times.
glo.timevar = [0 0]; % flag for time-varying exposure (second element is to tell read_scen which interval we need)
if isfield(glo,'int_scen') && ismember(c,glo.int_scen) % is c in the scenario range of the global?
    glo.timevar = [1 0]; % tell derivatives that we have a time-varying treatment
end

% Deal with fitting multiple data sets with common parameters
i_d     = 0;
E0_calc = glo.E0_calc(1,:); % use first row by default
if glo.timevar(1) == 1 % this requires scenarios to be used!
    if c >= 100 && ~isempty(names_sep) % then we have more data sets, and more separate parameters per set!
        i_d = floor(c/100); % extract the data set number from the treatment identifier
        for i_sep = 1:length(names_sep) % run through extra parameter names for separate sets
            par.(names_sep{i_sep}) = par.([names_sep{i_sep},num2str(i_d)]); % copy extra parameter level to par
        end
        if size(glo.E0_calc,1)>1            % if multiple rows have been defined
            E0_calc = glo.E0_calc(i_d+1,:); % use row for this study
        end
    end
end

XE = nan(4,length(X0)+2); % matrix to store the life-cycle events (birth, metamorphosis, puberty and starting length for the simulation)
% NOTE: good to initialise here, as there can be cases where a life-cycle
% event is not found, and then XE is undefined as output.

%% Prepare for calculations
% This part calls the ODE solver to calculate the output (the value of the
% state variables over time). There is generally no need to modify this
% part. The solver ode45 generally works well. For stiff problems, the
% solver might become very slow; you can try ode15s instead.

t     = t(:);   % force t to be a row vector (needed when useode=0)
t_rem = t;      % remember the original time vector (as we might add to it)

% The code below is meant for discontinous time-varying exposure. It will
% also work for constant exposure. Breaking up will not be as effective for
% piecewise polynomials (type 1), and slow, so better turn that off
% manually when you try splining. The code will run each interval between
% two exposure 'events' (in Tev) separately, stop the solver, and restart.
% This ensures that the discontinuities are no problem anymore.

TE  = 0; % dummy for time of events in the events function
Tev = [0 c]; % exposure profile events setting: without anything else, assume it is constant

InitialStep = max(t)/100; % specify initial stepsize
MaxStep     = max(t)/10;  % specify maximum stepsize
% For constant concentrations, and when breaking the time vector, we can
% use a default; Matlab uses as default the length of the time vector
% divided by 10.

if glo.timevar(1) == 1 % time-varying concentrations?
    % if it exists: use it to derive time vector for the events in the exposure scenario    
    Tev = read_scen(-2,c,-1,glo); % use read_scen to derive exposure concentration events
    % the -2 lets read_scen know we need events, the -1 that this is also needed for splines!
    
    min_t = max(min_t,length(Tev(Tev(:,1)<t(end),1))*2); 
    % For very long exposure profiles (e.g., FOCUS profiles), we now
    % automatically generate a larger time vector (twice the number of the
    % relevant points in the scenario). This is only used to set step size
    % when not using break_time, and for locating maximum length under
    % no-shrinking.
    
    if break_time == 0 && size(Tev,1) > 2
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
% 
% NOTE: this is a problem for stdDEB with the AmP parameters. The AmP
% entries do NOT consider a brood pouch. Adding a delay here makes the AmP
% values close to useless!
bp = 0; % by default, no brood-pouch delay
if isfield(glo,'Tbp') && glo.Tbp > 0 % if there is a global specifying a brood-pouch delay ...
    tbp = t(t>glo.Tbp)-glo.Tbp; % extra times needed to calculate brood-pounch delay
    t   = unique([t;tbp]);      % add the shifted time points to account for brood-pouch delay
    bp  = 1;                    % signal rest of code that we need brood-pouch delay
end

% When an animal cannot shrink in length, we need a long time vector as we
% need to catch the maximum length over time. Using an entry in the events
% function may be possible as alternative.
if glo.len == 2 && length(t) < min_t % make sure there are at least min_t points
    t = unique([t;(linspace(t(1),t(end),min_t))']);
end

% specify options for the ODE solver
options = odeset; % start with default options
% Note: since odeset takes considerable time, it is smart to set all
% options in ONE call below.

% This needs further study ... to find optimal settings for the solvers
optionsE = odeset(options,'Events',@eventsfun,'RelTol',1e-9,'AbsTol',1e-9); % specify MORE tightened tolerances for embryo
% Note: at this point, the embryo always uses ode45. I think the stiff
% solvers are only needed for specific cases with toxicant stress (e.g.,
% fast kinetics and toxicant-induced shrinking).
switch stiff(2)
    % separate one for embryos: ode45 is fine, with tight tolerances
    case 1 % for ODE15s, slightly tighter tolerances seem to suffice (for ODE113: not tested yet!)
        options = odeset(options,'Events',@eventsfun,'RelTol',1e-4,'AbsTol',1e-7,'InitialStep',InitialStep,'MaxStep',MaxStep); % specify tightened tolerances
    case 2 % somewhat tighter tolerances ...
        options = odeset(options,'Events',@eventsfun,'RelTol',1e-5,'AbsTol',1e-8,'InitialStep',InitialStep,'MaxStep',MaxStep); % specify tightened tolerances
    case 3 % for ODE45, very tight tolerances seem to be necessary
        options = odeset(options,'Events',@eventsfun,'RelTol',1e-9,'AbsTol',1e-9,'InitialStep',InitialStep,'MaxStep',MaxStep); % specify MORE tightened tolerances
end
% Note: setting tolerances is pretty tricky. For some cases, tighter
% tolerances are needed but not for others. For ODE45, tighter tolerances
% seem to work well, but not for ODE15s.

flag_skip = 0;
if isfield(glo,'stdDEB_start') && ~isempty(glo.stdDEB_start) % then we have already saved all the proper starting values etc.
    % So extract the saved values for THIS data set. Note that stdDEB_start
    % has separate rows for different data sets!
    glo.E0mat = glo.stdDEB_start{i_d+1}.E0mat;
    glo.Lmat  = glo.stdDEB_start{i_d+1}.Lmat;
    X0        = glo.stdDEB_start{i_d+1}.X0;
    flag_skip = 1; % flag to skip the initial things
end

spAm = par.spAm(1); % max. surface-specific assimilation rate (J/(cm2 d))
v    = par.v(1);    % energy conductance (cm/d)
dE   = glo.dV; % g/cm3; assume density of reserve equals that of structure (which seems to be true in AmP)
wE   = 23.9;   % g/C-mol; chemical indices times molecular weight
% based on general composition of reserve: 1*12+1.8*1+0.5*16+0.15*14
muE  = 550e3;  % J/C-mol; chemical potential of reserve
sEm  = spAm/v; % max reserve density
omegaV = sEm * wE/(dE*muE); % contribution of reserve to volume (or wwt)
% Note: for now, do this always. It is needed when we use wwt (or dwt) as
% output, but also when we use egg or body weight (dry or wet) as
% additional zero-variate output.

% disp(['Reserve wt as fraction of total body wt: ',num2str(omegaV/(1+omegaV))])
% error

%% Embryo stage
% First run through the embryo stage. I fill a table with results for
% various values of E0. This allows me to start the runs from birth with
% the correct state values. Furthermore, it provides E0 as function of e
% (or f) that is needed to calculate egg costs. Furthermore, it gives Lb
% that is needed for the abj model. Note that the table will now be filled
% with values that lead to roughly e=0.3-1.5 at birth. This seems wide
% enough for most relevant cases. Interpolation works fine as the
% relationship between E0 and e is quite linear in most cases
% (extrapolation may thus not be too bad either).

if flag_skip == 0

    % Unpack some of the parameters as we need them here
    TA   = par.TA(1);   % Arrhenius temperature
%     spAm = par.spAm(1); % max. surface-specific assimilation rate (J/(cm2 d))
    spM  = par.spM(1);  % volume-specific somatic maintenance costs (J/(cm3 d))
%     v    = par.v(1);    % energy conductance (cm/d)
    EG   = par.EG(1);   % volume-specific costs for growth (J/cm3)
    kap  = par.kap(1);  % allocation fraction to soma (-)
    EHb  = par.EHb(1);  % maturity level at birth
    f    = par.f(1);    % scaled food density (-)
    % NOTE: if f varies between treatments, this needs to be included here as
    % well (if we want f to change the start situation for the neonate).

    % Temperature corrections
    FT   = exp(TA/glo.Tref - TA/glo.T); % factor by which to modify rates/times
    spAm = spAm * FT;   % modify specific assimilation
    spM  = spM  * FT;   % modify volume-specific maintenance
    v    = v    * FT;   % modify energy conductance

    % -------------------------------------------------------------------------
    % Crude estimation for foetal development
    Tb   = (EHb * (kap/(1-kap)) * 27/(EG*v^3))^(1/3); % minimum estimate for time at birth
    Etst = EG * (v^3)/(27*kap) * Tb^3; % this is how much reserve is used
    te   = [0;max(1,Tb*4)]; % time vector for embryo; take quite a bit more as eggs develop slower than a foetus without maintenance
    X0e  = [0.75*Etst;0;1e-6]; % initial states for fresh egg (note that L0 cannot be zero)
    % Use 0.75*Etst as 0.25*Etst is added in while loop below

    glo.Lmat   = [NaN NaN NaN NaN]; % matrix with length at birth, metamorphosis, puberty and starting length (if L0 is used)
    glo.E0mat  = []; % initialise to prevent error
    cE         = 0; % for embryo, assume no stress ... this only runs scenario 0 for the embryo ...
    % CHECK! User, make sure that a scenario c=0 makes sense in derivatives!

    % -------------------------------------------------------------------------
    % First take big steps with E0 to find a point where there certainly is birth
    IE = 4; % just to get the while loop started
    i  = 0; % keep track of how many rounds we take
    while isempty(IE) || IE(end) ~= 1  
        X0e(1) = X0e(1) + 0.25*Etst;
        [~,~,~,YE,IE] = ode45(@derivatives_pre,te,X0e,optionsE,par,cE,glo); % simulate embryo
        i = i + 1;
        if i > 20 % this is just to avoid getting stuck
            if isempty(IE)
                error('Birth cannot be found. Check if time vector needs to be longer')
            else
                % Note: there are cases when kappa->1 and EHb->0 that lead
                % to weird things. With much higher egg costs, there will
                % be birth at some point, but I think this can never be a
                % meaningful result. It will then stop with IE=4.
                warning('(call_deri) Birth cannot be found; the parameter set is likely nonsense. Skipping parameter set.')
                return % the current parameterisation does not make sense, so we can return already
            end
        end
    end
    e = YE(end,1) * v / (spAm * (YE(end,glo.locL))^3); % reserve density at birth
    % Note that YE(end,1) is E and YE(end,glo.locL) is L

    % -------------------------------------------------------------------------
    % Estimate a step size that should lead to steps of e of 0.1
    X0e_tmp    = X0e; % copy start vector
    X0e_tmp(1) = X0e(1) + 0.01*Etst; % increase egg costs a little bit
    [~,~,~,YE1,~] = ode45(@derivatives_pre,te,X0e_tmp,optionsE,par,cE,glo); % simulate embryo
    e1 = YE1(end,1) * v / (spAm * (YE1(end,glo.locL))^3); % reserve density at birth

    stepE = Etst*0.01*0.08/(e1-e);  % step to reduce/increase E0 with for testing
    % this aims for steps of e of 0.08, assuming linearity
    
    E0mat_down      = nan(50,5); % predefine matrix with NaNs
    E0mat_up        = nan(50,5); % predefine matrix with NaNs
    E0mat_down(1,:) = [X0e(1) e YE(end,:)]; % collect the first results in the matrix
    E0mat_up(1,:)   = [X0e(1) e YE(end,:)]; % collect the first results in the matrix
    if e > 0.7
        X0e(1) = X0e(1) - stepE; % reduce E0 by x%
        stepEi = stepE; % make sure it is defined, before hitting while loop
    else % when e is already quite small, we need a smaller stepsize
        X0e(1) = X0e(1) - 0.5*stepE; % reduce E0 by x%
        stepEi = 0.5 * stepE; % make sure it is defined, before hitting while loop
    end

    % -------------------------------------------------------------------------
    % First go to lower values of e and E0
    i = 1; % keep track of how many entries we have
    j = 0; % keep track of how many times we have decrease the step size
    while e > 0.3 && j < 3 % continue till e is low enough or we have decreased step size enough
        [~,~,~,YE,IE] = ode45(@derivatives_pre,te,X0e,optionsE,par,cE,glo); % simulate embryo
        if isempty(IE) || IE(end) ~= 1
            % NOTE: assume that we have not found birth because the reserve
            % runs out before maturity at birth is reached. That will stop
            % the solver with an event IE=4.
            if IE(end) ~= 4 % this should not happen
                error('time vector may not be not long enough!')
            else % we can continue the while loop, with smaller step
                X0e(1) = X0e(1) + stepEi; % redo last step
                stepEi = 0.5 * stepEi;    % make step smaller
                X0e(1) = X0e(1) - stepEi; % take a new, smaller, step
                j      = j + 1; % count another decrease in step size
            end
        else
            e = YE(end,1) * v / (spAm * (YE(end,glo.locL))^3); % reserve density at birth
            % Note that YE(end,1) is E and YE(end,glo.locL) is L

            i = i + 1; % increase counter
            E0mat_down(i,:) = [X0e(1) e YE(end,:)]; % collect the results in the matrix
            % calculate new E0 based on previous result!
            stepEi = 0.08*(E0mat_down(i-1,1)-E0mat_down(i,1))/(E0mat_down(i-1,2)-E0mat_down(i,2));
            if e > 0.7 % reduce E0 by an amount aimed at a step of 0.1 in e
                X0e(1) = X0e(1) - stepEi; % take a new step
            else % when e is quite small, we need a smaller stepsize
                stepEi = 0.5 * stepEi;    % make step smaller
                X0e(1) = X0e(1) - stepEi; % take a new step
            end
        end
    end

    % NOTE: there is a critical value for E0 below which an animal cannot reach
    % birth. Therefore, there is also a lower boundary for e at birth. This may
    % be pretty high, intuitively, for accelerating species (since their Lm,
    % while they are embryo, is pretty low). We estimate it with e_crit below.
    % The algorithm above should have the lowest e close to this value.
    e_crit = E0mat_down(i,2+glo.locL) / (kap * spAm / spM); % e_crit is approx last l_crit (Lb/Lm)

    if min(E0mat_down(:,2)) > 0.4 && min(E0mat_down(:,2))/e_crit > 1.10
        error('Range for scaled reserve density (e) versus egg costs (E0) does not seem to go low enough!')
        % This should not happen anymore.
    end

    % -------------------------------------------------------------------------
    % Next, go to higher values of e and E0, if needed
    e      = E0mat_up(1,2); % start with scaled reserve density that started the previous run
    i      = 1; % keep track of how many entries we have
    if e > 0.7
        X0e(1) = E0mat_up(1,1) + stepE; % increase E0
    else % when e is quite small, we need a smaller stepsize
        X0e(1) = E0mat_up(1,1) + 0.5*stepE; % increase E0
    end
    while e < 1.5 || isnan(e) % continue till e is high enough
        [~,~,~,YE,IE] = ode45(@derivatives_pre,te,X0e,optionsE,par,cE,glo); % simulate embryo
        if isempty(IE) || IE(end) ~= 1
            error('find out why birth cannot be found!') % this should not happen since we are increasing E0
        else
            e = YE(end,1) * v / (spAm * (YE(end,glo.locL))^3); % reserve density at birth
            % Note that YE(end,1) is E and YE(end,glo.locL) is L

            i = i + 1; % increase counter
            E0mat_up(i,:) = [X0e(1) e YE(end,:)]; % collect the results in the matrix
            % calculate new E0 based on previous result!
            if e > 0.7 % reduce E0 by an amount aimed at a step of 0.1 in e
                X0e(1) = X0e(1) + 0.12*(E0mat_up(i-1,1)-E0mat_up(i,1))/(E0mat_up(i-1,2)-E0mat_up(i,2));
            else % when e is quite small, we need a smaller stepsize
                X0e(1) = X0e(1) + 0.06*(E0mat_up(i-1,1)-E0mat_up(i,1))/(E0mat_up(i-1,2)-E0mat_up(i,2));
            end

        end
    end
    % NOTE: this code assumes that the event we're stopping on is birth. If
    % there are other stopping events possible for the embryo, this may need to
    % be revised.

    if max(E0mat_up(:,2)) < 1.5
        error('Range for scaled reserve density (e) versus egg costs (E0) does not go high enough!')
        % I don't think this can ever happen
    end

    % -------------------------------------------------------------------------
    E0mat_down = E0mat_down(~isnan(E0mat_down(:,1)),:); % remove NaNs
    E0mat_up   = E0mat_up(~isnan(E0mat_up(:,1)),:);     % remove NaNs
    E0mat_tot  = [flip(E0mat_up(2:end,:),1) ; E0mat_down]; % combine matrices, in ordered version
    E0mat_tot(E0mat_tot(:,2)>2,:) = []; % remove the ones where e > 2

    % The tests below are tricky. Some parameter sets, especially with
    % acceleration, do not allow a wide range of egg costs. So the table
    % may be short. In some cases, this is problematic, and I try to catch
    % them below. However, this needs very careful reconsideration when f
    % is a function of time and/or set in derivatives. Note that warnings
    % are produced in several cases. When fitting or profiling, these
    % warnings may clogg up the screen ...
    if isempty(E0mat_tot)
        warning('(call_deri) The egg-costs table is empty! Skipping parameter set.')
        return % the current parameterisation does not make sense, so we can return already
    end
    if size(E0mat_tot,1) < 8 % the table is shorter than we expected ...
        if E0_calc(2) == 3 && any(glo.moa(1:2)==1) % then we calculate egg costs from e, and e is affected by the toxicant
            error('Range for scaled reserve density (e) versus egg costs (E0) is shorter than expected ...')
            % a warning, and skipping the parameter set, may be an option here as well.
            % Note: we don't know at this point which values of e will be
            % needed in derivatives, when fitting the treatments.
            % Therefore, there might still be considerable extrapolations
            % there.
        end
    end
    if all(E0_calc == 2) % if all are 2, then we always use f=1
        if min(E0mat_tot(:,2))>1 || max(E0mat_tot(:,2))<1 % then f=1 is not within the range of the table
            warning('(call_deri) The egg-costs table excludes f=1! Skipping parameter set.')
            return % the current parameterisation does not make sense, so we can return already
        end
    else % we use the actual f at some point
        if min(E0mat_tot(:,2))>f || max(E0mat_tot(:,2))<f % then f is not within the range of the table
            warning('(call_deri) The egg-costs table excludes the set f! Skipping parameter set.')
            return % the current parameterisation does not make sense, so we can return already
        end
        if E0_calc(1) == 1 && f < e_crit % error turned into a warning: we can safely skip this parameter set as it makes no sense
            warning('(call_deri) The requested f is too low to produce an egg that is able to finish its development! Skipping parameter set.')
            return % the current parameterisation does not make sense, so we can return already
            % Note: now, this check is only made when we use f to determine
            % starting values. We may consider also performing this test
            % when f is used to calculate egg costs (so E0_calc(2)==1).
            % This is not critical for running the model, but it is weird
            % to allow the mother to produce eggs that cannot hatch (and
            % egg costs are extrapolations).
        end
    end

    glo.E0mat = E0mat_tot(:,[1 2]); % this matrix is used as lookup table for egg costs in derivatives
    if E0_calc(1) == 1
        X0(1:3) = (interp1(E0mat_tot(:,2),E0mat_tot(:,3:5),f))'; % make new starting values for f in this treatment, starting at birth
    else
        X0(1:3) = (interp1(E0mat_tot(:,2),E0mat_tot(:,3:5),1))'; % make new starting values for f=1 for all treatment, starting at birth
    end
    XE(1,:)   = [0 1 X0']; % store the states at the event 'birth' at t=0

    X0(2)       = X0(2) * (1+1e-10); % increase maturity a tiny bit so the birth event is not caught again in the next step (might save some calculation time)
    glo.Lmat(1) = X0(glo.locL);      % make Lb available for acceleration in derivatives when needed, and let derivatives know we passed birth
    % E0        = (interp1(E0mat_tot(:,2),E0mat_tot(:,1),1)); % also have the new E0 (at f=1) as optional output!

end

%% Starting at a specified L0
% There may be cases where we don't want to start at birth, but rather at
% some length L0>Lb. This requires additional simulation to find the time
% at which that length is reached. Furthermore, we may pass metamorphosis
% and puberty along the way, so these would need to be caught as well.

if isfield(par,'Lw0') && flag_skip == 0 % then we want to start at a fixed length after birth, and need to calculate it
    if glo.len ~= 0 % only if the switch is set to 1 or 2!
        glo.Lmat(4) = -1 * glo.delM * par.Lw0(1); % physical length to volumetric length
    else % if switch is set to 0, assume parameter is initial wet weight!
        glo.Lmat(4) = -1 * (par.Lw0(1)/(1+omegaV*f))^(1/3); % wet weight to volumetric length
        % Note: I use f here, so I assume the time from birth to start of
        % experiment is spent at a food density f. It would be possible to
        % make it dwt as well, but then also consider changing the output
        % mapping below.
    end
    % Set starting length in glo.Lmat to minus value so the events function
    % knows we need to catch a specific length, and what that length is.
    
    if -1 * glo.Lmat(4) < glo.Lmat(1)
        error('The specified L0 should be larger than Lb.')
    end
    Lm = f*kap*spAm/spM; % maximum volumetric length at f
    if par.EHj(1) == 0 && -1 * glo.Lmat(4) > Lm
        warning('(call_deri) The L0 set cannot be reached at the current values of the basic model parameters!')
        return % the current parameterisation does not make sense, so we can return already
        % Note: when fitting, it may be better to simply return empty
        % matrices, rather than throw an error. However, this should not be
        % a common problem. When there is acceleration, we don't know what
        % Lm will be, at this point ... so the error will occur later.
    end

    kM   = spM/EG;
    Em   = spAm/v;
    g    = EG/(kap*Em);
    rB   = (1/3)*kM*g/(1+g);    
    Tend = log(1/(1-0.95))/rB; % if the juvenile would grow von B, this would be the time to reach 95% of the final size
    
    flag_L0 = 0; % flag to signal that specified length was found
    tout    = 0; % start with t=0
    while flag_L0 == 0 % loop until we flag that we're ready
        % Code in this loop is copied from the 'official' run below. There,
        % more comments are provided about what's going on.
        
        t_tmp = [tout;mean([tout,Tend]);Tend]; % take time vector that's left (3 elements to prevent the solver from returning all values)
        [tout_tmp,Xout_tmp,TE,YE,IE] = ode45(@derivatives_pre,t_tmp,X0(1:3),optionsE,par,cE,glo); % run the life cycle in the control
        
        % look if we have an stopping event (metamorphosis or puberty)
        if ~isempty(IE)
            if any(IE == 2) && isnan(glo.Lmat(2)) % we stopped at metamorphosis (for the first time!)
                glo.Lmat(2)     = YE(end,glo.locL); % remember the length at metamorphosis
                XE(2,:)         = [TE(end) 2 YE(end,:)]; % store the states at the event 'metamorphosis'
                Xout_tmp(end,2) = Xout_tmp(end,2) * (1+1e-10); % increase maturity a tiny bit so the event is not caught again in the next step
            end
            if any(IE == 3) && isnan(glo.Lmat(3)) % we stopped at puberty (for the first time!)
                glo.Lmat(3)     = YE(end,glo.locL); % remember the length at puberty
                XE(3,:)         = [TE(end) 3 YE(end,:)]; % store the states at the event 'puberty'
                Xout_tmp(end,2) = Xout_tmp(end,2) * (1+1e-10); % increase maturity a tiny bit so the event is not caught again in the next step
                if par.EHj(1)>0 && isnan(glo.Lmat(2)) % then we have puberty, but not had metamorphosis yet!
                    glo.Lmat(2) = glo.Lmat(3); % assume that puberty and metamorphis are at same point
                    par.EHj(1)  = par.EHp(1); % also modify the parameter so derivatives functions properly
                end
            end
            
            if any(IE == 4) % we stopped at specified length
                glo.Lmat(4)  = YE(end,glo.locL); % now set it to the positive value, so we won't stop again at this point
                XE(4,:)      = [0 4 YE(end,:) (X0(4:6))']; % store the states at the event 'specified length'
                X0(1:3)      = Xout_tmp(end,:)'; % update X0 for next round; the 'official' run
                X0(glo.locR) = 0; % reset any reproduction to zero
                flag_L0      = 1; % set flag to 1: we can stop here
            else
                tout = tout_tmp(end);
                X0(1:3)      = Xout_tmp(end,:)'; % update X0 for next round of the while loop
            end
        end
        
        if tout_tmp(end) == Tend && flag_L0 ~= 1
            error('Specified starting length was not found.')
        end
        
    end
end

if flag_skip == 0
    % Place the starting values in the output variable stdDEB_start. When we
    % fit the tox parameters, there is no need to recalculate them since they
    % only depend on the basic parameters!
    stdDEB_start.E0mat = glo.E0mat;
    stdDEB_start.Lmat  = glo.Lmat;
    stdDEB_start.X0    = X0;
else
    stdDEB_start = glo.stdDEB_start{i_d+1}; % just to make sure: make it (potential) output again
end

%% Regular run for juvenile/adult stages
% Run the ODE solver piece-wise across all stopping events (metamorphosis
% and puberty, and perhaps others as well). I initialise tout and Xout with
% NaNs. I think this will be faster than using 'cat', especially for long
% time vectors.

T = Tev(:,1); % time vector with events
if T(end) > t(end) % scenario may be longer than t(end)
    T(T>t(end)) = []; % remove all entries that are beyond the last time point
    % this may remove one point too many, but that will be added next
end
if T(end) < t(end) % scenario may (now) be shorter than we need
    T = cat(1,T,t(end)); % then add last point from t
end
% % If a lag time is used, it is a good idea to add it as an event as ODE45
% % does not like such a switch either.
% Tlag = par.Tlag(1);
% if Tlag > 0
%     T = unique([par.Tlag(1);T]);
% end

% Breaking up the time vector is faster and more robust when calibrating on
% data from time-varying exposure that contain discontinuities. However,
% when running through very detailed exposure profiles (e.g., those from
% FOCUS), it is probably best to use break_time=0. Testing indicates that
% breaking up will be much, much slower than simply running the ODE solver
% across the entire profile (and produces very similar output).

t = unique([T;t;(T(1:end-1)+T(2:end))/2]); % combine T, t, and halfway-T into new time vector
% this hopefully prevents the ODE solver from missing exposure pulses

tout    = nan(length(t)+10,1); % this vector will collect the time output from the solver (10 more than needed)
Xout    = nan(length(tout),length(X0)); % this matrix will collect the states output from the solver (10 more than needed)
ind_i   = 1; % index for where we are in tout and Xout
tout(1) = 0; % start tout at t=0, so the while loop starts properly
% Note: I initialise tout and Xout with more than the elements of t because
% the events will be added as well! The number is rather arbitrary but
% should be equal to, or more than, the number of events.

if break_time == 0

    while tout(ind_i) < t(end) % loop until we calculated the last time point in t

        t_tmp = [tout(ind_i);t(t>tout(ind_i))]; % take time vector that's left
        add_t = 0; % (re)set flag to zero; will be set to 1 when we add a point in t_tmp
        if length(t_tmp) == 2 % it t_tmp only has 2 elements ...
            t_tmp = [t_tmp(1);mean(t_tmp);t_tmp(2)]; % add the mean in the middle
            add_t = 1; % flag that we have added a time point
        end
        % NOTE: a short time vector may result as we're stopping on
        % life-cycle events with the events function.

        % call the correct ODE solver
        switch stiff(1)
            case 0
                [tout_tmp,Xout_tmp,TE,YE,IE] = ode45(@derivatives,t_tmp,X0,options,par,c,glo);
            case 1
                [tout_tmp,Xout_tmp,TE,YE,IE] = ode113(@derivatives,t_tmp,X0,options,par,c,glo);
            case 2
                [tout_tmp,Xout_tmp,TE,YE,IE] = ode15s(@derivatives,t_tmp,X0,options,par,c,glo);
        end

        % look if we have a stopping event (metamorphosis or puberty)
        if ~isempty(IE)
            % NOTE: there are two separate if statements, which is done to
            % allow the case where EHj=EHp. If they are exactly the same, the
            % events function returns two events!
            if any(IE == 2) && isnan(glo.Lmat(2)) % we stopped at metamorphosis (for the first time!)
                glo.Lmat(2)     = YE(end,glo.locL); % remember the length at metamorphosis
                XE(2,:)         = [TE(end) 2 YE(end,:)]; % store the states at the event 'metamorphosis'
                Xout_tmp(end,2) = Xout_tmp(end,2) * (1+1e-10); % increase maturity a tiny bit so the event is not caught again in the next step
            end
            if any(IE == 3) && isnan(glo.Lmat(3)) % we stopped at puberty (for the first time!)
                glo.Lmat(3)     = YE(end,glo.locL); % remember the length at puberty
                XE(3,:)         = [TE(end) 3 YE(end,:)]; % store the states at the event 'puberty'
                Xout_tmp(end,2) = Xout_tmp(end,2) * (1+1e-10); % increase maturity a tiny bit so the event is not caught again in the next step
                if par.EHj(1)>0 && isnan(glo.Lmat(2)) % then we have puberty, but not had metamorphosis yet!
                    glo.Lmat(2) = glo.Lmat(3); % assume that puberty and metamorphis are at same point
                    par.EHj(1)  = par.EHp(1); % also modify the parameter so derivatives functions properly
                end
                % NOTE: this last part is needed when fitting EHj. It may
                % coincide or even exceed EHp. This is problematic as, without
                % safeguard, acceleration will continue for ever (as maturity
                % does not increase beyond puberty). In the fit results, EHj
                % may still exceed EHp, but not in the model output.
            end
        end
        % Note: the last condition in the if statements is used to make sure
        % that we calculate Lj and Lp only once in the life time of an
        % individual. If we allow EH to decrease (rejuvenation), we may need to
        % reconsider this: the second Lj or Lp may be different from the first
        % one!

        % collect output in correct location
        if add_t == 1 % then we added a time point in between
            ind_rem = tout_tmp == t_tmp(2); % index for added time point
            tout_tmp(ind_rem)   = []; % looks a bit awkward, but we may have a stopping event in there as well ...
            Xout_tmp(ind_rem,:) = []; % ... so better look explicitly where tout_tmp == t_tmp(2)
        end
        % Note: this is not really needed as additional time points are weeded
        % out later anyway. However, if there are many of these cases in a
        % calculation, we may run out of pre-initialised Xout (and the
        % calculations will slow down).
        nt = length(tout_tmp); % length of the output time vector
        tout(ind_i:ind_i-1+nt)   = tout_tmp; % add tout_tmp to correct position in tout
        Xout(ind_i:ind_i-1+nt,:) = Xout_tmp; % add Xout_tmp to correct position in Xout

        ind_i = ind_i-1+nt;       % update ind_i for next round
        X0    = Xout_tmp(end,:)'; % update X0 for next round

        % Note: if an event is identified, the ODE solver will return an extra
        % time point with output, which is stored in tout and Xout. This is
        % fine as the code below will sieve out the requested time points of t.
        %
        % Note: the way ind_i is used, there will be no double time points in
        % tout and Xout.

    end

else % THIS COULD DO WITH MORE CHECKING, THOUGH IT LOOKS GOOD!

    % Transferring the interval to derivatives is a huge time saver! It is
    % not safe to use i as index to Tev since we may add Tlag to T (which
    % is not in Tev); so better calculate it explicitly. Using <find> may
    % be a bit slow, that is why this is done outside of the loop calling
    % the ODE solver.
    ind_Tev = zeros(length(T)-1,1);
    for i = 1:length(T)-1 % run through all intervals between events
        ind_Tev(i) = find(Tev(:,1) <= T(i),1,'last');
    end

    for i = 1:length(T)-1 % run through all intervals between events

        % create a new master time vector for this interval!
        T_tmp = [T(i);T(i+1)]; % start with start and end time for this period
        T_tmp = unique([T_tmp;t(t<T_tmp(2) & t>T_tmp(1))]); % and add the time points from t that fit in there
        glo.timevar(2) = ind_Tev(i); % tell derivatives which part of Tev we're in

        while tout(ind_i) < T_tmp(end) % loop until we calculated the last time point in T_tmp

            t_tmp = [tout(ind_i);T_tmp(T_tmp>tout(ind_i))]; % take time vector that's left
            add_t = 0; % (re)set flag to zero; will be set to 1 when we add a point in t_tmp
            if length(t_tmp) == 2 % if t_tmp only has 2 elements ...
                t_tmp = [t_tmp(1);mean(t_tmp);t_tmp(2)]; % add the mean in the middle
                add_t = 1; % flag that we have added a time point
            end
            % NOTE: a short time vector may result as we're stopping on
            % life-cycle events with the events function.

            % call the correct ODE solver
            switch stiff(1)
                case 0
                    [tout_tmp,Xout_tmp,TE,YE,IE] = ode45(@derivatives,t_tmp,X0,options,par,c,glo);
                case 1
                    [tout_tmp,Xout_tmp,TE,YE,IE] = ode113(@derivatives,t_tmp,X0,options,par,c,glo);
                case 2
                    [tout_tmp,Xout_tmp,TE,YE,IE] = ode15s(@derivatives,t_tmp,X0,options,par,c,glo);
            end

            % look if we have a stopping event (metamorphosis or puberty)
            if ~isempty(IE)
                % NOTE: there are two separate if statements, which is done to
                % allow the case where EHj=EHp. If they are exactly the same, the
                % events function returns two events!
                if any(IE == 2) && isnan(glo.Lmat(2)) % we stopped at metamorphosis (for the first time!)
                    glo.Lmat(2)     = YE(end,glo.locL); % remember the length at metamorphosis
                    XE(2,:)         = [TE(end) 2 YE(end,:)]; % store the states at the event 'metamorphosis'
                    Xout_tmp(end,2) = Xout_tmp(end,2) * (1+1e-10); % increase maturity a tiny bit so the event is not caught again in the next step
                end
                if any(IE == 3) && isnan(glo.Lmat(3)) % we stopped at puberty (for the first time!)
                    glo.Lmat(3)     = YE(end,glo.locL); % remember the length at puberty
                    XE(3,:)         = [TE(end) 3 YE(end,:)]; % store the states at the event 'puberty'
                    Xout_tmp(end,2) = Xout_tmp(end,2) * (1+1e-10); % increase maturity a tiny bit so the event is not caught again in the next step
                    if par.EHj(1)>0 && isnan(glo.Lmat(2)) % then we have puberty, but not had metamorphosis yet!
                        glo.Lmat(2) = glo.Lmat(3); % assume that puberty and metamorphis are at same point
                        par.EHj(1)  = par.EHp(1); % also modify the parameter so derivatives functions properly
                    end
                    % NOTE: this last part is needed when fitting EHj. It may
                    % coincide or even exceed EHp. This is problematic as, without
                    % safeguard, acceleration will continue for ever (as maturity
                    % does not increase beyond puberty). In the fit results, EHj
                    % may still exceed EHp, but not in the model output.
                end
            end
            % Note: the last condition in the if statements is used to make sure
            % that we calculate Lj and Lp only once in the life time of an
            % individual. If we allow EH to decrease (rejuvenation), we may need to
            % reconsider this: the second Lj or Lp may be different from the first
            % one!

            % collect output in correct location
            if add_t == 1 % then we added a time point in between
                ind_rem = tout_tmp == t_tmp(2); % index for added time point
                tout_tmp(ind_rem)   = []; % looks a bit awkward, but we may have a stopping event in there as well ...
                Xout_tmp(ind_rem,:) = []; % ... so better look explicitly where tout_tmp == t_tmp(2)
            end
            % Note: this is not really needed as additional time points are weeded
            % out later anyway. However, if there are many of these cases in a
            % calculation, we may run out of pre-initialised Xout (and the
            % calculations will slow down).
            nt = length(tout_tmp); % length of the output time vector
            tout(ind_i:ind_i-1+nt)   = tout_tmp; % add tout_tmp to correct position in tout
            Xout(ind_i:ind_i-1+nt,:) = Xout_tmp; % add Xout_tmp to correct position in Xout

            ind_i = ind_i-1+nt;       % update ind_i for next round
            X0    = Xout_tmp(end,:)'; % update X0 for next round

            % Note: if an event is identified, the ODE solver will return an extra
            % time point with output, which is stored in tout and Xout. This is
            % fine as the code below will sieve out the requested time points of t.
            %
            % Note: the way ind_i is used, there will be no double time points in
            % tout and Xout.

        end
    end
end

%% Output mapping
% _Xout_ contains a row for each state variable. It can be mapped to the
% data. If you need to transform the model values to match the data, do it
% here. 

% Here, translate from the state volumetric body weight to output in
% physical length, if needed. Furthermore, if no physical shrinking is
% possible for the species, the max over time is used.
if glo.len ~= 0 % only if the switch is set to 1 or 2!
    L  = Xout(:,glo.locL); % take body length from the correct location in Xout
    Lw = L/glo.delM;       % estimated physical body length    
    if glo.len == 2 % when animal cannot shrink in length (but does on volume!)
        Lw = cummax(Lw); % replace L by the maximum achieved over time
    end
    Xout(:,glo.locL) = Lw; % replace body weight by physical body length
else % then use volume or wet weight as measure
    L   = Xout(:,glo.locL); % take body length from the correct location in Xout
    E   = Xout(:,1); % take reserve in Joule from the correct location in Xout
    e   = E./(sEm*(L.^3)); % scaled reserve density
    wwt = (L.^3) .* (1+omegaV*e); % need to include reserve into the wwt!
    % FV  = (L.^3) ./ wwt; % fraction of wwt that is structure
    Xout(:,glo.locL) = wwt; % replace state variable body length by total wet weigth
end
% Note: we can easily get dwt as output by using dwt=glo.dV*wwt. But then
% also consider changing the use of par.Lw0 for the initial body size.

Xout(:,glo.locS) = max(0,Xout(:,glo.locS)); % make sure survival does not get negative
% In some case, it may become just a bit negative, which means zero.

if bp == 1 % if we use a brood-pouch delay ... 
    [~,loct] = ismember(tbp,tout); % find where the extra brood-pouch time points are in the long Xout
    Xbp      = Xout(loct,glo.locR); % only keep the brood-pouch ones we asked for
end
t = t_rem; % return t to the original input vector

% Select the correct time points to return to the calling function
[~,loct] = ismember(t,tout); % find where the requested time points are in the long Xout
Xout     = Xout(loct,:);     % only keep the ones we asked for

if bp == 1 % if we use a brood-pouch delay ... 
    [~,loct] = ismember(tbp+glo.Tbp,t); % find where the extra brood-pouch time points SHOULD BE in the long Xout
    Xout(:,glo.locR)    = 0;   % clear the reproduction state variable
    Xout(loct,glo.locR) = Xbp; % put in the brood-pouch ones we asked for
end

% % To obtain the output of the derivatives at each time point. The values in
% % dXout might be used to replace values in Xout, if the data to be fitted
% % are the changes (rates) instead of the state variable itself.
% % dXout = zeros(size(Xout)); % initialise with zeros
% for i = 1:length(t) % run through all time points
%     dXout(i,:) = derivatives(t(i),Xout(i,:),par,c,glo); 
%     % derivatives for each stage at each time
% end

%% Zero-variate data point(s)
% Code below can be used to calculate the model value for a zero-variate
% data point. The code now calculates maximum length (accounting for
% acceleration). It is active when zvd.* is specified in the main script
% as a 2-element row vector: value and SD. More/other zero-variate data
% points can be added to the zvd structure.
% 
% Note: it should not be a problem to specify more zero-variate output
% below than are used by the main script. Only the ones defined in the main
% script will be used in transfer. However, it might slow down the
% calculation a bit, so we should consider commenting out the parts that
% are not needed.

if ~isempty(glo.zvd) % if there are zero-variate data defined (see byom_bioconc_extra)
    zvd  = glo.zvd;     % copy zero-variate data structure to zvd
    spAm = par.spAm(1); % max. surface-specific assimilation rate (J/(cm2 d))
    spM  = par.spM(1);  % volume-specific somatic maintenance costs (J/(cm3 d))
    kap  = par.kap(1);  % allocation fraction to soma (-)
    EHj  = par.EHj(1);  % maturity level at metamorphosis (J)
    
    if EHj > 0 % when we have acceleration, need to include it
        Lb  = glo.Lmat(1); % length at birth
        Lj  = glo.Lmat(2); % length at metamorphosis
        del = Lj/Lb; % acceleration factor at metamorphosis and onwards
        spAm = spAm * del; % change surface-specific assimilation rate
    end

    Lm         = (kap*spAm/spM); % maximum volumetric length at f=1
    zvd.Lwm(3) = Lm/glo.delM; % add model prediction for max physical length as third value in zvd
    Vwm        = (Lm^3) * (1+omegaV*1); % maximum total wet weight at f=1
    % zvd.Vwm(3) = Vwm; % add model prediction for max. wet weight as third value in zvd
    zvd.Wdi(3) = glo.dV * Vwm; % add model prediction for maximum dry weight as third value in zvd
    % Note: no need for temperature correction, as maximum body size is not
    % affected. If you use a zero-variate data point for max body size,
    % watch out with fitting f.

    % egg size as potential zero-variate outputs
    E0 = (interp1(glo.E0mat(:,2),glo.E0mat(:,1),1)); % derive egg costs for f=1
    zvd.Vw0(3) = E0*omegaV/sEm; % egg wet weight (cm^3 or g_wwt)
    zvd.Wd0(3) = glo.dV * (E0*omegaV/sEm); % egg dry weight (g)

end

%% Events function
% This subfunction catches the 'events': in this case, it looks for the
% points where maturity exceeds the thresholds for birth, metamorphosis and
% puberty (and starting length). This function is specific for the standard
% DEB model.
%
% Note that the eventsfun has the same inputs, in the same sequence, as
% <derivatives.html derivatives.m>.

function [value,isterminal,direction] = eventsfun(t,X,par,c,glo)

EHb  = par.EHb(1);  % maturity level at birth (J)
EHj  = par.EHj(1);  % maturity level at metamorphosis (J)
EHp  = par.EHp(1);  % maturity level at puberty (J)

% NOTE: we could extend the events function to catch threshold as well, but
% I don't think this is really needed.

nevents     = 3;                % number of events that we try to catch
value       = zeros(nevents,1); % initialise with zeros
value(1)    = X(2) - EHb;       % thing to follow is maturity (state 2) minus threshold for birth
if EHj > 0                      % only if there is abj acceleration ...
    value(2) = X(2) - EHj;      % thing to follow is maturity (state 2) minus threshold for metamorphosis
end
value(3)    = X(2) - EHp;       % thing to follow is maturity (state 2) minus threshold for puberty

% Code below is here if we want to start at length L0>Lb
if ~isnan(glo.Lmat(4)) && glo.Lmat(4) < 0 % it is set negative when we need to find it, otherwise we don't need to stop
    nevents  = nevents + 1;
    value(4) = X(glo.locL) + glo.Lmat(4); % note the plus sign since the length we search for is set to its minus value!
end

isterminal  = ones(nevents,1);  % DO stop the solver at an event
% isterminal  = (isnan(glo.Lmat))'; % only where a Li is not found yet, we'll stop!
% % This code line may be needed when we allow rejuvenation; then there can
% % be multiple crossings of EHi.
direction   = zeros(nevents,1);   % catch ALL zero crossing when function is increasing or decreasing

if ~isnan(glo.Lmat(1)) % then we are NOT looking for egg costs and birth
    return % stop the events function
end

%% Special events section for embryos
% Catch cases where the embryo fails to be born because its energy is
% insufficient. No need to continue the run and possibly run into numerical
% issues.

% Unpack some of the parameters as we need them here
spAm = par.spAm(1); % max. surface-specific assimilation rate (J/(cm2 d))
spM  = par.spM(1);  % volume-specific somatic maintenance costs (J/(cm3 d))
kap  = par.kap(1);  % allocation fraction to soma (-)
v    = par.v(1);    % energy conductance (cm/d)
% Temperature corrections are not needed since the scalings for E and L are
% not affected (the correction cancels out).

e = X(1) * v / (spAm * (X(glo.locL))^3); % scaled reserve density
l = X(glo.locL) / (kap * spAm / spM);    % scaled length

value(nevents+1)      = e - l; % check whether e becomes lower than l
isterminal(nevents+1) = 1; % stop the run
direction(nevents+1)  = 0; % stop in both directions (though only one makes sense)
