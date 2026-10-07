%% BYOM function call_deri.m (calculates the model output)
%
%  Syntax: [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)
%
% This function calls the ODE solver to solve the system of differential
% equations specified in <derivatives.html derivatives.m>. As input, it
% gets:
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
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_debkiss.html>

%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)

% These outputs need to be defined, even if they are not used
Xout2    = []; % additional uni-variate output, not used in this example
zvd      = []; % additional zero-variate output, not used in this example

% NOTE: this file is modified so that data and output can be as body length
% whereas the state variable remains body weight. Also, there is a
% calculation to prevent shrinking in physical length (as shells, for
% example, do not shrink).

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

stiff = glo.stiff; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
min_t = 500; % minimum length of time vector (affects ODE stepsize only, when needed)

if length(stiff) == 1
    stiff(2) = 1; % by default: normally tightened tolerances
end

% Unpack the vector X0v, which is X0mat for one scenario
X0 = X0v(2:end); % these are the intitial states for a scenario

% % if needed, extract parameters from par that influence initial states in X0
% % start from specified initial size
% Lw0 = par.Lw0(1);  % initial body length (mm)
% X0(glo.locL) = glo.dV*((Lw0*glo.delM)^3); % recalculate to dwt, and put this parameter in the correct location of the initial vector

% start from the value provided in X0mat
if glo.len ~= 0 % only if the switch is set to 1 or 2! 
    % Assume that X0mat includes initial length; if it includes mass, no correction is needed here
    % The model is set up in masses, so we need to translate the initial length
    % (in mm body length) to body mass using delM. 
    X0(glo.locL) = glo.dV*((X0(glo.locL) * glo.delM)^3); % initial length to initial dwt
end

% % start at hatching, so calculate initial length from egg weight
% WB0   = par.WB0(1);     % initial dry weight of egg
% yVA   = par.yVA(1);     % yield of structure on assimilates (growth)
% kap   = par.kap(1);     % allocation fraction to soma
% WVb   = WB0 *yVA * kap; % dry body mass at birth
% X0(glo.locL) = WVb;     % put this estimate in the correct location of the initial vector

% % if needed, calculate model values for zero-variate data from parameter set
% if ~isempty(zvd)
%     zvd.ku(3) = par.Piw(1) * par.ke(1); % add model prediction as third value
% end

%% Calculations
% This part calls the ODE solver to calculate the output (the value of the
% state variables over time). There is generally no need to modify this
% part. The solver ode45 generally works well. For stiff problems, the
% solver might become very slow; you can try ode15s instead.

c     = X0v(1); % the concentration (or scenario number)
t     = t(:);   % force t to be a row vector (needed when useode=0)
t_rem = t;      % remember the original time vector (as we might add to it)

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

% When an animal cannot shrink in length, we need a long time vector as we
% need to catch the maximum length over time.
if glo.len == 2 && length(t) < min_t % make sure there are at least min_t points
    t = unique([t;(linspace(t(1),t(end),min_t))']);
end

% specify options for the ODE solver
options = odeset; % start with default options
% Note: since odeset takes considerable time, it is smart to set all
% options in ONE call below.

% This needs further study ...
switch stiff(2)
    case 1 % for ODE15s, slightly tighter tolerances seem to suffice (for ODE113: not tested yet!)
        options = odeset(options,'Events',@eventsfun,'RelTol',1e-4,'AbsTol',1e-7); % specify tightened tolerances
    case 2 % somewhat tighter tolerances ...
        options = odeset(options,'Events',@eventsfun,'RelTol',1e-5,'AbsTol',1e-8); % specify tightened tolerances
    case 3 % for ODE45, very tight tolerances seem to be necessary
        options = odeset(options,'Events',@eventsfun,'RelTol',1e-9,'AbsTol',1e-9); % specify MORE tightened tolerances
end
% Note: setting tolerances is pretty tricky. For some cases, tighter
% tolerances are needed but not for others. For ODE45, tighter tolerances
% seem to work well, but not for ODE15s.

% The code below uses an ODE solver to go through the time vector in one
% run. The DEBtox packages use code that can break up the time vector in
% sections. That would be especially useful if conditions change
% (non-smoothly) over time. But here, we keep things simple.

TE = 0; % dummy for time of events
% call the ODE solver with an events functions ... additional output arguments for events:
% TE catches the time of an event
switch stiff(1)
    case 0
        [tout,Xout,TE,~,~] = ode45(@derivatives,t,X0,options,par,c,glo);
    case 1
        [tout,Xout,TE,~,~] = ode113(@derivatives,t,X0,options,par,c,glo);
    case 2
        [tout,Xout,TE,~,~] = ode15s(@derivatives,t,X0,options,par,c,glo);
end

if isempty(TE) || all(TE == 0) % if there is no event caught
    TE = +inf; % return infinity
end

%% Output mapping
% _Xout_ contains a row for each state variable. It can be mapped to the
% data. If you need to transform the model values to match the data, do it
% here. 

% here, translate from the state body weight to output in length, if needed
if glo.len ~= 0 % only if the switch is set to 1 or 2!
    % Here, I translate the model output in body mass to estimated physical body length
    W = Xout(:,glo.locL); % take dry body weight from the correct location in Xout
    L = (W/glo.dV).^(1/3)/glo.delM; % estimated physical body length    
    
    if glo.len == 2 % when animal cannot shrink in length (but does on volume!)
        L = cummax(L); % replace L by the maximum achieved over time
    end
    Xout(:,glo.locL) = L; % replace body weight by physical body length
end

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

% % To obtain the output of the derivatives at each time point. The values in
% % dXout might be used to replace values in Xout, if the data to be fitted
% % are the changes (rates) instead of the state variable itself.
% % dXout = zeros(size(Xout)); % initialise with zeros
% for i = 1:length(t) % run through all time points
%     dXout(i,:) = derivatives(t(i),Xout(i,:),par,c,glo); 
%     % derivatives for each stage at each time
% end

%% Events function
% This subfunction catches the 'events': in this case, it looks for birth
% and the point where size exceeds the threshold for puberty. This function
% should be adapted to the problem you are modelling.
%
% Note that the eventsfun has the same inputs, in the same sequence, as
% <derivatives.html derivatives.m>.

function [value,isterminal,direction] = eventsfun(t,X,par,c,glo)

% Note: glo needs to be an input to the events function as well, rather
% than a global since the inputs for the events function must match those
% for derivatives.

Lwp      = par.Lwp(1); % shell length at puberty
if glo.len ~= 0 % only if the switch is set to 1 or 2!
    % Translate shell length at puberty to body weight
    WVp = glo.dV*(Lwp * glo.delM).^3; % initial length to initial dwt
end

nevents  = 1;          % number of events that we try to catch

value       = zeros(nevents,1);  % initialise with zeros
value(1)    = X(glo.locL) - WVp; % thing to follow is body mass (state 1) minus threshold

isterminal  = zeros(nevents,1); % do NOT stop the solver at an event
direction   = zeros(nevents,1); % catch ALL zero crossing when function is increasing or decreasing
