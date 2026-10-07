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
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_debkiss_paper.html>

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
% This part can be modified by the user to decide whether to calculate the
% model results using the ODEs in <derivatives.html derivatives.m>, or the
% analytical solution in <simplefun.html simplefun.m>. Further, the events
% function can be turned on here, initial values can be determined by a
% parameter, and zero-variate data can be calculated.

stiff = glo.stiff; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
min_t = 500; % minimum length of time vector (affects ODE stepsize only, when needed)

if length(stiff) == 1
    stiff(2) = 1; % by default: normally tightened tolerances
end

% Unpack the vector X0v, which is X0mat for one scenario
X0 = X0v(2:end); % these are the intitial states for a scenario

% if needed, extract parameters from par that influence initial states in X0

% % start from specified initial size
% Lw0 = par.Lw0(1);  % initial body length (mm)
% X0(2) = glo.dV*((Lw0*glo.delM)^3); % recalculate to dwt, and put this parameter in the correct location of the initial vector

% start from the value provided in X0mat
if glo.len ~= 0 % only if the switch is set to 1 or 2! 
    % Assume that X0mat includes initial length; if it includes mass, no correction is needed here
    % The model is set up in masses, so we need to translate the initial length
    % (in mm body length) to body mass using delM. Assume that size is first state.
    X0(1) = glo.dV*((X0(1) * glo.delM)^3); % initial length to initial dwt
end

% % start at hatching, so calculate initial length from egg weight
% WB0   = par.WB0(1);     % initial dry weight of egg
% yVA   = par.yVA(1);     % yield of structure on assimilates (growth)
% kap   = par.kap(1);     % allocation fraction to soma
% WVb   = WB0 *yVA * kap; % dry body mass at birth
% X0(2) = WVb;            % put this estimate in the correct location of the initial vector

% % if needed, calculate model values for zero-variate data from parameter set
% if ~isempty(zvd)
%     zvd.ku(3) = par.Piw(1) * par.ke(1); % add model prediction as third value
% end

%% Calculations
% This part calls the ODE solver (or the explicit model in <simplefun.html
% simplefun.m>) to calculate the output (the value of the state variables
% over time). There is generally no need to modify this part. The solver
% ode45 generally works well. For stiff problems, the solver might become
% very slow; you can try ode15s instead.

c     = X0v(1); % the concentration (or scenario number)
t     = t(:);   % force t to be a row vector (needed when useode=0)
t_rem = t;      % remember the original time vector (as we might add to it)

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

% % To obtain the output of the derivatives at each time point. The values in
% % dXout might be used to replace values in Xout, if the data to be fitted
% % are the changes (rates) instead of the state variable itself.
% % dXout = zeros(size(Xout)); % initialise with zeros
% for i = 1:length(t) % run through all time points
%     dXout(i,:) = derivatives(t(i),Xout(i,:),par,c); 
%     % derivatives for each stage at each time
% end

% the part below is highly modified for the specific calculations in the
% DEBkiss1 package, to reproduce the figures in the DEBkiss publication

if X0v(4) > 0 % there is an initial repro buffer, so we started with a fresh egg!
    Xout(t>TE,1) = NaN; % don't include body length after point of hatching
end

% next, calculate respiration from the body length data
Lw = Xout(:,glo.locL); % shell length from the output matrix
LV = Lw*glo.delM; % translate to volumetric length
JA = par.sJAm(1)*(LV.^2); % embryo assimilation rate for f=1
JA = JA * par.fB(1); % calculate with reduced f for the embryo
JM = par.sJM(1)*(LV.^3); % flux for maintenance
JV = par.yVA(1)*(par.kap(1)*JA-JM); % flux fixed in structural biomass
JD = JA - JV; % dissipation is assimilation minus what is fixed in growth

Resp = JD * 130; % proportionality is uL O2/mg dissipated assimilates (see paper)
Resp(t<glo.Tlag(1),:) = 0; % during the lag phase, assume no respiration
Resp(t>TE,:) = NaN; % don't include respiration after hatching

Xout2{1}(:,1) = t; % time is on the x-axis of the extra data set
Xout2{1}(:,2) = Resp; % repiration on the y-axis

[~,loct] = ismember(t_rem,tout);   % find where the requested time points are in the long Xout
Xout     = Xout(loct,:);           % only keep the ones we asked for


%% Events function
% Modify this part of the code if _eventson_=1. This subfunction catches
% the 'events': in this case, this one catches birth and puberty.
%
% Note that the eventsfun has the same inputs, in the same sequence, as
% <derivatives.html derivatives.m>.

function [value,isterminal,direction] = eventsfun(t,X,par,c,glo)

WVp      = par.WVp(1); % body mass at puberty
nevents  = 2;          % number of events that we try to catch

value       = zeros(nevents,1);  % initialise with zeros
value(1)    = X(glo.locL) - WVp; % thing to follow is body mass (state 1) minus threshold
value(2)    = X(3);              % follow the egg buffer (state 1) to spot hatching
isterminal  = zeros(nevents,1);  % do NOT stop the solver at an event
direction   = zeros(nevents,1);  % catch ALL zero crossing when function is increasing or decreasing
