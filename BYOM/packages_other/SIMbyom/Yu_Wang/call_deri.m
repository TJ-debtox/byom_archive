%% SIMbyom function call_deri.m (calculates the model output)
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
% columns. This function calls <derivatives.html derivatives.m>. 
%
% * Author: Tjalling Jager
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_simbyom.html>
% 
%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM.

%% Start

function [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)

% These outputs need to be defined, even if they are not used
Xout2    = []; % additional uni-variate output, not used in this example
zvd      = []; % additional zero-variate output, not used in this example

%% Initial settings
% For the simulations, the regular options make little sense, so they have
% been removed. Keep _useode=1_ as this SIMbyom is not so useful for
% explicit functions.

useode   = 1; % calculate model using ODE solver (1) or analytical solution (0)
eventson = 0; % events function on (1) or off (0)
stiff    = 0; % set to 1 to use a stiff solver instead of the standard one

% Unpack the vector X0v, which is X0mat for one scenario
X0 = X0v(2:end); % these are the intitial states for a scenario

%% Calculations
% This part calls the ODE solver to calculate the output (the value of the
% state variables over time). There is generally no need to modify this
% part. The solver ode45 generally works well. For stiff problems, the
% solver might become very slow; you can try ode15s instead.

c  = X0v(1);     % the concentration (or scenario number)
t  = t(:);       % force t to be a row vector (needed when useode=0)

% specify options for the ODE solver
options = odeset; % start with default options
if eventson == 1
    options = odeset(options, 'Events',@eventsfun); % add an events function
end
options = odeset(options, 'RelTol',1e-4,'AbsTol',1e-7); % specify tightened tolerances
% options = odeset(options,'InitialStep',max(t)/1000,'MaxStep',max(t)/100); % specify smaller stepsize

TE = 0; % dummy for time of events

if useode == 1 % use the ODE solver to calculate the solution
    % call the ODE solver (use ode15s for stiff problems, escpecially for pulses)
    if isempty(options.Events) % if no events function is specified ...
        switch stiff
            case 0
                [~,Xout] = ode45(@derivatives,t,X0,options,par,c,glo);
            case 1
                [~,Xout] = ode113(@derivatives,t,X0,options,par,c,glo);
            case 2
                [~,Xout] = ode15s(@derivatives,t,X0,options,par,c,glo);
        end
    else % with an events functions ... additional output arguments for events:
        % TE catches the time of an event, YE the states at the event, and IE the number of the event
        switch stiff
            case 0
                [~,Xout,TE,YE,IE] = ode45(@derivatives,t,X0,options,par,c,glo);
            case 1
                [~,Xout,TE,YE,IE] = ode113(@derivatives,t,X0,options,par,c,glo);
            case 2
                [~,Xout,TE,YE,IE] = ode15s(@derivatives,t,X0,options,par,c,glo);
        end
    end
else % alternatively, use an explicit function provided in simplefun!
    Xout = simplefun(t,X0,par,c);
end

if isempty(TE) || all(TE == 0) % if there is no event caught
    TE = +inf; % return infinity
end

%% Output mapping
% _Xout_ contains a row for each state variable. If you want to transform
% the model outputs, do it here.

% % To obtain the output of the derivatives at each time point. For loop can
% % be removed if derivatives accepts vectors for the states. The values in
% % dXout might be used to replace values in Xout, if the data to be fitted
% % are the changes (rates) instead of the state variable itself).
% dXout = zeros(size(Xout)); % initialise with zeros
% for i = 1:length(t) % run through all time points
%     dXout(i,:) = derivatives(t(i),Xout(i,:),par,c); 
%     % derivatives for each stage at each time
% end

%% Events function
% Modify this part of the code if _eventson_=1. This subfunction catches
% the 'events': in this case. This function should be adapted to the
% problem you are modelling (this one matches the byom_bioconc_... files,
% so it is of little use here). You can catch more events by making a
% vector out of _values_.
%
% Note that the eventsfun has the same inputs, in the same sequence, as
% <derivatives.html derivatives.m>.

function [value,isterminal,direction] = eventsfun(t,X,par,c,glo)

% EVENTS FUNCTION IS NOT USED FOR THIS EXAMPLE

Ct      = par.Ct(1); % threshold external concentration where degradation stops
nevents = 1;         % number of events that we try to catch

value       = zeros(nevents,1); % initialise with zeros
value(1)    = X(1) - Ct;        % thing to follow is external concentration (state 1) minus threshold
isterminal  = zeros(nevents,1); % do NOT stop the solver at an event
direction   = zeros(nevents,1); % catch ALL zero crossing when function is increasing or decreasing
