%% BYOM function call_deri.m (calculates the model output)
%
%  Syntax: [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)
%
% This function calls the explicit function(s) in <simplefun.html
% simplefun.m> to calculate damage and survival (for SD only). As input, it
% gets:
%
% * _t_   the time vector
% * _par_ the parameter structure
% * _X0v_   a vector with initial states and one concentration (scenario number)
% * _glo_ is the structure with information (normally global)
%
% The output _Xout_ provides a matrix with time in rows, and states in
% columns. This function calls <simplefun.html simplefun.m>. The optional
% output _TE_ is not used in this package. 
%
% Note: the files in this directory are for GUTS-RED-SD only!
%
% * Author: Tjalling Jager
% * Date: June 2023
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_guts.html>

%  Copyright (c) 2012-2023, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Start

function [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)

% These outputs need to be defined, even if they are not used
Xout2    = []; % additional uni-variate output, not used in this example
zvd      = []; % additional zero-variate output, not used in this example

%% Initial settings
% This part organises a few things. Initial values can be determined by a
% parameter (overwrite parts of _X0_), and zero-variate data can be
% calculated. See the example BYOM files for more information. Note that
% this version of call_deri only works with simplefun.

% Unpack the vector X0v, which is X0mat for one scenario
X0 = X0v(2:end); % these are the intitial states for a scenario

%% Calculations
% This part calls the explicit model in <simplefun.html simplefun.m>) to
% calculate the output (the value of the state variables over time). There
% is generally no need to modify this part. This version does NOT use an
% ODE solver.

c  = X0v(1);     % the concentration (or scenario number)
t  = t(:);       % force t to be a row vector (needed when useode=0)

if isfield(glo,'int_scen') && ismember(c,glo.int_scen) % is c in the scenario range global?
    % then we have a time-varying concentration
    error('The functions in this directory can be used with constant exposure only!')
end

TE = +inf; % dummy for time of events

% to calculate damage and survival, use an explicit function provided in simplefun
[Xout,haz] = simplefun(t,X0,par,c,glo);

% NEW: this is to accommodate data sets with time-to-death, which require a
% hazard rate at the time of observation of death (next to the survival
% over time for animals that survive or that are removed or go missing).
% Here, the hazard rate is calculated in simplefun as well.
if isfield(glo,'hazout') && glo.hazout == 1
    Xout2 = haz(:); % make sure it is a column vector
end

%% Output mapping
% _Xout_ contains a row for each state variable. It can be mapped to the
% data. If you need to transform the model values to match the data, do it
% here. 
%
% Stochastic death is calculated in simplefun.m directly. 

% Xout(:,1) = Xout(:,1).^3; % e.g., do something on first column, like cube it ...

