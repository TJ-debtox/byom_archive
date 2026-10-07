%% BYOM function call_deri.m (calculates the model output)
%
%  Syntax: [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)
%
% This package does not use the ODE solver, so this function is simplified.
% However, we keep some things in here to make this package in line with
% the general flow of BYOM. This function calls the explicit function(s) in
% <simplefun.html simplefun.m>. As input, it gets:
%
% * _t_   the time vector (in this package: concentration)
% * _par_ the parameter structure
% * _X0v_   a vector with initial states and one concentration (not used)
% * _glo_ the structure with various types of information (used to be global)
%
% The output _Xout_ provides a matrix with exposure concentration in rows,
% and states in columns.
%
% * Author: Tjalling Jager
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_doseresp.html>
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

% Unpack the vector X0v, which is X0mat for one scenario
X0 = X0v(2:end); % these are the intitial states for a scenario
% In the dose-response package, this is not used. However, _X_0 remains in
% this function to keep it in line with the rest of BYOM.

TE = +inf; % this is needed for compatibaility with rest of BYOM

%% Calculations
% This part calls the explicit model in <simplefun.html simplefun.m>) to
% calculate the output (the value of the state variables over time). There
% is generally no need to modify this part. Note that in the dose-response
% package, _t_ is used for the concentration!

c  = X0v(1);     % the concentration (or scenario number)
t  = t(:);       % force t to be a row vector (needed when useode=0)

% use an explicit function provided in simplefun!
Xout = simplefun(t,X0,par,c,glo);

%% Output mapping
% _Xout_ contains a row for each state variable. It can be mapped to the
% data. If you need to transform the model values to match the data, do it
% here. 

% Xout(:,1) = Xout(:,1).^3; % e.g., do something on first column, like cube it ...
