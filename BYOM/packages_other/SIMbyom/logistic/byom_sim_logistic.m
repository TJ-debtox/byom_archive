%% SIMbyom, byom_sim_logistic.m, logistic growth 
%
% * Author: Tjalling Jager
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_simbyom.html>
%
% This is the classical ODE for logistic growth with predation. More
% background in <http://en.wikipedia.org/wiki/Logistic_function>.
%
% Model system:
%
% $$ \frac{dN}{dt} = r N \left(1-\frac{N}{K} \right) -p $$
% 
%  Copyright (c) 2012-2021, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Initial things
% Make sure that this script is in a directory somewhere *below* the BYOM
% folder.

clear, clear global % clear the workspace and globals
global glo X0mat    % allow for global parameters in structure glo
diary off           % turn of the diary function (if it is accidentaly on)
% set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine % set path to the BYOM/engine directory
glo.basenm  = mfilename; % remember the filename for THIS file for the plots

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios. 

% scenarios specify different starting conditions
X0mat(1,:) = [1 2 3 4 5];
X0mat(2,:) = [10 12 50 100 150];

%% Values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 
  
% syntax: par.name = value
par.r   = 0.1;  % growth rate
par.K   = 100;  % carrying capacity
par.p   = 1;    % predation rate

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

glo.t = linspace(0,100,1000); % time vector for the model curves

% specify the y-axis labels for each state variable
glo.ylab{1} = 'N';
% specify the x-axis label (same for all states)
glo.xlab    = 'time';

glo.leglab1 = 'Scen. '; % legend label before the 'scenario' number
glo.leglab2 = ''; % legend label after the 'scenario' number
 
% If needed, a stiff solver or events function can be selected in call_deri

%% Optional settings for the simulations

opt_sim.plot_int = 30; % interval for plotting (how many time points to do in one go), default: 10
opt_sim.axrange  = [0 160]; % axis ranges for states, default: autoscale

opt_sim.plottype = 3; % default: 4
% 1) 3d plot
% 2) 2d plot
% 3) states in subplots
% 4) scenarios in subplots
% 5) dx/dt versus x
% 6) 2d plot for each scenario separate
% 7) 3d plot for each scenario separate

% equilibrium points (optional; comment out if not known)
opt_sim.Xeq(1,:) = [par.K*(0.5+sqrt(0.25-par.p/(par.r*par.K)))];
opt_sim.Xeq(2,:) = [par.K*(0.5-sqrt(0.25-par.p/(par.r*par.K)))];

%% Calculations and plotting
% Here, the script is called that will do the calculation and the plotting.
% The yellow points in the graph indicate the equilibrium points, and the
% squares the starting points. In the top of the file sim_and_plot.m, there
% are several more options to modify the plotting.

sim_and_plot(par,opt_sim); % call the script which calculates and plots
