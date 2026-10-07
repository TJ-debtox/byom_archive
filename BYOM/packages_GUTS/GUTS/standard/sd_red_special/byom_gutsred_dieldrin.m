%% BYOM, byom_gutsred_dieldrin.m
%
% * Author: Tjalling Jager
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_guts.html>
%
% BYOM is a General framework for simulating model systems. The files in
% this directory use an analytical solution only, and therefore
% <derivatives.html derivatives.m> will be missing.
%
% *The model:* fitting survival data with the
% <http://www.debtox.info/about_guts.html GUTS> special cases based on the
% reduced model (TK and damage dynamics lumped) and stochastic death (SD).
% 
% *This script:* Long acute toxicity test for guppy (_Poecilia reticulata_)
% exposed to the insecticide dieldrin. Data from: Bedaux and Kooijman
% (1994), <http://dx.doi.org/10.1007/BF00469427>. In this script, only few
% options for post-calculations are provided; these can be copied from the
% byom_guts_extra script if needed.
% 
%  Copyright (c) 2012-2023, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Initial things
% Make sure that this script is in a directory somewhere *below* the BYOM
% folder.

clear, clear global % clear the workspace and globals
global DATA W X0mat % make the data set and initial states global variables
global glo          % allow for global parameters in structure glo
diary off           % turn of the diary function (if it is accidentaly on)
% set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine(0) % set path to the BYOM/engine directory
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 0; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. 

% observed number of survivors, time in days, conc. in ug/L
DATA{1} = [ -1  0  3.2  5.6  10  18  32  56  100
		0     20   20   20  20  20  20  20   20
		1     20   20   20  20  18  18  17    5
		2     20   20   19  17  15   9   6    0
		3     20   20   19  15   9   2   1    0
		4     20   20   19  14   4   1   0    0
		5     20   20   18  12   4   0   0    0
		6     20   19   18   9   3   0   0    0
		7     20   18   18   8   2   0   0    0 ];

% scaled damage, can have no observations
DATA{2} = 0;   

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = DATA{1}(1,2:end); % scenarios (concentrations or identifiers)
X0mat(2,:) = 1;                % initial survival probability
X0mat(3,:) = 0;                % initial scaled damage

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% global parameters for GUTS purposes (in this directory: SD only)
glo.locD = 2; % location of scaled damage in the state variable list
glo.locS = 1; % location of survival probability in the state variable list

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.kd = [0.5  1 1e-3   20 0];   % dominant rate constant, d-1
par.mw = [5    1 0     1e6 1];   % median threshold for survival (ug/L)
par.hb = [0.01 1 1e-4    1 1];   % background hazard rate (1/d)
par.bw = [0.05 1 1e-6  1e6 0];   % killing rate (L/ug/d)

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival probability';
glo.ylab{2} = ['scaled damage (',char(181),'g/L)'];
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'conc. '; % legend label before the 'scenario' number
glo.leglab2 = [char(181),'g/L']; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the
% plotting. The files in this directory always apply the analytical
% solution for damage in simplefun.

% par = start_vals_guts(par); % experimental start-value finder; use at your own risk
% % Note: start_vals will now overwrite the parameter structure par!

opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_optim.it     = 1; % show iterations of the optimisation (1, default) or not (0)
opt_plot.bw      = 1; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
% opt_plot.statsup = [2]; % vector with states to suppress in plotting fits

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

opt_tktd.obspred = 0; % plot predicted-observed plots (1) or not (0)
plot_tktd(par_out,opt_tktd,[]); % leaving opt_conf empty suppresses all CIs for these plots

%% Profiling the likelihood
% By profiling you make robust confidence intervals for one or more of your
% parameters. 

opt_prof.detail   = 1; % detailed (1) or a coarse (2) calculation
opt_prof.subopt   = 10; % number of sub-optimisations to perform to increase robustness
opt_prof.brkprof  = 2; % when a better optimum is located, stop (1) or automatically refit (2)

calc_proflik(par_out,'all',opt_prof,opt_optim);  % calculate a profile
% Enter single parameter names, a cell array of names, or 'all' to profile
% all fitted parameters.

%% Other post-analysis can be used here as well
% Copy code from the other scripts in the 'standard' GUTS directory (e.g.,
% byom_guts_extra.m) to add various options for post-analysis here as well.
