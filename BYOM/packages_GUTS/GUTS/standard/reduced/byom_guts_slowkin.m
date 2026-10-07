%% BYOM, byom_guts_slowkin.m
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
% reduced model (TK and damage dynamics lumped): SD, IT and mixed (or GUTS
% proper). 
%
% *This script:* A simulated data set with 'slow kinetics', which means
% that the dominant rate constant _kd_ goes to zero. The data were created
% with SD. This illustrates how the 'slow kinetics' option can be used. The
% files in this directory use the analytical solution only!
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

pathdefine(0) % set path to the BYOM/engine directory (option 1 uses parallel toolbox)
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 0; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. 

% observed number of survivors, time in days, simulated data for slow kinetics
DATA{1} = [-1     0     2     4     8    16 
     0    20    20    20    20    20
     1    20    20    20    16    10
     2    20    20    14     8     0
     3    20    15     7     3     0
     4    20    11     1     0     0];
 
% scaled damage, can have no observations
DATA{2} = 0;   

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = DATA{1}(1,2:end); % scenarios (concentrations)
X0mat(2,:) = 1;                % initial survival probability
X0mat(3,:) = 0;                % initial scaled damage

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% global parameters for GUTS purposes
glo.sel  = 1; % select death mechanism: 1) SD 2) IT 3) mixed
glo.locD = 2; % location of scaled damage in the state variable list
glo.locS = 1; % location of survival probability in the state variable list

% start values for SD
% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.kd = [0.06  1 1e-3  10  0];  % dominant rate constant, d-1
par.mw = [0.18  1 1e-6 1e6  0];  % median threshold for survival (ug/L)
par.hb = [0.001 1 0      1  1];  % background hazard rate (1/d)
par.bw = [2     1 1e-6 1e6  0];  % killing rate (L/ug/d) (SD only)
par.Fs = [1     0 1    100  1];  % fraction spread of threshold distribution (IT only)

switch glo.sel % make sure that right parameters are fitted
    case 1 % for SD ...
        par.Fs(2) = 0; % never fit the threshold spread!
    case 2 % for IT ...
        par.bw(2) = 0; % never fit the killing rate!
    case 3 % mixed
        % do nothing: fit all parameters
end

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set

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

opt_optim.fit = 1; % fit the parameters (1), or don't (0)
opt_optim.it  = 0; % show iterations of the optimisation (1, default) or not (0)
 
glo.fastslow = 's'; % tell simplefun that we need slow (s) or fast (f) kinetics
% slow kinetics implies that mw is now used for mw/kd and bw becomes bw*kd!
% so we need to modify initial values for some parameters:
par.kd(2) = 0;    % don't fit the dominant rate constant anymore
par.mw(1) = 0.7;  % NOTE: this is now mw/kd!
par.bw(1) = 0.07; % NOTE: this is now bw*kd!

opt_optim.type = 1; % optimisation method 1) simplex, 2) simulated annealing (experimental), % 3) swarm (experimental)

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

% Below, several options for post analyses are provided.

return % stop here, you can run the next sections afterwards if you like

%% Likelihood region
% Make intervals on model predictions by using a sample of parameter sets
% taken from the joint likelihood-based conf. region. Sub-optimisations set
% to zero for this demo; in general, it is a good idea to take 10
% sub-optimisations for robustness.

opt_prof.subopt   = 0; % number of sub-optimisations to perform to increase robustness
opt_prof.detail   = 2; % detailed (1) or a coarse (2) calculation
opt_prof.brkprof  = 2; % when a better optimum is located, stop (1) or automatically refit (2)
opt_likreg.burst  = 2000; % number of random sets from parameter space evaluated every iteration

par_better = calc_likregion(par_out,500,opt_likreg,opt_prof,opt_optim); % second argument is target for number of samples in inner region (-1 to re-use saved sample from previous runs)
% Second entry is the number of accepted parameter sets to aim for. Use -1
% here to use a saved set.

if ~isempty(par_better) % if the profiling found a better optimum ...
    calc_and_plot(par_better,opt_plot); % calculate model lines and plot them
    par_out = par_better; % use the new parameter structure for further analyses below
end

%% Plot results with confidence intervals
% The following code can be used to make plots with confidence intervals.
% Options for confidence bounds on model curves can be set using opt_conf
% (see prelim_checks). The plot_tktd function makes multiplots for the
% survival data, which are more readable when plotting with various
% intervals.

opt_conf.type    = 2; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_tktd.max_exp = 0; % set to 1 to maximise exposure/damage plots on exposure rather than damage
% Scaling on damage is better here, though damage is now damage-time for
% slow kinetics!

plot_tktd(par_out,opt_tktd,opt_conf); 

%% Calculate LCx versus time
% Here, the LCx (by default the LC50) is calculated at several time points.
% LCx values are also printed on screen. If a sample from parameter space
% is available (e.g., from the slice sampler or the likelihood region), it
% can be used to calculate confidence bounds. 
% 
% Options for LCx (with confidence bounds) can be set using opt_lcx_lim
% (see prelim_checks). Note that opt_conf.type=-1 skips CIs.

opt_conf.type    = 2; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_lcx_lim.Feff = 0.50; % effect level (>0 en <1), x/100 in LCx

% We can directly use the fast method specific for this reduced GUTS folder
Tend = [2 4 10 20 30]; % times at which to calculate LCx, relative to control
calc_lcx_lim_guts_red(par_out,Tend,opt_lcx_lim,opt_conf); % calculates LCx values, CI requires that there is a mat file with sample

%% Calculate LPx
% We can use the function calc_lpx_lim to calculate an exposure
% multiplication factor: with which factor do we need to multiply an
% exposure scenario to obtain x% effect at the end of the scenario? Here,
% done for 10% effect at the end of the exposure scenario. Note that
% opt_conf.type = -1 skips CIs.
%
% This is the fast method, which is closely linked to this GUTS model
% (reduced model, standard directory). 

opt_conf.type     = 2; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set  = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_lcx_lim.Feff  = 0.10; % effect level (>0 en <1), x/100 in LCx/LPx

Tend  = []; % time at which to calculate LPx (empty is end of profile)
LPx_mon = calc_lpx_lim_guts_red(par_out,Tend,'profile_monit.txt',opt_lcx_lim,opt_conf);
