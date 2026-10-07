%% BYOM, byom_guts_calanus_marecotox.m
% This is the file to reproduce the case study in Section 3.5.3 of the
% book (Fig. 3.6).
%
% * Author: Tjalling Jager
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems. The files in
% this directory use an analytical solution only, and therefore
% <derivatives.html derivatives.m> will be missing.
%
% *The model:* GUTS hazard model following Jager et al (2011, DOI
% 10.1021/es103092a). This set of files is taken from the BYOM GUTS package
% that can be downloaded from <http://www.debtox.info/byom.html>. For this
% example, the files from the `sd_special' directory are used (stochastic
% death only). This allows for the analytical solution of the simple hazard
% model to be used (in the file simplefun.m).
%
% *This script:* This data set from Øverjordet et al (2014),
% <http://dx.doi.org/10.1016/j.aquatox.2014.06.019>. The raw data were
% obtained from the original author and consitutes survival data for
% _Calanus finmarchicus_ exposed to mercury over 4 days. Note that the
% resulting parameter estimates are different from the ones in the chapter:
% it turns out that that was the result of a local minimum (the profile for
% hb nicely shows that).
% 
%  Copyright (c) 2012-2021, Tjalling Jager.
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
% concentrations) in columns. First column are the exposure times, first
% row are the concentrations or scenario numbers. The number in the top
% left of the matrix indicates how to calculate the likelihood:
%
% * -1 for multinomial likelihood (for survival data)
% * 0  for log-transform the data, then normal likelihood
% * 0.5 for square-root transform the data, then normal likelihood
% * 1  for no transformation of the data, then normal likelihood

% observed number of survivors, time in days, conc in ug/L
DATA{1} = [ -1     0     7    12    19    32    54    90   150
             0    60    30    30    30    30    30    30    30
             1    60    30    30    30    30    29    26    23
             2    60    29    28    29    29    23     9     2
             3    60    28    27    25    23    14     4     1
             4    60    27    26    24    17     9     1     0];

% damage (no data) (previously: scaled internal concentration)
DATA{2} = 0;   

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = DATA{1}(1,2:end); % scenarios (concentrations)
X0mat(2,:) = 1; % initial survival probability
X0mat(3,:) = 0; % initial damage (previously: scaled internal concentration)

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% global parameters for GUTS purposes (in this directory: SD only)
glo.locD = 2; % location of scaled damage in the state variable list
glo.locS = 1; % location of survival probability in the state variable list

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.kd    = [0.4   1 0.001 10 0]; % elimination rate constant (d-1)
par.mw    = [10    1 0    1e6 1]; % NEC for survival (ug/L)
par.hb    = [0.01  1 0      1 1]; % background hazard rate (1/d)
par.bw    = [0.023 1 1e-6  10 0]; % killing rate (L/ug/d) (SD and mixed)

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival probability';
glo.ylab{2} = 'scaled damage (mg/L)'; % (previously: scaled internal concentration)
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = ''; % legend label before the 'scenario' number
glo.leglab2 = ' mg/L'; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. 

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

opt_prof.detail   = 2; % detailed (1) or a coarse (2) calculation
opt_prof.subopt   = 10; % number of sub-optimisations to perform to increase robustness
opt_prof.brkprof  = 2; % when a better optimum is located, stop (1) or automatically refit (2)

calc_proflik(par_out,'all',opt_prof,opt_optim);  % calculate a profile
% Enter single parameter names, a cell array of names, or 'all' to profile
% all fitted parameters.

%% Other files: simplefun
% To archive analyses, publishing them with Matlab is convenient. To keep
% track of what was done, the file derivatives.m can be included in the
% published result.
% 
% <include>simplefun.m</include>

%% Other files: call_deri
% To archive analyses, publishing them with Matlab is convenient. To keep
% track of what was done, the file call_deri.m can be included in the
% published result.
%
% <include>call_deri.m</include>