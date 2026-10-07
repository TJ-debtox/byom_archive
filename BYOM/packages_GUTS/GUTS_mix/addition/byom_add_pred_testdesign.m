%% BYOM, byom_add_pred_testdesign.m
%
% * Author: Tjalling Jager
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m> or <simplefun.html
% simplefun.m>, and <call_deri.html call_deri.m> may need to be modified to
% the particular problem as well. The files in the engine directory are
% needed for fitting and plotting. Results are shown on screen but also
% saved to a log file (|results.out|).
%
% *The model:* fitting binary mixtures for survival data with the
% <http://www.debtox.info/about_guts.html GUTS> special cases based on the
% reduced model (TK and damage dynamics lumped): SD or IT. All
% calculations performed in <simplefun.html simplefun.m>. This is for
% damage addition.
%
% *This script:* This script is set up to use a previously-made MAT file
% (from fitting on two single-exposure data sets, or combined with
% combine_mat_files.m) to make predictions for mixture test design.
% Additionally, it shows how to use glo.mix_fact=100 for if you need more
% than 9 treatments.
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

pathdefine(1) % set path to the BYOM/engine directory (option 1 uses parallel toolbox)
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

% Treatments are referred to by a two-digit code as in xy, which would mean
% concentration number x for chemical A and concentration number y for
% chemical B. Here, the first chemical is IMI (A) and second THI (B). The
% global variable tab_mix is a translation table for simplefun. In the
% future, I would like to implement this deeper into BYOM, probably into
% make_scen. Note that you can use any coding that you like, as long as
% each code refers to the correct concentration combination (which is laid
% down in tab_mix).
%
% The code below makes use of the option in BYOM to have multiple data sets
% per state, so the single exposures of IMI are a separate data set from
% the single exposures of THI. This means they are plotted in separate
% graphs, but has no consequences for the calculations. For mixture
% survival, you can add DATA{3,1}, or more data sets, with their own unique
% treatment coding.
% 
% The exposure definition is done with make_scen. All possible single
% exposures are defined; the mixtures are all a combination of single
% exposures. In simplefun, the scenario identifier is separated into a
% scenario for chemical A and a scenario for chemical B. So, scenario 42
% implies taking exposure to chemical A from scenario 40, and exposure to
% chemical B from scenario 2. By default, we use a factor of 10 to separate
% the scenario into chemical A and B. This can be modified to 100, by
% changing glo.mix_fact. In that case, IDs 1-99 are for chemical A, and
% 100-9900 (in steps of 100) for chemical B.

glo.mix_fact = 100; % factor for scenario IDs of the mixture analysis

% Use a test scenario with sequential exposures. 
CwA = [1 100 200 300
    0 20 30 50
    1 0 0 0
    4 0 0 0];

make_scen(2,CwA);

CwB = [1 1 2 3
    0  0 0 0
    2 20 30 50
    3 0 0 0
    4 0 0 0];

make_scen(2,CwB);
          
%% Create a table with nicer labels for the legends
% Creating a Matlab table in glo.LabelTable will replace the automatically
% generated labels (from glo.leglab1/2 and the scenario identifiers) with a
% dedicated text for each scenario.

% Start with chemical 
Scenario = [0 101 202 303]';
Label    = {'control';'A to B 20 mg/L';'A to B 30 mg/L';'A to B 50 mg/L'};

glo.LabelTable = table(Scenario,Label); % create a Matlab table for the labels

%% Call the Matlab GUI open-file element to load MAT file

[conf_type,~,par] = select_pred([1 1 0]); 

% Note that this function, with these settings, loads a pre-saved MAT file.
% The parameter structure <par> is loaded from the selected MAT file. This
% contains starting values for the parameter-space explorer, as generated
% by <combine_mat_files>. It also defines glo.sel.
% 
% Note that structure glo is not loaded from file. This is done to
% accommodate use of the MAT file generated by combine_mat_files from the
% single exposure fits.

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = glo.LabelTable.Scenario'; % scenarios to be run (treatment codes)
X0mat(2,:) = 1;      % initial survival probability
X0mat(3,:) = 0;      % initial scaled damage
X0mat(4,:) = 0;      % initial scaled damage

%% Initial values for the model parameters
% All parameters have been loaded from file using select_pred. All global

% global parameters for GUTS purposes
glo.locS = 1; % location of survival in the state-variable vector
glo.locD = [2 3]; % location of damage states in the state-variable vector

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

glo.t = linspace(0,4,100);

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival probability';
glo.ylab{2} = 'scaled damage PRO (mg/L)';
glo.ylab{3} = 'scaled damage TRI (mg/L)';
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'scen. '; % legend label before the 'scenario' number
glo.leglab2 = ''; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Plot results with confidence intervals
% The following code can be used to make plots with confidence intervals.
% Options for confidence bounds on model curves can be set using opt_conf
% (see prelim_checks). The plot_tktd function makes multiplots for the
% survival data, which are more readable when plotting with various
% intervals.
%
% Note that the MAT file selected with select_pred is used for loading the
% parameters and the sample for plot_tktd. This is present in glo.mat_nm.

opt_conf.type    = conf_type; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
% Note: using the type that relates to the MAT file you selected!

opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_conf.set_zero = 'hb'; % parameter name to set to zero (usually the background hazard)
% Last line forces hb=0 (both in the best-fit parameter set and in the
% sample).

% Calculate and plot dedicated TKTD plots. These plots are more readable when plotting CIs.
opt_tktd.preds   = 1; % set to 1 only plot predictions from X0mat without data
opt_tktd.addzero = 1; % set to 1 to always add a concentration zero to X0mat
opt_tktd.obspred = 2; % plot predicted-observed plots (1) or not (0), (2) makes 1 plot for multiple data sets

plot_tktd([],opt_tktd,opt_conf); % empty input is for par_out, which is read from mat file

