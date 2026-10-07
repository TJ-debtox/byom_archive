%% BYOM, byom_guts_parspace.m
%
% * Author     : Tjalling Jager
% * Date       : November 2021
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
% *This script:* This script uses the parameter-space explorer (parspace)
% for very simple, automated, data analysis. It loads input data sets and
% exposure profiles in <http://www.openguts.info openGUTS> format. All
% files are added in the sub-folders to perform the extended ring-test
% calculations. The files in this directory use the analytical solution
% only!
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

pathdefine(1) % set path to the BYOM/engine directory
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 0; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. First column are the exposure times, first
% row are the concentrations or scenario numbers. The number in the top
% left of the matrix indicates how to calculate the likelihood (-1 for
% for survival data). 
% 
% NOTE: time MUST be entered in DAYS for the estimation of starting values
% to provide proper search ranges!

% The code below uses the Matlab file-select GUI element to select one (or
% more) openGUTS input file(s). Next, it translates them into the BYOM
% format. The function load_data is adapted from openGUTS, and also sets
% glo.basenm to a string that includes the filename of the data file.
[filename,filepath] = uigetfile('input_data/*.txt','Select file(s) for calibration','MultiSelect','on'); % use Matlab GUI to select file(s)
if ~iscell(filename) && numel(filename) == 1 && filename(1) == 0 % if cancel is pressed ...
    return % simply stop
end
[data_all,scen_tot,unit,error_flag] = load_data_openguts(filename,filepath); % load and prepare the data set
if error_flag == 1 % see if there are errors spotted in the input file
    return % simply return (stop and do nothing)
end
DATA = cell(length(data_all),1); % empty cell array to catch the output correctly
DATA(:,1) = data_all(:,1); % place the output into the DATA structure
   
% scaled damage, can have no observations
DATA(:,2) = {0};   

% Use a dialog to let the user select whether to perform an SD or IT!
x = questdlg('Select death mechanism','User select','SD','IT','SD');
sel = 1; % by default, use SD
if strcmp(x,'IT') % unless user presses IT
    sel = 2;
end

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = scen_tot; % scenarios (concentrations or identifiers)
X0mat(2,:) = 1;        % initial survival probability
X0mat(3,:) = 0;        % initial scaled damage

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% global parameters for GUTS purposes
glo.sel  = sel; % select death mechanism: 1) SD 2) IT 3) mixed
glo.locD = 2; % location of scaled damage in the state variable list
glo.locS = 1; % location of survival probability in the state variable list

% start values and ranges are not used with parspace, but the fit flag and log setting are
% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.kd = [1    1 1e-3  100 0];   % dominant rate constant, d-1
par.mw = [1    1    0  1e6 1];   % median threshold for survival (ug/L)
par.hb = [0.01 1    0    1 1];   % background hazard rate (1/d)
par.bw = [1    1 1e-6  1e6 0];   % killing rate (L/ug/d) (SD only)
par.Fs = [2    1    1  100 1];   % fraction spread of threshold distribution (IT only)
% NOTE: we use parspace for optimisation here, and an automated generation
% of search ranges in startgrid_guts. Therefore, the only aspect of the
% parameter structure _par_ that is used is the fit column.

% After optimisation, copy-paste relevant lines from screen below. For
% these analyses, that could be done when _hb_ is estimated from the
% control treatment.

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
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival probability';
glo.ylab{2} = ['scaled damage (',unit,')'];
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'conc. '; % legend label before the 'scenario' number
glo.leglab2 = unit; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the
% plotting. Options for the plotting can be set using opt_plot (see
% prelim_checks.m). Options for the optimisation routine can be set using
% opt_optim. The files in this directory always apply the analytical
% solution for damage in simplefun.

opt_optim.type   = 4; % optimisation method 1) simplex, 4 parameter-space explorer
opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_plot.bw      = 1; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
% opt_plot.statsup = [2]; % vector with states to suppress in plotting fits
opt_optim.ps_slice = 0; % 0) standard openGUTS algorithm, 1) use slice sampler for main rounds

opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots to monitor progress of parameter-space explorer
opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
opt_optim.ps_rough = 0; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);
opt_optim.ps_dupl  = 0; % set to 1 to remove duplicates from sample in calc_parspace (slow for rough=0!)

% The switch fit_tox is used to select whether to fit the control treatment
% (0), the toxicity treatments (1; the control is shown as well), or both
% together (2). A good strategy is to fit the background hazard to the
% control data first (fit_tox=0), copy the best values into the parameter
% matrix above, and then fit the toxicity parameter to the complete data
% set (fit_tox=1). The code below automatically does that.

fit_tox = 2; % set by user
skip_sg = 0; % set to 1 to skip startgrid completely (use ranges in <par> structure)

par_out = automatic_runs_basic(fit_tox,par,-1,skip_sg,opt_optim,opt_plot);

% Dedicated TKTD plots; these plots are more readable, especially when plotting CIs
opt_tktd.repls   = 0; % plot individual replicates (1) or means (0)
opt_conf.type    = 3; % make intervals from parspace explorer
opt_conf.lim_set = 2; % use limited set of n_lim points (1) or outer hull (2) to create CIs
plot_tktd(par_out,opt_tktd,opt_conf);

% Below, several options for post analyses are provided. Note that the
% parameter-space explorer optimises AND provides CIs on model parameters
% AND saves a sample for CIs on predictions. So no need to profile the
% likelihood again or generate a sample for making CIs on predictions.

% return

%% Calculate LCx versus time
% Here, the LCx (by default the LC50) is calculated at several time points.
% LCx values are also printed on screen. If a sample from parameter space
% is available (e.g., from the slice sampler or the likelihood region), it
% can be used to calculate confidence bounds. 
% 
% Options for LCx (with confidence bounds) can be set using opt_lcx_lim
% (see prelim_checks). Note that opt_conf.type=-1 skips CIs.

% UNCOMMENT FOLLOWING LINE(S) TO CALCULATE
opt_conf.type    = 3; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_lcx_lim.Feff = 0.50; % effect level (>0 en <1), x/100 in LCx/LPx

% Tend = [4 7 14 21 28 30]; % times at which to calculate LCx, relative to control
% calc_lcx_lim_guts_red(par_out,Tend,opt_lcx_lim,opt_conf); % calculates LCx values, CI requires that there is a mat file with sample

% This is the slower general method
opt_ecx.Feff      = [0.10 0.50]; % effect levels (>0 en <1), x/100 in ECx
opt_ecx.plot_log  = 1; % set to 1 to plot ECx-time on log-log scale
opt_ecx.notitle   = 1; % set to 1 to suppress titles above ECx plots
Tend = [4 5 6 7 8 10 14 18 20];
calc_ecx(par_out,Tend,opt_ecx,opt_conf); % general method for ECx values

return

%% Calculate LPx
% We can use the function calc_lpx_lim to calculate an exposure
% multiplication factor: with which factor do we need to multiply an
% exposure scenario to obtain x% effect at the end of the scenario? Here,
% done for 10% effect at the end of the exposure scenario. Note that
% opt_conf.type = -1 skips CIs.
%
% This is the fast method, which is closely linked to this GUTS model
% (reduced model, standard directory). 

opt_conf.type     = 3; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set  = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_lcx_lim.Feff  = 0.10; % effect level (>0 en <1), x/100 in LCx/LPx
opt_lcx_lim.notitle = 1; % set to 1 to suppress titles above ECx plots

% use Matlab GUI to load exposure profile from text file for predictions
[fname_prof,filepath] = uigetfile('input_profile/*.txt','Select a text file with exposure profile for predictions','MultiSelect','off'); % use Matlab GUI to select file(s)
if ~iscell(fname_prof) && numel(fname_prof) == 1 && fname_prof(1) == 0 % if cancel is pressed ...
    return % simply stop
end

Tend = []; % time at which to calculate LPx (empty is end of profile)
LPx  = calc_lpx_lim_guts_red(par_out,Tend,[filepath,fname_prof],opt_lcx_lim,opt_conf);

% % This is the slower general method
% opt_ecx.Feff      = [0.10]; % effect levels (>0 en <1), x/100 in ECx
% opt_ecx.notitle   = 1; % set to 1 to suppress titles above ECx plots
% calc_epx(par_out,[filepath,fname_prof],[],opt_ecx,opt_conf,opt_tktd); % general method for EPx values
