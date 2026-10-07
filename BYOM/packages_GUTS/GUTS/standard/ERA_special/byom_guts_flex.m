%% BYOM, byom_guts_flex.m
%
% * Author: Tjalling Jager
% * Date: September 2023
% * Web support: <http://www.debtox.info/byom.html>
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
% *This script:* This loads an openGUTS file from the directory input_data
% and uses it. It can also load multiple files for a joint fit.
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

% NOTE: time MUST be entered in DAYS for the estimation of starting to
% provide proper search ranges!

% The code below uses the Matlab file-select GUI element to select one (or
% more) openGUTS input file(s). Next, it translates them into the BYOM
% format. The function load_data is adapted from openGUTS, and also sets
% glo.basenm to a string that includes the filename of the data file.
[filename,filepath] = uigetfile('input_data/*.txt','Select openGUTS input file(s) for calibration','MultiSelect','on'); % use Matlab GUI to select file(s)
if ~iscell(filename) && numel(filename) == 1 && filename(1) == 0 % if cancel is pressed ...
    return % simply stop
end
[data_all,scen_tot,unit,error_flag] = load_data_openguts(filename,filepath); % load and prepare the data set
if error_flag == 1 % see if there are errors spotted in the input file
    return % simply stop, errors will have been shown on screen
end
% Note that load_data_openguts also defines the exposure scenario for each
% data set, using linear extrapolation (as in openGUTS).

DATA      = cell(length(data_all),2); % empty cell array to catch the output correctly
DATA(:,1) = data_all(:,1);            % place the output into the DATA structure
   
% scaled damage, can have no observations
DATA(:,2) = {0};   

% Note: optionally, add some info to the MAT filename. The MAT filename
% will already include MoA and feedbacks, but if you want to try other
% things as well (changing opt, calibrating on the validation data, etc)
% it can be helpful to change the name to use this script but not 
% overwrite previous MAT files.
% 
% glo.basenm  = [mfilename,'_CAL1']; % remember the filename for THIS file for the plots

save([glo.basenm,'_DATA'],'DATA','W') % save MAT file with data set 
% Saving the data set is handy to allow for a simple reconstruction of the
% calibrations, without needing to define the data again, in the same way.
% Note that if you use DATAx and Wx, you need to save them as well!

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = scen_tot; % scenarios (concentrations)
X0mat(2,:) = 1;        % initial survival probability
X0mat(3,:) = 0;        % initial scaled damage

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% global parameters for GUTS purposes
glo.locS = 1; % location of survival probability in the state variable list
glo.locD = 2; % location of scaled damage in the state variable list

% NOTE: both SD and IT can be calculated automatically.
% 
% The start values and ranges below are not needed (filled later by
% startgrid_guts); only the fit/fix mark is relevant here. But make sure
% that the value in the first column is within the bounds (and not zero for
% log-scale parameters).

% start values and ranges are not used with parspace, but the fit flag and log setting are
% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.kd = [1    1 1e-3  100 0];   % dominant rate constant, d-1
par.mw = [1    1    0  1e6 1];   % median threshold for survival (ug/L)
par.hb = [0.01 1    0    1 1];   % background hazard rate (1/d)
par.bw = [1    1 1e-6  1e6 0];   % killing rate (L/ug/d) (SD only)
par.Fs = [2    1    1  100 1];   % fraction spread of threshold distribution (IT only)
% NOTE: we use parspace for optimisation here, and an automated generation
% of search ranges. Therefore, the only aspect of the parameter structure
% _par_ that is used is the fit column.

% After optimisation, copy-paste relevant lines from screen below. For
% these analyses, that would only be useful when _hb_ is estimated from the
% control treatment only.

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival probability';
glo.ylab{2} = ['scaled damage (',unit,')'];
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'conc. '; % legend label before the 'scenario' number
glo.leglab2 = [unit]; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the
% plotting. Options for the plotting can be set using opt_plot (see
% prelim_checks.m). Options for the optimisation routine can be set using
% opt_optim. The files in this directory always apply the analytical
% solution for damage in simplefun.

opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_optim.it     = 0; % show iterations of the optimisation (1, default) or not (0)
opt_plot.bw      = 1; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.statsup = []; % vector with states to suppress in plotting fits
basenm_rem       = glo.basenm; % remember basename as we will modify it!

% Select whether to fit the control treatment (0), the toxicity treatments
% (1; the control is shown as well), or both together (2). A good strategy
% is to fit the background hazard to the control data first (fit_tox=0),
% copy the best values into the parameter matrix above, and then fit the
% toxicity parameter to the complete data set (fit_tox=1). The code below
% automatically keeps the parameters fixed that need to be fixed.
% Select what to fit with fit_tox (this is a 2-element vector). First 
% element of fit_tox is which part of the data set to use:
%   fit_tox(1) = -1  control survival (c=0) only
%   fit_tox(1) = 0   not used for GUTS
%   fit_tox(1) = 1   fit tox parameters, but when fitting, keep hb fixed; run through all
%               elements in SEL sequentially and provide a table at the end
%   fit_tox(1) = 2   fit tox parameters, but when fitting, also fit hb; run through all
%               elements in SEL sequentially and provide a table at the end
% 
% Second element of fit_tox is whether to fit or only to plot:
%   fit_tox(2) = 0   don't fit; for standard optimisations, plot results for
%               parameter values in [par], for parspace optimisations, use saved mat file.
%               (there is now no difference between fit_tox(1) set to 1 or 2!)
%   fit_tox(2) = 1   fit parameters
% 
% A good strategy is to fit the background hazard to the control data first
% (fit_tox=[-1 1]), copy the best value into the parameter matrix above,
% and rerun this script to fit the toxicity parameter to the complete data
% set (fit_tox=[1 1]). The code below automatically keeps the parameters
% fixed that need to be fixed.
% 
% SEL is a vector with the death mechanisms to to run through, sequentially:
% 1) stochastic death
% 2) individual tolerance
% 3) both (GUTS proper)
% 
% For fit_tox(1)=-1, the setting of SEL has no impact; it must be defined
% to prevent errors. For fit_tox(1)~=-1, setting is relevant. Note that if
% you run multiple death mechanisms, the last one will remain in the
% memory. So, if you continue with the code below (e.g., plotting CIs on
% model curves), only the last setting in SEL will be used!

SEL     = [1 2]; % These are the death mechanisms that will be run automatically for fit_tox 1 and 2

% ===== FITTING CONTROLS ==================================================
% Simplex fitting works fine for control data
fit_tox = [-1 1]; % set by user (see table above)

opt_optim.type = 1; % optimisation method: 1) default simplex, 4) parspace explorer
par_out = automatic_runs_guts(fit_tox,par,[],SEL,[],opt_optim,opt_plot); % script to run the calculations and plot, automatically
par     = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par

% ===== FITTING TOX DATA ==================================================
% Use the parameter-space explorer for fitting the treatments. Note that,
% with fit_tox(1)=1, the tox parameters in par are replaced (when fitted) with
% estimates based on the data set using startgrid_guts (called in
% automatic_runs_guts).
fit_tox = [1 1]; % set by user (see table above)

opt_optim.type     = 4; % optimisation method 1) simplex, 4 parameter-space explorer
% Note: to use saved set for parameter-space explorer, use fit_tox(2)=0!
opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots of parameter space to monitor progress
opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS (2 for extra rough)
opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);

[par_out,best_sel] = automatic_runs_guts(fit_tox,par,[],SEL,[],opt_optim,opt_plot); % script to run the calculations and plot, automatically
% Note: automatic_runs will return the BEST parameter set in par_out.

% Code below plots the best death mechanism. It plots the fit without CIs.
% The automatic_runs returns the parameter set and indices for the best
% settings. 
glo.sel    = SEL(best_sel(1));   % change global for sel
glo.basenm = basenm_rem; % return the basenm to this filename 
% The plots saved will get their filename based on the name of THIS
% script (so without the specification of MoA and feedbacks). 

% Dedicated TKTD plots. These plots are more readable, especially when plotting CIs.
opt_tktd.repls   = 0; % plot individual replicates (1) or means (0)
opt_tktd.obspred = 0; % plot predicted-observed plots (1) or not (0)
opt_tktd.max_exp = 0; % set to 1 to maximise exposure/damage plots on exposure rather than damage

plot_tktd(par_out,opt_tktd,[]);
% Leave the options opt_conf empty to suppress all CIs for these plots.
% Plotting CIs requires a sample as saved by the various methods available
% in BYOM (see examples directory).

%% Plot results with confidence intervals
% The following code can be used to make plots with confidence intervals.
% Options for confidence bounds on model curves can be set using opt_conf
% (see prelim_checks). The plot_tktd function makes multiplots for the
% survival data, which are more readable when plotting with various
% intervals.
%
% Note that parspace has already calculated CIs on model parameters, and
% saved a .mat file with a sample. There is no need to run additional
% functions to generate CIs.

opt_conf.type     = 3; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set  = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_tktd.notitle  = 0; % set to 1 to suppress titles above plots (becomes messy for GUTS)
opt_tktd.lim_data = 1; % set to 1 to limit axes to data
opt_tktd.obspred  = 1; % plot predicted-observed plots (1) or not (0)

% Below some tricks to allow plotting results for SD and IT. The
% automatic_runs has modified glo.basenm so we will reconstruct it to load
% the correct MAT file.
for i = 1:length(SEL)
    glo.sel = SEL(i); % take the next death mechanism
    glo.basenm = [basenm_rem,'_sel',sprintf('%d',glo.sel)]; % create a new basenm for each death mechanism
    % Calculate and plot dedicated TKTD plots. These plots are more readable
    % when plotting CIs.
    plot_tktd([],opt_tktd,opt_conf);
    % Note that first input is left empty: this is for the parameter structure,
    % which is then obtained from the saved sample.
    % leave the options opt_conf empty to suppress all CIs for these plots
end

return

%% Calculate LCx versus time
% Here, the LCx (by default the LC50) is calculated at several time points.
% LCx values are also printed on screen. If a sample from parameter space
% is available (e.g., from the slice sampler or the likelihood region), it
% can be used to calculate confidence bounds. 
% 
% Options for LCx (with confidence bounds) can be set using opt_lcx_lim
% (see prelim_checks). Note that opt_conf.type=-1 skips CIs.
% 
% Note that in this script, running the code in this section will only plot
% the results from the LAST run you performed above. Same trick can be used
% as above to calculate LCx,t for both SD and IT.

opt_conf.type    = 3; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_lcx_lim.Feff = 0.50; % effect level (>0 en <1), x/100 in LCx/LPx

% This is for the slower general method
opt_ecx.Feff      = [0.10 0.50]; % effect levels (>0 en <1), x/100 in ECx
opt_ecx.notitle   = 1; % set to 1 to suppress titles above ECx plots

Tend = [1:8]; % times at which to calculate LCx, relative to control

% Below some tricks to allow calculating LCx for SD and IT. The
% automatic_runs has modified glo.basenm so we will reconstruct it to load
% the correct MAT file.
for i = 1:length(SEL)
    glo.sel = SEL(i); % take the next death mechanism
    glo.basenm = [basenm_rem,'_sel',sprintf('%d',glo.sel)]; % create a new basenm for each death mechanism
    calc_lcx_lim_guts_red([],Tend,opt_lcx_lim,opt_conf); % fast method
    % calc_ecx([],Tend,opt_ecx,opt_conf); % slow but general method for ECx values
end

return

%% Calculate LPx
% We can use the function calc_lpx_lim to calculate an exposure
% multiplication factor: with which factor do we need to multiply an
% exposure scenario to obtain x% effect at the end of the scenario? Here,
% done for 10% effect at the end of the exposure scenario. Note that
% opt_conf.type=-1 skips CIs.
%
% This is the fast method, which is closely linked to this GUTS model
% (reduced model, standard directory). 
% 
% Note that in this script, running the code in this section will only plot
% the results from the LAST run you performed above. Same trick can be used
% as above to calculate LCx,t for both SD and IT.

opt_conf.type     = 3; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set  = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_lcx_lim.Feff  = 0.10; % effect level (>0 en <1), x/100 in LCx/LPx
opt_lcx_lim.notitle = 0; % set to 1 to suppress titles above ECx plots

% This is for the slower general method
opt_ecx.Feff      = [0.10]; % effect levels (>0 en <1), x/100 in ECx
opt_ecx.notitle   = 0; % set to 1 to suppress titles above ECx plots

% use Matlab GUI to load exposure profile from text file for predictions
[fname_prof,filepath] = uigetfile('input_profile/*.txt','Select a text file with exposure profile for predictions','MultiSelect','off'); % use Matlab GUI to select file(s)
if ~iscell(fname_prof) && numel(fname_prof) == 1 && fname_prof(1) == 0 % if cancel is pressed ...
    return % simply stop
end

Tend    = []; % time at which to calculate LPx (empty is end of profile)

% Below some tricks to allow calculating LPx for SD and IT. The
% automatic_runs has modified glo.basenm so we will reconstruct it to load
% the correct MAT file.
for i = 1:length(SEL)
    glo.sel = SEL(i); % take the next death mechanism
    glo.basenm = [basenm_rem,'_sel',sprintf('%d',glo.sel)]; % create a new basenm for each death mechanism
    LPx = calc_lpx_lim_guts_red([],Tend,[filepath,fname_prof],opt_lcx_lim,opt_conf);
    % calc_epx([],[filepath,fname_prof],[],opt_ecx,opt_conf,opt_tktd); % general method for EPx values
end
