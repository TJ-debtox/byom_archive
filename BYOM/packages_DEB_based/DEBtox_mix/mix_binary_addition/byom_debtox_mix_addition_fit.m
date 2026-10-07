%% BYOM, binary mixture addition: byom_debtox_mix_addition_fit.m
%
% *Table of contents*

%% About
% * Author     : Tjalling Jager
% * Date       : September 2023
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* simple DEBtox model for toxicants, based on DEBkiss and
% formulated in compound parameters. The model includes flexible modules
% for toxicokinetics/damage dynamics and toxic effects. The DEBkiss e-book
% (see <http://www.debtox.info/book_debkiss.html>) provides a partial
% description of the model; the publication of Jager in Ecological
% Modelling contains the full details:
% <https://doi.org/10.1016/j.ecolmodel.2019.108904>. This is the model
% version for mixture analysis, this folder specifically for damage
% addition.
%
% *This script:* fitting two single exposures together under the assumption
% of damage addition. This script can also be used to fit singles together
% with the mixture. The saved MAT file from combine_mat_files is used to
% provide the ranges for the tox parameters and initial values for the
% basic parameters. Here, there is only one control as all experiments
% where run simultaneously. If there are more controls (and if they
% differ), more care is needed ... NOTE: we cannot use separate parameters
% for each data set in the same way as the 'regular' DEBtox2019 package.
% The glo.names_sep is based on a specific definition of scenario
% identifiers, but we now already use that for the two compounds!
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

% NOTE: time MUST be entered in DAYS for the estimation of starting values
% to provide proper search ranges! Controls must use identifiers 0 for the
% true control and 0.1 for the solvent control. If you want to use another
% identifier for the solvent, change parameter id_solvent in
% automatic_runs. Note that the entire data set is defined in the script
% data_all. Below, the correct parts of the data set are assigned to the
% correct state and data set.

TR           = 0.5; % transformations for continuous data
glo.mix_fact = 10;  % factor for scenario IDs of the mixture analysis (used in derivatives to derive proper exposure scenarios for each compound)
data_all;           % read all data sets from file (this implements TR as well)

% Code below selects the partial or complete data set. PYR is chemical A,
% data set 1. FLU is chemical B, data set 2. Make sure mat_combine is run
% with these settings as well! Mixture uses data set 3.

% scaled damage
DATA{1,1} = [0]; % there are never data for this state
DATA{2,1} = [0]; % there are never data for this state
% DATA{3,1} = [0]; % there are never data for this state

% body length (mm) on each observation time (d), concentrations in uM
DATA{1,2} = L(:,[1 i_LC+1 i_LA+1]);
W{1,2}    = Lw(:,[i_LC i_LA]);
DATA{2,2} = L(:,[1 i_LB+1]);
W{2,2}    = Lw(:,[i_LB]);
% DATA{3,2} = L(:,[1 i_LM+1]);
% W{3,2}    = Lw(:,[i_LM]);

% cumulative reproduction (nr. offspring per mother) (mm) on each
% observation time (d), concentrations in uM
DATA{1,3} = R(:,[1 i_RC+1 i_RA+1]);
W{1,3}    = Rw(:,[i_RC i_RA]);
DATA{2,3} = R(:,[1 i_RB+1]);
W{2,3}    = Rw(:,[i_RB]);
% DATA{3,3} = R(:,[1 i_RM+1]);
% W{3,3}    = Rw(:,[i_RM]);

glo.Tbp = 3; % shift model output for brood-pouch delay (rather than change the data set)

% survivors on each observation time (d), concentrations in uM
DATA{1,4} = S(:,[1 i_SC+1 i_SA+1]);
DATA{2,4} = S(:,[1 i_SB+1]);
% DATA{3,4} = S(:,[1 i_SM+1]);

% scaled damage
DATA{1,5} = [0]; % there are never data for this state
DATA{2,5} = [0]; % there are never data for this state
% DATA{3,5} = [0]; % there are never data for this state

% Note: optionally, add some info to the MAT filename. The MAT filename
% will already include MoA and feedbacks, but if you want to try other
% things as well (changing opt, calibrating on the validation data, etc)
% it can be helpful to change the name to use this script but not 
% overwrite previous MAT files.

% glo.basenm  = [mfilename,'_CAL1']; % remember the filename for THIS file for the plots
save([glo.basenm,'_DATA'],'DATA','W') % save MAT file with data set 
% Saving the data set is handy to allow for a simple reconstruction of the
% calibrations, without needing to define the data again, in the same way.
% Note that if you use DATAx and Wx, you need to save them as well!

%% Create a table with nicer labels for the legends
% Creating a Matlab table in glo.LabelTable will replace the automatically
% generated labels (from glo.leglab1/2 and the scenario identifiers) with a
% dedicated text for each scenario.

glo.LabelTable = table(Scenario,Label); % create a Matlab table for the labels

%% Call the Matlab GUI open-file element to load MAT file

[conf_type,~,par_saved] = select_pred([1 1 2]); 

% Note that this function, with these settings, loads a pre-saved MAT file.
% It loads the parameters in par_saved, as well as all settings used to
% generate the MAT file (in glo) apart from the exposure profile (which we
% have already defined above). This only works when the MAT file was
% generated by BYOM v.6.0 BETA 7 or newer as older versions won't save glo.
% When using an old MAT file (pre BYOM v6), make sure to set glo.Tbp to the
% correct value if it needs to be >0. Also, the data set is loaded, when it
% was saved in the calibration file with the same name (..._DATA.mat).

glo = rmfield(glo,'mat_nm'); % but remove the name of the MAT file as we'll make a new one

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios. This needs to go after select_pred, otherwise
% the saved one is always used.

X0mat      = [idC idA idB]; % the scenarios (here identifiers for BOTH SINGLE EXPOSURES!) 
% X0mat      = [idC idA idB idM]; % the scenarios (here identifiers for BOTH MIXTURE AND SINGLE EXPOSURES!) 
X0mat(2,:) = 0; % initial values state 1 (scaled damage A)
X0mat(3,:) = 0; % initial values state 2 (body length, initial value overwritten by L0)
X0mat(4,:) = 0; % initial values state 3 (cumulative reproduction)
X0mat(5,:) = 1; % initial values state 4 (survival probability)
X0mat(6,:) = 0; % initial values state 1 (scaled damage B)

% Put the position of the various states in globals, to make sure that the
% correct one is selected for extra things (e.g., for plotting in
% plot_tktd, in call_deri for accommodating 'no shrinking', for population
% growth rate).
glo.locD = [1 5]; % location of scaled damage in the state variable list
glo.locL = 2; % location of body size in the state variable list
glo.locR = 3; % location of cumulative reproduction in the state variable list
glo.locS = 4; % location of survival probability in the state variable list

glo.mix_fact  = 10; % factor for scenario IDs of the mixture analysis
% This needs to be redefined, as glo from the saved set was used by the
% call to select_pred.

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. All global
% parameters, as part of the structure glo, are also loaded by
% select_pred already. Feel free to overwrite them with other values.
  
par = par_saved; % the saved parameter structure is used as basis
% only the fit mark for the basic parameters is modified below; for the tox
% parameters, the ranges from the MAT file must be taken as starting point.

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.L0(2) = 1; % body length at start experiment (mm)
par.Lp(2) = 1; % body length at puberty (mm)
par.Lm(2) = 1; % maximum body length (mm)
par.rB(2) = 1; % von Bertalanffy growth rate constant (1/d)
par.Rm(2) = 1; % maximum reproduction rate (#/d)
par.f(2)  = 0; % scaled functional response
par.hb(2) = 1; % background hazard rate (d-1)
par.a(2)  = 0; % coefficient for Weibull backgound hazard (-)

% extra parameters for special situations
par.Lf(2)   = 0; % body length at half-saturation feeding (mm)
par.Lj(2)   = 0;  % body length at end acceleration (mm)
par.Tlag(2) = 0;  % lag time for start development

ind_tox = find(strcmp(fieldnames(par),'Tlag'))+1; % index where tox parameters start
% Toxicity parameters are those after Tlag. This has to be done in this way
% as we use the saved par to provide initial values and ranges.

% % When using the parameter-space explorer, start values and ranges are not
% % needed (filled later by startgrid_debtox); only the fit/fix mark is
% % relevant here. But make sure that the value in the first column is within
% % the bounds (and not zero for log-scale parameters).
% par.kdA  = [0.62   1 0.01  10 0]; % dominant rate constant (d-1)
% par.kdB  = [0.14   1 0.01  10 0]; % dominant rate constant (d-1)
% par.zb   = [0.17   1 0.01   1 1]; % effect threshold energy budget ([C])
% par.bb   = [45     1  0.1 1e3 0]; % effect strength energy-budget effect (1/[C])
% par.zs   = [0.32   1 0.01   1 1]; % effect threshold survival ([C])
% par.bs   = [0.26   1 0.01 1e3 0]; % effect strength survival (1/([C] d))
% par.WB   = [0.95   1  0.2   5 1]; % weight factor to translate B to A (-)
% par.IAB  = [1      0 -100 100 1]; % interaction factor on damage addition (-)

% After optimisation, copy-paste relevant lines (fitted parameters) from screen below!
 

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% These are already loaded from the saved MAT file, we have an additional
% state here (as we have two damage states). And we might want to make them
% a bit shorter.

% specify the y-axis labels for each state variable
glo.ylab{1} = ['damage A (',char(181),'M)'];
glo.ylab{2} = 'body length (mm)';
if isfield(glo,'Tbp') && glo.Tbp > 0
    glo.ylab{3} = ['cumul. repro. (',num2str(glo.Tbp),'d)'];
else
    glo.ylab{3} = 'cumul. repro.';
end
glo.ylab{4} = 'survival frac. (-)';
glo.ylab{5} = ['damage B (',char(181),'M)'];

% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
% Note: the legend labels will not be used when we make a glo.LabelTable

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. 
% 
% NOTE: for this package, the options useode and eventson in glo will not
% be functional: the ODE solver is always used, and the events function as
% well.
 
% -------------------------------------------------------------------------
% Configurations for the ODE solver
glo.stiff = [0 3]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Second argument is for default sloppy (0), normally tight (1), tighter
% (2), or very tight (3) tolerances. Use 1 for quick analyses, but check
% with 3 to see if there is a difference! Especially for time-varying
% exposure, there can be large differences between the settings!
glo.break_time = 0; % break time vector up for ODE solver (1) or don't (0)
% Note: for mixture analysis, glo.break_time is not implemented yet! This
% needs more thought as there could be two time-varying exposure scenarios
% (for chemical A and B) at the same time.
% -------------------------------------------------------------------------

opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_plot.bw      = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.repls   = 0; % set to 1 to plot replicates, 0 to plot mean responses
basenm_rem       = glo.basenm; % remember basename as automatic_runs may modify it!

% Note: using opt_optim.type = 4 (parspace explorer) would be served by
% setting a few more options. These settings can be ignored for other
% optimisation routines, but skip_sg must be defined before running
% automatic_runs.
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);
opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots of parameter space to monitor progress
opt_optim.ps_profs = 0; % when set to 1, makes profiles and additional sampling for parameter-space explorer
opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS
skip_sg            = 1; % set to 1 to skip startgrid completely (use ranges in <par> structure)

% The switch fit_tox is used to select whether to fit the control treatment
% (0), the toxicity treatments (1; the control is shown as well), or both
% together (2). A good strategy is to fit the basic parameters to the
% control data first (fit_tox=0), copy the best values into the parameter
% matrix above, and then fit the toxicity parameter to the complete data
% set (fit_tox=1). The code below automatically does that.

% ===== FITTING CONTROLS ==================================================
% Simplex fitting works fine for control data. Note that hb is fitted along
% with the other parameters.

% opt_optim.type = 1; % optimisation method: 1) default simplex, 4) parspace explorer
% fit_tox = [0 0 0.1]; % fit control parameters in data set
% par_out = automatic_runs_basic(fit_tox,par,ind_tox,skip_sg,opt_optim,opt_plot);
% % script to run the calculations and plot, automatically
% par = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par
% 
% % Note: there is no strict need to fit the controls again. We could also
% % directly use the values from the saved parameter set. 

% ===== FITTING TOX DATA ==================================================
% Here, we can also use the parameter-space explorer for fitting the
% treatments. For the explorer, note that, with fit_tox=1, the tox
% parameters in par are replaced (when fitted) with estimates based on the
% data set using startgrid_debtox (called in automatic_runs). This is still
% quite experimental, so you may need to restart with manually-adapted
% ranges! Furthermore, this is really slow ... (the parallel toolbox really
% helps here!)

opt_optim.type = 4; % optimisation method 1) simplex, 4 parameter-space explorer
fit_tox = [1]; % fit tox parameters in data set
par_out = automatic_runs_basic(fit_tox,par,ind_tox,skip_sg,opt_optim,opt_plot);
disp_settings_debtox2019 % display some information on the settings on screen

opt_tktd.repls   = 0; % plot individual replicates (1) or means (0)
opt_tktd.transf  = 1; % set to 1 to calculate means and SEs including transformations
opt_tktd.obspred = 1; % plot predicted-observed plots (1) or not (0)
opt_tktd.max_exp = 1; % set to 1 to maximise exposure/damage plots on exposure rather than damage

% dedicated TKTD plots; these plots are more readable, especially when plotting CIs
opt_tktd.repls = 0; % plot individual replicates (1) or means (0)
if opt_optim.type == 4    % when the parspace explorer was used, we can plot with CIs ...
    opt_conf.type    = 3; % make intervals from parspace explorer
    opt_conf.lim_set = 2; % use limited set of n_lim points (1) or outer hull (2) to create CIs
    plot_tktd(par_out,opt_tktd,opt_conf);
else
    plot_tktd(par_out,opt_tktd,[]); % leave the options opt_conf empty to suppress all CIs for these plots
    % Note: plotting CIs requires a sample as saved by the various
    % methods available in BYOM (see examples directory).
end
