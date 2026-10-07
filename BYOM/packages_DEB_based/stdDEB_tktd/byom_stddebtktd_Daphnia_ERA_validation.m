%% BYOM, validation with stdDEB model: byom_stddebtktd_Daphnia_ERA_validation.m
%
% *Table of contents*

%% About
% * Author     : Tjalling Jager
% * Date       : June 2022
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* standard DEB model in terms of powers (energy fluxes in J/d)
% and body length. In this package, toxicant stress is included. The model
% equations are those from the 'Tromso coffee mug' ;-). Only change is that
% the ODE for structure is phrased in volumetric length rather than volume.
% The TKTD model is lifted from DEBtox2019. The 2023 publication of Jager
% et al in Ecological Modelling contains the full details:
% <https://doi.org/10.1016/j.ecolmodel.2022.110187>.
%
% *This script:* The water flea _Daphnia magna_ exposed to fluoranthene;
% data set from Jager et al (2010),
% <http://dx.doi.org/10.1007/s10646-009-0417-z>. Published as case study in
% Jager & Zimmer (2012). Simultaneous fit on growth, reproduction and
% survival. Here the raw data on individual females is used, such that the
% special functions for censoring the reproduction data can be used. This
% script demonstrates validation. Making use of the parameter-space
% explorer. The controls of the validation data set will be fitted, but the
% tox parameters and all settings taken from the saved MAT file. This is an
% illustration and not a validation as the same data set is used as for
% calibration!
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
diary off           % turn off the diary function (if it is accidentaly on)
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
% automatic_runs.

TR   = 0.5; % transformations for continuous data (1 is no transformation)
opt  = 2; % select an option (opt=1 is generally recommended for Daphnia, but requires closer look at the data set)

% Options to deal with repro data set:
% 0) Check whether we can use a single intermoult period for the entire
%    data set. Screen output will show mean intermoult times and brood 
%    sizes across the replicates, as function of brood number and treatment.
% 1) Cumulate reproduction, but remove the time points with zero
%    reproduction. Good for clutch-wise reproduction.
% 2) Cumulate reproduction, but don't remove zeros. Good for continuous
%    reproduction or when animals are not followed individually.
% 3) Shift neonate release back to the previous moult. When this option is
%    used, don't shift the model predictions with <glo.Tbp>: the data now
%    represent egg production rather than neonate release.

% The data set is specified in a separate function. This is done so that we
% can easily combine or swap data sets in calibration and validation. Last
% entry is the number of the data set (if there is more than one, the
% scenario identifiers of the next one will have 100 added, the next one
% 200 added, etc.).

% Put the position of the various states in globals, to make sure correct
% one is selected for extra things (e.g., to prevent shrinking in
% call_deri, and for using plot_tktd). This is placed here as it is used in
% the data file that is called below, to put the data in the right position
% of the DATA structure.
glo.locL = 3; % location of body size in the state variable list
glo.locR = 4; % location of cumulative reproduction in the state variable list
glo.locD = 5; % location of damage in the state variable list 
glo.locS = 6; % location of survival probability in the state variable list

fnames = {'data_FLU_Dmagna'};
read_datafiles(fnames,TR,opt); % function defines DATA, W and glo.LabelTable

% Note: optionally, add some info to the MAT filename. 
glo.basenm  = [mfilename,'_FLU']; % remember the filename for THIS file for the plots

if opt == 0 % check intermoult duration (and mean brood size)
    return % we need to stop to check the results (no data are created)
end 
% Note: glo.Tbp will be read from MAT file below.

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

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios. Note that for stdDEB, the initial values cannot
% be randomly chosen. We need to simulate the individual up to the point
% where we start the model analysis (e.g., at birth). Therefore, most of
% the values in X0mat will not be used in the analysis.
% 
% Note that the complete X0mat for the calibration is also loaded with the
% MAT file. However, we need a new one for validation, so it is overwritten
% below.

X0mat = [glo.LabelTable.Scenario]'; % the scenarios (here identifiers) 
X0mat(2,:) = 0;   % initial reserve energy (J) NOT USED
X0mat(3,:) = 0;   % initial maturity level (J) NOT USED
X0mat(4,:) = 0;   % initial volumetric length (cm) NOT USED
X0mat(5,:) = 0;   % initial cumul. repro NOT USED
X0mat(6,:) = 0;   % initial values scaled damage
X0mat(7,:) = 1;   % initial values survival probability

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. All global
% parameters, as part of the structure glo, are also loaded by
% select_pred already. Feel free to overwrite them with other values.

% glo.T      = 20 + 273.15; % actual temperature to Kelvin
% glo.len    = 2;       % switch to fit physical length (0=wwt, 1=phys. length, 2=phys. length, no shrinking) (used in call_deri.m)

% =========================================================================
%      Parameters for Daphnia magna (Bas Kooijman, Andre Gergs. 2019. 
%      AmP Daphnia magna, version 2019/03/16.)
% =========================================================================

% par = par_saved; % take ALL parameters from the saved set (incl. basic parameters!)
% % This way, the calibration controls will define the basic parameters for validation as well.

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.TA   = [6400     0 1000  1e6 1]; % Arrhenius temperature (K)
par.spAm = [313.169  0    0  1e6 1]; % max. surface-specific assimilation rate (J/(cm2 d))
par.spM  = [1200     0    0  1e6 1]; % volume-specific somatic maintenance costs (J/(cm3 d))
par.spT  = [0        0    0  1e6 1]; % surface-specific maintenance costs (J/(cm2 d))
par.kJ   = [0.2537   0    0  1e6 1]; % maturity maintenance rate constant (1/d)
par.EG   = [4400     0    0  1e6 1]; % volume-specific costs for growth (J/cm3)
par.EHb  = [0.05464  0 1e-6  1e6 1]; % maturity level at birth (J)
par.EHj  = [0        0    0  1e6 1]; % maturity level at metamorphosis (J)
par.EHp  = [1.09     1 1e-6  1e6 1]; % maturity level at puberty (J)
par.v    = [0.1858   0    0  1e6 1]; % energy conductance (cm/d)
par.kap  = [0.5809   1 0.01 0.99 1]; % allocation fraction to soma (-)
par.kapR = [0.95     0 0.01 0.99 1]; % reproduction efficiency (-)
par.f    = [1        1    0    2 1];   % scaled food density (-)
par.hb   = [0.01     1 1e-3 0.07 0]; % background hazard rate (d-1)
par.a    = [1        0  0.1   10 0]; % coefficient for Weibull backgound hazard (-) (1 is off)
% par.Lw0  = [0.8      0    0  1e6 1]; % starting at a *physical* length L0>Lb (when glo.len=0, use wet weight)

glo.names_sep = {}; % names of parameters that can differ between data sets
% glo.names_sep = {'f';'L0'}; % names of parameters that can differ between data sets
% par.L01   = [0.9 1 0.5  1.5  1]; % body length at start experiment (mm)
% par.f1    = [1     0 0    2    1]; % scaled functional response
% par.f2    = [0.85  1 0    2    1]; % scaled functional response

% NOTE: some parameters for the AmP entry need to be refitted to obtain a
% closer correspondence to the controls of this data set.

% =========================================================================
%      TKTD Parameters for this data set
% =========================================================================
ind_tox = length(fieldnames(par))+1; % index where tox parameters start
% the parameters below this line are all treated as toxicity parameters!

% Values below are not used; they are overwritten by the values from the
% selected MAT file.
par.kd   = [0.08   1 0.01  10 0]; % dominant rate constant (d-1)
par.zb   = [0.10   1 0    1e6 1]; % effect threshold energy budget ([C])
par.bb   = [72     1 1e-6 1e6 0]; % effect strength energy-budget effect (1/[C])
par.zs   = [0.24   1 0    1e6 1]; % effect threshold survival ([C])
par.bs   = [0.64   1 1e-6 1e6 0]; % effect strength survival (1/([C] d))

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% All plot labels have been loaded from file using select_pred.

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
% well. Also note that glo.stiff and glo.break_time are loaded from the MAT
% file as well, but overwritten below. For pulsed exposure scenarios, we
% best use glo.break_time=1, and possibly a different solver than for
% constant exposure.

% -------------------------------------------------------------------------
glo.stiff = [0 3]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Second argument is for default sloppy (0), normally tight (1), tighter
% (2), or very tight (3) tolerances. Use 1 for quick analyses, but check
% with 3 to see if there is a difference! Especially for time-varying
% exposure, there can be large differences between the settings!
glo.break_time = 0; % break time vector up for ODE solver (1) or don't (0)
% Note: breaking the time vector is a good idea when the exposure scenario
% contains discontinuities. Don't use for continuous splines (type 1) as it
% will be much slower. For FOCUS scenarios (high time resolution), breaking
% up is not efficient and does not appear to be necessary.
% -------------------------------------------------------------------------

opt_optim.fit     = 1; % fit the parameters (1), or don't (0)
opt_plot.bw       = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot    = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.repls    = 0; % set to 1 to plot replicates, 0 to plot mean responses
basenm_rem        = glo.basenm; % remember basename as automatic_runs may modify it!

% Select what to fit with fit_tox (this is a 3-element vector). NOTE: use
% identifier c=0 for regular control, and c=0.1 for solvent control. When
% entering more than one data set, use 100 and 100.1 for the controls of
% the second data set (and 101, 102 ... for the treatments), 200 and 200.1
% for the controls of the third data set. etc. 
% 
% First element of fit_tox is which part of the data set to use:
%   fit_tox(1) = -2  comparison between control and solvent control (c=0 and c=0.1)
%   fit_tox(1) = -1  control survival (c=0) only
%   fit_tox(1) = 0   controls for growth/repro (c=0) only, but not for survival
%   fit_tox(1) = 1   all treatments, but, when fitting, keep all control parameters fixed; 
%               run through all elements in MOA and FEEDB sequentially and 
%               provide a table at the end (plots are made incl. control)
% 
% Second element of fit_tox is whether to fit or only to plot:
%   fit_tox(2) = 0   don't fit; for standard optimisations, plot results for
%               parameter values in [par], for parspace optimisations, use saved mat file.
%   fit_tox(2) = 1   fit parameters
% 
% Third element is what to use as control (fitted for fit_tox(1) = -1 or 0)
% (if this element is not present, only regular control will be used)
%   fit_tox(3) = 1   use regular control only (identifier 0)
%   fit_tox(3) = 2   use solvent control only (identifier 0.1)
%   fit_tox(3) = 3   use both regular and control (identifier 0 and 0.1)
% 
% The strategy in this script is to fit hb to the control data first. Next,
% fit the basic parameters to the control data. Finally, fit the toxicity
% parameter to the complete data set (keeping basic parameters fixed. The
% code below automatically keeps the parameters fixed that need to be
% fixed. 
% 
% MOA: Mode of action of toxicant as set of switches
% [assimilation/feeding, maintenance costs (somatic and maturity), growth costs, repro costs] 
% [1 0 0 0 0]   assimilation/feeding
% [0 1 0 0 0]   costs for maintenance 
% [0 0 1 1 0]   costs for growth and reproduction
% [0 0 0 1 0]   costs for reproduction
% [0 0 0 0 1]   hazard for reproduction
% 
% FEEDB: Feedbacks to use on damage dynamics as set of switches
% [surface:volume on uptake, surface:volume on elimination, growth dilution, losses with reproduction] 
% [1 1 1 1]     all feedbacks 
% [1 1 1 0]     classic DEBtox (no losses with repro)
% [0 0 1 0]     damage that is diluted by growth
% [0 0 0 0]     damage that is not diluted by growth

% These are the MoA's and feedback configurations that will be run
% automatically when fit_tox(1) = 1. For the controls, we can call
% automatic_runs with empty ones.
MOA   = [];
FEEDB = [];

% ===== FITTING CONTROLS ==================================================
% Simplex fitting works fine for control data
opt_optim.type = 1; % optimisation method: 1) default simplex, 4) parspace explorer

% NOTE: if we fit the controls from the validation data set, we need to
% remove the loaded glo.stdDEB_start from calibration! This element was
% also loaded with glo by the call to select_pred. Since the start values
% depend on the basic parameters, they need to be re-established for the
% validation data set.
glo.stdDEB_start = {}; % make sure it is empty and not defined!

% % Compare controls in data set
% fit_tox = [-2 1 3];
% automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot);
% % script to run the calculations and plot, automatically
% % NOTE: I think the likelihood-ratio test is often too strict, and that
% % using both controls should be the default situation.
% return

% Fit control parameters (not hb) in data set
glo.stdDEB_start = {}; % make sure it is empty and not defined!
fit_tox = [0 1 3]; % use both controls
par_out = automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot);
% script to run the calculations and plot, automatically
par = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par

% Next, we need to clear the zero-variate data points, so they are not used
% again for fitting the treatments.
glo.zvd = []; % this should be enough, since it makes the output zvd of call_deri empty as well

% Fit hb in data set
fit_tox = [-1 1 3]; % use both controls
par_out = automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot); 
% script to run the calculations and plot, automatically
par = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par

% NOTE: even when you don't want to fit any of the control parameters, it
% is good to run with fit_tox(1)=0. This will make sure that the start
% values are defined (that won't change in the analysis when fitting only
% the tox parameters), which will speed up things considerably.
% 
% NOTE: I start fitting the control parameters as that should lead to a
% proper definition of the initial values in glo.stdDEB_start!

%% Get par from SAVED set, and copy our par into that!

par_new = copy_par(par,par_saved); % copy the fitted parameters from the saved set into a new par!
% This is needed since the saved MAT can be made with more (or less)
% parameters, since separate parameters can be used for different data
% sets. This function also takes over the fit marks from the saved set.

% calc_and_plot(par_new,opt_plot); % then we can plot with the new par as well

%% Plot results with confidence intervals
% The following code can be used to make plots with confidence intervals.
% Options for confidence bounds on model curves can be set using opt_conf
% (see prelim_checks). The plot_tktd function makes multiplots for the
% effects data, which are more readable when plotting with various
% intervals.

opt_conf.type     = conf_type; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set  = 2; % use limited set of n_lim points (1), or outer hull (2, not for Bayes) to create CIs
% Note: using a limited set is good for test design, but for validation,
% using the outer hull is a better idea (but will be slower).
opt_tktd.repls    = 0; % plot individual replicates (1) or means (0)
opt_tktd.transf   = 1; % set to 1 to calculate means and SEs including transformations
opt_tktd.max_exp  = 0; % set to 1 to maximise exposure/damage plots on exposure rather than damage
opt_tktd.preds    = 0; % set to 1 to only plot predictions from X0mat without data
opt_tktd.addzero  = 1; % set to 1 to always add a concentration zero to X0mat
opt_tktd.obspred  = 0; % plot predicted-observed plots (1) or not (0)
opt_tktd.sppe     = 0; % set to 1 to calculate SPPEs (relative error at end of test)
% Note: this last option will produce something like the SPPE
% (opt_tktd.obspred must be set to 1 to calculate SPPE and other metrics).

fit_tox = [1 1 3]; % use both controls
% change X0mat to avoid plotting controls not used for fitting (assumes
% that regular control has identifier 0 and solvent control 0.1)
switch fit_tox(3)
    case 1 % remove solvent control
        X0mat(:,X0mat(1,:)==0.1) = [];
    case 2 % remove regular control
        X0mat(:,X0mat(1,:)==0) = [];
    case 3
        % leave all in
end

% % Option 1: use basic and tox parameters from the saved set.
% plot_tktd([],opt_tktd,opt_conf); 
% % Note that first input is empty: this is for the parameter structure,
% % which is then obtained from the saved sample.

% Option 2: use fitted (tox) parameters from saved set, but the fixed
% (basic) parameters as presented in the structure <par>. The option in
% opt_conf makes sure this is handled properly.
opt_conf.use_par_out = 1; % set to 1 to use par as entered into the plotting function for CIs, rather than from saved set
plot_tktd(par_new,opt_tktd,opt_conf); % note that par_new contains the tox parameters from par_saved (see use of copy_par above)

% % Alternative plot, using calc_and_plot with some trickery!
% opt_conf.use_par_out = 1; % set to 1 to use par as entered in this function for CIs, rather than from saved set
% % this option tells calc_conf to NOT use the fixed parameters from the saved set!
% [out_conf,par_mod] = calc_conf(par_new,opt_conf); % par_mod contains fixed parameters from par, but fitted parameters from saved set
% calc_and_plot(par_mod,opt_plot,out_conf);
