%% BYOM, ERA-special with simple-compound model: byom_debtox_folsomia_Cd.m
%
% *Table of contents*

%% About
% * Author     : Tjalling Jager
% * Date       : November 2022
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_debtox2019.html>
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
% (see <http://www.debtox.info/book_debkiss.html>) provides a description
% of the model, as well as the publication of Jager in Ecological
% Modelling: <https://doi.org/10.1016/j.ecolmodel.2019.108904>.
%
% *This script:* The springtail _Folsomia candida_ exposed to cadmium;
% data set from Jager et al (2004), <http://dx.doi.org/10.1021/es0352348>.
% Simultaneous fit on growth and reproduction. This script was used for the
% case study in the DEBkiss e-book (version 3.0).
% 
%  Copyright (c) 2012-2022, Tjalling Jager.
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
% automatic_runs.

% scaled damage
DATA{1} = [0]; % there are never data for this state
        
% body size (mg fwt) on each observation time (d), concentrations in mg/kg
% food
DATA{2} = [0.5	0	64	139	300	646	1392	3000
0	0.0018	0.0018	0.0018	0.0018	0.0018	0.0018	0.0018
16	0.072906667	0.077726667	0.074526667	0.06714	0.0716	0.072406667	0.041553333
23	0.165593333	0.172106667	0.153926667	0.14354	0.115153333	0.130886667	0.060193333
30	0.225493333	0.23206	0.193986667	0.18548	0.16162	0.15982	0.080446667
37	0.241893333	0.244306667	0.21086	0.199233333	0.187413333	0.160133333	0.093733333
44	0.24978	0.260906667	0.222533333	0.214573333	0.208013333	0.17944	0.09874
51	0.263126667	0.270673333	0.251053333	0.246033333	0.241233333	0.18244	0.088906667
58	0.277366667	0.266686667	0.25828	0.246066667	0.22772	0.180666667	0.094993333
65	0.293493333	0.2635	0.249033333	0.23784	0.237293333	0.184426667	0.096253333
72	0.271473333	0.278586667	0.264326667	0.253433333	0.244773333	0.186853333	0.102126667
86	0.309766667	0.257766667	0.259045455	0.277572727	0.252081818	0.208746667	0.119946154
100	0.298733333	0.258233333	0.268128571	0.282357143	0.233766667	0.227077778	0.11421
114	0.286475	0.2071	0.251925	0.25146	0.214933333	0.1973	0.109133333];

DATA{2}(2:end,2:end) = (DATA{2}(2:end,2:end)).^(1/3); % mg to mm volumetric length (assume density 1 mg_wet/mm3)

W{2} = [20	20	20	20	20	20	20
15	15	15	15	15	15	15
15	15	15	15	15	15	15
15	15	15	15	15	15	15
15	15	15	15	15	15	15
15	15	15	15	15	15	15
15	15	15	15	15	15	15
15	15	15	15	15	15	15
15	15	15	15	15	15	15
15	15	15	15	15	15	15
12	9	11	11	11	15	13
6	6	7	7	6	9	10
4	3	4	5	3	6	9];

% cumulative reproduction (nr. offspring) (mm) on each observation time
% (d), concentrations in mg/kg food
DATA{3} = [0.5	0	64	139	300	646	1392	3000
0	0	0	0	0	0	0	0
22	31.06666667	31.13333333	32.73333333	33.2	30.46666667	34.26666667	13.93333333
29	127.6666667	136.5619048	120.3333333	115.7714286	107.6	106.4666667	31.8
36	255.8	273.4849817	214.1333333	216.6175824	199.1333333	159	58.66666667
43	347.0142857	357.8849817	276.2	288.4747253	256.6666667	189.7857143	87.66666667
50	376.7285714	430.2849817	305.9142857	340.6747253	324.8095238	244.3190476	107.1333333
57	451.0922078	490.018315	365.9142857	430.9054945	400.5595238	277.5498168	116.6333333
64	490.7285714	526.018315	426.9142857	492.443956	419.0210623	318.8575092	132.6333333
% 71	543.0922078	592.6546787	458.0681319	538.5208791	471.790293	351.0113553	145.4904762
% 78	588.1922078	627.8769009	514.8181319	605.8542125	524.6993839	384.9344322	168.4904762
% 85	654.8922078	638.2102342	535.6181319	625.9653236	562.7993839	420.6267399	191.9520147
% 92	716.1422078	678.2102342	568.760989	662.7153236	598.0850982	466.8994672	213.1186813
% 99	758.5422078	735.2102342	577.260989	675.0010379	632.9184316	481.0105783	229.755045
% 106	789.7422078	767.2102342	600.660989	701.1677045	662.4184316	520.5105783	236.380045
% 113	840.9922078	796.7102342	600.660989	714.4177045	698.4184316	546.0105783	267.380045
% 120	888.3255411	805.7102342	620.3276557	734.9177045	698.4184316	562.0105783	276.6022672
];

W{3} = [15	15	15	15	15	15	15
15	15	15	15	15	15	15
15	14	15	14	15	15	15
15	13	15	13	15	15	15
14	15	15	14	15	14	15
14	15	14	15	14	15	15
11	15	13	13	12	13	14
11	12	13	13	13	13	14
% 11	11	13	13	13	13	14
% 10	9	12	12	11	13	13
% 10	6	10	9	10	13	13
% 8	6	7	8	7	11	12
% 5	4	6	7	6	9	11
% 5	4	5	6	4	6	8
% 4	4	3	4	2	6	9
% 3	1	3	4	1	4	9
];

DATA{4} = []; % survival not considered

% Put the position of the various states in globals, to make sure that the
% correct one is selected for extra things (e.g., for plotting in
% plot_tktd, in call_deri for accommodating 'no shrinking', for population
% growth rate).
glo.locD = 1; % location of scaled damage in the state variable list
glo.locL = 2; % location of body size in the state variable list
glo.locR = 3; % location of cumulative reproduction in the state variable list
glo.locS = 4; % location of survival probability in the state variable list

% Note: optionally, add some info to the MAT filename. The MAT filename
% will already include MoA and feedbacks, but if you want to try other
% things as well (changing opt, calibrating on the validation data, etc)
% it can be helpful to change the name to use this script but not 
% overwrite previous MAT files.

% glo.basenm  = [mfilename,'_CAL']; % remember the filename for THIS file for the plots
save([glo.basenm,'_DATA'],'DATA','W') % save MAT file with data set 
% Saving the data set is handy to allow for a simple reconstruction of the
% calibrations, without needing to define the data again, in the same way.
% Note that if you use DATAx and Wx, you need to save them as well!

glo.Tbp = 0; % no brood pouch for Folsomia candida

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat      = [0 64 139 300 646 1392 3000]; % the scenarios (here identifiers) 
X0mat(2,:) = 0; % initial values state 1 (scaled damage)
X0mat(3,:) = 0; % initial values state 2 (body length, initial value overwritten by L0)
X0mat(4,:) = 0; % initial values state 3 (cumulative reproduction)
X0mat(5,:) = 1; % initial values state 4 (survival probability)

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 
  
% global parameters as part of the structure glo
glo.FBV    = 0.008;    % dry weight egg as fraction of structural body weight (-) (for losses with repro; approx. for F. candida)
glo.KRV    = 1;       % part. coeff. repro buffer and structure (kg/kg) (for losses with reproduction)
glo.kap    = 0.8;     % approximation for kappa (for starvation response)
glo.yP     = 0.8*0.8; % product of yVA and yAV (for starvation response)
glo.Lm_ref = 0.7;    % reference max length for scaling rate constants
glo.len    = 1;       % switch to fit length 1) with shrinking, 2) without shrinking (used in call_deri.m)
% NOTE: the settings above are species specific! Make sure to use the same
% settings for validation and prediction as for calibration! 
% 
% NOTE: For arthropods, one would generally want to fit the model without
% shrinking (since the animals won't shrink in length). However, here we
% have the data as volumetric length (cubic root of fresh weight).

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.L0   = [0.12  1 0.05 0.5  1]; % body length at start experiment (mm)
par.Lp   = [0.37  1 0.1  1.0  1]; % body length at puberty (mm)
par.Lm   = [0.66  1 0.1  1.5  1]; % maximum body length (mm)
par.rB   = [0.1   1 0.01 0.2  1]; % von Bertalanffy growth rate constant (1/d)
par.Rm   = [14    1 5    30   1]; % maximum reproduction rate (#/d)
par.f    = [1     0 0    2    1]; % scaled functional response
par.hb   = [0.001 0 1e-3 0.07 0]; % background hazard rate (d-1)
par.a    = [1     0 0.1  10   0]; % coefficient for Weibull backgound hazard (-)
% - Note 1: it does not matter whether length measures are entered as actual
% length or as volumetric length (as long as the same measure is used
% consistently). 
% - Note 2: hb is fitted on log-scale. This is especially helpful for
% fit_tox(1)=-2 when hb is fitted along with the other parameters. The
% simplex fitting has trouble when one parameter is much smaller than
% others. 
% - Note 3: the min-max ranges in par are appropriate for Daphnia magna. For
% other species, these ranges and starting values need to be modified.

% % When using multiple data sets, we can let some parameters differ.
% glo.names_sep = {'f';'L0'}; % names of parameters that can differ between data sets
% par.L01   = [0.9 1 0.5  1.5  1]; % body length at start experiment (mm)
% par.f1    = [0.9 1 0    2    1]; % scaled functional response
% % Note: using separate parameters for separate data sets requires using
% % specific identifiers, and using exposure scenarios with make_scen (also
% % for constant exposure and for the controls)!

% extra parameters for special situations
par.Lf   = [0 0 0 1e6 1];    % body length at half-saturation feeding (mm)
par.Lj   = [0.37 1 0 0.6 1]; % body length at end acceleration (mm)
par.Tlag = [0 0 0 1e6 1];    % lag time for start development

ind_tox = length(fieldnames(par))+1; % index where tox parameters start
% the parameters below this line are all treated as toxicity parameters!

% When using the parameter-space explorer, start values and ranges are not
% needed (filled later by startgrid_debtox); only the fit/fix mark is
% relevant here. But make sure that the value in the first column is within
% the bounds (and not zero for log-scale parameters).
par.kd   = [0.08   1 0.01  1 0]; % dominant rate constant (d-1)
par.zb   = [0.1    1 0    80 1]; % effect threshold energy budget ([C])
par.bb   = [1e-4   1 1e-6 1e-3 0]; % effect strength energy-budget effect (1/[C])
par.zs   = [1e6   0 0    1e6 1]; % effect threshold survival ([C])
par.bs   = [1     0 1e-6 1e6 0]; % effect strength survival (1/([C] d))
 
% After optimisation, copy-paste relevant lines (fitted parameters) from screen below!

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = ['scaled damage (mg/kg)'];
glo.ylab{2} = 'volumetric length (mm)';
if isfield(glo,'Tbp') && glo.Tbp > 0
    glo.ylab{3} = ['cumul. repro. (shift ',num2str(glo.Tbp),'d)'];
else
    glo.ylab{3} = 'cumul. repro. (no shift)';
end
glo.ylab{4} = 'survival fraction (-)';

% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = ''; % legend label before the 'scenario' number
glo.leglab2 = 'mg/kg food'; % legend label after the 'scenario' number
% Note: these legend labels will not be used when we make a glo.LabelTable

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
% Note: breaking the time vector is a good idea when the exposure scenario
% is specified as type 2-4 (even for scenarios with constant exposure).
% Don't use for continuous splines (type 1) as it will be much slower. For
% FOCUS scenarios (high time resolution), breaking up is not efficient and
% does not appear to be necessary.
% -------------------------------------------------------------------------

opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_plot.bw      = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 0; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.repls   = 0; % set to 1 to plot replicates, 0 to plot mean responses
basenm_rem       = glo.basenm; % remember basename as automatic_runs may modify it!

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
% automatically when fit_tox(1) = 1. This needs to be defined here, even
% for control fits when it is not used.

MOA   = [1 1 0 0 0]; % assimilation and maintenance costs
FEEDB = [0 0 0 0]; % no feedbacks
% FEEDB = [0 0 0 0; 1 1 1 0];               % you may want to try more feedbacks ...
% FEEDB = allcomb([0 1],[0 1],[0 1],[0 1]); % or test ALL feedback configurations

% For fit_tox(1)~=1, the setting of MOA and FEEDB has no impact; they must
% be defined to prevent errors. For fit_tox(1)=1, setting is relevant. Note
% that if you run multiple configurations, the last one will remain in the
% memory, unless glo.moa and glo.feedb are redefined. 

% Note: when using the parspace explorer, the automatic_runs uses the start
% ranges as produced by startgrid_debtox. These are preliminary! It may be
% needed to tweak these, but that should then be done in startgrid_debtox
% or in automatic_runs. So watch out when the analysis runs into min-max
% bounds.

% ===== FITTING CONTROLS ==================================================
% Simplex fitting works fine for control data
opt_optim.type   = 1; % optimisation method: 1) default simplex, 4) parspace explorer
opt_plot.statsup = [glo.locD,glo.locS]; % vector with states to suppress in plotting fits
opt_plot.limax   = 1; % if set to 1, limit axes to the data set for each stage

% opt_optim.type     = 4; % optimisation method 1) simplex, 4 parameter-space explorer
% opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots of parameter space to monitor progress
% opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS (2 for extra rough)

% % Compare controls in data set
% fit_tox = [-2 1 3];
% automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot);
% % script to run the calculations and plot, automatically
% % NOTE: I think the likelihood-ratio test is often too strict, and that
% % using both controls should be the default situation.

% % Fit hb in data set
% fit_tox = [-1 1 1]; % use regular control
% par_out = automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot); 
% % par_out = automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot,opt_prof); 
% % script to run the calculations and plot, automatically
% par = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par

% Fit other control parameters (not hb) in data set
fit_tox = [0 1 1]; % use regular control
par_out = automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot);
% par_out = automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot,opt_prof);
% script to run the calculations and plot, automatically
par = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par

% return

% ===== FITTING TOX DATA ==================================================
% Use the parameter-space explorer for fitting the treatments. Note that,
% with fit_tox=1, the tox parameters in par are replaced (when fitted) with
% estimates based on the data set using startgrid_debtox (called in
% automatic_runs_debtox2019).
% 
% Note: you can use the rough settings to find better ranges and restart
% with refined settings. With opt_optim.ps_profs = 0, the algorithm will
% provide new search ranges on screen that can be directly copied-pasted
% into this script. You may comment out the fitting of the controls above.
% Make sure that skip_sg is then set to 1 to avoid startgrid_debtox to
% overwrite par again.

opt_optim.type     = 4; % optimisation method 1) simplex, 4 parameter-space explorer
opt_plot.statsup   = [glo.locS]; % vector with states to suppress in plotting fits
% Note: to use saved set for parameter-space explorer, use fit_tox(2)=0!
opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots of parameter space to monitor progress
opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS (2 for extra rough)
% -----------------------------------------------------------------------
% THIS BLOCK: SETTINGS FOR ROUGHLY GOING THROUGH MANY CONFIGURATIONS
% opt_optim.ps_rough = 2; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS (2 for extra rough)
opt_optim.ps_profs = 0; % when set to 1, makes profiles and additional sampling for parameter-space explorer
% glo.stiff          = [2 1]; % for quick exploration of many options, could try stiff solver with sloppy tolerances if default setting is slow/stuck
glo.stiff          = [0 3]; % use ode45 with very strict tolerances
skip_sg            = 0; % set to 1 to skip startgrid completely (use ranges in <par> structure)
% -----------------------------------------------------------------------
% % THIS BLOCK: SETTINGS TO REFINE FROM RANGES COPIES FROM SCREEN INTO THIS SCRIPT
% opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
% glo.stiff          = [0 3]; % use ode45 with very strict tolerances
% skip_sg            = 1; % set to 1 to skip startgrid completely (use ranges in <par> structure)
% -----------------------------------------------------------------------

fit_tox = [1 1 1]; % use regular control, and fit tox data
[par_out,best_MoaFb] = automatic_runs_debtox2019(fit_tox,par,ind_tox,skip_sg,MOA,FEEDB,opt_optim,opt_plot);
% script to run the calculations and plot, automatically
% Note: automatic_runs will return the BEST parameter set in par_out.
disp_settings_debtox2019 % display a bit more info on the settings on screen
print_par(par_out) % print the complete best parameter vector that can be copied-pasted
% This includes the control parameters (the parspace explorer itself will
% only plot the results for the fitted parameters).

% Note: these results cannot be translated to the DeEP tool since
% acceleration is included through par.Lj.

%% Plot results with confidence intervals
% The following code can be used to make plots with confidence intervals.
% Options for confidence bounds on model curves can be set using opt_conf
% (see prelim_checks). The plot_tktd function makes multiplots for the
% effects data, which are more readable when plotting with various
% intervals.
% 
% Note: the set up below plots the best fit only (when multiple combinations are tried).

opt_conf.type     = 3; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set  = 2; % use limited set of n_lim points (1) or outer hull (2, not for Bayes) to create CIs
opt_tktd.repls    = 0; % plot individual replicates (1) or means (0)
opt_tktd.transf   = 1; % set to 1 to calculate means and SEs including transformations
opt_tktd.obspred  = 1; % plot predicted-observed plots (1) or not (0)
opt_tktd.max_exp  = 0; % set to 1 to maximise exposure/damage plots on exposure rather than damage
opt_tktd.sppe     = 1; % set to 1 to calculate SPPEs (relative error at end of test)
opt_tktd.statsup  = [glo.locS]; % states to suppress from the plots (e.g., locS)
opt_tktd.lim_data = 1; % set to 1 to limit axes to data

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

% Below some tricks to allow plotting results for selected configurations.
glo.moa    = MOA(best_MoaFb(1),:); % change global for MoA
glo.feedb  = FEEDB(best_MoaFb(2),:); % change global for feedback configuration
% glo.moa    = MOA(1,:);   % change global for MoA to first one
% glo.feedb  = FEEDB(1,:); % change global for feedback configuration to first one

glo.mat_nm = [basenm_rem,'_moa',sprintf('%d',glo.moa),'_feedb',sprintf('%d',glo.feedb)];
% The glo.mat_nm specifies the filename to read the parameters and sample
% from, as needed for plot_tktd. 
glo.basenm = basenm_rem; % return the basenm to this filename 
% The plots saved will get their filename based on the name of THIS
% script (so without the specification of MoA and feedbacks). 

plot_tktd([],opt_tktd,opt_conf); 
% leave the options opt_conf empty to suppress all CIs for these plots
% leave par_out (first input) empty to read it from saved mat file. If we
% ran more configurations, we should watch out when using par_out
% (automatic_run returns the best-fitting one).
