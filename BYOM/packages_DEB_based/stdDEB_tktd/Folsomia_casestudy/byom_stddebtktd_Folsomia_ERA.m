%% BYOM, byom_stddebtktd_Folsomia_ERA.m
%
% * Author: Tjalling Jager
% * Date: June 2022
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
% *This script:* the colembolan Folsomia candida exposed to chlorpyrifos in
% food. This is a tricked-out version, using automatic_runs_debtox2019,
% which can fit different configurations sequentially. All treatments have
% a scenario prepared with make_scen. This script was used for calibration
% of the case study in the manuscript of Jager et al (2023).
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

% wet body weight (mg) on each observation time (d), concentrations
% in mg/kg food. 
DATA{3} = [0	0	0.1	1 2 3 4 5
    0   NaN NaN NaN NaN NaN NaN NaN
    11	0.034435714	0.031236364	0.031653846	0.030528571	0.032414286	0.03024	0.032857143
    16	0.071646667	0.066173333	0.077033333	0.063886667	0.06492	0.06834	0.06484
    23	0.141813333	0.120026667	0.12588	0.127233333	0.115426667	0.12238	0.119661538
    30	0.197253333	0.175433333	0.17492	0.19766	0.165713333	0.15356	0.166169231
    37	0.220106667	0.207553333	0.196426667	0.203186667	0.18576	0.192426667	0.191833333
    44	0.233993333	0.217633333	0.213446667	0.22014	0.201006667	0.21082	0.206666667
    51	0.260473333	0.223593333	0.228893333	0.23876	0.226526667	0.216653333	0.215209091
    58	0.26656	0.235986667	0.236578571	0.249026667	0.215546667	0.22774	0.223554545
    65	0.276306667	0.265233333	0.253357143	0.265046667	0.25124	0.248173333	0.249081818
    79	0.276706667	0.27322	0.270514286	0.28732	0.261893333	0.261938462	0.24556
    93	0.287323077	0.280353846	0.2852	0.309384615	0.284664286	0.278958333	0.2519
    107	0.273115385	0.26517	0.296146154	0.306445455	0.293938462	0.2788125	0.2693
    121	0.29285	0.28484	0.309818182	0.282125	0.290608333	0.2887	0.2779
    138	0.306171429	0.29865	0.32534	0.296733333	0.29668	0.30742	0.286257143];

W{3} = [0 0 0 0 0 0 0
    14	11	13	14	14	15	14
    15	15	15	15	15	15	15
    15	15	15	15	15	15	13
    15	15	15	15	15	15	13
    15	15	15	15	15	15	12
    15	15	15	15	15	15	12
    15	15	15	15	15	15	11
    15	15	14	15	15	15	11
    15	15	14	15	15	15	11
    15	15	14	15	15	13	10
    13	13	14	13	14	12	7
    13	10	13	11	13	8	7
    12	5	11	8	12	4	6
    7	4	10	6	10	5	7 ]; 

DATA{3}(2:end,2:end) = 1e-3 * DATA{3}(2:end,2:end); % mg wwt to g wwt

DATA{3}(2,2:end) = 1.8e-6; % add initial wet weight as data point for each treatment! (1.8 ug)
W{3}(1,1:end)    = 15; % also provide weights for that initial value

% cumulative reproduction (nr. eggs per individual) on each observation
% time (d), concentrations in mg/kg food
DATA{4} = [0.5	0	0.1	1 2 3 4 5
    0	0	0	0	0	0	0	0
    8	0	0	0	0	0	0	0
    10	0	0	0	0	0	0	0
    12	0	0	0	0	0	0	0
    15	0	0	0	0	0	0	0
    17	11.9	10.1	6.5	2.6	0	0	0
    19	29.8	20.6	18.4	21.3	14.3	2.875	0
    22	29.8	20.6	23.6	21.3	19	6	0
    24	55.8	61	46.2	53.8	25.7	29.57142857	1.714285714
    26	116.4	126.7	114.5	107	94.4	34.19642857	1.714285714
    29	119.4	129.6	128.5	117.7	107.1	51.82142857	1.714285714
    31	119.4	134.6	142.1	117.7	108.4333333	55.67857143	12.04761905
    34	235.3	241.3	236	208.3	177.5444444	90.39285714	12.04761905
    36	238.4	241.3	252.7	231	228.1	116.6785714	12.04761905
    38	241.4	243	252.7	231	230.4333333	116.6785714	12.04761905
    40	248.6	243	252.7	231	230.4333333	132.1071429	12.04761905
    43	312.4	306.3333333	301.5	276.8888889	292.9888889	152.3571429	12.04761905
    45	333.775	334.7777778	327.7	337.031746	329.2388889	188.6904762	12.04761905];

% weights: number of individuals from which reproduction was determined
W{4} = [10	10	10	10	10	10	10
    10	10	10	10	10	10	10
    10	10	10	10	10	10	10
    10	10	10	10	10	10	7
    10	10	10	10	10	8	8
    10	10	10	10	10	8	7
    10	10	10	10	10	8	8
    10	10	10	10	10	8	6
    10	10	10	10	10	7	7
    10	10	10	10	10	8	7
    10	10	10	10	10	8	7
    10	10	10	10	9	7	6
    10	10	10	10	9	7	6
    10	10	10	10	9	7	3
    10	10	10	10	9	6	2
    10	10	10	9	9	7	2
    10	9	10	9	9	8	2
    8	9	10	7	8	6	2];

% survivors on each observation time (d), concentrations in mg/kg food
DATA{6} = [-1	0	0.1	1 2 3 4 5
    0	30	30	30	30	30	30	30
    8	30	30	30	30	30	30	30
    10	30	30	30	30	30	30	30
    12	30	30	30	30	30	30	23
    15	30	30	30	30	30	24	19
    17	30	30	30	30	30	22	14
    19	30	30	30	30	30	22	14
    22	30	30	30	30	30	21	11
    24	30	30	30	30	30	20	11
    26	30	29	30	30	30	20	9
    29	30	28	30	30	30	18	8
    31	30	28	30	30	30	17	6
    34	30	28	30	30	30	17	6
    36	30	27	30	30	30	17	3
    38	30	27	30	30	29	17	2
    40	30	26	30	29	29	17	2
    43	30	25	29	28	29	17	2
    45	28	25	29	26	28	17	2];

% Define scenarios for the treatments (mg/kg food)
Cw = [0  0 0.1 1 2 3 4 5
      0 0 0 0.93 2 4.31 9.28 20
     140 0 0 0.93 2 4.31 9.28 20];

glo.scen_plot = 1; % make a plot for the scenarios
make_scen(2,Cw); % type 2 creates block pulses (fine for constant exposure)

% Create a table with nicer custom labels for the legends
Scenario = [0 0.1 1 2 3 4 5]'; % scenario identifiers that get a label
Label    = {'control';'solvent control';'0.93 mg/kg food';'2 mg/kg food';'4.31 mg/kg food';'9.28 mg/kg food';'20 mg/kg food'};
glo.LabelTable = table(Scenario,Label); % create a Matlab table for the labels

% Note: optionally, add some info to the MAT filename. The MAT filename
% will already include MoA and feedbacks, but if you want to try other
% things as well (changing opt, calibrating on the validation data, etc)
% it can be helpful to change the name to use this script but not 
% overwrite previous MAT files.

% glo.basenm  = [mfilename,'_maintxt']; % remember the filename for THIS file for the plots
save([glo.basenm,'_DATA'],'DATA','W') % save MAT file with data set 
% Saving the data set is handy to allow for a simple reconstruction of the
% calibrations, without needing to define the data again, in the same way.
% Note that if you use DATAx and Wx, you need to save them as well!

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios. Note that for stdDEB, the initial values cannot
% be randomly chosen. We need to simulate the individual up to the point
% where we start the model analysis (e.g., at birth). Therefore, most of
% the values in X0mat will not be used in the analysis.

X0mat = [glo.LabelTable.Scenario]'; % the scenarios (here identifiers) 
X0mat(2,:) = 0;   % initial reserve energy (J) NOT USED
X0mat(3,:) = 0;   % initial maturity level (J) NOT USED
X0mat(4,:) = 0;   % initial volumetric length (cm) NOT USED
X0mat(5,:) = 0;   % initial cumul. repro NOT USED
X0mat(6,:) = 0;   % initial values scaled damage
X0mat(7,:) = 1;   % initial values survival probability

% Put the position of the various states in globals, to make sure correct
% one is selected for extra things (e.g., to prevent shrinking in
% call_deri, and for using plot_tktd). This is placed here as it is used in
% the data file that is called below, to put the data in the right position
% of the DATA structure.
glo.locL = 3; % location of body size in the state variable list
glo.locR = 4; % location of cumulative reproduction in the state variable list
glo.locD = 5; % location of damage in the state variable list 
glo.locS = 6; % location of survival probability in the state variable list

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. Global
% parameters as part of the structure glo. Note that glo.dV is not used
% here.

% Global settings for the analysis; several of these values will be
% species-specific. These are not taken from AmP.
glo.yP      = 0.64;  % product of yVE and yEV (-)
glo.len     = 0;     % switch to fit physical length (0=wwt, 1=phys. length, 2=phys. length, no shrinking) (used in call_deri.m)
glo.Tref    = 20 + 273.15; % reference temperature to Kelvin
glo.T       = 20 + 273.15; % actual temperature to Kelvin
glo.Tbp     = 0; % time (d) that eggs spend in brood pouch (model output for repro will be shifted)
glo.FBV     = 0.007; % dry weight egg as fraction of structural body weight (-) (for losses with repro; approx. for Folsomia, see Jager, 2020)
glo.KRV     = 1;     % part. coeff. repro buffer and structure (kg/kg) (for losses with reproduction)
glo.Lwm_ref = 0.24;  % reference max *physical* length for scaling rate constants (cm) (observed value in AmP entry)
% Note: using a fixed reference length for scaling is helpful for comparing
% parameters between data sets, and absolutely needed when fitting on
% multiple data sets where animals reach a different maximum length.
% 
% Note: Crommentuijn et al 1997 specifies 20C for the experiments. FBV from
% SI of Jager (2020) is 0.008; for some reason, the value in the main text
% is set at 0.007 ...

zvd.Wd0 = [0.41e-6 0.02e-6]; % zero-variate data point for egg dry weight (g) at f=1, with normal s.d. (from thesis Hamda, page 198)

glo.E0_calc = [1 1]; % strategy for maternal effects rule
% first element, initial values from provided f in treatment (1) or always use f=1 (2)
% second element, for repro, egg costs from f (1), always use f=1 (2), or from actual reserve status of mother (3)

% =========================================================================
% Parameters for Folsomia candida (Natanael T. Hamda, Starrlight Augustine,
% Bas Kooijman. 2016. AmP Folsomia candida, version 2016/02/09.)
% =========================================================================
% global settings for conversions from AmP
glo.delM = 0.1766; % shape corrector (-)
glo.dV   = 0.17;   % dry weight density (g/cm3)
% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.TA   = [8000     0 1000  1e6 1]; % Arrhenius temperature (K)
par.spAm = [84.56    0    0  1e6 1]; % max. surface-specific assimilation rate (J/(cm2 d))
par.spM  = [845.6    1    0  1e6 1]; % volume-specific somatic maintenance costs (J/(cm3 d))
par.spT  = [0        0    0  1e6 1]; % surface-specific maintenance costs (J/(cm2 d))
par.kJ   = [0.002    0    0  1e6 1]; % maturity maintenance rate constant (1/d)
par.EG   = [4314     0    0  1e6 1]; % volume-specific costs for growth (J/cm3)
par.EHb  = [0.005193 1 1e-6  1e6 1]; % maturity level at birth (J)
par.EHj  = [0.05483  1    0  1e6 1]; % maturity level at metamorphosis (J)
par.EHp  = [0.1112   1 1e-6  1e6 1]; % maturity level at puberty (J)
par.v    = [0.003881 1    0  1e6 1]; % energy conductance (cm/d)
par.kap  = [0.2226   1 0.01 0.99 1]; % allocation fraction to soma (-)
par.kapR = [0.95     0 0.01 0.99 1]; % reproduction efficiency (-)
par.f    = [0.9      0    0    2 1]; % scaled food density (-)
par.hb   = [0.01     1 1e-3 0.07 0]; % background hazard rate (d-1)
par.a    = [1        0 0.1  10   0]; % coefficient for Weibull backgound hazard (-) (1 is off)
% par.Lw0  = [1.8e-6   0 0 1e6 1]; % starting at a *physical* length L0>Lb (when glo.len=0, use wet weight)

% NOTE: the AmP entry assumes f=0.9 for the data set of Crommentuijn et al.
% so we stick to that value here as well.
% 
% NOTE: Crommentuijn et al (1997) specify 18 ug for 1-day old animals, but
% this cannot be correct. Here, 1.8 ug (following Jager, 2020) is added to
% the data set.
% 
% NOTE: some parameters for the AmP entry need to be refitted to obtain a
% closer correspondence to the controls of this data set. This choice is
% discussed in the paper of Jager et al (2023).

% =========================================================================
%      TKTD Parameters for this data set
% =========================================================================
ind_tox = length(fieldnames(par))+1; % index where tox parameters start
% the parameters below this line are all treated as toxicity parameters!

% When using the parameter-space explorer, start values and ranges are not
% needed (filled later by startgrid_debtox); only the fit/fix mark is
% relevant here. But make sure that the value in the first column is within
% the bounds (and not zero for log-scale parameters).
par.kd   = [0.08   1 0.01  10 0]; % dominant rate constant (d-1)
par.zb   = [0.10   1 0    1e6 1]; % effect threshold energy budget ([C])
par.bb   = [72     1 1e-6 1e6 0]; % effect strength energy-budget effect (1/[C])
par.zs   = [0.24   1 0    1e6 1]; % effect threshold survival ([C])
par.bs   = [0.64   1 1e-6 1e6 0]; % effect strength survival (1/([C] d))

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

% specify the y-axis labels for each state variable
glo.ylab{1} = 'reserve (J)';
glo.ylab{2} = 'maturity (J)';
glo.ylab{3} = 'body wet weight (g)';
if isfield(glo,'Tbp') && glo.Tbp > 0
    glo.ylab{4} = ['cumul. repro (shift ',num2str(glo.Tbp),'d)'];
else
    glo.ylab{4} = 'cumul. reproduction';
end
glo.ylab{5} = ['damage (',char(181),'g/L)'];
glo.ylab{6} = 'survival fraction (-)';
% specify the x-axis label (same for all states)
% if isfield(par,'L0')
%     glo.xlab    = 'time since start (days)';
% else
    glo.xlab    = 'time (days)';
% end
glo.leglab1 = 'scen. '; % legend label before the 'scenario' number
glo.leglab2 = ''; % legend label after the 'scenario' number
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
% contains discontinuities. Don't use for continuous splines (type 1) as it
% will be much slower. For FOCUS scenarios (high time resolution), breaking
% up is not efficient and does not appear to be necessary.
% -------------------------------------------------------------------------

opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_plot.bw      = 1; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.repls   = 0; % set to 1 to plot replicates, 0 to plot mean responses
opt_plot.limax   = 1; % if set to 1, limit axes to the data set for each stage
opt_plot.notitle = 1; % set to 1 to suppress plotting titles on graphs with fits
opt_plot.statsup = [1 2 5]; % vector with states to suppress in plotting fits
opt_plot.y_zero  = 1;  % set to 1 to force y-axis to start at zero

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

MOA   = [0 0 0 0 1]; % hazard to the embryo
FEEDB = [0 0 0 0]; % no feedbacks
% FEEDB = [1 0 0 1]; % juveniles more sensitive than adults
% MOA   = [0 0 0 1 0]; % reproduction costs
% FEEDB = [0 0 0 0;1 1 1 1;0 0 1 0;1 1 1 0]; % four 'typical' feedbakcs
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

% Results from all-possible-feedbacks for different pMoAs, copied from screen
% 
% MoA  feedbacks   MLL     delta-AIC     prob.    best
% ====================================================================
% 00001  0000    1202.86      28.11       0.00 
% 00001  0001    1195.24      12.87       0.00 
% 00001  0010    1210.60      43.61       0.00 
% 00001  0011    1209.66      41.71       0.00 
% 00001  0100    1210.87      44.14       0.00 
% 00001  0101    1208.27      38.93       0.00 
% 00001  0110    1211.12      44.64       0.00 
% 00001  0111    1211.03      44.46       0.00 
% 00001  1000    1196.84      16.09       0.00 
% 00001  1001    1188.80       0.00       1.00      *
% 00001  1010    1219.40      61.21       0.00 
% 00001  1011    1214.17      50.74       0.00 
% 00001  1100    1205.56      33.52       0.00 
% 00001  1101    1192.00       6.40       0.04 
% 00001  1110    1212.43      47.26       0.00 
% 00001  1111    1212.37      47.14       0.00 
% ====================================================================
% 
% MoA  feedbacks   MLL     delta-AIC     prob.    best
% ====================================================================
% 00010  0000    1204.52      13.23       0.00 
% 00010  0001    1197.90       0.00       1.00      *
% 00010  0010    1212.51      29.22       0.00 
% 00010  0011    1212.32      28.83       0.00 
% 00010  0100    1219.44      43.07       0.00 
% 00010  0101    1217.18      38.56       0.00 
% 00010  0110    1219.72      43.63       0.00 
% 00010  0111    1219.65      43.48       0.00 
% 00010  1000    1199.07       2.33       0.31 
% 00010  1001    1198.09       0.36       0.83 
% 00010  1010    1217.42      39.03       0.00 
% 00010  1011    1214.40      33.00       0.00 
% 00010  1100    1207.57      19.32       0.00 
% 00010  1101    1198.59       1.38       0.50 
% 00010  1110    1212.75      29.70       0.00 
% 00010  1111    1212.53      29.25       0.00 
% ====================================================================

% MoA  feedbacks   MLL     delta-AIC     prob.    best
% ====================================================================
% 00110  0000    1227.97      18.60       0.00 
% 00110  0001    1218.67       0.00       1.00      *
% 00110  0010    1246.45      55.55       0.00 
% 00110  0011    1238.99      40.63       0.00 
% 00110  0100    1231.46      25.58       0.00 
% 00110  0101    1219.65       1.96       0.38 
% 00110  0110    1231.48      25.63       0.00 
% 00110  0111    1231.52      25.70       0.00 
% 00110  1000    1233.73      30.13       0.00 
% 00110  1001    1233.52      29.70       0.00 
% 00110  1010    1277.05     116.75       0.00 
% 00110  1011    1279.05     120.76       0.00 
% 00110  1100    1235.95      34.56       0.00 
% 00110  1101    1235.77      34.21       0.00 
% 00110  1110    1276.74     116.14       0.00 
% 00110  1111    1284.97     132.61       0.00 
% ====================================================================

% % ===== SIMULATING CONTROLS WITH AmP VALUES ===============================
% glo.stdDEB_start = {}; % make sure it is empty and not defined!
% fit_tox = [0 0 3]; % use both controls, don't fit
% automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot);
% 
% return

% Provide fitted parameters as starting values to aid optimisation ...
par.spM    = [    2030.9  1          0      1e+06 1]; 
par.EHb    = [  0.003147  1      1e-06      1e+06 1]; 
par.EHj    = [  0.075188  1          0      1e+06 1]; 
par.EHp    = [   0.29522  1      1e-06      1e+06 1]; 
par.v      = [ 0.0066135  1          0      1e+06 1]; 
par.kap    = [   0.44356  1       0.01       0.99 1]; 
par.hb     = [ 0.0026503  1      0.001       0.07 0]; 

% ===== FITTING CONTROLS ==================================================
% Simplex fitting works fine for control data
opt_optim.type = 1; % optimisation method: 1) default simplex, 4) parspace explorer

% opt_optim.type     = 4; % optimisation method 1) simplex, 4 parameter-space explorer
% opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots of parameter space to monitor progress
% opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS (2 for extra rough)

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
% par_out = automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot,opt_prof);
% script to run the calculations and plot, automatically
par = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par

% Next, we need to clear the zero-variate data points, so they are not used
% again for fitting the treatments.
glo.zvd = []; % this should be enough, since it makes the output zvd of call_deri empty as well

% NOTE: even when you don't want to fit any of the control parameters, it
% is good to run with fit_tox(1)=0. This will make sure that the start
% values are defined (that won't change in the analysis when fitting only
% the tox parameters), which will speed up things considerably.

% Fit hb in data set
fit_tox = [-1 1 3]; % use both controls
par_out = automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot); 
% par_out = automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot,opt_prof); 
% script to run the calculations and plot, automatically
par = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par

% NOTE: I start fitting the control parameters, rather than hb, as that
% should lead to a proper definition of the initial values in
% glo.stdDEB_start!

% return

% ===== FITTING TOX DATA ==================================================
% Use the parameter-space explorer for fitting the treatments. Note that,
% with fit_tox=1, the tox parameters in par are replaced (when fitted) with
% estimates based on the data set using startgrid_debtox (called in
% automatic_runs_debtox2019). This is still quite experimental, so you may
% need to restart with manually-adapted ranges! Furthermore, this is really
% slow ... (the parallel toolbox really helps here!).
% 
% Note: you can use the rough settings to find better ranges and restart
% with refined settings. With opt_optim.ps_profs = 0, the algorithm will
% provide new search ranges on screen that can be directly copied-pasted
% into this script. You may comment out the fitting of the controls above.
% Make sure that skip_sg is then set to 1 to avoid startgrid_debtox to
% overwrite par again.

opt_optim.type     = 4; % optimisation method 1) simplex, 4 parameter-space explorer
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
% % THIS BLOCK: SETTINGS TO REFINE FROM RANGES COPIED FROM SCREEN INTO THIS SCRIPT
% opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
% glo.stiff          = [0 3]; % use ode45 with very strict tolerances
% skip_sg            = 1; % set to 1 to skip startgrid completely (use ranges in <par> structure)
% -----------------------------------------------------------------------

opt_plot = []; % empty options prevents plotting from automatic_runs

fit_tox = [1 1 3]; % use both controls, and fit tox data
[par_out,best_MoaFb] = automatic_runs_debtox2019(fit_tox,par,ind_tox,skip_sg,MOA,FEEDB,opt_optim,opt_plot);
% script to run the calculations and plot, automatically
% Note: automatic_runs will return the BEST parameter set in par_out.
disp_settings_debtox2019 % display a bit more info on the settings on screen
print_par(par_out) % print the complete best parameter vector that can be copied-pasted
% This includes the control parameters (the parspace explorer itself will
% only plot the results for the fitted parameters).

%% Plot results with confidence intervals
% The following code can be used to make plots with confidence intervals.
% Options for confidence bounds on model curves can be set using opt_conf
% (see prelim_checks). The plot_tktd function makes multiplots for the
% effects data, which are more readable when plotting with various
% intervals.

opt_conf.type     = 3; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set  = 2; % use limited set of n_lim points (1) or outer hull (2, not for Bayes) to create CIs
opt_tktd.repls    = 0; % plot individual replicates (1) or means (0)
opt_tktd.transf   = 1; % set to 1 to calculate means and SEs including transformations
opt_tktd.lim_data = 1; % set to 1 to limit axes to data

opt_tktd.max_exp  = 0; % set to 1 to maximise exposure/damage plots on exposure rather than damage
opt_tktd.obspred  = 0; % plot predicted-observed plots (1) or not (0)
opt_tktd.sppe     = 0; % set to 1 to calculate SPPEs (relative error at end of test)
% Note: this last option will produce something like the SPPE
% (opt_tktd.obspred must be set to 1 to calculate SPPE and other metrics).

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

% Code below plots the best or a selected MoA/feedback configuration. It
% plots the fit with or without CIs (opt_conf.type=0 or empty opt_conf in
% the call to plot_tktd). The automatic_runs returns the parameter set and
% indices for the best settings.
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
% leave par_out (first input) empty to read it from saved mat file
