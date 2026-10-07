%% BYOM, byom_stddeb_Folsomia.m
%
% * Author: Tjalling Jager
% * Date: February 2023
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
% and body length. In this folder, no toxicant stress is included in the
% model. The model equations are those from the 'Tromso coffee mug' ;-).
% Only change is that the ODE for structure is phrased in volumetric length
% rather than volume. The 2023 publication of Jager et al in Ecological
% Modelling contains the full details:
% <https://doi.org/10.1016/j.ecolmodel.2022.110187>.
%
% *This script:* the colembolan _Folsomia candida_ exposed to chlorpyrifos
% in food. Only the controls are used here. Data set from Crommentuijn et
% al, analysed in Jager et al (2007)
% <http://dx.doi.org/10.1016/j.envpol.2006.04.028> and used as case study
% in Jager (2020) <https://doi.org/10.1016/j.ecolmodel.2019.108904>.
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
DATA{3} = [0	0	0	0.93	2	4.31	9.28	20
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

% % Remove last time points!
% ind_rem = find(DATA{3}(:,1) > 51); % indices for points beyond t=51
% DATA{3}(ind_rem,:) = []; % remove the last points from the data matrix
% W{3}(ind_rem-1,:) = []; % remove the last points from the weights matrix

% cumulative reproduction (nr. eggs per individual) on each observation
% time (d), concentrations in mg/kg food
DATA{4} = [0.5	0	0	0.93	2	4.31	9.28	20
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

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios. Note that for stdDEB, the initial values cannot
% be randomly chosen. We need to simulate the individual up to the point
% where we start the model analysis (e.g., at birth).

X0mat      = [0];   % the scenarios
X0mat(2,:) = 0;   % initial reserve energy (J) NOT USED
X0mat(3,:) = 0;   % initial maturity level (J) NOT USED
X0mat(4,:) = 0;   % initial volumetric length (cm) NOT USED
X0mat(5,:) = 0;   % initial cumul. repro NOT USED

% Put the position of the various states in globals, to make sure correct
% one is selected for extra things (e.g., to prevent shrinking in
% call_deri, and for using plot_tktd).
glo.locL = 3; % location of body size in the state variable list
glo.locR = 4; % location of cumulative reproduction in the state variable list

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. Global
% parameters as part of the structure glo. Note that glo.dV is not used
% here.

% Global settings for the analysis; several of these values will be
% species-specific. These are not taken from AmP.
glo.yP     = 0.64;  % product of yVE and yEV (-)
glo.len    = 0;     % switch to fit physical length (0=wwt, 1=phys. length, 2=phys. length, no shrinking) (used in call_deri.m)
glo.Tref   = 20 + 273.15; % reference temperature to Kelvin
glo.T      = 20 + 273.15; % actual temperature to Kelvin
glo.Tbp    = 0;     % brood pouch delay

zvd.Wd0 = [0.41e-6 0.02e-6]; % zero-variate data point for egg dry weight (g) at f=1, with normal s.d.
% zvd.Vwm = [300e-6 20e-6]; % zero-variate data point for maximum wet weight (g) at f=1, with normal s.d.

glo.E0_calc = [1 1]; % strategy for maternal effects rule (just 1 row allowed!)
% first element, initial values from provided f in treatment (1) or always use f=1 (2)
% second element, for repro, egg costs from f (1), always use f=1 (2), or from actual reserve status of mother (3)
% 
% It is possible to make glo.E0_calc with different rows, and to use
% separate parameters per study, as in the TKTD part of this package. Make
% sure to use identifiers 0-99 for set 1, 100-199 for set 2, etc. And,
% create exposure scenarios with make_scen (just zero, for example).

% =========================================================================
% Parameters for Folsomia candida (Natanael T. Hamda, Starrlight Augustine,
% Bas Kooijman. 2016. AmP Folsomia candida, version 2016/02/09.)
% =========================================================================
% global settings for conversions from AmP
glo.delM   = 0.1766; % shape corrector (-)
glo.dV     = 0.17;  % dry weight density (g/cm3)
% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.TA   = [8000     0 1000  1e6]; % Arrhenius temperature (K)
par.spAm = [84.56    0    0  1e6]; % max. surface-specific assimilation rate (J/(cm2 d))
par.spM  = [845.6    1    0  1e6]; % volume-specific somatic maintenance costs (J/(cm3 d))
par.spT  = [0        0    0  1e6]; % surface-specific maintenance costs (J/(cm2 d))
par.kJ   = [0.002    0    0  1e6]; % maturity maintenance rate constant (1/d)
par.EG   = [4314     0    0  1e6]; % volume-specific costs for growth (J/cm3)
par.EHb  = [0.005193 1 1e-6  1e6]; % maturity level at birth (J)
par.EHj  = [0.05483  1    0  1e6]; % maturity level at metamorphosis (J)
par.EHp  = [0.1112   1 1e-6  1e6]; % maturity level at puberty (J)
par.v    = [0.003881 1    0  1e6]; % energy conductance (cm/d)
par.kap  = [0.2226   1 0.01 0.99]; % allocation fraction to soma (-)
par.kapR = [0.95     0 0.01 0.99]; % reproduction efficiency (-)
par.f    = [0.9      0    0    2];   % scaled food density (-)
% par.Lw0  = [1.8e-6   0 0 1e6]; % starting at a *physical* length L0>Lb (when glo.len=0, use wet weight)

% NOTE: some parameters for the AmP entry need to be refitted to obtain a
% closer correspondence to the controls of this data set.

% % For testing the starvation response, define glo.starv, and also
% % uncomment the starvation section in derivatives.
% glo.starv = [100 200 0.3]; % FOR TESTING: temporary starvation between starv(1) and starv(2) at f=starv(3)
% glo.starv = [inf inf 0.2]; % FOR TESTING: temporary starvation between starv(1) and starv(2) at f=starv(3)

% Provide fitted parameters as starting values to aid optimisation ...
% since the optimisation starting with the AmP values gets stuck in a local
% minimum.
par.spM    = [    2030.9  1          0      1e+06 1]; 
par.EHb    = [  0.003147  1      1e-06      1e+06 1]; 
par.EHj    = [  0.075188  1          0      1e+06 1]; 
par.EHp    = [   0.29522  1      1e-06      1e+06 1]; 
par.v      = [ 0.0066135  1          0      1e+06 1]; 
par.kap    = [   0.44356  1       0.01       0.99 1]; 
par.hb     = [ 0.0026503  1      0.001       0.07 0]; 

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

% glo.t   = linspace(0,80,100); % time vector for the model curves in days

% specify the y-axis labels for each state variable
glo.ylab{1} = 'reserve (J)';
glo.ylab{2} = 'maturity (J)';
glo.ylab{3} = 'body wet weight (g)';
glo.ylab{4} = 'cumul. repro (eggs)';
% specify the x-axis label (same for all states)
if isfield(par,'L0')
    glo.xlab    = 'time since start (days)';
else
    glo.xlab    = 'time since birth (days)';
end
glo.leglab1 = 'scen. '; % legend label before the 'scenario' number
glo.leglab2 = ''; % legend label after the 'scenario' number

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

glo.stiff = [0 3]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Second argument is for normally tight (1), tighter (2), or very tight (3)
% tolerances. Use 1 for quick analyses, but check with 3 to see if there is
% a difference!

opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_optim.it     = 1; % show iterations of the optimisation (1, default) or not (0)
opt_plot.bw      = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 0; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.limax   = 1; % if set to 1, limit axes to the data set for each stage
% opt_plot.statsup = [1 2]; % vector with states to suppress in plotting fits

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them
