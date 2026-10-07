%% BYOM, byom_stddeb_test.m
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
% *This script:* Testing the stdDEB model with the AmP entry for _Gammarus
% pulex_. This is just a small part of the data that are used in the AmP
% entry.
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

% =========================================================================
%      Growth curve for Gammarus pulex
% =========================================================================
% State is for body length (mm -> cm)
DATA{3} = [1 1
    7.82806921963  2.24083499405;
    25.1401453784  3.04268161551;
    47.8644855123  4.11177221179;
    67.2052645916  4.00716856156;
    83.9080130181  8.27561685472;
    110.847467347  8.70515706318;
    171.070552983  8.76491419232;
    220.784527377  10.6368836829;
    473.56234498   11.3038466501];

DATA{3}(2:end,2:end) = DATA{3}(2:end,2:end) / 10; % from mm to cm

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios. Note that for stdDEB, the initial values cannot
% be randomly chosen. We need to simulate the individual up to the point
% where we start the model analysis (e.g., at birth).

X0mat      = 1;   % the scenarios
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
glo.yP    = 0.64;  % product of yVE and yEV (-)
glo.len   = 1;     % switch to fit physical length (0=wwt, 1=phys. length, 2=phys. length, no shrinking) (used in call_deri.m)
glo.Tref  = 20 + 273.15; % translate reference temperature to temperature in Kelvin
glo.T     = 13 + 273.15; % translate actual temperature to temperature in Kelvin
% glo.Tbp   = 30; % testing effect of brood pouch delay

glo.E0_calc = [1 1]; % strategy for maternal effects rule
% first element, initial values from provided f in treatment (1) or always use f=1 (2)
% second element, for repro, egg costs from f (1), always use f=1 (2), or from actual reserve status of mother (3)
% 
% It is possible to make glo.E0_calc with different rows, and to use
% separate parameters per study, as in the TKTD part of this package. Make
% sure to use identifiers 0-99 for set 1, 100-199 for set 2, etc. And,
% create exposure scenarios with make_scen (just zero, for example).

% =========================================================================
%   Parameters for Gammarus pulex (Elke Zimmer, Bas Kooijman, 
%   Annika Mangold-Doering. 2021. AmP Gammarus pulex, version 2021/07/30.)
% =========================================================================
% global settings for conversions from AmP
glo.delM  = 0.23453; % shape corrector (-)
glo.dV    = 0.17;  % dry weight density (g/cm3)
% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.TA   = [10000    0 1000  1e6]; % Arrhenius temperature (K)
par.spAm = [43.0629  0    0  1e6]; % max. surface-specific assimilation rate (J/(cm2 d))
par.spM  = [338.417  0    0  1e6]; % volume-specific somatic maintenance costs (J/(cm3 d))
par.spT  = [0        0    0  1e6]; % surface-specific maintenance costs (J/(cm2 d))
par.kJ   = [0.002    0    0  1e6]; % maturity maintenance rate constant (1/d)
par.EG   = [4446.22  0    0  1e6]; % volume-specific costs for growth (J/cm3)
par.EHb  = [0.06315  0 1e-6  1e6]; % maturity level at birth (J)
par.EHj  = [2.108    0    0  1e6]; % maturity level at metamorphosis (J)
par.EHp  = [2.109    0 1e-6  1e6]; % maturity level at puberty (J)
par.v    = [0.018066 0    0  1e6]; % energy conductance (cm/d)
par.kap  = [0.88685  0 0.01 0.99]; % allocation fraction to soma (-)
par.kapR = [0.95     0 0.01 0.99]; % reproduction efficiency (-)
par.f    = [1        1    0    2]; % scaled food density (-)

% % For starting at a length L0>Lb, uncomment line below and set value
% par.Lw0 = [0.4        0 0 1e6]; % starting at a *physical* length L0>Lb (when glo.len=0, use wet weight)

% NOTE: only f is refitted to obtain a closer correspondence to the data
% set.

% % For testing the starvation response, define glo.starv, and also
% % uncomment the starvation section in derivatives.
% glo.starv = [100 200 0.3]; % FOR TESTING: temporary startvation between starv(1) and starv(2) at f=starv(3)
% glo.starv = [inf inf 0.2]; % FOR TESTING: temporary startvation between starv(1) and starv(2) at f=starv(3)

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

% glo.t   = linspace(0,100,100); % time vector for the model curves in days

% specify the y-axis labels for each state variable
glo.ylab{1} = 'reserve (J)';
glo.ylab{2} = 'maturity (J)';
glo.ylab{3} = 'body length (cm)';
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

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them
