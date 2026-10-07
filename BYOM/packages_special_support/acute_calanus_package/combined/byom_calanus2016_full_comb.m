%% BYOM, byom_calanus2016_full_comb.m GUTS
%
% * Author: Tjalling Jager 
% * Date: December 2021
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* fitting survival data and body residues simultaneously with
% the <http://www.debtox.info/about_guts.html GUTS> special cases based on
% the internal concentration: SD, IT and mixed. The internal concentration
% in structure is directly linked to the effect, without damage stage.
%
% *This script:* C2-naphthalene in _Calanus finmarchicus_. Fitting both the
% survival data and the two-compartment TK model on the body residues.
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

% observed number of survivors, time in days, conc in uM
DATA{1} = [-1	0	1.03	1.97	4.83	10.3	29.3
    0	56	56	28	28	28	28
    1	56	56	28	28	25	0
    2	56	56	28	28	10	0
    3	56	56	28	26	1	0
    4	56	55	26	21	0	0
    5	56	55	26	18	0	0
    6	56	55	25	6	0	0];
        
% internal concentrations in umol/kg, time in days
DATA{4} = [0.5	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2
    0	2.634920501	1.128226067	2.067572648	1.239467263	2.133021081	3.968671089	2.232029925	1.71563171
    0.25	346.903685	326.4266006	330.6922294	NaN	NaN	NaN	NaN	NaN
    0.5	616.6576112	592.2330082	677.762585	NaN	NaN	NaN	NaN	NaN
    0.75	762.6419853	750.075518	854.8605579	840.4228865	NaN	NaN	NaN	NaN
    1	1187.964104	1224.36273	1075.951767	1121.508449	NaN	NaN	NaN	NaN
    2	2027.432656	2037.744655	1638.995127	NaN	NaN	NaN	NaN	NaN
    3	2524.564775	2471.508177	2675.329348	NaN	NaN	NaN	NaN	NaN
    4	2685.458367	2782.061779	2345.599143	3347.038305	NaN	NaN	NaN	NaN
    4.25	2194.05948	2119.389788	2501.947601	3182.207062	NaN	NaN	NaN	NaN
    4.5	2309.995772	2353.777302	2290.190177	2226.329017	NaN	NaN	NaN	NaN
    4.75	2462.489785	2316.686876	2141.947838	1862.544099	NaN	NaN	NaN	NaN
    5	1772.724311	2170.229973	2019.317227	1786.46631	NaN	NaN	NaN	NaN
    6	2076.314829	1390.663864	742.5848544	NaN	NaN	NaN	NaN	NaN
    7	1800.811532	1862.208204	1223.001179	1050.560266	NaN	NaN	NaN	NaN
    8	1920.780515	1490.778938	1228.435888	1602.794556	NaN	NaN	NaN	NaN];

W{4} = 21 * ones(size(DATA{4})-1); % each point is a pooled sample of 21 animals

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = [DATA{1}(1,2:end) 1.2]; % scenarios (concentrations)
X0mat(2,:) = 1; % initial survival probability
X0mat(3,:) = 0; % initial internal concentration in structure
X0mat(4,:) = 0; % initial internal concentration in lipids
X0mat(5,:) = 0; % initial internal concentration in total animal

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% global parameter
glo.FLS  = 0.2; % VL/VS
glo.Tdep = 4; % depuration time for fitting
glo.sel  = 3; % select death mechanism: 1) SD 2) IT 3) mixed
% globals used for GUTS purposes
glo.locD = 2; % location of the dose metric in the state variable list
glo.locS = 1; % location of survival probability in the state variable list

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.keS    = [   0.57192  1       0.01         10 1]; % elimination rate constant for structure (d-1)
par.kuS    = [    1460.8  1       0.01      1e+06 1]; % uptake rate constant for structure (L/L/d)
par.keL    = [  0.017739  1       0.01         10 1]; % elimination rate constant for lipids (d-1)
par.kuL    = [   0.88245  1       0.01      1e+06 1]; % uptake rate constant for lipids (L/L/d)
par.mS     = [    7803.2  1          0      1e+06 1]; % median threshold for survival (umol/L)
par.hb     = [  0.004495  1          0      1e+06 1]; % background hazard rate (1/d)
par.bS     = [0.00025833  1          0       1000 1]; % killing rate (L/umol/d)
par.Fs     = [    1.0017  1          1        100 1]; % fraction spread of NEC distribution (-) (IT and mixed)
% Note: these starting values result from better optima found by profiling
% the likelihood function.

switch glo.sel % make sure that right parameters are fitted
    case 1 % SD
        par.Fs(2) = 0; % never fit the NEC spread!
    case 2 % IT
        par.bS(2) = 0; % never fit the killing rate!
end

%% Time vector and labels for plots
% Specify what to calculate and what to plot.

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival probability';
glo.ylab{2} = ['internal conc., structure (',char(181),'mol/L)'];
glo.ylab{3} = ['internal conc., lipid sac (',char(181),'mol/L)'];
glo.ylab{4} = ['internal conc., total (',char(181),'mol/kg)'];

% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = ''; % legend label before the 'scenario' number
glo.leglab2 = [char(181),'M']; % legend label after the 'scenario' number

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

glo.stiff    = [0 1]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Seems that standard solver performs good enough, with tight tolerances.
% Second argument is for normally tight (1), tighter (2), or very tight (3)
% tolerances.
glo.break_time = 0; % break time vector up for ODE solver (1) or don't (0)

opt_optim.fit  = 1; % fit the parameters (1), or don't (0)
opt_optim.it   = 1; % show iterations of the simplex optimisation (1, default) or not (0)
opt_plot.sho   = 0; % set to 1 to show all scenarios, 0 to only plot model for scenarios with data
opt_plot.bw    = 1; % plot in black and white
opt_plot.cn    = 0; % if set to 1, connect model line to points (only for bw=1)
opt_plot.limax = 1; % if set to 1, limit axes to the data set for each stage

opt_plot.statsup = [2 3]; % vector with states to suppress in plotting fits
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

%% Profiling the likelihood
% By profiling you make robust confidence intervals for one or more of your
% parameters. 
 
opt_prof.detail   = 2; % detailed (1) or a coarse (2) calculation
opt_prof.subopt   = 30; % number of sub-optimisations to perform to increase robustness
opt_prof.brkprof  = 2; % when a better optimum is located, stop (1) or automatically refit (2)

% % UNCOMMENT FOLLOWING LINE(S) TO CALCULATE
% calc_proflik(par_out,'all',opt_prof,opt_optim);  % calculate a profile
% % Enter single parameter names, a cell array of names, or 'all' to profile
% % all fitted parameters.

% CALCULATION OF INTERVALS COMMENTED OUT AS THIS IS VERY TIME CONSUMING
