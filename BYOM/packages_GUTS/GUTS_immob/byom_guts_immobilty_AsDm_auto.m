%% BYOM, byom_guts_immobility_AsDm_auto.m
%
% * Author     : Tjalling Jager
% * Date       : April 2024
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> is modified to the particular problem as well. The files in
% the engine directory are needed for fitting and plotting. Results are
% shown on screen but also saved to a log file (|results.out|).
%
% *The model:* fitting survival data with an additional state of
% immobility. The files in this folder are based on the full GUTS
% (toxicokinetics and damage dynamics as separate compartments). However,
% there are switches to make either TK or damage fast (or both). Death and
% immobility can both be determined by damage, or immobility can be linked
% to the internal concentration (which is changing faster than damage).
% Note that internal concentration can double as a damage compartment by
% fixing Kiw=1. A publication presenting the model is underway. Using
% GUTS-immobility requires BYOM version 6.9 or newer!
%
% *This script:* demonstrating fitting with data set for _Asellus aquaticus_
% and deltamethrin. The data set is from Bayer study M-643156-01-1, and was
% published in the paper by Bauer et al. (2024) in Environ. Toxicol. Chem.:
% <https://doi.org/10.1002/etc.5761>. Three entries corrected 8/4/2024.
% Here, automatic_runs is used to sequentially run through various model
% configurations.
% 
%  Copyright (c) 2012-2024, Tjalling Jager.
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

pathdefine(1) % set path to the BYOM/engine directory (use 1 for parallel-version)
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 0; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. First column are the exposure times, first
% row are the concentrations or scenario numbers. The number in the top
% left of the matrix indicates how to calculate the likelihood:
%
% * -3 for unconditional multinomial likelihood (for immobility data; at this moment, SSQ is used)
% * -1 for multinomial likelihood (for survival data)
% * 0  for log-transform the data, then normal likelihood
% * 0.5 for square-root transform the data, then normal likelihood
% * 1  for no transformation of the data, then normal likelihood
  
% NOTE: use DATA{3} for active, DATA{4} for immobile, DATA{5} for dead
% animals, and DATA{6} for dead OR immobile animals (when no distinction is
% made). Each individual needs to be in ONE of these states. Use NaN where
% no observation is made. The data at t=0 should be the TOTAL number of
% individuals at the start of the experiment. This should be done for EACH
% state!

% internal concentration, has no observations
DATA{1} = 0;  
% scaled damage, can have no observations
DATA{2} = 0;   

% use DATA{3} for number of certainly active animals
DATA{3} = [-3	0 0 0 0 0 0	1 2 3
0.00	50	50	50	50	50	50	50	50	50
0.04	49	49	50	50	50	50	20	26	33
0.17	50	48	50	50	50	50	2	1	4
0.33	50	47	47	50	50	50	1	1	3
1.00	49	44	46	48	50	49	1	2	2
3.00	48	44	45	48	50	49	37	33	36
4.00	48	44	NaN	NaN	NaN	NaN	38	NaN	NaN
4.04	48	44	NaN	NaN	NaN	NaN	5	NaN	NaN
4.17	48	44	NaN	NaN	NaN	NaN	0	NaN	NaN
4.33	48	44	NaN	NaN	NaN	NaN	0	NaN	NaN
5.00	46	44	45	47	50	49	0	36	36
6.00	NaN	NaN	45	47	NaN	NaN	NaN	37	NaN
6.04	NaN	NaN	45	47	NaN	NaN	NaN	13	NaN
6.17	NaN	NaN	45	47	NaN	NaN	NaN	1	NaN
6.33	NaN	NaN	45	47	NaN	NaN	NaN	0	NaN
7.00	46	44	45	47	50	49	16	0	36
8.00	46	44	NaN	NaN	NaN	NaN	23	NaN	NaN
8.04	45	44	NaN	NaN	NaN	NaN	7	NaN	NaN
8.17	46	44	NaN	NaN	NaN	NaN	0	NaN	NaN
8.33	46	43	NaN	NaN	NaN	NaN	0	NaN	NaN
9.00	46	43	44	47	50	47	1	12	36
11.00	45	43	44	47	49	46	5	13	35
11.04	NaN	NaN	NaN	NaN	49	46	NaN	NaN	10
11.17	NaN	NaN	NaN	NaN	49	46	NaN	NaN	0
11.33	NaN	NaN	NaN	NaN	49	45	NaN	NaN	0
12.00	45	43-1	44	46	49	45	3	14	1 % correction made by Gaiac 8/4/2024
12.04	NaN	NaN	44	46	NaN	NaN	NaN	0	NaN
12.17	NaN	NaN	44	46	NaN	NaN	NaN	0	NaN
12.33	NaN	NaN	44	46	NaN	NaN	NaN	0	NaN
13.00	NaN	NaN	44	46	NaN	NaN	NaN	0	NaN
14.00	NaN	NaN	NaN	NaN	48	44	NaN	NaN	17
15.00	NaN	NaN	44	46	NaN	NaN	NaN	2	NaN
16.00	NaN	NaN	NaN	NaN	48	44	NaN	NaN	18
17.00	NaN	NaN	44	44	NaN	NaN	NaN	2	NaN
18.00	NaN	NaN	44	44	47	44	NaN	2	19
20.00	NaN	NaN	NaN	NaN	47	43	NaN	NaN	18
22.00	NaN	NaN	NaN	NaN	47	43	NaN	NaN	20];

% use DATA{4} for number for certainly immobile animals
DATA{4} = [-3	0 0 0 0 0 0	1 2 3
    0.00	50	50	50	50	50	50	50	50	50
0.04	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
0.17	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
0.33	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
1.00	0	1	0	0	0	0	47	47	47
3.00	0	0	0	0	0	0	8	15	8
4.00	0	0	NaN	NaN	NaN	NaN	7	NaN	NaN
4.04	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
4.17	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
4.33	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
5.00	0	0	0	0	0	0	41	7	2
6.00	NaN	NaN	0	0	NaN	NaN	NaN	5	NaN
6.04	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
6.17	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
6.33	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
7.00	0	0	0	0	0	0	24	42	1
8.00	0	0	NaN	NaN	NaN	NaN	15	NaN	NaN
8.04	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
8.17	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
8.33	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
9.00	0	0	0	0	0	0	37	28	0
11.00	0	0	0	0	0	0	30	24	0
11.04	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
11.17	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
11.33	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
12.00	0	0	0	0	0	0	23	21	32
12.04	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
12.17	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
12.33	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
13.00	NaN	NaN	0	0	NaN	NaN	NaN	30	NaN
14.00	NaN	NaN	NaN	NaN	0	0	NaN	NaN	13
15.00	NaN	NaN	0	0	NaN	NaN	NaN	26	NaN
16.00	NaN	NaN	NaN	NaN	0	0	NaN	NaN	8
17.00	NaN	NaN	0	0	NaN	NaN	NaN	19	NaN
18.00	NaN	NaN	0	0	0	0	NaN	17	7
20.00	NaN	NaN	NaN	NaN	0	0	NaN	NaN	5
22.00	NaN	NaN	NaN	NaN	0	0	NaN	NaN	2];
  
% use DATA{5} for number for certainly dead animals
DATA{5} = [-3	0 0 0 0 0 0	1 2 3
    0.00	50	50	50	50	50	50	50	50	50
0.04	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
0.17	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
0.33	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
1.00	1	5	4	2	0	1	2	1	1
3.00	2	6	5	2	0	1	5	2	6
4.00	2	6	NaN	NaN	NaN	NaN	5	NaN	NaN
4.04	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
4.17	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
4.33	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
5.00	4	6	5	3	0	1	9	7	12
6.00	NaN	NaN	5	3	NaN	NaN	NaN	8	NaN
6.04	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
6.17	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
6.33	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
7.00	4	6	5	3	0	1	10	8	13
8.00	4	6	NaN	NaN	NaN	NaN	12	NaN	NaN
8.04	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
8.17	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
8.33	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
9.00	4	7	6	3	0	3	12	10	14
11.00	5	7	6	3	1	4	15	12+1	15 % correction made by Gaiac 8/4/2024
11.04	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
11.17	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
11.33	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
12.00	5	5+3	6	4	1	5	24	15	17 % correction made by Gaiac 8/4/2024
12.04	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
12.17	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
12.33	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
13.00	NaN	NaN	6	4	NaN	NaN	NaN	20	NaN
14.00	NaN	NaN	NaN	NaN	2	6	NaN	NaN	20
15.00	NaN	NaN	6	4	NaN	NaN	NaN	22	NaN
16.00	NaN	NaN	NaN	NaN	2	6	NaN	NaN	24
17.00	NaN	NaN	6	6	NaN	NaN	NaN	29	NaN
18.00	NaN	NaN	6	6	3	6	NaN	31	24
20.00	NaN	NaN	NaN	NaN	3	7	NaN	NaN	27
22.00	NaN	NaN	NaN	NaN	3	7	NaN	NaN	28];

% use DATA{6} for number for immobile OR dead animals (when no distinction is made)
DATA{6} = [-3	0 0 0 0 0 0	1 2 3
0.00	50	50	50	50	50	50	50	50	50
0.04	1	1	0	0	0	0	30	24	17
0.17	0	2	0	0	0	0	48	49	46
0.33	0	3	3	0	0	0	49	49	47
1.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
3.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
4.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
4.04	2	6	NaN	NaN	NaN	NaN	45	NaN	NaN
4.17	2	6	NaN	NaN	NaN	NaN	50	NaN	NaN
4.33	2	6	NaN	NaN	NaN	NaN	50	NaN	NaN
5.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
6.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
6.04	NaN	NaN	5	3	NaN	NaN	NaN	37	NaN
6.17	NaN	NaN	5	3	NaN	NaN	NaN	49	NaN
6.33	NaN	NaN	5	3	NaN	NaN	NaN	50	NaN
7.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
8.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
8.04	5	6	NaN	NaN	NaN	NaN	43	NaN	NaN
8.17	4	6	NaN	NaN	NaN	NaN	50	NaN	NaN
8.33	4	7	NaN	NaN	NaN	NaN	50	NaN	NaN
9.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
11.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
11.04	NaN	NaN	NaN	NaN	1	4	NaN	NaN	40
11.17	NaN	NaN	NaN	NaN	1	4	NaN	NaN	50
11.33	NaN	NaN	NaN	NaN	1	5	NaN	NaN	50
12.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
12.04	NaN	NaN	6	4	NaN	NaN	NaN	50	NaN
12.17	NaN	NaN	6	4	NaN	NaN	NaN	50	NaN
12.33	NaN	NaN	6	4	NaN	NaN	NaN	50	NaN
13.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
14.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
15.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
16.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
17.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
18.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
20.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
22.00	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN];

check_data_guts_immob % quick check for obvious errors in the data set

% Define the exposure scenarios for each treatment. Note: in this case, it
% is easiest to define the scenarios as 3 separate matrices. The recent
% version of make_scen allows input of a cell array, and to use double time
% points for instant changes in exposure (i.e., pulses).

Cw{1} = [1.0000    1.0000
         0   16.6600
    1.0000   16.9200
    1.0000         0
    4.0000         0
    4.0000   13.7500
    5.0000   10.2700
    5.0000         0
    8.0000         0
    8.0000   15.7800
    9.0000   13.7800
    9.0000         0
   15.0000         0];

Cw{2} = [1.0000    2.0000
         0   15.0600
    1.0000   16.7800
    1.0000         0
    6.0000         0
    6.0000   12.1330
    7.0000   12.5000
    7.0000         0
   12.0000         0
   12.0000   17.8000
   13.0000   13.2330
   13.0000         0
   15.0000         0];

Cw{3} = [1.0000    3.0000
         0   15.6400
    1.0000   13.1060
    1.0000         0
   11.0000         0
   11.0000   16.2600
   12.0000   11.6640
   12.0000         0
   15.0000         0];
        
make_scen(4,Cw); % create the exposure scenarios with linear interpolation

% Make nice labels for all scenarios
Scenario = [0 1 2 3]';
Label    = {'control';'3x 3d 20ng/L';'3x 5d 20ng/L';'2x 10d 20ng/L'};

glo.LabelTable = table(Scenario,Label); % create a Matlab table for the labels

% NOTE: saving the DATA array here is possible, but will not automatically
% work when using the showcal file. The automatic_runs below will save each
% calibration in a MAT file with a different name. So, when running showcal
% you might be asked to manually select the correct data file.
save([glo.basenm,'_DATA'],'DATA','W') % save MAT file with data set 
% Saving the data set is handy to allow for a simple reconstruction of the
% calibrations, without needing to define the data again, in the same way.
% Note that if you use DATAx and Wx, you need to save them as well!

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = Scenario'; % scenarios to run (taken from the Scenario label table)
X0mat(2,:) = 0;       % initial internal concentration
X0mat(3,:) = 0;       % initial scaled damage
X0mat(4,:) = 1;       % initial active probability
X0mat(5,:) = 0;       % initial immobile probability
X0mat(6,:) = 0;       % initial dead probability
X0mat(7,:) = 0;       % initial immobile+dead probability

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% -------------------------------------------------------------------------
% Configurations for the model structure
% Death mechanism (glo.sel) and configuration (glo.damconfig) are part of
% the automated runs: in the block "Calculations and plotting" you can set
% SEL and CONFIG to run through multiple options and compare the results.
glo.fastslow  = [0 0]; % whether to use infinitely fast steady state (set to 1) for TK and damage, resp.
glo.onehit    = 1; % set to 1 (default) to use independent mechanisms for immobility and death; set to 2 to use 2-hit (under SD only)
% Note: slow kinetics is not (yet) implemented in this package
% Note: the default settings for the above options are [0 0], 1, 1
% -------------------------------------------------------------------------

% global parameters for GUTS purposes (don't change these)
glo.locC   = 1; % location of body residues in the state variable list
glo.locD   = 2; % location of scaled damage in the state variable list
glo.loc_a  = 3; % location of active probability in the state variable list
glo.loc_i  = 4; % location of immobile probability in the state variable list
glo.loc_d  = 5; % location of dead probability in the state variable list
glo.loc_id = 6; % location of immobile+dead probability in the state variable list

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.ke  = [130    1   1e-3 1e3 0];  % elimination rate constant, d-1
par.kr  = [1.46   1   1e-3 100 0];  % damage repair rate constant, d-1
par.Kiw = [1      0      0 1e6 1];  % bioconcentration factor (L/kg)
par.hb  = [0.001  1  0.001   1 1];  % background hazard rate (1/d)
% NOTE: when there is no mortality in the control, hb can be fixed to a
% low, reasonable, value. This may be better than setting it to zero. Kiw
% cannot be fitted, unless there are data on body residues as well. So,
% keep it to 1 to make internal concentration 'scaled'.

% par.mi  = [0.01  1 1e-3 1e6 1];  % joint threshold for immobility and death (ng/L)
par.mii = [0.001  1 1e-3  1e6 1];  % threshold for immobility (ng/L)
par.mid = [0.001  1 1e-3  1e6 1];  % threshold for death (ng/L)
% NOTE: this code allows to use EITHER 1 threshold for immobility and
% death, or to use separate ones. Use either mi OR mii and mid. When using
% separate states to drive immobility and death (damconfig 1 or 3), it
% makes most sense to have two thresholds as they are scaled in a different
% way. For IT, it will always be needed to have different thresholds.

% To run through SD and IT sequentially, both sets of specific parameters
% need to be defined.
par.bii = [1.15   1 1e-3 1e+4 0];  % killing rate immobility (L/ng/d) (SD only)
par.bid = [0.0142 1 1e-3 1e+4 0];  % killing rate death (L/ng/d) (SD only)
par.bir = [2850   1 1e-3 1e+4 0];  % killing rate recovery (L/ng/d) (SD only)

% par.Fs  = [24.5  1    1 100 1]; % fraction spread of threshold distribution (IT only)
par.Fsi = [20.9    1    1 100 1];  % fraction spread of immobility threshold distribution (IT only)
par.Fsd = [26.2    1    1 100 1];  % fraction spread of death threshold distribution (IT only)
% NOTE: this code allows to use EITHER one spread factor for the
% thresholds, or to use separate ones. Use either Fs OR Fsi and Fsd. It
% is advisable to always use two thresholds for IT.

% When using fast TK or fast damage repair, we don't need to fit all rate constants.
if glo.fastslow(1) == 1 % fast TK, so ke is not used
    par.ke([1 2 4])  = [100 0 100];  % elimination rate constant, d-1
end
if glo.fastslow(2) == 1 % fast repair, so kr is not used
    par.kr([1 2 4])  = [100 0 100];  % damage repair rate constant, d-1
end

% After optimisation, you can copy-paste relevant lines (fitted basic parameters or updated search ranges) from screen below!


%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

% specify the y-axis labels for each state variable
glo.ylab{glo.locC}   = 'scaled int conc';
glo.ylab{glo.locD}   = 'scaled damage';
glo.ylab{glo.loc_a}  = 'active';
glo.ylab{glo.loc_i}  = 'immobile';
glo.ylab{glo.loc_d}  = 'dead';
glo.ylab{glo.loc_id} = 'immobile+dead';
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'conc. '; % legend label before the 'scenario' number
glo.leglab2 = 'ng/L'; % legend label after the 'scenario' number
% Note: legend labels here will be ignored since a label-table was defined
% above.

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo.

% -------------------------------------------------------------------------
% Configurations for the ODE solver
glo.stiff       = [0 3]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Second argument for glo.stiff is for default sloppy (0), normally tight
% (1), tighter (2), or very tight (3) tolerances. Use 0 or 1 for quick
% analyses, but check with higher values to see if there is a difference!
% There can be large differences between the settings, depending on the
% system and the parameter values.
glo.break_time  = 1; % break time vector up for ODE solver (1) or don't (0)
% Note: breaking the time vector is a good idea when the exposure scenario
% contains discontinuities. Don't use for continuous splines (type 1) as it
% will be much slower. For FOCUS scenarios (high time resolution), breaking
% up is not efficient and does not appear to be necessary.
% -------------------------------------------------------------------------

% Options specific for simplex optimisation
opt_optim.fit   = 1; % fit the parameters (1), or don't (0)
opt_optim.it    = 1; % show iterations of the optimisation (1, default) or not (0)
opt_optim.simno = 2; % for simplex: number of runs of fminsearch (starting again from previous best)

% Options for plotting with calc_and_plot
opt_plot.bw     = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot  = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend

basenm_rem       = glo.basenm; % remember basename as we will modify it!

% Select whether to fit the control treatment (0), the toxicity treatments
% (1; the control is shown as well), or both together (2). A good strategy
% is to fit the background hazard to the control data first (fit_tox=0),
% copy the best values into the parameter matrix above, and then fit the
% toxicity parameter to the complete data set (fit_tox=1). The code below
% does that automatically, and keeps the parameters fixed that need to be
% fixed.
% 
% Select what to fit with fit_tox (this is a 2-element vector). First 
% element of fit_tox is which part of the data set to use:
%   fit_tox(1) = -1  control survival (c=0) only
%   fit_tox(1) = 0   not used for GUTS
%   fit_tox(1) = 1   fit tox parameters, but when fitting, keep hb fixed; run through all
%               elements in SEL and CONFIG sequentially and provide a table at the end
%   fit_tox(1) = 2   fit tox parameters, but when fitting, also fit hb; run through all
%               elements in SEL and CONFIG sequentially and provide a table at the end
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
% CONFIG is a vector with the damage configurations to to run through, sequentially:
% 1) TK->D, TK->immob, D->death
% 2) TK->D, D->immob/death
% 3) two separate reduced GUTS, D1->immob, D2->death
 
% NOTE: for the immobility package, fit_tox = 0 will fit hb only, on the
% control treatment using the same likelihood function as for the effects
% (so least squares at the moment). Further, it uses only the deaths data
% for this; if there are immobiles, we need to think a bit.
% 
% Note: the code below assumes that all controls use the identifier 0
% (zero). 

% % These are the death mechanisms and feedback configurations that will be run automatically
% % Code below will fit all 6 permutations.
% SEL    = [1;2]; % death mechanisms
% CONFIG = [1;2;3]; % model configurations

% To demonstrate, better keep it simpler to speed up calculations
SEL    = [1]; % death mechanisms
CONFIG = [1 3]; % model configurations

% ===== FITTING CONTROLS ==================================================
% Simplex fitting works fine for background hazard
fit_tox = [-1 1]; % set by user (see table above)
opt_optim.type = 1; % optimisation method: 1) default simplex, 4) parspace explorer
opt_plot.statsup = []; % vector with states to suppress in plotting fits
par_out = automatic_runs_guts_immob(fit_tox,par,[],SEL,CONFIG,opt_optim,opt_plot); % script to run the calculations and plot, automatically
% script to run the calculations and plot, automatically
par = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par

% ===== FITTING TOX DATA ==================================================
% Here, we use the parameter-space explorer for fitting the treatments.
% Note that, with fit_tox(1)=1, the tox parameters in par are replaced
% (when fitted) with estimates based on the data set using
% startgrid_guts_immob (called in automatic_runs_guts_immob). You can set
% skip_sg=0 to use the ranges in par (as defined above) instead.

fit_tox = [1 1]; % set by user

opt_optim.type     = 4; % optimisation method 1) simplex, 4 parameter-space explorer
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);
opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots to monitor progress of parameter-space explorer

% -------------------------------------------------------------------------
% THIS BLOCK: SETTINGS FOR SOMEWHAT ROUGH BUT QUICKER EXPLORATION 
% In general, it is a good idea to first perform an exploratory
% optimisation. The settings below skip the profiling step, but better
% search ranges will be printed on screen after optimisation, formatted to
% be copy-pasted in this script above. (for a refined run, make sure that
% skip_sg=1, otherwise startgrid_guts_immob will overwrite the par
% structure again)
opt_optim.ps_profs = 0; % when set to 1, makes profiles and additional sampling for parameter-space explorer
glo.stiff          = [2 1]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
opt_optim.ps_rough = 1; % parameter-space explorer, 0) settings as in openGUTS, 1) rough settings, 2) extra rough
% Note: extra rough may be good enough to run through large numbers of
% model options, but the true optimum may be missed!

% The following options change the settings of startgrid for automatically
% extracting search ranges for tox parameters. These settings are only used
% for the parspace explorer.
skip_sg    = 0; % set to 1 to skip startgrid completely (use ranges in <par> structure)
glo.lim_k  = 0; % set to 1 to limit both rate constants ke and kr to 0.01-10 d-1
glo.mw_log = 1; % set to 1 to use log-scale for all thresholds

% % You can change the basenm to keep the saved plots and MAT files separate.
% glo.basenm = [glo.basenm,'_prelim']; % add 'prelim' so results of initial run are kept separate
% % This would also require you to use the definition of glo.mat_nm with
% % _prelim in the plotting section.

% -------------------------------------------------------------------------
% % THIS BLOCK: SETTINGS TO REFINE FROM RANGES COPIED FROM SCREEN INTO THIS SCRIPT ABOVE
% % After the exploratory run, take the parameter results from screen (where
% % it says that it can be copy-pasted into your script). Rerun with the
% % settings below. 
% glo.stiff          = [2 2]; % ODE solver ode15s with somewhat tighter tolerances if ode45 is too slow
% % glo.stiff          = [0 3]; % ODE solver ode45 with tighter tolerances
% opt_optim.ps_profs = 1; % make profiles now
% skip_sg            = 1; % set to 1 to skip startgrid completely (use ranges in par structure instead)
% -------------------------------------------------------------------------

opt_plot.statsup = []; % vector with states to suppress in plotting fits
[par_out,best_sel] = automatic_runs_guts_immob(fit_tox,par,skip_sg,SEL,CONFIG,opt_optim,opt_plot);
% Note: automatic_runs will return the BEST parameter set in par_out.
disp_settings_guts_immob % display settings on screen for archiving
print_par(par_out) % print the complete best parameter vector that can be copied-pasted
% This includes the control parameters (the parspace explorer itself will
% only plot the results for the fitted parameters).

%% Plot results with or without confidence intervals
% The following code can be used to make plots with confidence intervals.
% Options for confidence bounds on model curves can be set using opt_conf
% (see prelim_checks). The plot_tktd function makes multiplots for the
% data, which are more readable when plotting with various intervals.

% Here, CIs are only plotted when the parspace explorer was used
if opt_optim.type == 4 % for the parspace explorer, we can make CIs already
    opt_conf.type = 3; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
else % otherwise, plot now without CIs
    opt_conf.type = 0; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
end

opt_conf.lim_set = 2; % use limited set of n_lim points (1) or outer hull (2, not for Bayes) to create CIs
opt_tktd.repls   = 1; % plot individual replicates (1) or means (0)
opt_tktd.obspred = 0; % plot predicted-observed plots (1) or not (0)
opt_tktd.max_exp = 0; % set to 1 to maximise exposure/damage plots on exposure rather than on damage
opt_tktd.flip    = 1; % set to 1 to flip row 1 and 2 of the plot around (so int. conc. is on row 1)
opt_tktd.plotexp = 1; % set to 1 to plot exposure profile as area in damage plots

% Below some tricks to allow plotting results for the best configuration.
% The automatic_runs has modified glo.basenm so we will reconstruct it to
% load the correct MAT file.
glo.sel       = SEL(best_sel(1));    % change global for sel
glo.damconfig = CONFIG(best_sel(2)); % change global for feedback configuration
glo.mat_nm    = [basenm_rem,'_sel',sprintf('%d',glo.sel),'_config',sprintf('%d',glo.damconfig)];
% glo.mat_nm    = [basenm_rem,'_prelim_sel',sprintf('%d',glo.sel),'_config',sprintf('%d',glo.damconfig)];
% %  Use the latter one if you added _prelim to the basenm above

plot_tktd([],opt_tktd,opt_conf); % calculate and plot dedicated TKTD plots.
% Note that first input is left empty: this is for the parameter structure,
% which is then obtained from the saved sample.

% Note: it is also possible to modify glo.basenm instead of glo.mat_nm. The
% plots will be the same, but saved plots will always get a name based on
% glo.basenm.

% % Below some tricks to allow plotting results for specific configurations.
% % The automatic_runs has modified glo.basenm so we will reconstruct it to
% % load the correct MAT file.
% glo.sel       = SEL(1); % change global for MoA
% glo.damconfig = CONFIG(1); % change global for feedback configuration
% glo.mat_nm    = [basenm_rem,'_sel',sprintf('%d',glo.sel),'_config',sprintf('%d',glo.damconfig)];
% plot_tktd([],opt_tktd,opt_conf); % calculate and plot dedicated TKTD plots.

%% Calculate ECx versus time
% Here, the ECx is calculated at several time points. When sufficient
% points are specified, a smooth line for ECx versus time will be produced.
% ECx values are also printed on screen. If a sample from parameter space
% is available (e.g., from the slice sampler or the likelihood region), it
% can be used to calculate confidence bounds. Note that opt_conf.type=-1
% skips CIs.
% 
% With the GUTS immobility package, for now, ECx is calculated for the
% category 'active' only. For that trait, it is easy to calculate ECx
% relative to the control (for death and immobile, the control is zero).

% Here, CIs are only plotted when the parspace explorer was used
if opt_optim.type == 4 % for the parspace explorer, we can make CIs already
    opt_conf.type = 3; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
else % otherwise, plot now without CIs
    opt_conf.type = 0; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
end

opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
Tend = [1:4]; % times at which to calculate LCx, relative to control

% Below some tricks to allow plotting results for the best configuration.
% The automatic_runs has modified glo.basenm so we will reconstruct it to
% load the correct MAT file.
glo.sel       = SEL(best_sel(1));    % change global for sel
glo.damconfig = CONFIG(best_sel(1)); % change global for feedback configuration
glo.mat_nm    = [basenm_rem,'_sel',sprintf('%d',glo.sel),'_config',sprintf('%d',glo.damconfig)];
% Note: it is also possible to modify glo.basenm instead of glo.mat_nm. The
% plots will be the same, but saved plots will always get a name based on
% glo.basenm.

% % UNCOMMENT FOLLOWING LINE(S) TO CALCULATE
% opt_ecx.Feff      = [0.10 0.50]; % effect levels (>0 en <1), x/100 in ECx
% opt_ecx.notitle   = 1; % set to 1 to suppress titles above ECx plots
% calc_ecx([],Tend,opt_ecx,opt_conf); % general method for ECx values
% % Note: the general calc_ecx can also work for immobility data