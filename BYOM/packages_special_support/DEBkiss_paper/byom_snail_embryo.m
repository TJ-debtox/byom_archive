%% BYOM, byom_snail_embryo.m, DEBkiss with Lymnaea 
%
% * Author: Tjalling Jager 
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_debkiss_paper.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* DEBkiss following Jager et al (2013),
% <http://dx.doi.org/10.1016/j.jtbi.2013.03.011>.
%
% *This script:* The case study in the original DEBkiss paper for the pond
% snail _Lymnaea stagnalis_. This is the embryonic development. Data from
% Horstmann (1958). There is no fitting in this script.
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
global DATAx Wx     % optional global for extra data set with different x-axis
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

% embryo body weight (ug) over time (days)
% note: for some time point, 0.01 d is added to ensure that all times are
% unique, and do not overlap with the next data set (needed for
% mat_combine later)
WtData=[2.01	0.45
    2.71	0.91
    3.21	1.36
    3.61	2.42
    3.8	2.72
    4.0	3.02
    4.01	4.23
    4.02	3.33
    4.2	4.69
    4.3	4.69
    4.4	5.74
    4.5	6.50
    4.9	9.52
    5.0	9.52
    5.2	11.34
    6.0	15.42
    6.1	17.23
    6.3	20.26
    6.5	23.28
    7.0	32.20
    7.3	34.62
    7.9	53.51
    7.91	51.70
    8.6	64.70
    9.0	86.16
    9.8	103.09];

WtData(:,2) = WtData(:,2)/1000; % from ug to mg

% a second data set from the same experiment
% note: for some time point, 0.01 d is added to ensure that all times are
% unique, and do not overlap with the next data set (needed for
% mat_combine later)
WtData2 = [0.0	0.235
    0.7	0.263
    0.9	0.238
    1.0	0.239
    1.1	0.242
    1.2	0.234
    1.21	0.330
    1.3	0.244
    1.6	0.344
    1.7	0.336
    1.9	0.339
    2.0	0.376
    2.3	0.476
    2.5	0.391
    2.51	0.452
    2.52	0.530
    2.7	0.620
    3.0	0.972
    3.01	1.033
    3.3	1.428
    3.6	2.416];

WtData2(:,2) = WtData2(:,2)/1000; % from ug to mg

% converting mg dwt to mm volumetric length, assuming dV=0.1
WtData(:,2)  = (WtData(:,2)/0.1).^(1/3);
WtData2(:,2) = (WtData2(:,2)/0.1).^(1/3);

% combine both data sets into a single one with one time vector
DATA{1} = mat_combine(0,[1 2;WtData],[1 2;WtData2]);

DATA{2}=0; % no data for reproduction
DATA{3}=0; % no data for egg buffer

% Additional data set for respiration in nL O2/hr
% the model predictions for respiration are made in call_deri.m
RespData = [1 2
    0.09139761	0.31700015
    0.3208776	0.48724174
    0.46427758	0.65748333
    0.70807088	0.91284572
    0.94465747	1.2533289
    1.1669975	1.3384497
    1.37495748	1.50869129
    1.60443748	1.67893288
    1.76208403	2.10453686
    2.041744	2.35989925
    2.22097728	2.61526163
    2.49326374	3.38134879
    2.68654341	4.40279834
    2.80832329	4.82840232
    2.98001623	6.02009346
    3.12301583	7.2117846
    3.20175569	7.63738858
    3.63706036	13.42560268
    3.75850659	14.70241461
    3.97831092	21.25671588
    3.98585126	20.32038713
    4.22130348	23.55497736
    4.34932249	26.36396362
    4.47717467	29.59855385
    4.96466118	30.36464101
    4.96979929	35.55700955
    4.97650552	36.74870069
    5.15547189	37.68502944
    5.17815965	34.70580159
    5.46903005	42.96251877
    5.80434176	47.64416253
    5.97086322	43.72860593
    6.08446883	46.70783378
    6.29813413	50.62339037
    6.55410541	56.41160448
    6.9631857	55.9008797
    7.30707205	75.30842111
    7.96555085	79.05373612
    7.9719568	81.01151442
    8.63223728	98.46127753
    9.02420165	105.0155788
    9.77012842	123.5719122];

RespData(2:end,2)=RespData(2:end,2)*24/1000; % from nL/h to uL/d

% another respiration set from the same experiment
RespData2 = [1 2
    0.1	0.605
0.3	0.684
0.5	0.745
0.7	0.924
1.0	1.129
1.2	1.375
1.4	1.629
1.6	2.044
1.8	2.180
2.0	2.771
2.2	2.966
2.5	3.733
2.7	4.424
2.8	5.114
3.0	6.064
3.2	7.896];
RespData2(2:end,2)=RespData2(2:end,2)*24/1000; % from nL/h to uL/d

% combine both data sets into a single additional data set
DATAx{1} = mat_combine(0,RespData,RespData2);

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.
X0mat = [   1  2 % the scenarios (here two embryo feeding levels) 
            0.01 (0.25e-3/0.1)^(1/3)   % initial volumetric length (mm)
            0    0    % initial cumul repro
            0.15   0.15 ]; % initial weight of buffer in egg (mg)
 
% Put the position of the various states in globals, to make sure correct
% one is selected for extra things (e.g., to prevent shrinking in
% call_deri, and for using plot_tktd).
glo.locL = 1; % location of body size in the state variable list
glo.locR = 2; % location of cumulative reproduction in the state variable list

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 
% global parameters as part of the structure glo
  
glo.delM  = 1;     % shape corrector (used in call_deri.m)
% shape corr. set to 1 as data are now on volumetric length
glo.dV    = 0.1;   % dry weight density (used in call_deri.m)
glo.len   = 1;     % switch to fit physical length (0=off, 1=on, 2=on and no shrinking) (used in call_deri.m)
glo.Tlag  = 2.5;   % delay for start development (only used for embryo)
glo.mat   = 0;     % include maturity maint. (0=off, 1=include)
% Note: if glo.len > 0 than the initial state for size in X0mat is length
% too (see call_deri.m)

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
% these are the estimates including the male function
par.sJAm = [0.119   0 0 1e6]; % specific assimilation rate 
par.sJM  = [0.00795 0 0 1e6]; % specific maintenance costs 
par.WB0  = [0.15    0 0 1e6]; % initial dry weight of egg
par.WVp  = [70.3    0 0 1e6]; % body mass at puberty
par.yAV  = [0.8     0 0 1];   % yield of assimilates on volume (starvation)
par.yBA  = [0.55    0 0 1];   % yield of egg buffer on assimilates
par.yVA  = [0.8     0 0 1];   % yield of structure on assimilates (growth)
par.kap  = [0.828   0 0 1];   % allocation fraction to soma
par.fB   = [0.5     0 0 2];   % scaled food level, embryo (adapted model)
par.f1   = [1       0 0 2];   % scaled food level, ad libitum
par.f2   = [1       0 0 2];   % scaled food level, regime 2
par.f3   = [1       0 0 2];   % scaled food level, regime 3
par.WVf  = [0       0 0 1e6]; % half-saturation body weight
    
%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

glo.t = linspace(0,10,100); % time vector for the model simulation in days

% specify the y-axis labels for each state variable
glo.ylab{1} = 'volumetric length (mm)';
glo.ylab{2} = 'cumulative reproduction (eggs)';
glo.ylab{3} = 'egg buffer (mg)';

% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'scen. '; % legend label before the 'scenario' number
glo.leglab2 = ''; % legend label after the 'scenario' number

% additional labels for extra data sets
glo.xlab2{1} = 'time (days)';
glo.ylab2{1} = ['respiration rate (',char(181),'L O_2/d)'];

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

opt_optim.fit    = 0; % fit the parameters (1), or don't (0)
opt_optim.it     = 0; % show iterations of the optimisation (1, default) or not (0)
opt_plot.bw      = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 0; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.statsup = [2]; % vector with states to suppress in plotting fits

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

%% Other files: derivatives
% To archive analyses, publishing them with Matlab is convenient. To keep
% track of what was done, the file derivatives.m can be included in the
% published result.
% 
% <include>derivatives.m</include>

%% Other files: call_deri
% To archive analyses, publishing them with Matlab is convenient. To keep
% track of what was done, the file call_deri.m can be included in the
% published result.
%
% <include>call_deri.m</include>

