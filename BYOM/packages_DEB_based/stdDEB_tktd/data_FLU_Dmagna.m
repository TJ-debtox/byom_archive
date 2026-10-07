function [DATA,W,LabelTable] = data_FLU_Dmagna(TR,opt,study)

% Handy function to load the data sets from standard text files, and
% prepare them for BYOM analyses. The files in this folder assume that the
% data sets are set up with identifiers for the treatments. Use 0 for the
% control, 0.1 for the solvent control, and identifiers 1:n for the
% treatments. For a second study (study=1) this function will automatically
% modify the data matrices, exposure scenarios and labels to create
% identifiers 100-199. This way, it is easy to work with multiple data
% sets, and combine them in different ways.
% 
% The water flea _Daphnia magna_ exposed to fluoranthene; data set from
% Jager et al (2010), <http://dx.doi.org/10.1007/s10646-009-0417-z>.
% Published as case study in Jager & Zimmer (2012). Here the raw data on
% individual females is used, such that the special functions for censoring
% the reproduction data can be used. When removing zeros from the repro data
% set, consider adding some -1's to indicate moults without reproduction. 
% Otherwise, all zeros will be removed, while they may be a 'true zero'.
%
% * Author: Tjalling Jager
% * Date: January 2022
% * Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2023, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

global glo

% There are two controls, half of the replicates at concentration zero is
% for a solvent control (the second half, with identifier 0.1). Body length
% is determined desctructively on a separate group of daphnids.
% Reproduction and survival are from the same animals.

% Body length in mm, concentrations in uM, time in days.
L = [TR	0	0	0	0	0	0.1	0.1	0.1	0.1	0.1	1	1	1	1	1	2	2	2	2	2	3	3	3	3	3
0	0.9	0.9	0.8	0.9	0.9	0.9	0.9	0.8	0.9	0.9	0.9	0.9	0.8	0.9	0.9	0.9	0.9	0.8	0.9	0.9	0.9	0.9	0.8	0.9	0.9
2	1.4	1.5	1.3	1.3	1.4	1.3	1.4	1.3	1.4	1.3	1.4	1.5	1.4	1.2	1.5	1.4	1.5	1.5	1.5	1.3	1.5	1.3	1.4	1.2	1.2
4	1.8	2	2.1	2	1.9	1.6	1.5	1.7	1.8	2	1.8	1.9	2	2	1.8	2	1.9	1.7	1.9	1.9	1.7	1.7	1.8	1.6	1.9
6	2.2	2.5	2.2	2.3	2.2	2.4	2.2	2.3	2.6	2.4	2.1	2.1	2	2.2	2.3	2.2	2	2.1	2	2.4	2	2.1	2.3	1.9	2.1
8	2.4	2.3	2.4	2.4	2.2	2.3	2.4	2.5	2.7	2.7	2.5	2.5	2.5	2.7	2.6	2.5	2.4	2.3	2.3	2.3	1.8	2.5	2.1	2.1	2.3
10	2.6	2.7	2.6	2.6	2.7	2.4	2.5	2.8	2.6	2.5	2.5	2.6	2.3	2.4	2.5	2.6	2.3	2.3	2.5	2.6	2.7	2.6	2.3	2.3	2.5
12	2.6	2.7	2.6	2.8	2.6	2.5	2.5	2.4	2.7	2.7	2.5	2.6	2.8	2.6	2.5	2.5	2.6	2.6	2.6	2.6	2.5	3	2.8	2.7	2.5
14	2.7	2.6	2.6	2.8	2.7	2.7	2.7	2.7	2.7	2.7	2.6	2.7	2.9	2.8	2.8	2.8	2.7	2.8	2.8	2.8	2.9	2.9	2.5	3	2.7
16	2.8	2.8	2.9	3	2.9	2.7	2.8	2.6	2.9	2.9	2.8	2.9	3	2.7	2.7	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
18	2.7	3.1	2.9	2.8	3.2	2.8	2.8	2.8	2.8	3	3	3	2.9	2.8	2.9	3.2	3.2	3	2.9	3.1	3	2.9	2.8	2.8	3
20	3	3.1	3.3	3	2.9	3	3	3.1	3.1	2.9	3	3.1	2.9	2.9	3.2	3.2	3	3.4	3.1	3.4	3	2.8	2.5	2.8	2.7
21	3.3	3.1	3.2	3.4	3.3	2.9	3.1	3.2	3	3	3.2	3.1	3	3	3	3	3.1	3	3.1	2.8	2.9	2.7	2.7	3.1	2.8
];

L(2:end,2:end) = L(2:end,2:end) * 0.10; % mm to cm for stdDEB
LW = ones(size(L)-1);

% Raw data on reproduction: clutch size of offspring per female (not
% cumulated yet!). Concentrations in uM, time in days.
R = [TR	0	0	0	0	0	0	0	0	0	0	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	1	1	1	1	1	1	1	1	1	1	2	2	2	2	2	2	2	2	2	2	3	3	3	3	3	3	3	3	3	3
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
6	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
8	0	0	0	0	0	0	0	0	0	12	0	10	0	0	7	0	0	0	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
10	0	2	7	0	0	0	14	12	17	1	12	0	0	16	1	10	0	6	0	11	0	16	21	19	0	14	8	11	11	0	2	0	0	5	0	0	0	0	1	0	0	NaN	0	0	12	0	0	4	0	0
12	18	2	22	1	21	0	0	22	4	24	29	20	0	0	31	34	0	0	33	30	12	0	0	0	27	22	0	0	2	0	0	0	0	0	0	0	0	0	0	12	0	NaN	0	0	0	0	0	0	NaN	NaN
14	33	6	0	0	37	22	27	11	3	0	0	27	0	4	0	0	0	10	0	0	5	24	0	25	0	29	27	27	33	0	0	0	12	0	0	0	0	0	3	29	NaN	NaN	0	NaN	0	0	0	0	NaN	NaN
16	0	18	26	0	0	26	44	0	41	33	30	0	0	32	32	28	0	0	36	35	4	29	23	24	32	0	53	35	1	28	0	8	0	10	0	0	0	3	6	3	NaN	NaN	0	NaN	0	0	0	0	NaN	NaN
18	23	0	0	0	34	0	0	8	0	42	43	36	0	0	36	0	NaN	15	0	0	0	0	0	0	0	41	0	0	22	0	0	0	12	0	0	0	0	0	0	32	NaN	NaN	0	NaN	NaN	NaN	0	0	NaN	NaN
20	0	8	24	33	0	29	35	28	40	0	0	0	0	44	0	24	NaN	24	11	22	10	0	0	0	6	0	0	0	0	0	0	14	0	0	0	NaN	11	0	0	0	NaN	NaN	0	NaN	NaN	NaN	0	0	NaN	NaN
21	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	NaN	0	0	NaN	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	NaN	0	0	0	0	NaN	NaN	0	NaN	NaN	NaN	0	NaN	NaN	NaN
];

% Create the data set with cumulatives and weights matrix from R
[Rc,RW] = makerepro_ind(R,opt); % we have individuals here, so use makerepro_ind

S = [-1 0 0.1 1 2 3
0	10	10	10	10	10
2	10	10	10	10	10
4	10	10	10	10	10
6	10	10	10	10	10
8	10	10	10	10	10
10	10	10	10	10	9
12	10	10	10	10	7
14	10	10	10	10	5
16	10	10	10	10	5
18	10	9	10	10	3
20	10	9	10	9	3
21	10	8	10	9	2
];

%% Define exposure scenarios

% Define a scenario for the exposure treatments (uM)
Cw = [0   1     2    3
      0 0.213 0.426 0.853
     21 0.213 0.426 0.853];

% Define exposure scenarios for controls (for id=0.1 this is no exposure,
% as it is the solvent control). This is handy to do separately as there
% are usually less events for the controls.
Cw0 = [0   0 0.1 
       0   0 0 
       21  0 0 ];

% Create a table with nicer custom labels for the legends
Scenario = [0 0.1 1 2 3]'; % scenario identifiers that get a label
Label    = {'control';'solvent control';'T1';'T2';'T3'};

%% Modify the data sets for different studies
% The data have to use identifiers. The first data set uses identifiers
% 0-99 as given. For second study (study=2), the matrices will be modified
% to make the identifiers go from 100-199, the third from 200-299, etc. The
% controls are then 0, 100, 200, etc. and the solvent controls 0.1, 100.1,
% 200.1 etc.

% Modify the scenario identifiers using the study number provided
S(1,2:end)   = S(1,2:end)   + (study-1)*100;
L(1,2:end)   = L(1,2:end)   + (study-1)*100;
Rc(1,2:end)  = Rc(1,2:end)  + (study-1)*100;
Cw0(1,2:end) = Cw0(1,2:end) + (study-1)*100;
Cw(1,2:end)  = Cw(1,2:end)  + (study-1)*100;
Scenario     = Scenario     + (study-1)*100;

glo.scen_plot = 0; % don't make a plot for the solvent scenario
make_scen(2,Cw0); % type 2 creates block pulses (fine for controls and constant exposure)
glo.scen_plot = 1; % but do make a plot for the exposure scenarios 
make_scen(2,Cw); % type 2 creates block pulses (fine for controls and constant exposure)

% Create a table with nice custom labels for the legends
LabelTable     = table(Scenario,Label); % create a Matlab table for the labels

DATA{1,glo.locD} = 0;   % there are never data for scaled damage (state 1)
DATA{1,glo.locL} = L;   % length data
W{1,glo.locL}    = LW;  % length weights data
DATA{1,glo.locR} = Rc;  % reproduction data
W{1,glo.locR}    = RW;  % reproduction weights data
DATA{1,glo.locS} = S;   % survival data

