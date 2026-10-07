% Daphnia magna exposed to PAHs, data set from Jager et al (2010),
% <http://dx.doi.org/10.1007/s10646-009-0417-z>. Published as case study in
% Jager & Zimmer (2012).
% 
% These are the mean responses for each treatments. Treatments are coded as
% XY, where X is the number of the pyrene exposure, and Y the number of the
% fluoranthene exposure. Single pyrene treatments are thus 10, 20, etc.
% Single fluoranthene exposures are 1, 2, etc. Concentrations in uM.
% Mixture 12 would thus be exposure scenario 1 for pyrene and scenario 2
% for fluoranthene.

% These identifiers will be used for the mixture analyses to select the
% correct part of the data set for each analysis. 
idC = [0 0.1];    % identifiers for control treatment(s)
idA = [10 20 40]; % identifiers for compound A (pyrene)
idB = [1 2 4];    % identifiers for compound B (fluoranthene)
idM = [11 22 33 13 31 44]; % identifiers for mixture

% Cumulative repro (live neonates)
R = [TR	0	0.1	10	20	40	1	2	4	11	22	33	13	31	44
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0
2	0	0	0	0	0	0	0	0	0	0	0	0	0	0
3	0	0	0	0	0	0	0	0	0	0	0	0	0	0
4	0	0	0	0	0	0	0	0	0	0	0	0	0	0
6	0	0	0	0	0	0	0	0	0	0	0	0	0	0
7	0	0	0	0	0	0	0	0	0	0	0	0	0	0
8	1.2	2.1	0	0	0	0	0	0	0	0	0	0	0	0
10	6.5	7.7	3.333333333	8.1	0	10	0.8	1.777777778	2.777777778	0	1	0	0	NaN
12	17.9	25.4	8.222222222	12.1	4.222222222	16.3	2	1.777777778	5.333333333	0	1	0	0	NaN
14	31.8	29.5	25.55555556	29.3	5.347222222	33.3	6.4	1.777777778	11.55555556	0	1	0	0	NaN
16	50.6	48.8	35.88888889	39.52222222	5.347222222	56.2	9.4	1.777777778	18	0	NaN	0	0	NaN
18	61.3	63.24444444	55.66666667	56.39722222	5.347222222	62.5	13.8	1.777777778	19.375	0	NaN	0	0	NaN
20	81	77.13333333	62.29166667	62.14722222	5.347222222	64.1	16.57777778	1.777777778	19.375	0	NaN	0	0	NaN
21	81	77.13333333	87.89166667	62.14722222	5.347222222	64.1	16.57777778	1.777777778	19.375	0	NaN	0	0	NaN];

% The indices i_** are used to find where the specific treatments are in
% the data set, for each endpoint.
i_RC = find(ismember(R(1,2:end),idC)==1);
i_RA = find(ismember(R(1,2:end),idA)==1);
i_RB = find(ismember(R(1,2:end),idB)==1);
i_RM = find(ismember(R(1,2:end),idM)==1);

% Survival
S = [-1	0	0.1	10	20	40	1	2	4	11	22	33	13	31	44
0	10	10	10	10	10	10	10	10	10	10	10	10	10	10
1	10	10	10	10	10	10	10	10	10	10	10	10	10	6
2	10	10	10	10	10	10	10	10	9	10	10	10	10	5
3	10	10	10	10	10	10	10	10	9	10	10	10	10	4
4	10	10	10	10	10	10	10	10	9	10	10	10	10	4
6	10	10	10	10	10	10	10	10	9	10	10	10	10	4
7	10	10	10	10	10	10	10	10	9	10	10	10	10	4
8	10	10	9	10	10	10	10	10	9	10	10	10	10	0
10	10	10	9	10	10	10	10	9	9	8	3	8	6	0
12	10	10	9	10	9	10	10	7	9	7	2	7	4	0
14	10	10	9	10	8	10	10	5	9	7	0	5	4	0
16	10	10	9	9	8	10	10	5	9	7	0	5	4	0
18	10	9	9	8	8	10	10	3	8	7	0	2	2	0
20	10	9	8	8	7	10	9	3	7	4	0	2	2	0
21	10	8	5	8	7	10	9	2	7	4	0	2	2	0];

% The indices i_** are used to find where the specific treatments are in
% the data set, for each endpoint.
i_SC = find(ismember(S(1,2:end),idC)==1);
i_SA = find(ismember(S(1,2:end),idA)==1);
i_SB = find(ismember(S(1,2:end),idB)==1);
i_SM = find(ismember(S(1,2:end),idM)==1);

% Repro weights
Rw = S(2:end,2:end);

% Body length (mm)
L = [TR	0	0.1	10	20	40	1	2	4	11	22	33	13	31	44
0	0.88	0.88	0.88	0.88	0.88	0.88	0.88	0.88	0.88	0.88	0.88	0.88	0.88	0.88
2	1.38	1.34	1.28	1.22	1.26	1.4	1.44	1.32	1.18	1.24	1.18	1.2	1.4	1.14
4	1.96	1.72	1.9	1.84	1.86	1.9	1.88	1.74	1.76	1.78	1.58	1.74	1.74	1.533333333
6	2.28	2.38	2.18	2.2	2.12	2.14	2.14	2.08	2.04	2.14	1.875	1.9	2.1	NaN
8	2.34	2.52	2.44	2.38	2.32	2.56	2.36	2.16	2.16	2.26	2.22	2.16	2.24	NaN
10	2.64	2.56	2.28	2.38	2.58	2.46	2.46	2.48	2.34	2.28	2.4	2.32	2.4	NaN
12	2.66	2.56	2.6	2.66	2.7	2.6	2.58	2.7	2.72	2.46	2.4	2.36	2.62	NaN
14	2.68	2.7	2.58	2.68	2.82	2.76	2.78	2.8	2.78	2.62	2.54	2.6	2.64	NaN
16	2.88	2.78	2.64	2.98	3	2.82	NaN	NaN	3	2.84	NaN	2.8	2.72	NaN
18	2.94	2.84	2.86	3.16	3.04	2.92	3.08	2.9	3.08	2.8	NaN	2.74	2.96	NaN
20	3.06	3.02	2.96	3.1	3.14	3.02	3.22	2.76	3.1	2.98	NaN	2.7	2.98	NaN
21	3.26	3.04	2.94	3.12	3.1	3.06	3	2.84	3	3.08	NaN	2.875	2.86	NaN];

% The indices i_** are used to find where the specific treatments are in
% the data set, for each endpoint.
i_LC = find(ismember(L(1,2:end),idC)==1);
i_LA = find(ismember(L(1,2:end),idA)==1);
i_LB = find(ismember(L(1,2:end),idB)==1);
i_LM = find(ismember(L(1,2:end),idM)==1);

% Body length weights
Lw = [5	5	5	5	5	5	5	5	5	5	5	5	5	5
5	5	5	5	5	5	5	5	5	5	5	5	5	5
5	5	5	5	5	5	5	5	5	5	5	5	5	3
5	5	5	5	5	5	5	5	5	5	4	4	5	0
5	5	5	5	5	5	5	5	5	5	5	5	5	0
5	5	5	5	5	5	5	5	5	5	5	5	5	0
5	5	5	5	5	5	5	5	5	5	5	5	5	0
5	5	5	5	5	5	5	5	5	5	5	5	5	0
5	5	5	5	5	5	0	0	5	5	0	5	5	0
5	5	5	5	5	5	5	5	5	5	0	5	5	0
5	5	5	5	5	5	5	5	5	5	0	4	5	0
5	5	5	5	5	5	5	5	5	5	0	4	5	0];

% Pyr (uM)
Cw1 = [0	0	0.1	10	20	30	40
0	0	0	0.0865	0.173	0.26	0.346
21	0	0	0.0865	0.173	0.26	0.346];

% Flu (uM)
Cw2 = [0	1	2	3	4
0	0.213	0.426	0.64	0.853
21	0.213	0.426	0.64	0.853];

glo.scen_plot = 0; % suppress plots for make_scen
make_scen(2,Cw1,Cw2); % define the exposure scenarios

% Create a label-table to provide more useful legends. For each treatment
% identifier a text string is defined.

% Pyrene single exposure (and controls)
Scenario = [0 0.1 10 20 30 40]';
Label    = {'control';'solvent';'PYR 0.0865';'PYR 0.173';'PYR 0.260';'PYR 0.346'};
% Fluoranthene single exposure (don't add control again)
Scenario = [Scenario; [1 2 3 4]'];
Label    = [Label;{'FLU 0.213';'FLU 0.426';'FLU 0.640';'FLU 0.853'}];
% And add the mixture exposure (don't add control again)
Scenario = [Scenario;[11 22 33 13 31 44]'];
Label    = [Label;{'MIX 11';'MIX 22';'MIX 33';'MIX 13';'MIX 31';'MIX 44'}];

