%% BYOM, script_guts.m, create data sets for testing GUTS.
%
% * Author: Tjalling Jager 
% * Date: Jan 2018
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_guts.html>
%
% This script creates random data sets assuming SD, IT or mixed death
% mechanisms. IT and mixed with log-normal or log-logistic distribution of
% thresholds.
%
% At this moment, all stochasticity is in the death mechanism, but it is
% easy to add stochasticity in TK as well. Damage is not used here explicitly.
% 
% At this moment, only constant exposure is accommodated, though code from
% the regular GUTS analyses can be used to apply time-varying exposure.
% 
%  Copyright (c) 2012-2023, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

%% Initial things

clear, clear global % clear the workspace and globals
format short g

disp(' ')
disp('===============================================================')
disp('Starting simulations')
set(0,'DefaultFigureWindowStyle','docked') ; 
% use this option to collect all figure into one window with tab controls

%% Choose parameter values and experimental set-up

kd  = 1;    % dominant rate constant
mw  = 3;    % median of threshold distribution
bw  = 1;    % killing rate (set high for IT)
hb  = 0.01; % background hazard rate
Fs  = 1;    % factor spread of threshold distribution (set 1 for SD)

cd   = [0 2 4 8 16];  % external concentrations
t    = [0 1 2 3 4]'; % observation times
n    = 10;   % number of simulated individuals
n2   = 100; % many individuals for underlying model curves (0 to skip)
dis  = 2; % 1) log-normal, 2) log-logistic (without statistics toolbox)
nDAT = 3; % number of data sets to produce

disp(' ')
disp('Parameter set used to create the simulated data')
fprintf('par.kd = [%5.5g 1 1e-3 10   1]; \n',kd)
fprintf('par.mw = [%5.5g 1 0    1e6  1]; \n',mw)
fprintf('par.hb = [%5.5g 1 1e-4 1    1]; \n',hb)
fprintf('par.bw = [%5.5g 1 1e-4 1e6  1]; \n',bw)
fprintf('par.Fs = [%5.5g 1 1    100  1]; \n',Fs)
disp(' ')

%% Simulate data sets

X0 = [0;0]; % initial values body residues and integrated hazard rate
options = [];
options = odeset(options,'RelTol',1e-4,'AbsTol',1e-7); % specify tightened tolerances
options = odeset(options,'InitialStep',max(t)/1000,'MaxStep',max(t)/100); % stepsize adaptation

Fs     = max(Fs,1+1e-10); % make sure Fs is always slightly larger than 1
sd_NEC = log10(Fs)/1.96;  % sd on 10log scale, calculated from Fs
beta   = log(39)/log(Fs); % shape parameter for log-logistic from Fs

% mean_nec_dist = 10^(log10(z)+(sd_NEC^2)/2); % this should be the mean NEC from the lognormal

for ndat = 1:nDAT % run through data sets to produce
    
    S_coll  = zeros(length(t),length(cd)); % initialise survival data matrix
    Dw_coll = zeros(length(t),length(cd)); % initialise damage data matrix
    
    for j = 1:length(cd) % run through exposure concentrations
        
        if dis == 1 % lognormal distribution
            mw_range = 10.^(normrnd(log10(mw),sd_NEC,n,1)); % random samples for the NEC
        else
            % % log-logistic distribution with statistics toolbox
            % z_range = random('loglogistic',log(z),1/beta,n,1); % samples for the NEC
            % log-logistic without toolboxes!!
            mw_range = mw*((1./rand(n,1))-1).^(-1/beta); % random samples for the NEC
        end
        
        coll_surv = zeros(length(t),n); % collect the simulated individuals in one matrix
        coll_dam  = zeros(length(t),n); % also collect the damage level of each individual
        
        for i = 1:n % run through individuals
            
            par = [kd;mw_range(i);bw;hb]; % give them their own parameter set
            % feel free to change other parameters per individual here too
            % par = [normrnd(ke,0.1*ke);z_range(i);kk;h0;PVd]; % normal noise on ke with CV=10%
            
            [tn1, Xout] = ode15s(@deri, t, X0, options, cd(j), par); % calculate ODEs
            
            Dw_out = Xout(:,1);   % damage level
            H_out  = Xout(:,2);   % integrated hazard over time
            S_out  = exp(-H_out); % survival probabilities from hazard rates
            
            p = [-diff(S_out);1+sum(diff(S_out))]; % multinomial probabilities
            p = max(0,p); % sometimes it gets slightly negative ...
            
            deaths = mnrnd(1,p); % take a sample from multinomial distribution
            S_sim = [1 1-cumsum(deaths)]; % and translate to survival data
            S_sim(end) = []; % last interval is after end of test, so remove
            
            coll_surv(:,i) = S_sim;  % store the sample for each individual
            coll_dam(:,i)  = Dw_out; % store the damage as well
        end
        
        S_coll(:,j)  = sum(coll_surv,2); % sum individuals to survivors in cohort
        Dw_coll(:,j) = mean(Dw_out,2);   % also store the mean damage
        
    end
    
    % put the data set on-screen
    disp(['Simulated data set ',num2str(ndat),' of ',num2str(nDAT)])
    DATA{ndat} = [-1 cd;t S_coll]; % data set!
    disp(DATA{ndat})
    
    % and make a plot of it
    h1{ndat} = figure;
    subplot(1,2,1)
    hold on
    set(gca,'LineWidth',1,'FontSize',12) % adapt axis formatting
    plot(tn1,Dw_coll,'ko','LineWidth',1,'MarkerFaceColor','y')
    xlabel('time (days)','FontSize',12)
    ylabel('scaled damage','FontSize',12)
    
    subplot(1,2,2)
    hold on
    set(gca,'LineWidth',1,'FontSize',12) % adapt axis formatting
    plot(tn1,S_coll,'ko:','LineWidth',1,'MarkerFaceColor','y')
    xlabel('time (days)','FontSize',12)
    ylabel('survivors','FontSize',12)
    ylim([0 n])
    
    
%     % You can calculate the LC50 using trimmed Spearman-Karber (end of test)
%     figure
%     hold on
%     set(gca,'LineWidth',1,'FontSize',12) % adapt axis formatting
%     
%     [mu,gsd,left,right] = TSK(cd,n-S_coll(end,:),n*ones(size(cd)),0.025,0.975);
%     
%     disp('===============================================================')
%     fprintf('Trimmed Spearman-Karber LC50: %g (%g - %g) \n',mu,left,right)
%     % disp('===============================================================')
%     plot(cd,S_coll(end,:)/(S_coll(end,1)),'ko:')
%     xlabel('concentration','FontSize',12)
%     ylabel('survival probability','FontSize',12)
%     ylim([0 1.02])

end
disp('===============================================================')

%% Now for the model lines (for large number of animals)

disp(' ')
% return % uncomment to NOT plot model lines (faster!)

drawnow

if n2 > 0 
    disp('Calculating the underlying model to plot model curve ...')
    t = linspace(0,max(t),100); % many time points
    
    S_coll2  = zeros(length(t),length(cd));
    Dw_coll2 = zeros(length(t),length(cd));
    
    if dis == 1 % log-normal distribution
        logmw_range = linspace(log10(mw)-3*sd_NEC,log10(mw)+3*sd_NEC,n2); % regular range of NECs
        prob_range = normpdf(logmw_range,log10(mw),sd_NEC); % with their prob densities
        mw_range   = 10.^logmw_range; % NECs on normal scale
    else % log-logistic
        Fs2 = 999^(1/beta); % This is the fraction spread for 99.9% of the distribution
        mw_range = linspace(mw/(1.5*Fs2),mw*Fs2,n2); % range of NECs
        prob_range = ((beta/mw)*(mw_range/mw).^(beta-1)) ./ ((1+(mw_range/mw).^beta).^2); % pdf for the log-logistic (Wikipedia)
    end
    prob_range   = prob_range / sum(prob_range); % normalise the densities to exactly one
    
    for j = 1:length(cd)
        
        coll_surv2 = zeros(length(t),n2);
        coll_dam2 = zeros(length(t),n2);
        
        for i = 1:n2
            
            par = [kd;mw_range(i);bw;hb];
            [tn2, Xout] = ode15s(@deri, t, X0, options, cd(j), par);
            
            Dw_out = Xout(:,1); % damage
            H_out  = Xout(:,2); % integrated hazard rate
            S_out  = exp(-H_out); % survival probabilities
            
            coll_surv2(:,i) = S_out*prob_range(i); % survival prob times the prob of this NEC
            coll_dam2(:,i) = Dw_out; % store thedamage as well
        end
        
        S_coll2(:,j) = sum(coll_surv2,2); % sum over all NECs in the range
        Dw_coll2(:,j) = mean(coll_dam2,2); % mean damage
    end
        
    for ndat = 1:nDAT % run through data sets to plot model line in there
        
        set(0,'CurrentFigure',h1{ndat})
        subplot(1,2,1)
        hold on
        plot(tn2,Dw_coll2,'b-','LineWidth',2)
        plot(tn1,Dw_coll,'ko','LineWidth',1,'MarkerFaceColor','y') % replot, so markers are on top
        
        subplot(1,2,2)
        plot(tn2,n*S_coll2,'b-','LineWidth',2)
        S_coll = DATA{ndat}(2:end,2:end);
        plot(tn1,S_coll,'ko:','LineWidth',1,'MarkerFaceColor','y')
        
    end
end
