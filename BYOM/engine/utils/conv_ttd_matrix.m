function [Dout,Wout] = conv_ttd_matrix(Din,Nstep)

% Usage: [Dout,Wout] = conv_ttd_matrix(Din,Nstep)
% 
% This function converts a time-to-death matrix into a regular survival
% matrix, for plotting purposes. The new data matrix gets a minus 1 in the
% upper-left corner, which signals that the data are now a regular survival
% counts matrix (in Din it would have been -4). The time vector in the new
% data format will run from zero to the largest absolute value in the data
% set, with Nstep points. For daily steps, take 1+length of test for Nstep.
% 
% The matrix Din should have a first row with a -4 and the treatment IDs.
% The first column is arbitrary and is not used. Use positive values for
% time-to-death when an individual has died. Use a negative value for the
% time-to-death when an individual has not died (use the last time point
% that it was seen alive, e.g., at the end of the test). If not all columns
% have the same length, use NaN for padding.
%
% Author     : Tjalling Jager 
% Date       : June 2023
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2024, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

D = Din(2:end,2:end);
if any(D(:)==0)
    error('The time-to-death matrix contains one or more zeros, which is impossible!')
end
c = Din(1,2:end); % these are the treatment indicators
t = linspace(0,max(abs(D(:))),Nstep); % new time vector for output

Dnew = nan(length(t),size(D,2));   % this will be the new survival matrix
Wnew = zeros(length(t),size(D,2)); % this will be the new weight (missing animals) matrix

for i=1:size(D,2)
    Dtmp = D(:,i); % results for treatment i
    Dnew(1,i) = sum(~isnan(Dtmp)); % number of individuals at t=0
    for j = 2:length(t)
        Dnew(j,i) = sum(Dtmp(Dtmp>0) > t(j)); % first count the ones that have not died yet
        Dnew(j,i) = Dnew(j,i) + sum(abs(Dtmp(Dtmp<0)) >= t(j)); % and add the ones that have not disappeared yet

        if j < length(t) % count missing/removed individuals
            % Note: this needs a number of animals at the time point where
            % they were still known to be alive.
            Wnew(j,i) = sum(abs(Dtmp(Dtmp<0 & Dtmp>-t(end))) < t(j+1) & abs(Dtmp(Dtmp<0 & Dtmp>-t(end))) > t(j));
        end
    end
end

Dout = [-1 c;t(:) Dnew]; % output matrix with time vector and treatments added
% Wout = [-1 c;t(:) Wnew]; % normally, time and scenario are not included in W, but prelim_checks will trim them off (I added them here for testing)
Wout = Wnew; % just return Wnew as 'weights' matrix, so without times and scenarios (better when called by a plotting function)