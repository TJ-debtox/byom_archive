% This script plots all exposure profiles from text files in a selected
% directory.
% 
% Author     : Tjalling Jager 
% Date       : November 2022
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2024, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

oldpath  = pwd;      % remember directory in which we started
selpath  = uigetdir(oldpath,'Select folder from which to plot exposure profiles'); % use Matlab GUI element to select folder!
cd(selpath)          % change directory to the one with the profiles
dirData  = dir;      % get the data for the current directory
dirIndex = [dirData.isdir];  % find the index for directories
fileList = {dirData(~dirIndex).name}'; % get a list of the files

% restrict list to .txt files
coll_rem = [];
for i = 1:length(fileList)
    if ~strcmp(fileList{i}(end-3:end),'.txt') % if this does not have a .txt extension ...
        coll_rem = cat(1,coll_rem,i); % add it to the remove list
    end
    
end
fileList(coll_rem) = []; % remove all non text files

if isempty(fileList)
    error('No .txt files found in the selected directory!')
end

Nmax = 16; % maximum number of panels per multiplot

% and plot them, in portions of Nmax
for j = 1:ceil(length(fileList)/Nmax)

    fileList_tmp = fileList((j-1)*Nmax+1:min(length(fileList),(j-1)*Nmax+Nmax));

    % make one plot with subplots for each state variable
    n = ceil(sqrt(length(fileList_tmp)));
    m = ceil(length(fileList_tmp)/n);

    figure
    for i = 1:length(fileList_tmp)
        subplot(m,n,i)
        hold on
        A = load(fileList_tmp{i});
        plot(A(:,1),A(:,2),'k-')
        xlabel('time')
        ylabel('concentration')
        title(fileList_tmp{i}, 'Interpreter', 'none')
    end

end

cd(oldpath) % change directory back to where we came from