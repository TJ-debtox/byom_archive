% A small script to check the data sets entered in the immobility script.
% Only the most obvious errors are spotted! (and some errors may be missed)
%
% Author     : Tjalling Jager 
% Date       : January 2024
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2024, Tjalling Jager.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE-MIT.txt file in the root directory of BYOM. 

data_coll = zeros(size(DATA{3}(3:end,2:end))); % intialise matrix with zeros
data_nan  = zeros(size(DATA{3}(3:end,2:end))); % intialise matrix with zeros
nr_data   = 0; % count data sets that have actual data (more than 1 element in the matrix)
t_vect    = DATA{3}(2:end,1);
for i = 3:6 % run through data sets
    data    = DATA{i}(3:end,2:end); % take out relevant data (not first row, not second row, not first column)
    ind_nan = isnan(data); % remember where the NaNs are
    data(ind_nan) = 0; % replace NaNs by zeros
    
    if numel(DATA{i}) < 2
        if i == 3
            error('There must be data for the number of active individuals')
        end
    else
        if ~isequal(t_vect,DATA{i}(2:end,1))
            error('The time vector must be the same for all data sets')
        end
        if ~isequal(DATA{i}(2,2:end),DATA{3}(2,2:end))
            error('The data format is not correct: the second rows in each data set must match (start nr. of individuals in each replicate)')
        end
        if ~isequal(size(DATA{i}),size(DATA{3}))
            error('The data format is not correct: all matrices must be the same size')
        end

        data_coll = data_coll + data; % add the data set to the previous one
        data_nan  = data_nan + ind_nan; % add the data set to the previous one
        nr_data   = nr_data + 1; % count number of data sets with data
    end
end

error_flag = 0;
for j = 1:size(data_coll,1) % run through rows (time points)
    ind_ok   = data_nan(j,:) ~= nr_data; % only check the treatments where not ALL entries across all types is NaN!
    data_all = DATA{3}(2,2:end); % sum of all individuals
    if ~isequal(data_coll(j,ind_ok),data_all(ind_ok))
        if error_flag == 0
            disp(' ')
        end
        disp(['Error for t = ',num2str(t_vect(j+1))])
        error_flag = 1;

        ind_check = find(ind_ok); % turn logicals into sequential indices
        for k = 1:length(ind_check)
            if data_coll(j,ind_check(k)) ~= data_all(ind_check(k))
                diff_check = data_all(ind_check(k)) - data_coll(j,ind_check(k));
                if diff_check > 0
                    disp(['   in treatment column ',num2str(ind_check(k)),' missing ',num2str(diff_check),' individuals'])
                else
                    disp(['   in treatment column ',num2str(ind_check(k)),' extra ',num2str(-diff_check),' individuals found'])
                end
            end
        end
    end
end

if error_flag == 1
    disp(' ')
    error('The data format is not correct: the sum across all categories must match the starting number of animals. See details above.')
end