% function [ResultsCell, ResultsTable, SumStats] = runall2(remove_outliers, plots)
% 
% runall2() is an overhaul of the run_all.m analysis script and represents
% a significant improvement in functionality and efficiency. One is now able to 
% select all four .xdf files for a given subject or four .xdf files
% representing four different subjects. In fact, you could compare two
% subjects on one or two conditions. The function was designed to have this
% flexibility, though emphasis is placed on processing the four .xdf files
% for a single participant, and full testing for bugs related to other
% work-flows has not been completed.
%
% FIXED most issues with preallocation and loop-based increases in variable size. 
% 
% 
% INPUT: All OPTIONAL (Defaults explained below)
% remove_outliers <BOOL> 'true' (1) or 'false' (0) indicating
%       whether or not outlier removal is performed. Defaults to false.
%
% plots <BOOL> Logical 'true' (1) or 'false' (0) indicating whether the
%       bar plots for the extracted saccade features are shown.
%
% OUTPUT: ALL OPTIONAL 
% ResultsTable <TABLE> A table with magnitude, accuracy, peak velocity, and
%       latency variables for each trial and condition if requested. 
%
% SumStats <STR ARRAY> A string array with information on the calibration
%           factors and summary stats for the Elink and EOG accuracy data.
%
% xdf_fname <CELL STRING> I added this simply to keep tack of which files I
% was dealing with. Completely optional. 
%
% DEPENDENCIES: get_features.m, linear_calibration.m, and xdf4runall.m
%
% CODE LIMITATIONS: Full bug testing for the case in which .XDF files of the same
% condition are requested has not be completed. This could prove useful if you 
% wanted to compare subjects of the same condition, but obviously there are
% work arounds to this. 
%
% Importantly, the results table will not be available if you have .XDF files
% from the same condition. It will be easy enough to reshape the
% ResultsCell and rename the headers to account for duplicate variable
% names.
%
% Overhaul by Ray MacNeil (UBC, Psychology), May 1, 2020
% GitHub Repository: https://github.com/Speleo4Life/occulo
% Branch: improved_eog_workflow
function [ResultsCell, ResultsTable, SumStats, xdf_fname] = runall2(remove_outliers, plots)

%% Exception Handling
if nargin < 1 || ~exist('remove_outliers', 'var') || isempty(remove_outliers) ||...
        ~islogical(remove_outliers)
    warning(['remove_outliers: no argument or an invalid argument. Defaulting'...
        ' to ''false''.']);
    remove_outliers = false;
    
end 
   
if nargin < 2 || ~exist('plots', 'var') || isempty(plots) || ~islogical(plots)
    warning('plots: no argument or an invalid argument. Defaulting to ''false''');
    plots = false;
end    

%% Get Data for Import

[xdf_fname, xdf_path] = uigetfile('*.xdf', 'Select a File','MultiSelect','on');

% xdf_fname = {'P021_C1.xdf','P021_C2.xdf','P021_C3.xdf', 'P021_C4.xdf'};
% xdf_path = '/Users/ray.macneil/Nextcloud/ubc_occulo/Analysis/Real_Participant_Data/XDFS/';


opbci_threshold = [15 100]; % Threshold for filtering EOG saccade
elink_threshold = [-inf inf];


% Exception handling for file selection
% Handle the case of a single selection which imports as a character vector
if ~iscell(xdf_fname)
    xdf_fname = cellstr(xdf_fname);
end

numFilesXdf = numel(xdf_fname);

if numFilesXdf > 4
    error(['Sorry, but runall2() currently only supports the selection of four'...
        'XDF files.']);
end


% Constants
points = ['A' 'B' 'C' 'D'];
numStat = 2; % Standard Deviation and Mean
numTrials = 10; % Ten trials per point, PER CONDITION
numVars = 8; % Number of variables we are declaring for EACH CONDITION
numCols = numVars * numel(points); % PER EACH CONDITION
numRows = numTrials + numStat;


% Preallocate Key Variables
SumStats = strings(6, numFilesXdf);
Headers = strings(1, numCols .* numFilesXdf);
ResultsTable = table('Size', [numRows, numCols .* numFilesXdf], 'VariableTypes',... 
    repelem({'single'}, numCols .* numFilesXdf)); 
Idx = reshape((1:(numCols * numFilesXdf)), numVars,...
        numel(points) * numFilesXdf)';
ResultsCell = cell(numRows+1, numCols, numFilesXdf);

%% Main Loop where II is the index corresponding to the .xdf file
for II = 1:numel(xdf_fname)
    
    if ~contains(xdf_fname{II}, 'C4')
        [~, ~, ~, el_dat_x, el_time_axis, ob_dat_x, ob_time_axis,...
            ev_dat, ev_ts] = xdf4runall('false', xdf_path, xdf_fname{II});
    else
        [~, ~, ~, ~, ~, ob_dat_x, ob_time_axis,...
            ev_dat, ev_ts] = xdf4runall('false', xdf_path, xdf_fname{II});
    end
    
    if isnan(ob_dat_x(1))
        ob_dat_x(1) = 0;
    end
    
    if isnan(el_dat_x(1))
        el_dat_x(1) = 0;
    end
    
    % Deal with NaN values
    ob_dat_x = fillmissing(ob_dat_x, "previous");
    el_dat_x = fillmissing(el_dat_x, "previous");
    event_struct = struct("time_stamps", ev_ts);
    event_struct.time_series = ev_dat;
    opbci_struct = struct("time_stamps", ob_time_axis, "time_series", ob_dat_x);
    elink_struct = struct("time_stamps", el_time_axis, "time_series", el_dat_x-960);
    
    if contains(xdf_fname{II}, 'C1')
        calibration_factor_eog = linear_calibration(event_struct,...
            opbci_struct, opbci_threshold, true);
        calibration_factor_elink = linear_calibration(event_struct,...
            elink_struct, elink_threshold, true);
    end
    
    ob_features = get_features(event_struct, opbci_struct,...
        calibration_factor_eog, "EOG", remove_outliers);
    
    if ~contains(xdf_fname{II}, 'C4')
        elink_features = get_features(event_struct, elink_struct,...
            calibration_factor_elink, "EyeLink", remove_outliers);
    else
        elink_features = ob_features;
    end
    
    SumStats(:,II) = ["Calibration factor - EOG: " + calibration_factor_eog;
        "EOG Mean Accuracies: " + num2str(ob_features.accuracy);
        "EOG Max Accuracy: " + ob_features.max_accuracy;
        "Calibration factor - Elink: " + calibration_factor_elink;
        "Elink Mean Accuracies: " + num2str(elink_features.accuracy);
        "Elink Max Accuracy: " + elink_features.max_accuracy];
   
    copy_table = cell(numRows, numCols);
    errors = NaN(numStat,numCols);
    condTag = string(regexp(xdf_fname{II}, 'C\d', 'match'));
    subjID = string(regexp(xdf_fname, 'P\d*', 'match'));
    
    for ii = 1:numel(points)
        point = points(ii);
        Headers(1,Idx(II*4-4+ii,:)) = [strcat(condTag,...
            "_EOG_MAG_", point),...
            strcat(condTag, "_EOG_ACC_", point),...
            strcat(condTag, "_EOG_PVEL_", point),...
            strcat(condTag, "_EOG_LAT_", point),...
            strcat(condTag, "_ELINK_MAG_", point),...
            strcat(condTag, "_ELINK_ACC_", point),...
            strcat(condTag, "_ELINK_PVEL_", point),...
            strcat(condTag, "_ELINK_LAT_", point)];
        
        to_append = [ob_features.all_values.(point)', ...
            ob_features.all_accuracy.(point)', ...
            ob_features.all_peak_vel.(point)', ...
            ob_features.all_latency.(point)', ...
            elink_features.all_values.(point)', ...
            elink_features.all_accuracy.(point)', ...
            elink_features.all_peak_vel.(point)', ...
            elink_features.all_latency.(point)'];
        
        if remove_outliers
            % Remove outliers again based on Elink, EOG difference
            elink_col = elink_features.all_values.(point);
            eog_col = ob_features.all_accuracy.(point);
            diff_col = abs(elink_col-eog_col);
            diff_iqr = iqr(diff_col);
            iqr_filt = abs(eog_col - elink_col) > mean(diff_col) + 1.5*diff_iqr;
            rows_to_remove = [];
            
            for jj = 1:length(ob_features.all_values.(point))
                
                if isnan(eog_col(jj)) || isnan(elink_col(jj)) || iqr_filt(jj)
                    rows_to_remove = [rows_to_remove jj]; %#ok<AGROW>
                end
                
            end
            
            % Remove appropriate rows in to_append
            for jj=1:length(rows_to_remove)
                to_append(rows_to_remove(jj), :) = NaN;
            end
        end
        
        for jj = 1:numVars
            
            new_col = [nanmean(to_append(:,jj)); nanstd(to_append(:,jj))];
            errors(:,ii*numVars-numVars+jj) = new_col;
        end
        
        copy_table(:,Idx(ii,:)) = vertcat(num2cell(to_append),...
            num2cell(errors(:,Idx(ii,:))));
    end
    
Z = numel(points);
CopyHead = reshape(Headers(Idx((II*Z-Z+1):(II*Z),:)).',1, numCols);


% Format copy_table for Excel
% Mag Acc PV Lat for A B C D for EOG, then Elink
copy_table = [cellstr(CopyHead); copy_table]; %#ok<AGROW>
el_half_idx = cellfun(@(x) contains(x, 'ELINK'), copy_table(1,:));
ob_half_idx = ~el_half_idx;
ResultsCell(:,:,II) = [copy_table(:, ob_half_idx), copy_table(:, el_half_idx)];
    
    
if plots


    % Bar plot of accuracies
    figure()
    bins = categorical({'A','B','C','D','Overall'});
    bins = reordercats(bins,{'A','B','C','D','Overall'});
    bars = bar(bins, [ob_features.accuracy; elink_features.accuracy]);
    set(bars, {'DisplayName'}, {'OpenBCI','Eyelink'}')
    ylabel("Angular Difference");
    title("Accuracy");
    legend();
    hold off


    % Bar plot of Peak Velocities
    figure()
    bins = categorical({'A','B','C','D','Overall'});
    bins = reordercats(bins,{'A','B','C','D','Overall'});
    bars = bar(bins, [ob_features.peak_vel; elink_features.peak_vel]);
    set(bars, {'DisplayName'}, {'OpenBCI','Eyelink'}')
    ylabel("Peak Velocity (degrees/s)");
    title("Peak Velocity");
    legend();
    hold off

    % Bar plot of latencies
    figure()
    bins = categorical({'A','B','C','D','Overall'});
    bins = reordercats(bins,{'A','B','C','D','Overall'});
    bars = bar(bins, [ob_features.latency*1000; elink_features.latency*1000]);
    set(bars, {'DisplayName'}, {'OpenBCI','Eyelink'}')
    ylabel("Latency (ms)");
    title("Latency");
    legend();
    hold off

end
    
end
Headers = reshape(ResultsCell(1,:,:), 1, numCols .* numFilesXdf);
ResultsTable(:,:) = reshape(ResultsCell(2:end,:,:), numTrials+numStat,...
    numCols .* numFilesXdf);
Cols2Del = cellfun(@(x) contains(x, 'C4_ELINK'), Headers);
FinalHeaders = Headers(~Cols2Del);
ResultsTable(:,Cols2Del) = [];
try
ResultsTable.Properties.VariableNames = FinalHeaders;
catch
fprintf('Variable Name Conflicts. Format the table from the ''ResultsCell'' output.')
end
end