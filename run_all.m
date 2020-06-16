

[xdf_fname, xdf_path] = uigetfile('*.xdf', 'Select a File','MultiSelect','on');


% Handle the case of a single selection which imports as a character vector
if ~iscell(xdf_fname)
    xdf_fname = cellstr(xdf_fname);
end

if MechEngTbl && numel(xdf_fname) > 1 
    error(['Sorry, runall2() does not support the MechEng data output type if'...
        'more than one .xdf file has been selected!']);
end

if numel(xdf_fname) > 4
    error('Sorry, but runall() currently only deals with four XDF files.');
end

    
    
[~, ~, ~, el_dat_x, el_time_axis, ob_dat_x, ob_time_axis,...
    ev_dat, ev_ts] = xdf4runall('false', xdf_path, xdf_fname{II});

% Parameters
% calibration_factor_eog = 0.405666328;
% calibration_factor_elink = 0.075804808;
remove_outliers = false;
opbci_threshold = [15 100]; % Threshold for filtering EOG saccade
elink_threshold = [-inf inf];

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
elink_struct = struct("time_stamps", el_time_axis, "time_series", el_dat_x-1000);

calibration_factor_eog = linear_calibration(event_struct, opbci_struct,...
    opbci_threshold);
calibration_factor_elink = linear_calibration(event_struct, elink_struct,...
    elink_threshold);



ob_features = get_features(event_struct, opbci_struct, calibration_factor_eog,... 
    "EOG", remove_outliers);
elink_features = get_features(event_struct, elink_struct, calibration_factor_elink,...
    "EyeLink", remove_outliers);

sumstats = ["Calibration factor - EOG: " + calibration_factor_eog;
    "EOG Mean Accuracies: " + num2str(ob_features.accuracy);
    "EOG Max Accuracy: " + ob_features.max_accuracy;
    "Calibration factor - Elink: " + calibration_factor_elink;
    "Elink Mean Accuracies: " + num2str(elink_features.accuracy);
    "Elink Max Accuracy: " + elink_features.max_accuracy];

points = ['A' 'B' 'C' 'D'];
nStat = 2; % Standard Deviation and Mean
nTrials = 10; % Ten trials per point
nVars = 8; % Number of variables we are declaring

% Now we only need to change these parameters, and our code will always
% work without having to debug/edit extensively. 
nCols = nVars * numel(points); 
nRows = nTrials;
copy_table = cell(nRows, nCols);
errors = NaN(2,nCols);
headers = strings(1,nCols);
idx = reshape((1:nCols), nVars, numel(points))';
condTag = string(regexp(xdf_fname, 'C\d', 'match'));
subjID = string(regexp(xdf_fname, 'P\d*', 'match'));

for ii = 1:numel(points)
    point = points(ii);
    headers(idx(ii,:)) = [strcat(condTag, "_EOG_MAG_", point),...
        strcat(condTag, "_EOG_ACC_", point),...
        strcat(condTag, "_EOG_PVEL_", point),...
        strcat(condTag, "_EOG_LAT_", point),...
        strcat(condTag, "_ELINK_MAG_", point),...
        strcat(condTag, "_ELINK_ACC_", point),...
        strcat(condTag, "_ELINK_PVEL_", point),...
        strcat(condTag, "_ELINK_LAT_", point)
        ];
    
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
        elink_col = elink_features.all_values.(point); %#ok<UNRCH>
        eog_col = ob_features.all_accuracy.(point);
        diff_col = abs(elink_col-eog_col);
        diff_iqr = iqr(diff_col);
        iqr_filt = abs(eog_col - elink_col) > mean(diff_col) + 1.5*diff_iqr;
        rows_to_remove = [];
       
        for jj = 1:length(ob_features.all_values.(point))
            
            if isnan(eog_col(jj)) || isnan(elink_col(jj)) || iqr_filt(jj)
                rows_to_remove = [rows_to_remove jj];
            end
            
        end
          
        % Remove appropriate rows in to_append
        for jj=1:length(rows_to_remove)
            to_append(rows_to_remove(jj), :) = NaN;
        end
    end
    
    for jj = 1:nVars
        new_col = [nanmean(to_append(:,jj)); nanstd(to_append(:,jj))];
        errors(:,ii*nVars-nVars+jj) = new_col;
    end
    
    copy_table(:,idx(ii,:)) = num2cell(to_append);
end


% Format copy_table for Excel 
% Mag Acc PV Lat for A B C D for EOG, then Elink
copy_table = [cellstr(headers); copy_table; num2cell(errors)];
el_half_idx = cellfun(@(x) contains(x, 'ELINK'), copy_table(1,:)); 
ob_half_idx = ~el_half_idx;
write2table = [copy_table(2:end, ob_half_idx), copy_table(2:end, el_half_idx)];
M = length(write2table);

if ~MechEngTbl

    if II > 1
        ResultsTable((II*M-M+1):II*M) = cell2table(write2table, 'VariableNames',...
            cellstr(headers));
    elseif II == 1
        ResultsTable = cell2table(write2table, 'VariableNames', cellstr(headers));
    end

else
    empty_col = cellstr(strings(size(copy_table,1),1));
    copy_table = [copy_table(:,ob_half_idx), empty_col,...
        copy_table(:,el_half_idx)];
    ReusultsTable = copy_table;
end

N1 = cell2table(num2cell(randn(5,3)),'VariableNames', {'C1_X', 'C1_Y', 'C1_Z'})
N2 = cell2table(num2cell(randn(5,3)),'VariableNames', {'C2_X', 'C2_Y', 'C2_Z'})

II*-32+1


% MechEng Table



close all;

% if CSV == false
%     fSpec1 = strcat(repmat('%-13s',1,32), '\n');
%     fSpec2 = strcat(repmat('%-13.4f ',1,32), '\n');
%     fileID = fopen(subjID + condTag + ".txt", 'w')
%     fprintf(fileID, fSpec1, headers);
%     fprintf(fileID, fSpec2, write_table);
% elseif (CSV == true) && (II == 1)
%     csvwrite


% % Bar plot of accuracies
% figure()
% bins = categorical({'A','B','C','D','Overall'});
% bins = reordercats(bins,{'A','B','C','D','Overall'});
% bars = bar(bins, [ob_features.accuracy; elink_features.accuracy]);
% set(bars, {'DisplayName'}, {'OpenBCI','Eyelink'}')
% ylabel("Angular Difference");
% title("Accuracy");
% legend();
% hold off
% 
% 
% % Bar plot of Peak Velocities
% figure()
% bins = categorical({'A','B','C','D','Overall'});
% bins = reordercats(bins,{'A','B','C','D','Overall'});
% bars = bar(bins, [ob_features.peak_vel; elink_features.peak_vel]);
% set(bars, {'DisplayName'}, {'OpenBCI','Eyelink'}')
% ylabel("Peak Velocity (degrees/s)");
% title("Peak Velocity");
% legend();
% hold off
% 
% % Bar plot of latencies
% figure()
% bins = categorical({'A','B','C','D','Overall'});
% bins = reordercats(bins,{'A','B','C','D','Overall'});
% bars = bar(bins, [ob_features.latency*1000; elink_features.latency*1000]);
% set(bars, {'DisplayName'}, {'OpenBCI','Eyelink'}')
% ylabel("Latency (ms)");
% title("Latency");
% legend();
% hold off
% 
