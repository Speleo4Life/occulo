function [copy_table, copy_removed_table,...
    copy_table_no_outliers, copy_removed_table_no_outliers, ...
    calibration_factor_eog, calibration_factor_elink, R_value] = ...
    run_single_condition(xdf_fname, xdf_path, calibration_factor_eog, calibration_factor_elink, R_value, iqr_scale, filter)

AnalyzeXDF_oct2019(xdf_fname, xdf_path);

opbci_threshold = evalin('base','opbci_threshold');
elink_threshold = evalin('base','elink_threshold');
Q = evalin('base', 'Q');

condition = extractAfter(xdf_fname(1:length(xdf_fname) - 4),"_");

% Deal with NaN values
% Calculate number of corrupted data
% elink_filled_percent = sum(isnan(el_dat_x))/length(el_dat_x);
% eog_filled_percent = sum(isnan(ob_dat_x))/length(ob_dat_x);

if isnan(ob_dat_x(1))
    ob_dat_x(1) = 0;
end
if isnan(el_dat_x(1))
    el_dat_x(1) = 0;
end
ob_dat_x = fillmissing(ob_dat_x, "previous");
el_dat_x = fillmissing(el_dat_x, "previous");

event_struct = struct("time_stamps", ev_ts);
event_struct.time_series = ev_dat;
opbci_struct = struct("time_stamps", ob_time_axis, "time_series", ob_dat_x);
elink_struct = struct("time_stamps", el_time_axis, "time_series", el_dat_x-1000); % Confirm this?

if ~strcmp(condition, "C1")
    % Apply filter to data
    R_value = find_R_on_all(event_struct, opbci_struct, opbci_threshold, iqr_scale);
    opbci_struct = struct("time_stamps", ob_time_axis, "time_series", ...
        filter(R_value, Q, opbci_struct, false, event_struct));
else
    % If C1, recalculate and assign correct R, linear calibration values
    R_value = find_R(event_struct, opbci_struct, opbci_threshold, iqr_scale);
    opbci_struct = struct("time_stamps", ob_time_axis, "time_series", filter(R_value, Q, opbci_struct, false, event_struct));

    calibration_factor_eog = linear_calibration(event_struct, opbci_struct, opbci_threshold, iqr_scale);
    calibration_factor_elink = linear_calibration(event_struct, elink_struct, elink_threshold, iqr_scale);
end

ob_features = get_features(event_struct, opbci_struct, calibration_factor_eog, "EOG", false, iqr_scale);
ob_features_no_outliers = get_features(event_struct, opbci_struct, calibration_factor_eog, "EOG", true, iqr_scale);
% Only run if not condition 4
if ~strcmp(condition, "C4")
    elink_features = get_features(event_struct, elink_struct, calibration_factor_elink, "EyeLink", false, iqr_scale);
    elink_features_no_outliers = get_features(event_struct, elink_struct, calibration_factor_elink, "EyeLink", true, iqr_scale);
else
    elink_features = ob_features;
    elink_features_no_outliers = ob_features_no_outliers;
end

copy_table = [];
copy_removed_table = [];
errors = [];
copy_table_no_outliers = [];
copy_removed_table_no_outliers = [];
errors_no_outliers = [];

headers = strings([1,0]);
points = ['A' 'B' 'C' 'D'];
diff_cols = []; % Matrix to keep track of difference Elink - EOG
diff_cols_no_outliers = [];
for i = 1:length(points)
    point = points(i);
    headers = [headers, strcat('EOG_Mag_', point), ...
        strcat('EOG_Acc_', point), ...
        strcat('EOG_PV_', point), ...
        strcat('EOG_Lat_', point), ...
        strcat('Elink_Mag_', point), ...
        strcat('Elink_Acc_', point), ...
        strcat('Elink_PV_', point), ...
        strcat('Elink_Lat_', point)
        ];
    
    to_append = [ob_features.all_values.(point)', ...
        ob_features.all_accuracy.(point)', ...
        ob_features.all_peak_vel.(point)', ...
        ob_features.all_latency.(point)', ...
        elink_features.all_values.(point)', ...
        elink_features.all_accuracy.(point)', ...
        elink_features.all_peak_vel.(point)', ...
        elink_features.all_latency.(point)'];
    
    to_append_no_outliers = [ob_features_no_outliers.all_values.(point)', ...
        ob_features_no_outliers.all_accuracy.(point)', ...
        ob_features_no_outliers.all_peak_vel.(point)', ...
        ob_features_no_outliers.all_latency.(point)', ...
        elink_features_no_outliers.all_values.(point)', ...
        elink_features_no_outliers.all_accuracy.(point)', ...
        elink_features_no_outliers.all_peak_vel.(point)', ...
        elink_features_no_outliers.all_latency.(point)'];
    
    % Remove outliers again based on Elink, EOG difference
    elink_col = elink_features_no_outliers.all_values.(point);
    eog_col = ob_features_no_outliers.all_values.(point);
    diff_col = abs(elink_col-eog_col);
    diff_q1q3 = quantile(diff_col, [0.25, 0.75]);
    diff_iqr = iqr(diff_col);
    iqr_filt = ~(diff_col <= diff_q1q3(2) + iqr_scale*diff_iqr);
    
    diff_cols = [diff_cols diff_col'];
    diff_col_no_outliers = diff_col;
    diff_col_no_outliers(iqr_filt) = NaN;
    diff_cols_no_outliers = [diff_cols_no_outliers diff_col_no_outliers'];

    full_rows_to_remove = [];
    eog_rows_to_remove = [];
    for j = 1:length(ob_features_no_outliers.all_values.(point))
        if isnan(elink_col(j))
            full_rows_to_remove = [full_rows_to_remove j];
        elseif isnan(eog_col(j)) || iqr_filt(j)
            eog_rows_to_remove = [eog_rows_to_remove j];
        end
    end     

    % Remove appropriate rows in to_append_no_outliers
    for j=1:length(eog_rows_to_remove)
        to_append_no_outliers(eog_rows_to_remove(j), 1:4) = NaN; % Remove EOG outliers
    end
    for j=1:length(full_rows_to_remove)
        to_append_no_outliers(full_rows_to_remove(j), :) = NaN;
    end
        
    for j = 1:size(to_append,2)
        new_col = [nanmean(to_append(:,j));nanstd(to_append(:,j))];
        errors = [errors, new_col];
        
        new_col_no_outliers = [nanmean(to_append_no_outliers(:,j));nanstd(to_append_no_outliers(:,j))];
        errors_no_outliers = [errors_no_outliers, new_col_no_outliers];
    end
    
    copy_table = [copy_table, to_append];
    copy_table_no_outliers = [copy_table_no_outliers, to_append_no_outliers];
    
    % Find number of removed points for each feature
    % Total missing
    % Total removed
    % Total
    for j = 1:8
        if j == 4
            num_missing = ob_features.num_removed_latency.(point);
            num_missing_no_outliers = ob_features_no_outliers.num_removed_latency.(point);
        elseif j == 8
            num_missing = elink_features.num_removed_latency.(point);
            num_missing_no_outliers = elink_features_no_outliers.num_removed_latency.(point);
        else
            num_missing = 0;
            num_missing_no_outliers = 0;
        end
        copy_removed_table = [copy_removed_table [num_missing;
            sum(isnan(to_append(:,j))) - num_missing;
            size(to_append,1) - num_missing]];
        copy_removed_table_no_outliers = [copy_removed_table_no_outliers [num_missing_no_outliers;
            sum(isnan(to_append_no_outliers(:,j))) - num_missing_no_outliers;
            size(to_append_no_outliers,1) - num_missing_no_outliers]];
    end
end

% Format copy_table for Excel 
% Mag Acc PV Lat for A B C D for EOG, then Elink
copy_table = [headers; copy_table; errors];
copy_removed_table = [headers; copy_removed_table];
eog_half = [];
removed_eog_half = [];
elink_half = [];
removed_elink_half = [];

copy_table_no_outliers = [headers; copy_table_no_outliers; errors_no_outliers];
copy_removed_table_no_outliers = [headers; copy_removed_table_no_outliers];
eog_half_no_outliers = [];
removed_eog_half_no_outliers = [];
elink_half_no_outliers = [];
removed_elink_half_no_outliers = [];

for i = 1:size(copy_table, 2)/8
    eog_half = [eog_half copy_table(:, (8*(i-1) + 1):(8*i-4))];
    elink_half = [elink_half copy_table(:, (8*(i-1)+ 5):(8*i))];
    removed_eog_half = [removed_eog_half copy_removed_table(:, (8*(i-1) + 1):(8*i-4))];
    removed_elink_half = [removed_elink_half copy_removed_table(:, (8*(i-1)+ 5):(8*i))];
    
    eog_half_no_outliers = [eog_half_no_outliers copy_table_no_outliers(:, (8*(i-1) + 1):(8*i-4))];
    elink_half_no_outliers = [elink_half_no_outliers copy_table_no_outliers(:, (8*(i-1)+ 5):(8*i))];
    removed_eog_half_no_outliers = [removed_eog_half_no_outliers copy_removed_table_no_outliers(:, (8*(i-1) + 1):(8*i-4))];
    removed_elink_half_no_outliers = [removed_elink_half_no_outliers copy_removed_table_no_outliers(:, (8*(i-1)+ 5):(8*i))];
end

diff_headers = ["Diff_Mag_A", "Diff_Mag_B", "Diff_Mag_C", "Diff_Mag_D"];
if ~strcmp(condition,"C4")
    empty_col = string(zeros(size(copy_table,1), 1));
    empty_col(:) = "";
    copy_table = [eog_half empty_col elink_half empty_col ... 
        [diff_headers; diff_cols; nanmean(diff_cols); nanstd(diff_cols)]];
    copy_table_no_outliers = [eog_half_no_outliers empty_col elink_half_no_outliers ... 
        empty_col [diff_headers; diff_cols_no_outliers; nanmean(diff_cols_no_outliers); nanstd(diff_cols_no_outliers)]];

    empty_col = string(zeros(size(copy_removed_table,1), 1));
    empty_col(:) = "";
    copy_removed_table = [removed_eog_half empty_col removed_elink_half];
    copy_removed_table_no_outliers = [removed_eog_half_no_outliers empty_col removed_elink_half_no_outliers];
else
    copy_table = [eog_half];
    copy_removed_table = [removed_eog_half];
    copy_table_no_outliers = [eog_half_no_outliers];
    copy_removed_table_no_outliers = [removed_eog_half_no_outliers];
end

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
end
