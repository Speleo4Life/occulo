format long;
% AnalyzeXDF_oct2019;

% Parameters
calibration_factor_eog = 0.405666328;
calibration_factor_elink = 0.075804808;
remove_outliers = true;

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

% calibration_factor_eog = linear_calibration(event_struct, opbci_struct, opbci_threshold)
% calibration_factor_elink = linear_calibration(event_struct, elink_struct, elink_threshold)

ob_features = get_features(event_struct, opbci_struct, calibration_factor_eog, "EOG", remove_outliers);
elink_features = get_features(event_struct, elink_struct, calibration_factor_elink, "EyeLink", remove_outliers);
% elink_features = ob_features;


"Calibration factor - EOG: " + calibration_factor_eog
"EOG Mean Accuracies: " + num2str(ob_features.accuracy)
"EOG Max Accuracy: " + ob_features.max_accuracy

"Calibration factor - Elink: " + calibration_factor_elink
"Elink Mean Accuracies: " + num2str(elink_features.accuracy)
"Elink Max Accuracy: " + elink_features.max_accuracy

copy_table = [];
errors = [];
headers = strings([1,0]);
points = ['A' 'B' 'C' 'D'];
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
    
    if remove_outliers
        % Remove outliers again based on Elink, EOG difference
        elink_col = elink_features.all_values.(point);
        eog_col = ob_features.all_accuracy.(point);
        diff_col = abs(elink_col-eog_col);
        diff_iqr = iqr(diff_col);
        iqr_filt = abs(eog_col - elink_col) > mean(diff_col) + 1.5*diff_iqr;
        rows_to_remove = [];
        for j = 1:length(ob_features.all_values.(point))
            if isnan(eog_col(j)) || isnan(elink_col(j)) || iqr_filt(j)
                rows_to_remove = [rows_to_remove j];
            end
        end     
          
        % Remove appropriate rows in to_append
        for j=1:length(rows_to_remove)
            to_append(rows_to_remove(j), :) = NaN;
        end
    end
    
    for j = 1:size(to_append,2)
        new_col = [nanmean(to_append(:,j));nanstd(to_append(:,j))];
        errors = [errors, new_col];
    end
    copy_table = [copy_table, to_append];
end



% Format copy_table for Excel 
% Mag Acc PV Lat for A B C D for EOG, then Elink
copy_table = [headers; copy_table; errors];
eog_half = [];
elink_half = [];
for i = 1:size(copy_table, 2)/8
    eog_half = [eog_half copy_table(:, (8*(i-1) + 1):(8*i-4))];
    elink_half = [elink_half copy_table(:, (8*(i-1)+ 5):(8*i))];
end
empty_col = string(zeros(size(copy_table,1), 1));
empty_col(:) = "";
copy_table = [eog_half empty_col elink_half];
close all;





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
