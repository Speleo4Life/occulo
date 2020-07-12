clear;
addpath('Filters','-end');

%-------------------------------
% Parameters
opbci_threshold = [-inf inf]; % Threshold for filtering saccade
elink_threshold = [-inf inf];
outlier_iqr_scale = 1.5;

% Blinks
remove_blinks = false;

% For WH filter Q = 0.0005
% LR filter Q = 0.000005
% Else Q = 0.012
filter_types = {@Bandpass};
filter = filter_types{1};
Q = 0.000005;

excel_fname = "test.xlsx";
plot_folder_name = "test";
write_excel = true;
save_plots = true;
save_csv = false;

run_all_participants = false;
participants_to_run = ["P19", "P23"];
%-------------------------------
main_path = uigetdir;
xdf_paths = dir(main_path);
xdf_paths = xdf_paths([xdf_paths(:).isdir]);
xdf_paths = xdf_paths(~ismember({xdf_paths(:).name},{'.','..'}));

conditions = ["C1", "C2", "C3", "C4"];

% Format headers
points = ["A" "B" "C" "D"];
if write_excel
    headers = strings(3, 33);
    headers(1, 1) = "EOG";
    headers(1, 18) = "Elink";
    for i = 1:length(points)
        headers(2, (i-1)*4 + 1) = points(i);
        headers(2, (i-1)*4 + 18) = points(i);
        headers(3, (i-1)*4 + 1: (i-1)*4 + 4) = ["Mag" "Acc" "Pvel" "Lat"];
        headers(3, (i-1)*4 + 18: (i-1)*4 + 21) = ["Mag" "Acc" "Pvel" "Lat"];
    end
    headers = [headers strings(3, 1) [strings(2, 4); "Diff_A" "Diff_B" "Diff_C" "Diff_D"]];
    for i = 1:length(conditions)
        if strcmp(conditions(i), "C4")
            xlswrite(excel_fname, headers(:, 1:16), conditions(i), "B1");
            xlswrite(excel_fname, headers(:, 1:16), conditions(i) + "_out", "B1");
        else
            xlswrite(excel_fname, headers, conditions(i), "B1");
            xlswrite(excel_fname, headers, conditions(i) + "_out", "B1");
        end
    end
end

average_matrices = struct("C1", [], "C2", [], "C3", [], "C4", []);
std_matrices = struct("C1", [], "C2", [], "C3", [], "C4", []);
average_matrices_out = struct("C1", [], "C2", [], "C3", [], "C4", []);
std_matrices_out = struct("C1", [], "C2", [], "C3", [], "C4", []);
participant_array = [];

total_to_write = struct("C1", [], "C2", [], "C3", [], "C4", []);
total_to_write_out = struct("C1", [], "C2", [], "C3", [], "C4", []);
calibration_table = strings(length(participant_array) + 2, 6);
calibration_table(1, 2) = "EOG";
calibration_table(1, 5) = "Elink";
calibration_table(2, :) = ["Participant", "C. factor", "Outlier %", "R", "C. factor", "Outlier %"];

row = 3;

for i = 1:length(xdf_paths)
    xdf_path = [main_path, '\', xdf_paths(i).name];
    [~, participant] = fileparts(xdf_path);
    
    if ~run_all_participants && ~ismember(participant, participants_to_run)
        continue;
    end
    
    participant_array = [participant_array convertCharsToStrings(participant)];
    
    [calibration_factor_eog, calibration_factor_elink, R_value, ...
        to_copy, to_copy_removed, ...
        to_copy_outliers_removed, to_copy_removed_outliers_removed] = run_participant(xdf_path, outlier_iqr_scale, filter);
        
    eog_removed_total = 0;
    eog_total_data = 0;
    elink_removed_total = 0;
    elink_total_data = 0;
    % Loop through conditions
    for j = 1:length(conditions)
        condition = conditions(j);
        
        % Keep track of averages and stds
        % Second last row
        average_matrices.(condition) = [average_matrices.(condition); to_copy.(condition)(end - 1, :)];
        average_matrices_out.(condition) = [average_matrices_out.(condition); to_copy_outliers_removed.(condition)(end - 1, :)];
        % Final row
        std_matrices.(condition) = [std_matrices.(condition); to_copy.(condition)(end, :)];
        std_matrices_out.(condition) = [std_matrices_out.(condition); to_copy_outliers_removed.(condition)(end, :)];

        to_write = [strings(size(to_copy.(condition), 1) - 1, 1), to_copy.(condition)(2:end, :)];
        to_write_out = [strings(size(to_copy_outliers_removed.(condition), 1) - 1, 1), ...
            to_copy_outliers_removed.(condition)(2:end, :)];
        to_write(1) = participant;
        to_write(size(to_copy.(condition), 1) - 2, 1) = "Avg";
        to_write(size(to_copy.(condition), 1) -1, 1) = "Std";
        to_write_out(1) = participant;
        to_write_out(size(to_copy_outliers_removed.(condition), 1) - 2, 1) = "Avg";
        to_write_out(size(to_copy_outliers_removed.(condition), 1) -1, 1) = "Std";

        total_to_write.(condition) = [total_to_write.(condition); to_write; ...
            [["Missing"; "Removed"; "Total"] to_copy_removed.(condition)(2:end, :) ...
            strings(3, size(to_write,2)- 1 - size(to_copy_removed.(condition),2))]; strings(1, length(to_write))];
        total_to_write_out.(condition) = [total_to_write_out.(condition); to_write_out; ...
            [["Missing"; "Removed"; "Total"] to_copy_removed_outliers_removed.(condition)(2:end, :) ... 
             strings(3, size(to_write_out,2)- 1 - size(to_copy_removed_outliers_removed.(condition),2))]; strings(1, length(to_write_out))];
        
        eog_removed_total = eog_removed_total + nansum(double(to_copy_removed_outliers_removed.(condition)(end - 1, 1:16)));
        eog_total_data = eog_total_data + nansum(double(to_copy_removed_outliers_removed.(condition)(end, 1:16)));
        if ~strcmp(condition, "C4")
            elink_removed_total= elink_removed_total + nansum(double(to_copy_removed_outliers_removed.(condition)(end - 1, 18:end)));
            elink_total_data = elink_total_data + nansum(double(to_copy_removed_outliers_removed.(condition)(end, 18:end)));
        end
    end
    % Calibration table
    calibration_table(row, :) = [participant_array(end), calibration_factor_eog, eog_removed_total/eog_total_data*100, ...
        R_value, calibration_factor_elink, elink_removed_total/elink_total_data*100];
    row = row + 1;
end

close all;
if write_excel
    num_participants = length(participant_array);
    % Write averages and std's to Excel
    averages_table = [strings(3, 2) headers];
    averages_table_out = [strings(3, 2) headers];
    
    correlation_headers = ["Mag_EOG" "Mag_EL" "Acc_EOG" "Acc_EL" "Pvel_EOG" "Pvel_EL" "Lat_EOG" "Lat_EL"];
    observation_matrix = [];
    observation_matrix_out = [];
    
    for i = 1:length(conditions)
       condition = conditions(i);
       
       current_mean = mean(double(average_matrices.(condition)), 1);
       current_averages = [[condition; strings(length(participant_array) + 3, 1)] ... 
           [[participant_array', average_matrices.(condition)];
            ["Avg" current_mean];
            ["Std" vecnorm(double(std_matrices.(condition)), 2, 1)/length(participant_array)];
            strings(2, size(average_matrices.(condition), 2) + 1)]];
        
        current_mean_out = mean(double(average_matrices_out.(condition)), 1);
        current_averages_out = [[condition; strings(length(participant_array) + 3, 1)] ...
            [[participant_array', average_matrices_out.(condition)];
            ["Avg" current_mean_out];
            ["Std" vecnorm(double(std_matrices_out.(condition)), 2, 1)/length(participant_array)];
            strings(2, size(average_matrices_out.(condition), 2) + 1)]];
        
        if(strcmp(condition, "C4"))
            % Take half of average matrix
            current_averages = [current_averages strings(size(current_averages, 1), size(averages_table, 2) - size(current_averages, 2))];
            current_averages_out = [current_averages_out strings(size(current_averages_out, 1), size(averages_table, 2) - size(current_averages_out, 2))];
        else
            % Observation matrix
            for j = 1:4
                observation_matrix = [observation_matrix; current_mean([1 18 2 19 3 20 4 21]+4*(j-1))];
                observation_matrix_out = [observation_matrix_out; current_mean_out([1 18 2 19 3 20 4 21]+4*(j-1))];
            end
        end
        
        averages_table = [averages_table; current_averages];
        averages_table_out = [averages_table_out; current_averages_out];
        
%         start_row = 4 + (i-1)*(length(participant_array) + 4);
%         averages_table(start_row, 1) = condition;
%         averages_table_out(start_row, 1) = condition;
% 
%         averages_table(start_row:start_row + length(participant_array) - 1, 2:2 + size(average_matrices.(condition), 2)) = ...
%             [participant_array', average_matrices.(condition)];
%         averages_table(start_row + length(participant_array), 2:2 + size(average_matrices.(condition), 2)) = ...
%             ["Avg" mean(double(average_matrices.(condition)), 1)];
%         averages_table(start_row + length(participant_array) + 1, 2:2 + size(std_matrices.(condition), 2)) = ...
%             ["Std" vecnorm(double(std_matrices.(condition)), 2, 1)/length(participant_array)];
% 
%         averages_table_out(start_row:start_row + length(participant_array) - 1, 2:2 + size(average_matrices_out.(condition), 2)) = ...
%             [participant_array', average_matrices_out.(condition)];
%         averages_table_out(start_row + length(participant_array), 2:2 + size(average_matrices_out.(condition), 2)) = ...
%             ["Avg" mean(double(average_matrices_out.(condition)), 1)];
%         averages_table_out(start_row + length(participant_array) + 1, 2:2 + size(std_matrices_out.(condition), 2)) = ...
%             ["Std" vecnorm(double(std_matrices_out.(condition)), 2, 1)/length(participant_array)];
    end

    % Create correlation tables
    [R, P] = corrcoef(observation_matrix);
    [R_out, P_out] = corrcoef(observation_matrix_out);
    
    correlations_sheet = [["Correlation Coefficients" correlation_headers;
        correlation_headers' R] strings(length(R)+1, 1) ... 
        ["P-values" correlation_headers;
        correlation_headers' P]];
    correlations_sheet_out = [["Correlation Coefficients" correlation_headers;
        correlation_headers' R_out] strings(length(R)+1, 1) ... 
        ["P-values" correlation_headers;
        correlation_headers' P_out]];
    
    
    % Write all to excel
    for i = 1:length(conditions)
        condition = conditions(i);
        xlswrite(excel_fname, total_to_write.(condition), condition, "A4");
        xlswrite(excel_fname, total_to_write_out.(condition), condition + "_out", "A4");
    end
    xlswrite(excel_fname, averages_table, "Averages", "A1");
    xlswrite(excel_fname, averages_table_out, "Averages_out", "A1");
    xlswrite(excel_fname, correlations_sheet, "Correlations", "A1");
    xlswrite(excel_fname, correlations_sheet_out, "Correlations_out", "A1");
    xlswrite(excel_fname, calibration_table, "Calibration", "A1");

    % Plot bar graphs
    if save_plots
        generateplots;
    end
end


