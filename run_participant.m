% Run on C1 - C4 of a single particpant
% Return tables with and without outliers
function [calibration_factor_eog, calibration_factor_elink, R_value, ...
    to_copy, to_copy_removed, ...
    to_copy_outliers_removed, to_copy_removed_outliers_removed] = run_participant(xdf_path, iqr_scale, filter)

calibration_factor_eog = 0; % Intialise values - these are set when running C1
calibration_factor_elink = 0;
R_value = 0;

xdf_path = [xdf_path, '\'];
trials = dir(fullfile(xdf_path, '*.xdf'));

to_copy = struct();
to_copy_removed = struct();
to_copy_outliers_removed = struct();
to_copy_removed_outliers_removed = struct();

for i = 1:length(trials)
    xdf_fname = trials(i).name;
    condition = extractAfter(xdf_fname(1:length(xdf_fname) - 4),"_");

    [copy_table, copy_removed_table,copy_table_no_outliers, copy_removed_table_no_outliers, ...
        calibration_factor_eog, calibration_factor_elink, R_value] = run_single_condition(xdf_fname, ...
        xdf_path, calibration_factor_eog, calibration_factor_elink, R_value, iqr_scale, filter);
    
    to_copy.(condition) = copy_table;
    to_copy_removed.(condition) = copy_removed_table;
    to_copy_outliers_removed.(condition) = copy_table_no_outliers;
    to_copy_removed_outliers_removed.(condition) = copy_removed_table_no_outliers;
end

end