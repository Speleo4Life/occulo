function [calibration_factor] = linear_calibration(event_struct, data, threshold, iqr_scale)

event_struct = event_struct;  % struct containing time_series (labels) and time_stamps (times)
data_struct = data; % Same as above, for EOG/Elink data

calibration_angles = [-22 -11 11 22];

% Struct with letter key, matrix value corresponding to each letter 
% [start ping end] for each row
times = struct('A', [], 'B', [], 'C', [], 'D', []);
points = fieldnames(times);
 
% Get start and stop times
% Extract [start ping end] times from time series and time stamps
% Assume ping and stop always follow start, in that order
for ind = 1:length(event_struct.time_series)
    if strcmp(event_struct.time_series(ind), "exp_start")
        exp_start_ind = ind;
        break
    end
    
    % Extract tag name, eg, t1_calib_H_A_start
    tag = strsplit(event_struct.time_series{ind},'_');  
    
    % If calibration and horizontal
    if startsWith(tag{1}, 't') ...
      && strcmp(tag{2}, 'calib') ...
      && strcmp(tag{3}, 'H') ...
      && strcmp(tag{5}, 'start')
        % Validate
        if ~ismember(tag{4}, points)
            error(tag{4} + " not a valid point")
        else
            % Add to matrix
            toAdd = [event_struct.time_stamps(ind) event_struct.time_stamps(ind+1) event_struct.time_stamps(ind+2)];
            times.(tag{4}) = [times.(tag{4}); toAdd];
        end
    end
end

% Values of each saccade
values = struct('A', [], 'B', [], 'C', [], 'D', []);

% For each point, add corresponding saccade voltage value to array - taking
% points between (start + (ping-start)/2) and ping
dat_x = data_struct.time_series;
% opbci_dat_x = data_struct.time_series(1, :);
% opbci_dat_x = detrend(opbci_dat_x);
% opbci_dat_x = bandpass(opbci_dat_x,[0.6,35],250);

for i = 1:length(points)
    point = points{i};
    
    % Loop through each row of times [start ping stop]
    for row = 1:size(times.(point), 1)
        start = times.(point)(row, 1);
        ping = times.(point)(row, 2);
        stop = times.(point)(row, 3);

        saccade = dat_x(data_struct.time_stamps >= start & data_struct.time_stamps <= stop); 
        cutoffs = findchangepts(saccade, "MaxNumChanges", 2); % Isolate saccade using changepoints
        
        % If insufficient changepoints found
        if length(cutoffs) < 2 || cutoffs(2) - cutoffs(1) < 3
%               continue
              cutoffs = [1 length(saccade)];
        end
        
        %Plot individual saccades
%         figure()
%         plot(saccade)
%         xline(cutoffs(1));
%         xline(cutoffs(2));

        % Only take peaks of saccade
        peaks = findpeaks(abs(saccade(cutoffs(1): cutoffs(2))));
        % If there are no peaks in isolated region, then take entire length
        if(isempty(peaks))
            peaks = findpeaks(abs(saccade));
        end
        peaks = peaks(peaks > threshold(1) & peaks < threshold(2));
        
        if strcmp(point, 'A') || strcmp(point, 'B')
            peaks = -peaks;
        end
        values.(point) = [values.(point) peaks]; % Append to all values
    end
end

% Remove outliers - 1.5*IQR
for i = 1:length(points)
    point = points{i};
    
    val_q1q3 = quantile(values.(point), [0.25, 0.75]);
    val_iqr = iqr(values.(point));
    val_filter = values.(point) >= val_q1q3(1) - iqr_scale*val_iqr & ...
        values.(point) <= val_q1q3(2) + iqr_scale*val_iqr;
    
    
    % Increase IQR until more than 80% of data points remain
    new_iqr_scale = iqr_scale;
    while nnz(val_filter)/numel(val_filter) < 0.8
        new_iqr_scale = new_iqr_scale + 0.1;
        val_filter = values.(point) >= val_q1q3(1) - new_iqr_scale*val_iqr & ...
            values.(point) <= val_q1q3(2) + new_iqr_scale*val_iqr;
    end
    
    values.(point) = values.(point)(val_filter);
end

% Calculate mean and error of each point
average_voltages = zeros(1, length(points));
errors = zeros(1, length(points));

for i = 1:length(points)
    average_voltages(i) = mean(values.(points{i}));
    errors(i) = std(values.(points{i}));
end

% Linear fit - m is calibration factor: angle = m*voltage
calibration_factor = calibration_angles/average_voltages;

% Plots
figure()
errorbar(calibration_angles,average_voltages, errors, ...
   's','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','black')
title('Calibration Curve : Voltage vs. Angle Plot');
hold on 
plot(calibration_angles, calibration_angles/calibration_factor);
xlabel('Saccade size - degrees'); 
ylabel('Voltage -  milli volts'); 
legend('Calibration experimental data','Fitted curve');

end
