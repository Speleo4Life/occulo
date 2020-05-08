function [calibration_factor] = linear_calibration(event_struct,...
    data_struct, threshold, plots) 

% event_struct: struct containing time_series (labels) and time_stamps (times)
% data_struct: Same as above, for EOG/Elink data

calibration_angles = [-22 -11 11 22];

% Struct with letter key, matrix value corresponding to each letter 
% [start ping end] for each row
times = struct('A', [], 'B', [], 'C', [], 'D', []);
points = fieldnames(times);
numEvt2Find = 3; % Per trial: Start, Ping, Stop
numCalTrial = 5; % Per point/marker
 
% Get start and stop times
% Extract [start ping end] indices from time series event labels and use
% these to assign the relevant timestamp values to the times struct
for ii = 1:numel(points)
    exp = strcat('t\d*_calib_H_', points{ii});
    idxHzCal = ~cellfun('isempty', regexp(event_struct.time_series, exp)); 
    times.(points{ii}) = reshape(event_struct.time_stamps(idxHzCal),numEvt2Find,...
        nnz(idxHzCal)/numEvt2Find)';
end

% Values of each saccade
values = struct('A', [], 'B', [], 'C', [], 'D', []);

% For each point, add corresponding saccade voltage value to array - taking
% points between (start + (ping-start)/2) and ping
dat_x = data_struct.time_series;


for i = 1:length(points)
    point = points{i};
    
    % Loop through each row of times [start ping stop]
    for row = 1:size(times.(point), 1)
        start = times.(point)(row, 1);
        ping = times.(point)(row, 2);
        stop = times.(point)(row, 3);
        saccade = dat_x((data_struct.time_stamps >= start) &...
              (data_struct.time_stamps < stop));

        
       % Loop through data
%         tic
%         saccade = [];
%         for ind = 1:length(data_struct.time_stamps)
%             
%             if data_struct.time_stamps(ind) > stop
%                 break
%             end
%             
%             if data_struct.time_stamps(ind) >= start
%                saccade = [saccade dat_x(ind)]; % Append to array
%             end
%             
%         end
%         toc
        % Isolate saccade using changepoints
        cutoffs = findchangepts(saccade, "MaxNumChanges", 2); 
        
        % If insufficient changepoints found
        if length(cutoffs) < 2
            cutoffs = [1 length(saccade)];
        end
        
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

% Remove outliers - 1.5*IQR from mean
for i = 1:length(points)
    point = points{i};
    
    val_avg = mean(values.(point));
    val_iqr = iqr(values.(point));
    
    val_filter = abs(values.(point)-val_avg) < 1.5*val_iqr;
    
    % Increase IQR until more than 80% of data points remain
    iqr_scale = 1.5;
    while nnz(val_filter)/numel(val_filter) < 0.8
        iqr_scale = iqr_scale + 0.1;
        val_filter = abs(values.(point)-val_avg) < iqr_scale*val_iqr;
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

if plots
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

end
