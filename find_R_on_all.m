function [R] = find_R_on_all(event_struct, data, threshold, iqr_scale)

event_struct = event_struct;  % struct containing time_series (labels) and time_stamps (times)
data_struct = data; % Same as above, for EOG/Elink data


% Struct with letter key, matrix value corresponding to each letter 
% [start ping end] for each row
times = struct('A', [], 'B', [], 'C', [], 'D', []);
points = fieldnames(times);
 
% Values of each saccade
% For each point, add corresponding saccade voltage value to array - taking
% points between (start + (ping-start)/2) and ping
dat_x = data_struct.time_series;
% opbci_dat_x = data_struct.time_series(1, :);
% opbci_dat_x = detrend(opbci_dat_x);
% opbci_dat_x = bandpass(opbci_dat_x,[0.6,35],250);
values = struct('A', [], 'B', [], 'C', [], 'D', []);

% Get start and stop times
% Extract [start ping end] times from time series and time stamps, extract
% saccade
% Assume ping and stop always follow start, in that order
for ind = 1:length(event_struct.time_series)
    % Extract tag name, eg, t20_exp_D_end
    tag = strsplit(event_struct.time_series{ind},'_');  
    
    extract = false;
    % Experimental values
    if startsWith(tag{1}, 't') && strcmp(tag{2}, 'exp') ...
      && strcmp(tag{4}, 'start')&& ismember(tag{3}, points)
            % Start, end
            point = tag{3};
            extract = true;
    % If calibration and horizontal
    elseif startsWith(tag{1}, 't') && strcmp(tag{2}, 'calib') ...
      && strcmp(tag{3}, 'H') && strcmp(tag{5}, 'start') && ismember(tag{4}, points)
            point = tag{4};
            extract = true;
    end
    
    % Extract single saccade
    if extract
        start = event_struct.time_stamps(ind);
        ping = event_struct.time_stamps(ind + 1);
        stop = event_struct.time_stamps(ind+2);
            
        saccade = dat_x(data_struct.time_stamps >= start & data_struct.time_stamps <= stop); 
        cutoffs = findchangepts(saccade, "MaxNumChanges", 2); % Isolate saccade using changepoints
        
        % If insufficient changepoints found
        if length(cutoffs) < 2 || cutoffs(2) - cutoffs(1) < 3
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

R = mean([std(abs(values.A/max(values.A))) std(abs(values.B/max(values.B))) ...
    std(abs(values.C/max(values.C))) std(abs(values.D/max(values.D)))])/1000;
end
