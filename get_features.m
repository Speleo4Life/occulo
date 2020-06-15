function features_struct = get_features(event_struct, data, calibration_factor, label, remove_outliers, iqr_scale)

event = event_struct;
data_struct = data;
if ~exist('label','var')
    label = "";
end

% Struct of experiment times, arranged in matrix with [start ping end] on
% each row
exp_times = struct('A', [], 'B', [], 'C', [], 'D', []);
points = fieldnames(exp_times);
time_step = data_struct.time_stamps(2)-data_struct.time_stamps(1); % Assume time step is uniform

% Extract times
for ind = 1:length(event.time_series)
    % Extract tag name, eg, t20_exp_D_end
    tag = strsplit(event.time_series{ind},'_');  
    
    if startsWith(tag{1}, 't') ...
      && strcmp(tag{2}, 'exp') ...
      && strcmp(tag{4}, 'start')
        % Validate
        if ~ismember(tag{3}, points)
            error(tag{3} + " not a valid point")
        else
            % Add to matrix
            toAdd = [event.time_stamps(ind) event.time_stamps(ind+1) event.time_stamps(ind+2)];
            exp_times.(tag{3}) = [exp_times.(tag{3}); toAdd];
        end
    end
end

% Values and errors of each saccade
% all_ prefix indicates the struct containing values with outliers marked
% NaN, while without the prefix the outliers are completely removed for
% calculations
exp_values = struct('A', [], 'B', [], 'C', [], 'D', []);
all_exp_values = struct('A', [], 'B', [], 'C', [], 'D', []);
exp_errors = struct('A', [], 'B', [], 'C', [], 'D', []);
all_exp_errors = struct('A', [], 'B', [], 'C', [], 'D', []);
accuracy_values = struct('A', [], 'B', [], 'C', [], 'D', []);
all_accuracy_values = struct('A', [], 'B', [], 'C', [], 'D', []);
peak_vel_values = struct('A', [], 'B', [], 'C', [], 'D', []);
all_peak_vel_values = struct('A', [], 'B', [], 'C', [], 'D', []);
latencies = struct('A', [], 'B', [], 'C', [], 'D', []);
all_latencies = struct('A', [], 'B', [], 'C', [], 'D', []);

% Extract angles from data
data_x = data_struct.time_series(1, :);

% if (is_opbci) % Process data if opbci
%     data_x = detrend(data_x);
%     data_x = bandpass(data_x,[0.6,35],250);
% end

angle_data = data_x*calibration_factor; % Convert to angles

% Loop through data
for i = 1:length(points)
    point = points{i};
    
    % Loop through each row of exp_times [start ping stop]
    for row = 1:size(exp_times.(point), 1)
        start = exp_times.(point)(row, 1);
        ping = exp_times.(point)(row, 2);
        stop = exp_times.(point)(row, 3);
        
        % Extract saccade
        saccade = angle_data(data_struct.time_stamps >= start & data_struct.time_stamps <= stop); 
        
        % Isolate saccade using changepoints
        cutoffs = findchangepts(saccade, "MaxNumChanges", 2); 
        
        get_latency = true; % Flag for considering latency of current saccade
        % If insufficient changepoints found
        if length(cutoffs) < 2 || cutoffs(2) - cutoffs(1) < 3
            get_latency = false; % Disregard this saccade for latency
            cutoffs = [1 length(saccade)];
            
            % Another way to deal with insufficient changepoints - improve
%             if length(cutoffs) == 0
%                 cutoffs = [1 length(saccade)];
%             elseif length(cutoffs) == 1
%                 if opbci.time_stamps(cutoffs(1)) < ping
%                     cutoffs = [cutoffs length(saccade)];
%                 else
%                     cutoffs = [1 cutoffs];
%                 end
%             end     
        end      
        
        %Plot individual saccades
%         figure()
%         plot((1:length(saccade))*time_step,saccade)
%         xline(cutoffs(1)*time_step);
%         xline(cutoffs(2)*time_step);
%         title(['Trial ' point]);

        peaks = findpeaks(abs(saccade(cutoffs(1): cutoffs(2))));
%         peaks = peaks(peaks > threshold(1) & peaks < threshold(2))

        % If there are no peaks in isolated region, then take entire length
        if(isempty(peaks))
            peaks = findpeaks(abs(saccade));
        end
        if(isempty(peaks))
            peaks = abs(saccade);
        end
        if strcmp(point, 'A') || strcmp(point, 'B')
            peaks = -peaks;
        end

        % Append average and std of trial to array
        exp_values.(point) = [exp_values.(point) mean(peaks)];
        exp_errors.(point) = [exp_errors.(point) std(peaks)];
        % Find peak velocity (magnitude)
        peak_vel_values.(point) = [peak_vel_values.(point) max(abs(gradient(saccade, time_step)))];
        % Find latency - NaN for undeterminable latency
        all_latencies.(point) = [all_latencies.(point) NaN];
        if get_latency
            all_latencies.(point)(end) = cutoffs(1)*time_step;
            latencies.(point) = [latencies.(point) cutoffs(1)*time_step];
        end
    end
end

% Get number of removed latency points
num_removed_latencies = struct();
for i = 1:length(points)
    point = points{i};
    num_removed_latencies.(point) = sum(isnan(all_latencies.(point)));
end

% Remove outliers - iqr_scale*IQR 
for i = 1:length(points)
    point = points{i};
    
    acc_q1q3 = quantile(exp_values.(point), [0.25, 0.75]);
    acc_iqr = iqr(exp_values.(point));
    vel_q1q3 = quantile(peak_vel_values.(point), [0.25, 0.75]);
    vel_iqr = iqr(peak_vel_values.(point));
    lat_q1q3 = quantile(latencies.(point), [0.25, 0.75]);
    lat_iqr = iqr(latencies.(point));
    
    % Increase iqr range until 80% data is kept
    % NO 80 percent limit for now
    acc_filter = exp_values.(point) >= acc_q1q3(1) - iqr_scale*acc_iqr & ...
        exp_values.(point) <= acc_q1q3(2) + iqr_scale*acc_iqr;
%     new_iqr_scale = iqr_scale;
%     while nnz(acc_filter)/numel(acc_filter) < 0.8
%         new_iqr_scale = new_iqr_scale + 0.1;
%         acc_filter = exp_values.(point) >= acc_q1q3(1) - new_iqr_scale*acc_iqr & ...
%           exp_values.(point) < acc_q1q3(2) + new_iqr_scale*acc_iqr;
%     end
    
    vel_filter = peak_vel_values.(point) >= vel_q1q3(1) - iqr_scale*vel_iqr & ...
        peak_vel_values.(point) <= vel_q1q3(2) + iqr_scale*vel_iqr;
%     new_iqr_scale = iqr_scale;
%     while nnz(vel_filter)/numel(vel_filter) < 0.8
%         new_iqr_scale = new_iqr_scale + 0.1;
%         vel_filter = peak_vel_values.(point) >= vel_q1q3(1) - new_iqr_scale*vel_iqr & ...
%           peak_vel_values.(point) <= vel_q1q3(2) + new_iqr_scale*vel_iqr;
%     end
    
    lat_filter = all_latencies.(point) >= lat_q1q3(1) - iqr_scale*lat_iqr & ...
        all_latencies.(point) <= lat_q1q3(2) + iqr_scale*lat_iqr;
%     new_iqr_scale = iqr_scale;
%     while nnz(lat_filter)/numel(lat_filter) < 0.8
%         % Deal with existing NaN values
%         if iqr_scale > 10
%             break;
%         end
%         new_iqr_scale = new_iqr_scale + 0.1;
%         lat_filter = all_latencies.(point) >= lat_q1q3(1) - new_iqr_scale*lat_iqr & ...
%           all_latencies.(point) <= lat_q1q3(2) + new_iqr_scale*lat_iqr;
%     end
    
    % Replace outlier values with NaN
    all_exp_values.(point) = exp_values.(point);
    all_exp_errors.(point) = exp_errors.(point);
    all_peak_vel_values.(point) = peak_vel_values.(point);
    
    if remove_outliers
        all_exp_values.(point)(~acc_filter) = NaN;
        all_exp_errors.(point)(~acc_filter) = NaN;
        all_peak_vel_values.(point)(~vel_filter) = NaN;
        all_latencies.(point)(~lat_filter) = NaN;

        % Apply filters
        exp_values.(point) = exp_values.(point)(acc_filter);
        exp_errors.(point) = exp_errors.(point)(acc_filter);
        peak_vel_values.(point)= peak_vel_values.(point)(vel_filter);
        latencies.(point) = all_latencies.(point)(~isnan(all_latencies.(point)));
    end
end

% Average accuracy (mean deviation and std)
%[A B C D Total]
values = [];
values_error = [];
accuracy = [];
accuracy_error = [];
peak_vel = [];
peak_vel_error = [];
latency = [];
latency_error = [];

% Extract features
calibration_angles = [-22 -11 11 22];
for i = 1:length(points) 
    point = points{i};
    accuracy_values.(point) = abs((exp_values.(point) - calibration_angles(i)));
    all_accuracy_values.(point) = abs((all_exp_values.(point) - calibration_angles(i)));
    point_latency = latencies.(point)(~isnan(latencies.(point)));
    
    values = [values mean(exp_values.(point))];
    values_error = [values_error std(exp_values.(point))];
    accuracy = [accuracy mean(accuracy_values.(point))];
    accuracy_error = [accuracy_error std(accuracy_values.(point))];
    peak_vel = [peak_vel mean(peak_vel_values.(point))];
    peak_vel_error = [peak_vel_error std(peak_vel_values.(point))];
    latency = [latency mean(point_latency)];
    latency_error = [latency_error std(point_latency)];
end

% Overall values and error
accuracy = [accuracy mean(accuracy)];
accuracy_error = [accuracy_error norm(accuracy_error)/length(points)];
peak_vel = [peak_vel mean(peak_vel)];
peak_vel_error = [peak_vel_error norm(peak_vel_error)/length(points)];
latency = [latency mean(latency)];
latency_error = [latency_error norm(latency_error)/length(points)];

% Plots

% % Accuracy of each point
% for i = 1:length(points)
%     point = points{i};
%     
%     figure()
%     errorbar(1:length(exp_values.(point)),exp_values.(point), exp_errors.(point), ...
%    's','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','black');
%     xlim([0 length(exp_values.(point))+1])
%     hold on 
%     yline(calibration_angles(i)); % Reference
%     title([label ': Accuracy - point ' num2str(point)]);
%     xlabel('Trial number'); 
%     ylabel('Angle'); 
%     legend('Experimentation data','Calibration value');
% end

% Accuracy of all points
% figure()
% colours = ['r' 'g' 'b' 'y'];
% hold on 
% for i = 1:length(points)
%     point = points{i};
%     scatter(exp_values.(point),linspace(0, 1, length(exp_values.(point))), 20, colours(i), "filled");
%     xline(calibration_angles(i)); % Reference
% end
% hold off

% % Bar plot of accuracies
% figure()
% bins = categorical({'A','B','C','D','Overall'});
% bins = reordercats(bins,{'A','B','C','D','Overall'});
% bar(bins, accuracy);
% ylabel("Relative Difference (x100%)");
% title("Accuracy");
% hold on 
% er = errorbar(bins, accuracy,accuracy_error,accuracy_error);
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% hold off
% 
% % Bar plot of peak velocities
% figure()
% bar(bins,peak_vel);
% ylabel("Peak Velocity (degrees/s)");
% title("Peak Velocities");
% hold on 
% er = errorbar(bins, peak_vel,peak_vel_error,peak_vel_error);
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% hold off
% 
% % Bar plot of latencies
% figure()
% bar(bins,latency*1000);
% title("Latency");
% ylabel("Latency (ms)");
% hold on 
% er = errorbar(bins, latency*1000, latency_error*1000, latency_error*1000);
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% hold off

features_struct = struct( ...
    "all_values", all_exp_values, "values", values, "values_error", values_error, ...
    "all_accuracy", all_accuracy_values, "accuracy", accuracy, ...
    "accuracy_error", accuracy_error,  "max_accuracy", max(accuracy),...
    "all_peak_vel", all_peak_vel_values, "peak_vel", peak_vel, "peak_vel_error", peak_vel_error, ...
    "all_latency", all_latencies, "latency", latency, "latency_error", latency_error,...
    "num_removed_latency", num_removed_latencies);
end
