function output = KFWesthInputAll(Q, R, inp, show_plot, event)

points = ['A', 'B', 'C', 'D'];
data_struct = inp;
output = inp.time_series;

% Extract times
for ind = 1:length(event.time_series)
    % Extract tag name, eg, t20_exp_D_end
    tag = strsplit(event.time_series{ind},'_');  
    
    extract = false;
    % If calibration and horizontal
    if startsWith(tag{1}, 't') && strcmp(tag{2}, 'exp') ...
      && strcmp(tag{4}, 'start')&& ismember(tag{3}, points)
            % Start, end
            point = tag{3};
            extract = true;
    % Experimental values
    elseif startsWith(tag{1}, 't') && strcmp(tag{2}, 'calib') ...
      && strcmp(tag{3}, 'H') && strcmp(tag{5}, 'start') && ismember(tag{4}, points)
            point = tag{4};
            extract = true;
    end
    
    if extract
        start = event.time_stamps(ind);
        stop = event.time_stamps(ind+2);
            
        % Single saccade
        saccade_position = data_struct.time_stamps >= start & data_struct.time_stamps <= stop;
        saccade = inp.time_series(saccade_position); 
        start_ind = find(saccade_position, 1, 'first');
        % Replace saccade with filtered values
        output(start_ind:start_ind + length(saccade) - 1) = KFWesthInput(Q, R, saccade, false, point);  
    end

end

if show_plot
    figure();
    plot(inp.time_series,'b');
    hold on; 
    plot(output, 'r'); 
    xlabel('Sample No.');
    ylabel('Signal Magnitude'); 
    title('Filtered Signal'); 
    legend('raw EOG','KF-Westheimer - EOG');
end

end