% Returns signal with blinks removed
function blinks_removed = flatten_blinks(data, data_time, first_event_time, xdf_fname)
    load Blink_Times.mat Blink_Times;
    file_info = strsplit(xdf_fname(1:end-4), "_");
    participant = file_info{1};
    condition = file_info{2};

    % Add offset by first event time
    blink_times = Blink_Times.(participant).(condition) + first_event_time;

    blinks_removed = data;
    % Flatten signal between blink times
    for row = 1:length(blink_times)
        to_remove = data_time >= blink_times(row, 1) & data_time <= blink_times(row, 2); 
        to_remove_data = data(to_remove);
        if ~isempty(to_remove_data)
            blinks_removed(to_remove) = to_remove_data(1); % Replace with first value in range
        end
    end
    
    figure()
    hold on;
    plot(data_time, data);
    for row = 1:length(blink_times)
        xline(blink_times(row, 1), "r");
        xline(blink_times(row, 2), "b");
    end
end