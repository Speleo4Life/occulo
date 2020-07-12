clear;
close all;
[xdf_fname, xdf_path] = uigetfile('..\Data\*.xdf', 'Select a File');
[event_struct, opbci_struct, elink_struct] = get_raw_data(xdf_fname, xdf_path);
csvs = dir(xdf_fname(1:end-4) + "*.csv");
csv_fname = csvs.name;
fname = strsplit(csv_fname(1: end-4),'_');
calibration_factor_elink = str2double(fname{3});
data = readmatrix(csv_fname, 'OutputType', 'string');
filternames = struct("LinearRecipoAll", "Linear Reciprocal", "KFWesthInputAll", "Westheimer", ...
    "Bandpass", "Bandpass");
participant = fname{1};

calib_trials = [];
exp_trials = [40 34 20 11 5];

raw_eog = opbci_struct.time_series*0.01;
elink = elink_struct.time_series*calibration_factor_elink;

if isnan(raw_eog(1))
    raw_eog(1) = 0;
end
if isnan(elink(1))
    elink(1) = 0;
end
raw_eog = fillmissing(raw_eog, "previous");
elink = fillmissing(elink, "previous");

% Create struct of filtered signals for each saccade
filtered_eog = struct();
for row = 1:size(data, 1)
    % Format as follows [name calibration_factor [data]]
    filtered_eog.(data(row, 1)) = double(data(row, 3:end))*double(data(row, 2));
end
filters = fieldnames(filtered_eog);

points = ["A" "B" "C" "D"];
% Assume ping and stop always follow start, in that order
for ind = 1:length(event_struct.time_series)
    % Extract tag name, eg, t20_exp_D_end
    tag = strsplit(event_struct.time_series{ind},'_');  
    exp_number =  tag{1}(2:end);

    extract = false;
    % Experimental values
    
    if ismember(str2double(exp_number), exp_trials) && ...
       startsWith(tag{1}, 't') && strcmp(tag{2}, 'exp')&& ...
       strcmp(tag{4}, 'start')&& ismember(tag{3}, points)
            point = tag{3};
            extract = true;
    % If calibration and horizontal
    elseif ismember(str2double(exp_number), calib_trials) && ...
       startsWith(tag{1}, 't') && strcmp(tag{2}, 'calib')&& ...
       strcmp(tag{3}, 'H') && strcmp(tag{5}, 'start') && ismember(tag{4}, points)
            point = tag{4};
            extract = true;
    end
    
    % Extract single saccade
    if extract
        start = event_struct.time_stamps(ind);
        ping = event_struct.time_stamps(ind + 1);
        stop = event_struct.time_stamps(ind+2);
        
        eog_inds = opbci_struct.time_stamps >= start & opbci_struct.time_stamps <= stop;
        elink_inds = elink_struct.time_stamps >= start & elink_struct.time_stamps <= stop;
        eog_saccade = raw_eog(eog_inds);
        elink_saccade = elink(elink_inds);
        eog_time = opbci_struct.time_stamps(eog_inds);
        elink_time = elink_struct.time_stamps(elink_inds);
        
        cutoffs = findchangepts(elink_saccade, "MaxNumChanges", 2); % Isolate saccade using changepoints
        % If insufficient changepoints found
        if length(cutoffs) < 2 || cutoffs(2) - cutoffs(1) < 3
              cutoffs = [1 length(elink_saccade)];
        end
        [peaks, peak_inds] = get_peaks(elink_saccade, point, cutoffs);
        peak_avg = mean(peaks);
        peak_std = std(peaks);
        
        %Plot Elink
        figure();
        hold on;
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
        yyaxis left
        set(gca,'ycolor','g');
        ylabel("Raw EOG Signal: millivolts");
        plot(eog_time, eog_saccade, 'LineWidth', 1,'Color','g', 'DisplayName','Raw EOG'); 
        yyaxis right
        set(gca,'ycolor','b') 
        ylabel("Saccade Magnitude (degrees)");
        plot(elink_time, elink_saccade, 'DisplayName','EyeLink');
        % Plot cutoffs
        line1 = xline(elink_time(cutoffs(1)), "-.");
        line2 = xline(elink_time(cutoffs(2)), "-.");
        line1.Annotation.LegendInformation.IconDisplayStyle = 'off';
        line2.Annotation.LegendInformation.IconDisplayStyle = 'off';
        % Plot peaks
        plot(elink_time(peak_inds), peaks, '-o', 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Peaks');
        yline(peak_avg, 'DisplayName', 'Average Peak');
        yline(peak_avg - peak_std, "--", 'DisplayName', 'Std Peak');
        y2 = yline(peak_avg + peak_std, "--", 'DisplayName', 'Std Peak');
        y2.Annotation.LegendInformation.IconDisplayStyle = 'off';
        xlabel("Time: seconds");
        title(strcat("Participant ", participant,", Trial ", exp_number, " (", point, "): ", "EyeLink Signal"));
        legend("Location", "southeast");
        
        % Plot each filter
        for i = 1:numel(filters)
            filtered_data = filtered_eog.(filters{i});
            filtered_saccade = filtered_data(eog_inds);
            
            % Find cutoffs again
            cutoffs = findchangepts(filtered_saccade, "MaxNumChanges", 2); % Isolate saccade using changepoints
            if length(cutoffs) < 2 || cutoffs(2) - cutoffs(1) < 3
              cutoffs = [1 length(eog_saccade)];
            end
            [peaks, peak_inds] = get_peaks(filtered_saccade, point, cutoffs);
            peak_avg = mean(peaks);
            peak_std = std(peaks);
            
            figure()          
            hold on;
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
            % Plot Raw EOG
            yyaxis left
            set(gca,'ycolor','g');
            ylabel("Raw EOG Signal: millivolts");
            plot(eog_time, eog_saccade, 'LineWidth', 1,'Color','g', 'DisplayName','Raw EOG'); 
            % Plot filtered signal
            yyaxis right
            set(gca,'ycolor','b') 
            ylabel("Saccade Magnitude: degrees");
            plot(eog_time, filtered_saccade, 'LineWidth', 1.5,'Color', 'b','DisplayName', filternames.(filters{i}));
            % Plot Cutoffs
            y_s = ylim;
            line1 = xline(eog_time(cutoffs(1)), "-.");
            line2 = xline(eog_time(cutoffs(2)), "-.");
            line1.Annotation.LegendInformation.IconDisplayStyle = 'off';
            line2.Annotation.LegendInformation.IconDisplayStyle = 'off';
            % Plot Peaks
            plot(eog_time(peak_inds), filtered_saccade(peak_inds), '-o', 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Peaks');
            yline(peak_avg, 'DisplayName', 'Avg Peak');
            yline(peak_avg - peak_std, "--", 'DisplayName', 'Std Peak');
            y2 = yline(peak_avg + peak_std, "--", 'DisplayName', 'Std Peak');
            y2.Annotation.LegendInformation.IconDisplayStyle = 'off';
            xlabel("Time: seconds");
            title(strcat("Participant ", participant,", Trial ", exp_number, " (", point, "): ", filternames.(filters{i})));
            legend("Location", "southeast");
        end
    end
end

function [peaks, inds] = get_peaks(saccade, point, cutoffs)
    [peaks, inds] = findpeaks(abs(saccade));
    % Isolated region
    peaks = peaks(inds >= cutoffs(1) & inds <= cutoffs(2)); 
    inds = inds(inds >= cutoffs(1) & inds <= cutoffs(2));
    % If there are no peaks in isolated region, then take entire length
    if(isempty(peaks))
        [peaks, inds] = findpeaks(abs(saccade));
    end
    if(isempty(peaks)) % If still empty, take entire saccade
        peaks = abs(saccade);
        inds = 1:length(saccade);
    end
    if strcmp(point, 'A') || strcmp(point, 'B')
        peaks = -peaks;
    end
end

function [event_struct, opbci_struct, elink_struct] = get_raw_data(xdf_fname, xdf_path)
SRATE = 250;
ISI = 1/SRATE;
% try
% cd C:/Users/VisionLab/Desktop/SaccadeStudy_newLSL/
% catch
% end
% [xdf_fname, xdf_path] = uigetfile('..\Data\*.xdf', 'Select a File');

xdf = load_xdf([xdf_path, xdf_fname]);

% Initialize variables
event_idx = NaN;
elink_idx = NaN;
opbci_idx = NaN;
[~, strm_idx_count] = size(xdf);
event_strm_exist = false;
elink_strm_exist = false;
opbci_strm_exist = false;

% Search XDF struct to identify streams
for i=1:strm_idx_count
    if strcmp(xdf{1,i}.info.name, 'EventMarkers')
        event_idx = i;   event_strm_exist = true;
    elseif strcmp(xdf{1,i}.info.name, 'EyeLink')
        elink_idx = i;   elink_strm_exist = true;
    elseif strcmp(xdf{1,i}.info.name, 'OpenBCI_EOG')
        opbci_idx = i;   opbci_strm_exist = true;
    end
end

% Check for errors in XDF streams
if (event_strm_exist && isnan(event_idx)) || (elink_strm_exist &&...
        isnan(elink_idx)) || (opbci_strm_exist && isnan(opbci_idx))
    error('The data seems to be corrupted. Please check the ''xdf'' stucture to determine the problem.')
end

%%% Read EyeLink Data from XDF %%%
if elink_strm_exist
   elink = xdf{1,elink_idx};
   el_time_axis = elink.time_stamps;
   el_dat_x = elink.time_series(1,:);
   el_dat_y = elink.time_series(2,:);
   el_dat_x(el_dat_x < 0) = NaN;
   el_dat_x(el_dat_x > 1920) = NaN;
   el_dat_y(el_dat_y < 0) = NaN;
   el_dat_y(el_dat_y > 1080) = NaN;
elseif ~elink_strm_exist
   elink = NaN; el_time_axis = NaN; el_dat_x = NaN; el_dat_y = NaN;
   el_dat_x(el_dat_x < 0) = NaN; el_dat_x(el_dat_x > 1920) = NaN;
   el_dat_y(el_dat_y < 0) = NaN; el_dat_y(el_dat_y > 1080) = NaN;
end

%%% Read OpenBCI Data from XDF %%%

if opbci_strm_exist
opbci = xdf{1,opbci_idx};
ob_time_axis = opbci.time_stamps;
ob_dat_x = opbci.time_series(1,:);
ob_dat_y = opbci.time_series(2,:);
elseif ~opbci_strm_exist
opbci = NaN;     ob_time_axis = NaN;
ob_dat_x = NaN;  ob_dat_y = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% EOG Filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ob_dat_x = detrend(ob_dat_x);
 ob_dat_y = detrend(ob_dat_y);
%  ob_dat_x = bandpass(ob_dat_x,[0.6,35],250);
%  ob_dat_y = bandpass(ob_dat_y,[0.6,35],250);
     
%%% Instatiate Event Data %%%

if event_strm_exist
event = xdf{1,event_idx};
ev_ts = event.time_stamps;
ev_dat = event.time_series;
elseif ~event_strm_exist 
event = NaN;
ev_ts = NaN;
ev_dat = NaN;
end

%%% Plot Data %%%%

% Create time axis
t1 = min([el_time_axis(1), ob_time_axis(1), ev_ts(1)]);


el_time_axis = el_time_axis - t1;
ob_time_axis = ob_time_axis - t1;
ev_ts = ev_ts - t1;

% fig=figure;
% 
% %%% Create Subplot for Gaze X Data %%%
% subplot(2,1,1); 
% hold on;
% plot(el_time_axis, el_dat_x, 'Color', [0, 0.447, 0.741]);
% plot(ob_time_axis, detrend(ob_dat_x),'Color', [0.635 0.078 0.184]);
% title('Horizontal Gaze');
% xlabel('Time (secs)');
% ylabel('Gaze Position (pixels)');
% hold off
% legend({'EyeLink X', 'EOG X'}, 'AutoUpdate', 'off');
% 
% %%% Create Subplot for Gaze Y Data %%%
% subplot(2,1,2);
% hold on;
% plot(el_time_axis, el_dat_y, 'Color', [0.929, 0.694, 0.125]);
% plot(ob_time_axis, detrend(ob_dat_y), 'Color', [0.494, 0.184, 0.556]);
% title('Vertical Gaze')
% xlabel('Time (secs)')
% ylabel('Gaze Position (pixels)')
% hold off
% legend({'EyeLink Y', 'EOG Y'}, 'AutoUpdate', 'off');
% 
% for jj = 1:2
%     subplot(2,1,jj)
%     h = gca;
%     
% for ii=1:length(ev_ts)
%     event_tag = event.time_series{ii};
%     event_tag_split = strsplit(event_tag,'_');
%     if startsWith(event_tag_split{1}, 't')
%         if ~strcmp(event_tag_split{2}, 'cal') && strcmp(event_tag_split{2}, 'exp')
%             if strcmp(event_tag_split{4}, 'start')
%                 line([ev_ts(ii) ev_ts(ii)], h.YLim,'Color',[0.1 0.8 0.1])
%             elseif strcmp(event_tag_split{4}, 'ping')
%                 line([ev_ts(ii) ev_ts(ii)], h.YLim,'Color',[0 1 1])
%             elseif strcmp(event_tag_split{4}, 'end')
%                 line([ev_ts(ii) ev_ts(ii)], h.YLim,'Color',[1 0 0])
%             end
%         elseif strcmp(event_tag_split{2}, 'calib') % && strcmp(event_tag_split{2}, 'exp')
%             if strcmp(event_tag_split{length(event_tag_split)}, 'start')
%                 line([ev_ts(ii) ev_ts(ii)], h.YLim,'Color',[0.1 0.8 0.1])
%             elseif strcmp(event_tag_split{length(event_tag_split)}, 'ping')
%                 line([ev_ts(ii) ev_ts(ii)], h.YLim,'Color',[0 1 1])
%             elseif strcmp(event_tag_split{length(event_tag_split)}, 'end')
%                 line([ev_ts(ii) ev_ts(ii)], h.YLim, 'Color',[1 0 0])
%             end
%         end
%     end
% end
% end

event_struct = struct("time_stamps", ev_ts);
event_struct.time_series = ev_dat;
opbci_struct = struct("time_stamps", ob_time_axis, "time_series", ob_dat_x);
elink_struct = struct("time_stamps", el_time_axis, "time_series", el_dat_x-1000);
end