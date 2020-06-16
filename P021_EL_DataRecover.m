%% Load and prepare data
[file, path] = uigetfile('*.edf');
cd(path)
FileName = Edf2Mat(file);
cd '/Users/ray.macneil/Nextcloud/ubc_occulo/Analysis'
PINFO = strsplit(file,{'_', '.'});  
PID = PINFO{1};
ConStr = PINFO{2};

if contains(ConStr, '1')
    VisOrMem = 'Visually';
else
    VisOrMem = 'Memory';
end

if contains(ConStr, '3') || contains(ConStr, '4')
    lights = 'Off';
else
    lights = 'On';
end

% Horizontal Calibration Time Series
EventIdxs =  ~cellfun('isempty', regexp(FileName.Events.Messages.info,...
    '^\[\''t.*\]$', 'match')) & ~contains(FileName.Events.Messages.info,...
    'calib_V');
EventMsgs = FileName.Events.Messages.info(EventIdxs)';
EventTime = FileName.Events.Messages.time(EventIdxs)';
EventMsgTimeVerify = [EventMsgs(1:2:end-1), num2cell(EventTime(1:2:end-1)),...
    EventMsgs(2:2:end), num2cell(EventTime(2:2:end))];
EventTime = cell2mat(EventMsgTimeVerify(:,[2,4]));
EvtTimeRngs = cell(numel(EventTime)/2,1);

for ii = 1:size(EvtTimeRngs,1)
    EvtTimeRngs(ii) = {EventTime(ii,1):4:EventTime(ii,2)};
end

CalEvtTimeIndices = cellfun(@(x) contains(x, 'calib'), EventMsgTimeVerify(:,1));
CalTimeSeries = EvtTimeRngs(CalEvtTimeIndices);
ExpTimeSeries = EvtTimeRngs(~CalEvtTimeIndices);
CalTimeSeries = [CalTimeSeries, cell(size(CalTimeSeries,1),1)];
ExpTimeSeries = [ExpTimeSeries, cell(size(ExpTimeSeries,1),1)];

ADJ = [0, 1, -1, 2, -2];

for ii = 1:size(CalTimeSeries,1)
    
    jj = 0;
    
    while isempty(CalTimeSeries{ii,2})
        jj = jj+1;
        
        if jj > 5
            break
        end
        
        [~,GzXIdx] = ismember(CalTimeSeries{ii,1}+ADJ(jj),...
            FileName.Samples.time(:,1));
        CalTimeSeries(ii,2) = {FileName.Samples.gx(nonzeros(GzXIdx),2)};
    end
    
end

for ii = 1:size(ExpTimeSeries,1)
    
    jj = 0;
    
    while isempty(ExpTimeSeries{ii,2})
        jj = jj+1;
        
        if jj > 5
            break
        end
        
        [~,GzXIdx] = ismember(ExpTimeSeries{ii,1}+ADJ(jj),...
            FileName.Samples.time(:,1));
        ExpTimeSeries(ii,2) = {FileName.Samples.gx(nonzeros(GzXIdx),2)};
    end
    
end

MarkerExpIndices = cell(4,1); 
Markers = {'A', 'B', 'C', 'D'};

for ii = 1:numel(Markers)
  MarkerExpIndices(ii,1) = {find(cellfun(@(x) contains(x, strcat('exp_',...
      Markers{ii})), EventMsgTimeVerify(:,1)))};
end

smXoffset = 161.40;
lgXoffset = 326.40;
home = 960;
TargetLocations = [home-lgXoffset, home-smXoffset, home+smXoffset, home+lgXoffset];

%% Get Spaghetti Plots

for jj = 1:size(MarkerExpIndices,1)
    
    Indices = MarkerExpIndices{jj}-size(CalTimeSeries,1);
    fig = figure();
    hold on
    
    for ii = 1:numel(MarkerExpIndices{1})
        data = ExpTimeSeries{Indices(ii),2};
        plot(data)
    end
    
    xTitle = xlabel('Sample Count (250 Hz)', 'FontSize', 14, 'FontWeight', 'bold');
    yTitle = ylabel('Gaze Position (Pixels)', 'FontSize', 14, 'FontWeight', 'bold');
    myYline = yline(TargetLocations(jj), '--', 'Target Position', 'FontSize', 14,...
        'FontWeight', 'bold', 'LineWidth', 2.5, 'Color', [0,0,50/255])
    myTitle = title(sprintf('%s: %s-Guided Saccades | Lights %s | %s Trials',...
        PID, VisOrMem, lights, Markers{jj}), 'FontSize', 14); 
    
    print(gcf,strcat(PID, ConStr, Markers{jj}, '.png'),'-dpng','-r400');
    hold off
end
   
%% Get Revised Table
GzXIdx = nonzeros(GzXIdx);
AvgRows = 11:12:191;
SDevRows = 12:12:192;

for jj = 1:NumUniquePIDs
 ResultsLongForm(AvgRows(jj),:) = num2cell(nanmean(table2array(ResultsLongForm((...
     AvgRows(jj)-10):(AvgRows(jj)-1),:)),1));
 ResultsLongForm(SDevRows(jj),:) = num2cell(nanstd(table2array(ResultsLongForm((...
     SDevRows(jj)-11):(SDevRows(jj)-2),:)),1));
end

LongFormVarNames = ResultsLongForm.Properties.VariableNames;
Col2DelIdx = contains(LongFormVarNames,'MAX_AMP_X');
ResultsLongForm(:,Col2DelIdx) = [];
Cols2Convert = contains(MyVarNames, 'MAX_GZX');
A = table2array(ResultsLongForm);

 