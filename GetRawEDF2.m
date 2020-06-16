%% Partition Raw EDF Data
% Raymond MacNeil
% UBC Vision Lab
% June 2020
function ExpTimeSeries = GetRawEDF2()
    
currentPath = pwd;
[FileNames, path] = uigetfile('*.edf','MultiSelect', 'on');
cd(path)

if ischar(FileNames)
    FileNames = cellstr(FileNames);
end

EdfCell = cell(numel(FileNames),1);
ExpTimeSeries = cell(40,2,numel(FileNames));

for ii = 1:length(EdfCell)
    EdfCell{ii,1} = Edf2Mat(FileNames{ii});
end
cd(currentPath);

for ii = 1:numel(FileNames)
    PINFO = strsplit(FileNames{1,ii},{'_', '.'});
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

    % Time Series
    EventIdxs =  ~cellfun('isempty', regexp(EdfCell{ii,1}.Events.Messages.info,...
        '^\[\''t.*\]$', 'match')) & ~contains(EdfCell{ii,1}.Events.Messages.info,...
        {'calib_V', 'calib_H'});
    EventMsgs = EdfCell{ii,1}.Events.Messages.info(EventIdxs)';
    EventTime = EdfCell{ii,1}.Events.Messages.time(EventIdxs)';
    EventMsgTimeVerify = [EventMsgs(1:2:end-1), num2cell(EventTime(1:2:end-1)),...
        EventMsgs(2:2:end), num2cell(EventTime(2:2:end))];
    EventTime = cell2mat(EventMsgTimeVerify(:,[2,4]));
    EvtTimeRngs = cell(numel(EventTime)/2,1);

    for jj = 1:size(EvtTimeRngs,1)
        EvtTimeRngs(jj) = {EventTime(jj,1):4:EventTime(jj,2)};
    end

    ExpTimeSeries(:,1,ii) = EvtTimeRngs;
    ADJ = [0, 1, -1, 2, -2];


    for kk = 1:size(ExpTimeSeries,1)

        jj = 0;

        while isempty(ExpTimeSeries{kk,2,ii})
            jj = jj+1;

            if jj > 5
                break
            end

            [~,GzXIdx] = ismember(ExpTimeSeries{kk,1,ii}+ADJ(jj),...
                EdfCell{ii,1}.Samples.time(:,1));
            ExpTimeSeries(kk,2,ii) = {EdfCell{ii,1}.Samples.gx(nonzeros(GzXIdx),2)};
        end
        
    end    
end

return