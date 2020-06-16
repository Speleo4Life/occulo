%% MODIFIED Closed-Eyes EyeLink Data Analysis FOR GETTING BLINK DATA
% Raymond MacNeil, 2020-June
% Vision Lab, Department of Psychology, University of British Columbia
% Updated June 15th, 2020
% GitHib Repo: https://github.com/Speleo4Life/occulo2.git

% Dependencies: varmaker2.m

%% DATA IMPORT, CONSTANTS, AND PREALLOCATION OF RESULTS MATRICES/CELLS/TABLES/ETC 

%----------------------------------------------------------------------------------%
%                        IMPORT DATA VIA GUI AND STORE IN STRUCT                   %
%----------------------------------------------------------------------------------%                            


What2Do = char;

while ~ismember(What2Do, {'Y', 'N'})
    What2Do = input('Load in preconfigured EventData4LeoHiro? Y/N [N]:','s');
    
    if isempty(What2Do)
        Wnat2Do = 'N'
        break
    end
end
        
if strcmp(What2Do, 'N')

    %clearvars
    currentPath = pwd;
    % Get .EDF files to be analyzed
    [fileIDs, fileDir]=uigetfile('MultiSelect','On','*.edf');
    fileIDs = transpose(sort(string(fileIDs)));
    N = numel(unique(extractBefore(fileIDs, '_'))); 


    % Import files batchwise and store in struct array
    cd(fileDir);
    for ii = 1:length(fileIDs)

       try
           subj = extractBefore(fileIDs(ii), '_');
           cond = extractBetween(fileIDs(ii), '_', '.edf');
           edfm = Edf2Mat(char(fileIDs(ii)));
           edfevt.(subj).(cond) = edfm.RawEdf.FEVENT;
       catch
           warning('There was a problem with file: %s', fileIDs(ii));
       end

    end
    cd(currentPath)

elseif strcmp(What2Do, 'Y')
    
 % Don't no the consequence of excluding these will have (normally stored when
 % loading files in via the GUI), but will  add in case to prevent breaking.
 fileIDs = ["P019_C1.edf"; "P019_C2.edf"; "P019_C3.edf"; "P020_C1.edf";...
     "P020_C2.edf"; "P020_C3.edf"; "P022_C1.edf"; "P022_C2.edf"; "P022_C3.edf";...
     "P023_C1.edf"; "P023_C2.edf"; "P023_C3.edf"; "P024_C1.edf"; "P024_C2.edf";...
     "P024_C3.edf"; "P025_C1.edf"; "P025_C2.edf"; "P025_C3.edf"; "P026_C1.edf";...
     "P026_C2.edf"; "P026_C3.edf"; "P027_C1.edf"; "P027_C2.edf"; "P027_C3.edf";...
     "P028_C1.edf"; "P028_C2.edf"; "P028_C3.edf"; "P029_C1.edf"; "P029_C2.edf";...
     "P029_C3.edf"; "P030_C1.edf"; "P030_C2.edf"; "P030_C3.edf"; "P032_C1.edf";...
     "P032_C2.edf"; "P032_C3.edf"; "P033_C1.edf"; "P033_C2.edf"; "P033_C3.edf";...
     "P034_C1.edf"; "P034_C2.edf"; "P034_C3.edf"; "P035_C1.edf"; "P035_C2.edf";...
     "P035_C3.edf"]
     N = numel(unique(extractBefore(fileIDs, '_')));     
    [fName, pName]=uigetfile('*.mat', 'Load the EventData File!');
    
    try
        load(fullfile(pName, fName))
    catch
        error('Hmmm...there seems to be a problem. Please try again.')
    end
    
end
%%
% Get Calibration Data For each Subject and Condition        %
T = LocalGetCalDat(N, edfevt, fileIDs);
pidFLDnames = char(fieldnames(edfevt));


%----------------------------------------------------------------------------------%
%%           Declare Constants and Preallocate Result Matrices/Cells                %
%----------------------------------------------------------------------------------%

letters = {'A'; 'B'; 'C'; 'D'}; 
seq = {'start'; 'end'};
numSaccParam = 15; % Number of Saccade Parameters
numTrialReps = 10; % Number of Trial Repetitions for each Marker Position
numTrialType = 4; % Number of Trial Types / Marker Positions
numConds = 3; % Number of Conditons Relevant for EyeLink Analysis
numTrialTot = numTrialReps .* numTrialType .* numConds; % Total Number of Trials
ResultsMTRX1 = zeros(N .* numSaccParam, numTrialTot);
ResultsAverages = zeros(N, numSaccParam .* (numTrialType .* numConds));
ResultsMTRX3 = zeros(N,(numTrialType .* numConds));
numColsTrack = 14;
TrackingMatrix = num2cell(NaN(numTrialReps .* numTrialType .* numConds,...
    numColsTrack, N));
TrackingMatrix = LocalFillTrackingMatrix(fileIDs, TrackingMatrix, numTrialReps);
LetterTrack = repmat(string(letters)',1,3);
% Create an array with the 'conversion' coefficients corresponding 
% to each trial type across the three conditions. 
CFT = repmat([1, 1, -1, -1],1,3); 

% Proximal marker (11°) center offset in PIXELS calculated using PsychoPy FUN
smX_eccentric_offset = 161.40665321197295; % 

% Distal marker (22°) center offset in PIXELS calculated using PsychoPy FUN
lgX_eccentric_offset = 326.4192270718748;  

% Degrees visual angle per pixel at a viewing distance of 53-cm
% This value is definitely accurate for our data/setup. 
% Calculated using PsychoPy FUN
pix2deg = 0.06537277453582728; 

% Homebase marker is at 960 pixels along x-axis
xHOME = 960; 

% CAUTION: THIS MAY HAVE VARIED FOR SUBJECTS RUN IN MARCH
% IT WILL BE NECESSARY TO FURTHER INVESTIGATE
yHOME = 588; % Homebase marker is at 588 pixels along y-axis for Ss P019 to P031

% Define marker positions (in pixels) given the above constants
MkrPos = repmat([xHOME-lgX_eccentric_offset, xHOME-smX_eccentric_offset,...
    xHOME+smX_eccentric_offset, xHOME+lgX_eccentric_offset],1,3);

%----------------------------------------------------------------------------------%
%%                              ISOLATE/SUBSET TRIALS                               %
%----------------------------------------------------------------------------------%

% Get trial START and END timestamp indices 
evt_indices = LocalGetTrialStartEndIndices(N, letters, seq, pidFLDnames, edfevt);

% Using the above info, determine the index range for each trial's timestamp range 
% Intialize and preallocate a table for storing the data
var_names = {'C1_A_IDX', 'C1_B_IDX', 'C1_C_IDX', 'C1_D_IDX',...
    'C2_A_IDX', 'C2_B_IDX', 'C2_C_IDX', 'C2_D_IDX',...
    'C3_A_IDX', 'C3_B_IDX', 'C3_C_IDX', 'C3_D_IDX'};
var_types = cellstr(repmat("cell", 1, numel(var_names)));
TEVT = table('Size', [N, numel(var_names)], 'VariableTypes', var_types,...
    'VariableNames', var_names);

% Extract the information using the data in evt_indices table
for ii = 1:N
    d = 1;
    for jj = 1:2:width(evt_indices)-1
        holder1 = table2array(evt_indices(ii,jj));
        holder2 = table2array(evt_indices(ii,jj+1));
        loopsz = numel(holder1{:});
        temp = cell(loopsz,1);
        
        for kk = 1:loopsz
            temp{kk}= holder1{:}(kk):holder2{:}(kk);  
        end
        TEVT{ii,d} = {temp};
        d=d+1;
    end
end

edfevt = LocalGetAmps(edfevt, pidFLDnames, letters, N);

%----------------------------------------------------------------------------------%
%%          DECLARE CONSTANTS FOR GAZE SWEEP PARAMETER EXTRACTION ALGORITHM         %
%----------------------------------------------------------------------------------%

% Create an anonymous function that shows corresponding indices of the 
% ENDSACC and ENDFIX gaze positions. For use when the algorithm fails to 
% determine saccade start and return index values and the user opts for
% manual specification.
FHM = @(x) table(transpose(1:height(x)), 'VariableNames', {'Index'});  

% Create an array of anonymous functions to dynmaically get min or max
% amplitude depending on trial type. NB. This functionality of anonymous functions 
% will likely be of particular value elsewhere in the future. 

FH1 = @(x) min(x);
FH2 = @(x) max(x);
FHT = repmat(repelem({FH1; FH2},2,1),numConds,1);

% Threshold parameters for saccade start and return search algorithm
SEARCH_THRESH_PIX = -25;
CHNGTOL = 100;
DIFF_TOLERANCE = -1400;
DIFF_ST_THRESH_PIX = -76;
limxUp = 2500;

%% ALGORITHM EXTRACTING GAZE PARAMETERS

% THE LOGICAL HIERARCHY OF THE NESTED LOOPS IS AS FOLLOWS: 
% ii: Looping at the subject-level
% jj: Looping at the condition- AND trial-type level
%     It follows a C1A, C1B, C1C, C1D, C2A, ..., C3D cycle
%     of twelve for each ii (subject). 
% kk: Looping at the trial level. There are ten trials 
%     per trial type ('A', 'B', 'C', and 'D' i.e. marker position).

for ii = 1:N
    FLname = pidFLDnames(ii,:);
    fprintf('Processing Participant %s \n', FLname)
    condFLDnames = fieldnames(edfevt.(pidFLDnames(ii,:)));
    CondLoop = repelem(string(condFLDnames), numel(letters));
    
    for jj = 1:numel(CondLoop)
        
        % Trigger update in data source depending on condition when jj is
        % equal to 1-4, it is C1, jj == 5:8 is C2, and jj == 9:12 is C3
        if (jj == 1) || CondLoop{jj}(2) ~= CondLoop{jj-1}(2)
            get = struct2table(edfevt.(FLname).(CondLoop{jj}));
        end
        
        idx = TEVT{ii,jj};
        SaccParams = zeros(numSaccParam,numel(idx{:}));
        CurrentCondNum = ceil(jj/numel(letters));

        for kk = 1:numel(idx{:})
            dud_trial = false;
            AnticipatoryStartRightDir = false;
            AnticipatoryStartWrongDir = false;
            evtidx = cell2mat(idx{:}(kk));
            evtsubTb = get(evtidx,:);
            
            %% MODIFIED THE SUBSET INDEX TO CONTAIN MORE EVENT INFO 
             % MODIFY AS YOU SEE FIT
            subset = strcmp(evtsubTb.codestring, 'STARTSACC') |...
                strcmp(evtsubTb.codestring, 'ENDSACC') |...
                strcmp(evtsubTb.codestring, 'STARTFIX') |...
                strcmp(evtsubTb.codestring, 'ENDFIX') |...  
                strcmp(evtsubTb.codestring, 'STARTBLINK ') |...
                strcmp(evtsubTb.codestring, 'ENDBLINK') |...
                evtsubTb.parsedby==48;
                evtsubT = evtsubTb(subset,:); 
       
%% MODIFIED SECTION THAT SKIPS EVENT PROCESSING & JUST STORES THE TRIAL EVENT TABLE 
            
            % Get the event info restricted to the current subject,
            % condition, and trial. It is the evtidx variable that is key
            % in allowing us to index into the relevant portion of the 
            TrialTbl = evtsubTb(subset,:);
            
            % Get the proper row index for indexing into the TrackingMatrix
            % Note the hiearchy described above, we're always keeping track
            % of the jj, kk, and ii indices. By substracting one from jj,
            % which is tracking the condition number, we can multiply by ten
            % (i.e. the number of trials per Target) and add the KKth trial
            % index to determine the proper row index. We simplify life by
            % distributing the subject indices (ii) across 'pages' (dimension 3)
            % of the cell matrix. We reshape later, however, for ease of
            % navigating the data structure.
           
            M1Idx = (jj-1) .* numTrialReps + kk + 1;
            TrackingMatrix(M1Idx, 12, ii) = {TrialTbl};
            dud_trial = true;
            SaccParams(1:end, kk) = NaN;
            
            
            continue
%% WHAT NORMALLY HAPPENS            
            
            inf_trial = evtsubT.gstx > limxUp; %#ok<*UNRCH>
            evtsubT(inf_trial, :) = [];      
            inf_end = find(evtsubT.genx > limxUp);
            
            if ~isempty(inf_end) && (evtsubT.gstx(inf_end) <= limxUp) 
                evtsubT.genx(inf_end) = evtsubT.gstx(inf_end);
            end
            WSevtsubT = evtsubT;            
            evtsubT(:,[5,6]) = num2cell(CFT(jj) .* table2array(evtsubT(:,[5,6])));
            CFT_SHIFT = circshift(CFT, 2);
            WSevtsubT(:,[5,6]) = num2cell(CFT_SHIFT(jj) .*... 
                table2array(WSevtsubT(:,[5,6])));
            Y = 0;
            GzStart = [];
            GzReturn = [];
            GzStThresh = false;
            GzRtThresh = false;
          TrialTbl = [FHM(evtsubT), evtsubT(:,[1,4,5,6])];

            
            Implement algorithm to obtain the gaze start index.
    
            while isempty(GzStart) || isnan(GzStart)
                
                % If iterative threshold reduction yields no result prompt
                % option for manual specification. 
                if (DIFF_ST_THRESH_PIX + Y) == SEARCH_THRESH_PIX
                    GzStThresh = true;
                    StartEndString = 'START';
                    [GzStart, NOgzst, dud_trial] = LocalManIdxSpec(FLname,...
                        CurrentCondNum, LetterTrack, jj, kk, TrialTbl,...
                        StartEndString)  

                    break
                end
                
                
                
% The algorithm for determining when gaze departs from home base works as
% follows: The events of a trial are first subsetted to ENDSACC and ENDFIX events.
% Then, it identifies indices of the events classfied as saccades that have
% (1) a horizontal gaze end position of 125 pixels (if necesssary, iteratively
% decreasing to 25 pixels at minimum) less ('A' and 'B' trials) or greater
% ('C' and 'D' trials) than the previous event; (2) the absolute difference
% is no greater than 1400 pixels; and (3) the gaze position is less or
% greater (AB and CD, respectively) than 100 pixels greater than home base
% (1060; 'A' and 'B' trials) or 100 pixels less than home base (860; 'C'
% and 'D' trials).
               
                GzStartAll = find([0; diff(evtsubT.genx) < (DIFF_ST_THRESH_PIX + Y)]...
                    & [0; diff(evtsubT.genx) > DIFF_TOLERANCE] & (evtsubT.genx...
                    < (CFT(jj) * (xHOME + (CFT(jj)* CHNGTOL)))) &...
                    strcmp(evtsubT.codestring, 'ENDSACC'));
                
%                  GzStartAll = find([0; diff(test.genx) < (DIFF_ST_THRESH_PIX + Y)]...
%                     & [0; diff(test.genx) > DIFF_TOLERANCE] & (test.genx...
%                     < (CFT(jj) * (xHOME + (CFT(jj)* CHNGTOL)))) &...
%                     strcmp(evtsubT.codestring, 'ENDSACC'));
                
                if ~isempty(GzStartAll)
                    GzStart = GzStartAll(1);
                end
                              
                Y = Y + 1;
            end
            
            M1Idx = (jj-1) .* numTrialReps + kk + 1;
            
            if GzStThresh
                
                TrackingMatrix(M1Idx, 3, ii) = {1};
                TrackingMatrix(M1Idx, 4, ii) = {NOgzst};
                TrackingMatrix(M1Idx, 5, ii) = {dud_trial};
                
                if ~dud_trial
                    TrackingMatrix(M1Idx, 6, ii) = {GzStart};
                end
                            
            else
                
                TrackingMatrix(M1Idx, 3, ii) = {0};
                TrackingMatrix(M1Idx, 5, ii) = {double(dud_trial)};
                TrackingMatrix(M1Idx, 6, ii) = {GzStart};   
            end
            
% Implement algorithm to obtain gaze return index
% This is very similar to the above loop. 
           
            Y = 0;
            
            while isempty(GzReturn) || isnan(GzReturn)
                
              
                
                
                if (DIFF_ST_THRESH_PIX + Y) == SEARCH_THRESH_PIX
      
                    GzRtThresh = true;
                    
                    try
                        
                        if ~isempty(evtsubT.genx) && ((evtsubT.genx(end) -...
                                min(evtsubT.genx)) < (std(evtsubT.genx) * 0.30))
                            NOgzre = 3;
                            GzReturn = length(evtsubT.genx);
                            break
                            
                        else
                            
                            StartEndString = 'RETURN';
                            [GzReturn, NOgzre, dud_trial] = LocalManIdxSpec(FLname,...
                                CurrentCondNum, LetterTrack, jj, kk, TrialTbl,...
                                StartEndString);
                            break
                        end
                    
                    catch 
                               
                       StartEndString = 'RETURN';
                        [GzReturn, NOgzre, dud_trial] = LocalManIdxSpec(FLname,...
                            CurrentCondNum, LetterTrack, jj, kk, TrialTbl,...
                            StartEndString);
                    end
                         
                end        
                             
                GzReturnAll = find([0; diff(evtsubT.genx(GzStart:end))...
                    > (abs(DIFF_ST_THRESH_PIX) - Y)] &...
                    [0; diff(evtsubT.genx(GzStart:end))...
                    < abs(DIFF_TOLERANCE)] &...
                    evtsubT.genx(GzStart:end) > (CFT(jj)*(MkrPos(jj) +...
                    (CFT(jj)*(100-Y)))) &...
                    strcmp(evtsubT.codestring(GzStart:end),...
                    'ENDSACC')) + (GzStart-1);
                if ~isempty(GzReturnAll)
                    GzReturn = GzReturnAll(1);
                end
                
           
                Y = Y + 1;  
            end
            
            if ~dud_trial && exist('GzStartAll', 'var')
            
                if ~issorted([GzStartAll; GzReturnAll])
                    [GzStart, GzReturn, DidResetStartReturn] = LocalGzMultiSweep(FLname,...
                        CurrentCondNum, LetterTrack, jj, kk, TrialTbl, GzStart,...
                        GzReturn);
                    
                    if DidResetStartReturn
                        AnticipatoryStartRightDir = true;
                    end
                    
                end
            
            end
                        

            if ~dud_trial && ~isempty(GzStart)
                Y = 0;    
                [~, AnticipatoryStartWrongDir] = DetectWrongStart(WSevtsubT,...
                    CFT_SHIFT, FLname, CurrentCondNum, LetterTrack, Y, jj, kk, TrialTbl,...
                    DIFF_TOLERANCE, xHOME, CHNGTOL, MkrPos, GzStart, GzReturn);
                
                if ~issorted([GzStart; GzReturn])
                    [GzStart, GzReturn] = LocalIndexAnomaly(FLname, CurrentCondNum,...
                        LetterTrack, jj, kk, TrialTbl, GzStart, GzReturn);
                end
                
            end
            
                           
            if GzRtThresh
                
                TrackingMatrix(M1Idx,  8, ii) = {1};
                TrackingMatrix(M1Idx,  9, ii) = {NOgzre};
                TrackingMatrix(M1Idx, 10, ii) = {double(dud_trial)};
                
                if ~dud_trial
                    TrackingMatrix(M1Idx, 11, ii) = {GzReturn};
                end
                
                TrackingMatrix(M1Idx, 12, ii) = {evt};
                TrackingMatrix(M1Idx, 13, ii) = {AnticipatoryStartRightDir};
                TrackingMatrix(M1Idx, 14, ii) = {AnticipatoryStartWrongDir};
                
            else
                TrackingMatrix(M1Idx,  8, ii) = {0};
                TrackingMatrix(M1Idx, 10, ii) = {double(dud_trial)};
                TrackingMatrix(M1Idx, 11, ii) = {GzReturn};
                TrackingMatrix(M1Idx, 12, ii) = {TrialTbl};
                TrackingMatrix(M1Idx, 13, ii) = {AnticipatoryStartRightDir};
                TrackingMatrix(M1Idx, 14, ii) = {AnticipatoryStartWrongDir};
            end
                         
            
            if isempty(GzStart) || isempty(GzReturn) || isnan(GzStart)...
                    || isnan(GzReturn)
                dud_trial = true;
            end
            
                    
            if ~dud_trial
                
                % Get saccadic latency
                if ~isempty(GzStartAll)
                    try
                        latency = evtsubT.sttime(GzStartAll(1))-evtsubT.sttime(1);
                    catch ME
                        warning(['Could not process latency for this trial.'...
                            ' Seems as if the CR and Pupil was lost.'])
                    end
                end
                
                % Get minnimum (really: maximum) X- and Y-gaze coordinates
                [min_gzx, Imng] = min(evtsubT.genx(1:end));
                min_gzx = CFT(jj) .* min_gzx;
                min_gzy = evtsubT.geny(Imng);
                
                % Get index range for finding the saccade end-point
                bound4fix = GzStart+1:GzReturn-1;
                [max_fix, land_idx] = max(evtsubT.evtdur(bound4fix));
                land_idx = bound4fix(land_idx);
                
                if ~isempty(max_fix)
                    land_drift_ampx = evtsubT.ampX(land_idx);
                    land_drift_ampxy = evtsubT.ampXY(land_idx);
                    land_drift_avel = evtsubT.avel(land_idx);
                end
                
                       
                try
                    pix_land_gzx = evtsubT.gavx(land_idx);
                    pix_land_gzy = evtsubT.gavy(land_idx);
                catch
                    pix_land_gzx = NaN;
                    pix_land_gzy = NaN;
                    dud_trial = true;
                end
                           
                peak_vel = max(evtsubT.pvel(GzStart:GzReturn-1));
                num_sacc = sum(strcmp('ENDSACC', evtsubT.codestring(GzStart:GzReturn-1)));                
                [max_fix, min_gzx, bound4fix, peak_vel, land_idx, pix_land_gzx,... 
                    pix_land_gzy,latency] = LocalCheckEmpty(max_fix, min_gzx,...
                    bound4fix, peak_vel, land_idx, pix_land_gzx, pix_land_gzy,...
                    latency); 
                             
                try
                    
                    % Get distance of fixation landing point from home base 
                    % in degrees of visual angle:

                    deg_land_gzx = (pix_land_gzx - xHOME) * pix2deg;
                    deg_land_gzy = (pix_land_gzy - yHOME) * pix2deg;
                    
                    % Calculate accuracy, the absolute difference between
                    % the landing point and the marker position WRT to both
                    % X- and Y-gaze coordinates in pixels and degrees of
                    % visual angle:

                    xyPixAcc = sqrt((pix_land_gzx - (MkrPos(jj))).^2 +...
                        (pix_land_gzy-yHOME).^2);
                    xyDegAcc = xyPixAcc * pix2deg;
                    
                    % Perform this calculation once more, but only WRT to
                    % the X-gaze coordinates at the fixation landing point: 
                
                    xPixAcc = sqrt((pix_land_gzx - (MkrPos(jj))).^2);
                    xDegAcc = xPixAcc * pix2deg;
                    
                catch
                    warning(['Something went wrong with accuracy and landing'... 
                        'point offset']);
                    xyPixAcc = NaN;
                    xyDegAcc = NaN;
                    xPixAcc = NaN;
                    xDegAcc = NaN;
                    deg_land_gzx = NaN;
                    deg_land_gzy = NaN;
                    dud_trial = true;
                end
            end
            
            try
                
                % Collect trial-level saccade metrics for batch input into
                % the results matrix 
                SaccParams(1,kk) = min_gzx;        SaccParams(2,kk) = num_sacc;
                SaccParams(3,kk) = land_drift_ampx; SaccParams(4,kk) = land_drift_ampxy;
                SaccParams(5,kk) = land_drift_avel; SaccParams(6,kk) = peak_vel;
                SaccParams(7,kk) = pix_land_gzx;   SaccParams(8,kk) = pix_land_gzy; 
                SaccParams(9,kk) = deg_land_gzx;   SaccParams(10,kk) = deg_land_gzy;        
                SaccParams(11,kk) = xyDegAcc;      SaccParams(12,kk) = xDegAcc;         
                SaccParams(13,kk) = latency;       SaccParams(14,kk) = AnticipatoryStartRightDir; 
                SaccParams(15,kk) = AnticipatoryStartWrongDir; 
            catch               
                warning('There was a problem parsing this trial...skipping to next.')
                dud_trial = true;              
            end
            
                    
            if dud_trial
                SaccParams(1:end, kk) = NaN;
            end
            
        end
        
        Midx = ((ii * numSaccParam) - numSaccParam + 1):(ii * numSaccParam);
        Nidx = ((jj * numTrialReps) - numTrialReps + 1):(jj * numTrialReps);
        ResultsMTRX1(Midx, Nidx) = SaccParams;
        Midx2 = ii;
        Nidx2 = ((jj * numSaccParam) - numSaccParam + 1):(jj * numSaccParam);
        ResultsAverages(Midx2, Nidx2) = transpose(nanmean(SaccParams,2));
        ResultsMTRX3(ii,jj) = sum(isnan(SaccParams(1,:)));
    end
end
    

%% FINALIZE DATA OUTPUT
%----------------------------------------------------------------------------------%
%                             RESULTS TABLE LONG FORMAT                            %
%----------------------------------------------------------------------------------%

rshpColIdxs = reshape(1:1:numTrialTot, numTrialReps, numTrialType * numConds)'; %#ok<*NOPTS,BDSCI>
rshpRowIdxs = reshape(1:1:(numSaccParam * N), numSaccParam, N)';
numColsLongForm = numSaccParam * numTrialType * numConds;
numRowsLongForm = numTrialReps * N;
ResultsLongForm = zeros(numRowsLongForm, numColsLongForm);

for jj = 1:N

    for kk = 1:size(rshpColIdxs,1)
        Midx = rshpColIdxs(1,:) + (numTrialReps * (jj-1));
        Nidx = rshpRowIdxs(1,:) + (numSaccParam * (kk-1));
        ResultsLongForm(Midx,Nidx) = ResultsMTRX1(rshpRowIdxs(jj,:),...
            rshpColIdxs(kk,:))';
    end
    
end

% Get Mean and Standard Deviation without Having to Overhaul Code

NumStats = 2; 
NumRows4HoldingCell = (N .* NumStats) + size(ResultsLongForm,1);
HoldingCell = cell(NumRows4HoldingCell, size(ResultsLongForm,2));
RowCountPerPID = numTrialReps + NumStats;
StartIDX = 1:RowCountPerPID:(NumRows4HoldingCell + 1 - RowCountPerPID);
IndexHolder = NaN(N, numTrialReps);
    for ii = 1:numel(StartIDX)
        IndexHolder(ii,:) = StartIDX(ii):1:(StartIDX(ii)+(numTrialReps-1));
    end
    
IndexHolder = reshape(IndexHolder', [], 1);
HoldingCell(IndexHolder,:) = num2cell(ResultsLongForm);
StatIdx = find(cellfun('isempty', HoldingCell(:,1)));

    for jj = 1:2:(numel(StatIdx)-1)
       Start = StatIdx(jj)-numTrialReps;
       End = StatIdx(jj)-1;
       HoldingCell(StatIdx(jj),:) = num2cell(nanmean(cell2mat(HoldingCell(Start:End,...
           :)),1));
       HoldingCell(StatIdx(jj+1),:) = num2cell(nanstd(cell2mat(HoldingCell(Start:End,...
           :))));
    end
    
                
VarsLongFormat = varmaker2(["C1", "C2", "C3"],["A", "B", "C", "D"],'false','false',...
    [], ["MAX_GZX", "NUM_SACC", "LAND_DRIFT_AMP_X", "LAND_DRIFT_AMP_XY",...
    "LAND_DRIFT_AVEL", "SACC_PVEL", "LAND_GZX_PIX", "LAND_GZY_PIX", "LAND_GZX_DEG",...
    "LAND_GZY_DEG", "ACC_XY_DEG", "ACC_X_DEG", "LATENCY", "ANT_ST_RIGHT",...
    "ANT_ST_WRONG"], 'true', [])';

PIDRowLabels = reshape(unique(string(regexp(fileIDs, 'P\d*', 'match'))), 1, []);
PIDRowLabels = varmaker2(PIDRowLabels,[], 'false', 'true', 12, [], 'false', []);
PIDRowLabels = strrep(PIDRowLabels, '11', 'AVG');
PIDRowLabels = strrep(PIDRowLabels, '12', 'STDEV');
ResultsLongForm = cell2table(HoldingCell,'VariableNames', VarsLongFormat,...
    'RowNames', PIDRowLabels);


%----------------------------------------------------------------------------------%
%                              RESULTS TABLE WIDE FORMAT                           %
%----------------------------------------------------------------------------------%

ResultsWideForm = NaN(N,numSaccParam * numTrialTot);
VarsWideFormat = varmaker2(["C1", "C2", "C3"],["A", "B", "C", "D"],'false','true',... 
    10, ["MAX_GZX", "NUM_SACC", "LAND_DRIFT_AMP_X", "LAND_DRIFT_AMP_XY",...
    "LAND_DRIFT_AVEL", "SACC_PVEL", "LAND_GZX_PIX", "LAND_GZY_PIX", "LAND_GZX_DEG",...
    "LAND_GZY_DEG", "ACC_XY_DEG", "ACC_X_DEG", "LATENCY", "ANT_ST_RIGHT",...
    "ANT_ST_WRONG"], 'true', [])';

for jj = 1:N
    RowIdx = (1:numSaccParam) + (numSaccParam * (jj-1));
    ResultsWideForm(jj,:) = reshape(ResultsMTRX1(RowIdx,:),1,[]);
end

ResultsWideForm = [T, array2table(ResultsWideForm, 'VariableNames', VarsWideFormat)];

%----------------------------------------------------------------------------------%
%%                              Reshape Tracking Matrix                             %
%----------------------------------------------------------------------------------%
TrackingMatrix = reshape(permute(TrackingMatrix, [1,3,2]), [],...
    size(TrackingMatrix,2)); 
ColsToKeep = ~strcmp('GzSt_TrackTbl', TrackingMatrix(1,:));
TrackingMatrix = TrackingMatrix(:,ColsToKeep);

%% LOCAL FUNCTIONS 

function T = LocalGetCalDat(N, edfevt, fileIDs)

% Preallocate a table for storing the data
T = table('Size', [N 4], 'VariableTypes',...
     {'cellstr', 'string', 'string', 'string'}, 'VariableNames',...
     {'PID', 'C1_CAL', 'C2_CAL', 'C3_CAL'});
T.Properties.Description = 'Closed-Eyes EyeLink Data | MacNeil | Vision Lab | 2020';
pidFLDnames = char(fieldnames(edfevt));
T.PID = categorical(unique(extractBefore(fileIDs, '_')));

% Extract calibration information and store it in our table
    for ii = 1:N
        FLname = pidFLDnames(ii,:);
        condFLDnames = fieldnames(edfevt.(pidFLDnames(ii,:)));

        for jj = 1:numel(condFLDnames)
            idx = find(([edfevt.(FLname).(condFLDnames{jj}).parsedby] == 188),1, 'last');
            if ~isempty(idx)
                T(ii,1+jj) = cellstr(edfevt.(FLname).(condFLDnames{jj})(idx).message);
            end
        end
    end
end

function edfevt = LocalGetAmps(edfevt, pidFLDnames, letters, N)  

    for ii = 1:N
        
        FLname = pidFLDnames(ii,:);
        condFLDnames = fieldnames(edfevt.(pidFLDnames(ii,:)));
        CondLoop = repelem(string(condFLDnames), numel(letters));

        for jj = 1:numel(CondLoop)

            if (jj == 1) || CondLoop{jj}(2) ~= CondLoop{jj-1}(2) 

                get = struct2table(edfevt.(FLname).(CondLoop{jj}));
                ampX = NaN(height(get),1);
                ampY = NaN(height(get),1);
                ampXY = NaN(height(get),1);

                for kk = 1:height(get)

                    if strcmp('ENDSACC', get.codestring(kk)) || strcmp('ENDFIX',...
                            get.codestring(kk))
                        ampX(kk) = (get.genx(kk) - get.gstx(kk))...
                            / ((get.eupd_x(kk) + get.supd_x(kk))/2);
                        ampY(kk) = (get.geny(kk) - get.gsty(kk))...
                            / ((get.eupd_y(kk) + get.supd_y(kk))/2);
                        ampXY(kk) = sqrt(ampX(kk).^2 + ampY(kk).^2);
                    end

                end

                get = addvars(get, ampX, 'After', 'genx');
                get = addvars(get, ampY, ampXY, 'After', 'eupd_y');
                evtdur = get.entime - get.sttime;
                get = addvars(get, evtdur, 'After', 'entime');
                get = removevars(get,{'hstx', 'hsty', 'henx', 'heny', 'sta',...
                   'ena', 'havx', 'havy', 'ava', 'status', 'flags', 'input',... 
                   'buttons', 'time', 'type', 'read', 'eye'});
                get = movevars(get, {'genx', 'gavx', 'ampX', 'ampXY'},... 
                    'After', 'gstx');
                get = movevars(get, 'ampY', 'After', 'gavy');
                get = movevars(get, 'codestring', 'Before', 'sttime');
                edfevt.(FLname).(CondLoop{jj}) = table2struct(get);
            
            end
            
        end
        
    end
    
end

function evt_indices = LocalGetTrialStartEndIndices(N, letters, seq, pidFLDnames,...
    edfevt)

var_names = varmaker2(["C1", "C2", "C3"],["A", "B", "C", "D"], 'false',...
    'false', [], ["ST", "EN"], 'true', "IDX");
var_types = cellstr(repmat("cell",1,numel(var_names))); 
evt_indices = table('Size', [N numel(var_names)],... 
    'VariableTypes', var_types, 'VariableNames', var_names); 

for ii = 1:N
    
    FLname = pidFLDnames(ii,:);
    condFLDnames = fieldnames(edfevt.(FLname));
    
    for jj = 1:numel(condFLDnames)
        get = struct2table(edfevt.(FLname).(condFLDnames{jj}));
        check = cellfun(@(x) isempty(x), get.message);
        get.message(check) = {""};
        extract = string(get.message);
        x = 8*jj-8;
        
        for kk = 1:numel(letters)
            
            for ll = 1:numel(seq)
                
                if ll == 1 
                    evt_indices(ii,x+(kk*2-1)) = {find(contains(extract,...
                        '_exp_') & contains(extract, letters(kk))...
                        & contains(extract, seq(ll)))}';
         
                elseif ll == 2
                     evt_indices(ii,x+(kk+kk)) = {find(contains(extract,...
                        '_exp_') & contains(extract, letters(kk))...
                        & contains(extract, seq(ll)))};        
                end
            end
        end  
    end        
end

end

function [GzStartorEnd, NOgzstORgzend, dud_trial] = LocalManIdxSpec(FLname,...
    CurrentCondNum, LetterTrack, jj, kk, TrialTbl, StartEndString)  

fprintf(['Participant: %s \n'...
'Condition: %d \n'...
'Marker: %s \n'...
'Trial: %d \n'], FLname, CurrentCondNum, LetterTrack(jj), kk)
disp(TrialTbl);
NOgzstORgzend = char('');
    
    while ~ismember(NOgzstORgzend, {'Y', 'N', 'WD'})
        
        NOgzstORgzend = input(sprintf(['Could not identify the gaze %s index.'...
            ' Manually specify? Y/N [N]: '], StartEndString), 's');
        if isempty(NOgzstORgzend)
            dud_trial = true;
            NOgzstORgzend = 0;
            GzStartorEnd = [];
            break
        end
    end
    
    if isempty(NOgzstORgzend) || strcmp(NOgzstORgzend, 'N')
        NOgzstORgzend = 0;
        dud_trial = true;
        GzStartorEnd = [];
        return
        
    elseif strcmp(NOgzstORgzend, 'Y')
        
        GzStartorEnd = 0;
        
        while ~ismember(GzStartorEnd, TrialTbl.Index)
            GzStartorEnd = input(['Enter a valid numeric index of your'...
                ' choice: ']);           
        end
        
        
            if GzStartorEnd == TrialTbl.Index(end) &&...
                    strcmp(StartEndString, 'RETURN')
                FlagNoReturn = char('');
                
                while ~ismember(FlagNoReturn, {'Y', 'N'})
                    FlagNoReturn = input('Flag trial as no return? Y/N [N]: ',...
                        's');
                    
                    if strcmp(FlagNoReturn, 'N') || isempty(FlagNoReturn)
                        NOgzstORgzend = 1;
                        
                    elseif strcmp(FlagNoReturn, 'Y')
                        NOgzstORgzend = 3;
                    end
                end
            else
                NOgzstORgzend = 1;
            end
                                
        dud_trial = false;
        return
    
    elseif strcmp(NOgzstORgzend, 'WD')
        dud_trial = true; 
        NOgzstORgzend = 2;
        GzStartorEnd = [];
        return
         
    else
        
        warning('Invalid input. Marking trial as dud.')
        NOgzstORgzend = 4;
        GzStartorEnd = [];
        dud_trial = true;
        return
        
    end 
end

function [GzStartWrong, AnticipatoryStartWrongDir] = DetectWrongStart(WSevtsubT,...
    CFT_SHIFT, FLname, CurrentCondNum, LetterTrack,Y, jj, kk, TrialTbl,...
    DIFF_TOLERANCE, xHOME, CHNGTOL, MkrPos, GzStart, GzReturn)  %#ok<INUSD>

MkrPosShift = circshift(MkrPos, 2); 
FALSE_ST_THRESH_PIX = -76;
AnticipatoryStartWrongDir = false;

GzStartWrong = find([0; diff(WSevtsubT.genx) < (FALSE_ST_THRESH_PIX + Y)] &...
    [0; diff(WSevtsubT.genx) > DIFF_TOLERANCE] & (WSevtsubT.genx < (CFT_SHIFT(jj) *...
    (xHOME + (CFT_SHIFT(jj)* CHNGTOL)))) & strcmp(WSevtsubT.codestring, 'ENDSACC'));


if ~isempty(GzStartWrong) && ~issorted([GzStart; GzStartWrong(1)])
    
    fprintf(['Participant: %s \n'...
        'Condition: %d \n'...
        'Marker: %s \n'...
        'Trial: %d \n\n'], FLname, CurrentCondNum, LetterTrack(jj), kk)
    disp(TrialTbl);
    
    
    fprintf('A directionally INCORRECT GAZE START was detected at INDEX %d. \n',...
        GzStartWrong(1))
    ConfirmWrongDirection = char('');
    AnticipatoryStartWrongDir = true;
    
    
%         while ~ismember(ConfirmWrongDirection, {'Y', 'N'})
%             ConfirmWrongDirection = input(['Flag trial as an aticipatory start'...
%                 ' in the wrong direction? Y/N [N]: '], 's');
% 
%             if strcmp(ConfirmWrongDirection, 'N') || isempty(ConfirmWrongDirection)
%                 AnticipatoryStartWrongDir = false;
%                 break
% 
%             elseif strcmp(ConfirmWrongDirection, 'Y')
%                 AnticipatoryStartWrongDir = true;
%                 break
%             end
%         end
    
%     GzReturnWrong = find([0; diff(WSevtsubT.genx(GzStartWrong:end))...
%     > (abs(FALSE_ST_THRESH_PIX) - Y)] & [0; diff(WSevtsubT.genx(GzStartWrong:end))...
%     < abs(DIFF_TOLERANCE)] & WSevtsubT.genx(GzStartWrong:end) > (CFT_SHIFT(jj) *...
%     (MkrPosShift(jj) + (CFT_SHIFT(jj)*(100-Y)))) &...
%     strcmp(WSevtsubT.codestring(GzStartWrong:end), 'ENDSACC')) + (GzStartWrong-1);    
%         
%     if ~isempty(GzReturnWrong) && ~isempty(GzReturn) && ~issorted([GzReturn;...
%             GzReturnWrong(1)])
%         fprintf('A directionally INCORRECT GAZE RETURN was detected at INDEX %d',...
%             GzReturnWrong(1))
%         pause()
%     end   
end

end

function [GzStart, GzReturn, DidResetStartReturn] = LocalGzMultiSweep(FLname,...
    CurrentCondNum, LetterTrack, jj, kk, TrialTbl, GzStart, GzReturn)

ResetGzStEnd = [NaN, NaN];
fprintf(['Participant: %s \n'...
'Condition: %d \n'...
'Marker: %s \n'...
'Trial: %d \n\n'], FLname, CurrentCondNum, LetterTrack(jj), kk)
disp(TrialTbl);
fprintf([' \n'...
    'Current Gaze Start Index: %d \n'...
    'Current Gaze Rreturn Index: %d \n'], GzStart, GzReturn)

ResetStartReturn = input(['Multiple gaze sweeps detected.'...
            ' Reset index values? Y/N [N]: '], 's');
        
    if isempty(ResetStartReturn) || strcmp(ResetStartReturn, 'N')
        DidResetStartReturn = false;
        return

    elseif strcmp(ResetStartReturn, 'Y')
        DidResetStartReturn = true;
       
        while nnz(ismember(ResetGzStEnd, TrialTbl.Index)) ~=2 
            ResetGzStEnd = input(['Enter two valid numeric indices of your'...
                ' choice corresponding to (1) the GzStart index; \n'...
                'and (2) the GzReturn index: ']);
        end

        GzStart = ResetGzStEnd(1);
        GzReturn = ResetGzStEnd(2);
        return
    end       
end


function [GzStart, GzReturn] = LocalIndexAnomaly(FLname,...
    CurrentCondNum, LetterTrack, jj, kk, TrialTbl, GzStart, GzReturn)

ResetGzStEnd = [NaN, NaN];
fprintf(['Participant: %s \n'...
'Condition: %d \n'...
'Marker: %s \n'...
'Trial: %d \n\n'], FLname, CurrentCondNum, LetterTrack(jj), kk)
disp(TrialTbl);
fprintf([' \n'...
    'Current Gaze Start Index: %d \n'...
    'Current Gaze Rreturn Index: %d \n'], GzStart, GzReturn)

ResetStartReturn = input(['Fixation bounds are nonsensical... \n'...
            ' Reset index values? Y/N [N]: '], 's');
        
    if isempty(ResetStartReturn) || strcmp(ResetStartReturn, 'N')
        return

    elseif strcmp(ResetStartReturn, 'Y')
       
       
        while nnz(ismember(ResetGzStEnd, TrialTbl.Index)) ~=2 
            ResetGzStEnd = input(['Enter two valid numeric indices of your'...
                ' choice corresponding to (1) the GzStart index; \n'...
                'and (2) the GzReturn index: ']);
        end

        GzStart = ResetGzStEnd(1);
        GzReturn = ResetGzStEnd(2);
        return
    end       
end

function [max_fix, min_gazex, bound4fix, peak_vel, land_idx, pix_land_gzx,... 
    pix_land_gzy, latency] = LocalCheckEmpty(max_fix, min_gazex, bound4fix, peak_vel,...
    land_idx, pix_land_gzx, pix_land_gzy, latency)

if isempty(max_fix)
    max_fix = NaN;
end

if isempty(min_gazex)
    min_gazex = NaN;
end

if isempty(bound4fix)
    bound4fix = NaN;
end

if isempty(peak_vel)
    peak_vel = NaN;
end

if isempty(land_idx)
    land_idx = NaN;
end

if isempty(pix_land_gzx)
    pix_land_gzx = NaN;
end

if isempty(pix_land_gzy)
    pix_land_gzy = NaN;
end

if isempty(latency)
    latency = NaN;
end

end

function TrackingMatrix = LocalFillTrackingMatrix(fileIDs, TrackingMatrix,...
    numTrialReps)

PIDs = unique(string(regexp(fileIDs, 'P\d*', 'match')));
M3 = size(TrackingMatrix, 3);
Header4Tracker = repmat({'PID', 'TrialID', 'GzSt_Thresh', 'GzSt_ManSelect',...
    'GzSt_DudTrial', 'GzSt_Index', 'GzSt_TrackTbl', 'GzRt_Thresh',... 
    'GzRt_ManSelect', 'GzRt_DudTrial','GzRt_Index', 'GzRt_TrackTbl',...
    'ANTCP_ST_RIGHT', 'ANTCP_ST_WRONG'},1,1,M3);

VarNames = varmaker2(["C1", "C2", "C3"],["A", "B", "C", "D"], 'false',...
    'true', numTrialReps, [], 'false', []);
VarNames = repmat(transpose(VarNames),1,1,M3);
PIDs = repmat(reshape(cellstr(PIDs),1,1,M3),size(TrackingMatrix,1),...
    1,1);

TrackingMatrix(1:end,1,:)=PIDs; 
TrackingMatrix(1:end,2,:)=VarNames;
TrackingMatrix = vertcat(Header4Tracker, TrackingMatrix);
  
end


