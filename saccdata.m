% Closed-Eyes EyeLink Data Analysis
% Raymond MacNeil, 2020-February
% Vision Lab, Department of Psychology, University of British Columbia
% Updated April 27th, 2020


%% Data Import 

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
%      edfall.(subj).(cond) = edfm;
       edfevt.(subj).(cond) = edfm.RawEdf.FEVENT;
   catch
       warning('There was a problem with file: %s', fileIDs(ii));
   end

end

%% Get calibration information for each subject and condition 

% Preallocate a table for storing the data
T = table('Size', [N 4], 'VariableTypes',...
     {'cellstr', 'string', 'string', 'string'}, 'VariableNames',...
     {'PID', 'C1_CAL', 'C2_CAL', 'C3_CAL'});
T.Properties.Description = 'Closed-Eyes EyeLink Data | MacNeil | Vision Lab | 2020';
pidFLDnames = char(fieldnames(edfevt));
T.PID = categorical(unique(extractBefore(fileIDs, '_')));

% Extract calibration information and store it in our table
for ii = 1:N
    fname = pidFLDnames(ii,:);
    condFLDnames = fieldnames(edfevt.(pidFLDnames(ii,:)));
    
    for jj = 1:numel(condFLDnames)
        idx = find(([edfevt.(fname).(condFLDnames{jj}).parsedby] == 188),1, 'last');
        if ~isempty(idx)
            T(ii,1+jj) = cellstr(edfevt.(fname).(condFLDnames{jj})(idx).message); 
        end
    end
end

clear cond subj ii jj edfm fileDir idx

%% Preallocate table to store trial START and END timestamp indices

% Initialize variable names and types for our table that will store store  
var_names = varmaker(["C1", "C2", "C3"],["A", "B", "C", "D"], 'false',...
    'false', [], ["ST", "EN"], 'true', "IDX");
var_types = cellstr(repmat("cell",1,numel(var_names))); 
evt_indices = table('Size', [N numel(var_names)],... 
    'VariableTypes', var_types, 'VariableNames', var_names); 

% Other key information for subsetting trial events
letters = {'A'; 'B'; 'C'; 'D'}; 
seq = {'start'; 'end'};

clear var_names var_types

%% Get START and END timestamp indices for each unique trial type
for ii = 1:N
    fname = pidFLDnames(ii,:);
    condFLDnames = fieldnames(edfevt.(fname));
    for jj = 1:numel(condFLDnames)
        get = struct2table(edfevt.(fname).(condFLDnames{jj}));
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

%% Using the above info, determine the index range for each trial's timestamp range

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

clear holder1 holder2 var_names var_types x ii jj kk loopsz d
clear temp new_var_names extract check 

%% Preallocate variables for main loop extracting event information
MATRIX = zeros(N,120);
MATRIX2 = zeros(N,12);

%% Declare variables for accuracy calculations
ABidxVals = [1 2 5 6 9 10];
CDidxVals = [3 4 7 8 11 12];
ABCDidxVals = [ABidxVals, CDidxVals];  
smX_eccentric_offset = 161.40665321197295; % Units are Pixels
lgX_eccentric_offset = 326.4192270718748; % Units are Pixels
degpix = 0.06537277453582728; % Degrees visual angle per Pixel
xHOME = 980; % Homebase marker is at 980 pixels along x-axis
ABS = [xHOME-lgX_eccentric_offset, xHOME-smX_eccentric_offset,...
    xHOME+smX_eccentric_offset, xHOME+lgX_eccentric_offset]';


%% Loop to calculate amplitude variables
for ii = 1:N
    fname = pidFLDnames(ii,:);
    condFLDnames = fieldnames(edfevt.(pidFLDnames(ii,:)));
    CondLoop = repelem(string(condFLDnames), numel(letters));


    for jj = 1:numel(CondLoop)
        
        
        if (jj == 1) || CondLoop{jj}(2) ~= CondLoop{jj-1}(2) 
            
            get = struct2table(edfevt.(fname).(CondLoop{jj}));
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
            get = removevars(get,{'hstx', 'hsty', 'henx', 'heny', 'sta', 'ena',...
               'havx', 'havy', 'ava', 'status', 'flags', 'input', 'buttons',...
               'time', 'type', 'read', 'eye'});
            get = movevars(get, {'genx', 'gavx', 'ampX', 'ampXY'}, 'After', 'gstx');
            get = movevars(get, 'ampY', 'After', 'gavy');
            get = movevars(get, 'codestring', 'Before', 'sttime');
            edfevt.(fname).(CondLoop{jj}) = table2struct(get);
        end
    end
end


%% Main Loop %%

for ii = 1:N
    fname = pidFLDnames(ii,:)
    condFLDnames = fieldnames(edfevt.(pidFLDnames(ii,:)))
    CondLoop = repelem(string(condFLDnames), numel(letters))
    
    for jj = 1:numel(CondLoop)
        
        dud_trial = false;
        idx = TEVT{ii,jj};
        key_evts = zeros(10,numel(idx{:}));
        
        for kk = 1:numel(idx{:})
            dud_trial = false;
            evtidx = cell2mat(idx{:}(kk));
            evtsubTb = get(evtidx,:);
            subset = strcmp(evtsubTb.codestring, 'ENDSACC') |... 
                strcmp(evtsubTb.codestring, 'ENDFIX'); 
            evtsubT = evtsubTb(subset,:);
               
            if any(jj == ABidxVals) % A OR B
                
                
                % Determine appropriate ABS index for accuracy computations
                if any(jj == ABidxVals(1:2:end-1)) 
                   abs_idx = 1;
                else
                   abs_idx = 2;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Key variables for analysis %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                inf_trial = find(evtsubT.gstx==100000000 | evtsubT.gstx > 2000);
                evtsubT(inf_trial, :) = [];
                gz_start =  find([0; diff(evtsubT.genx) < -100] & [0; diff(evtsubT.genx)...
                    > -800] & evtsubT.genx < (980-100),1, 'first');
                gz_return = find([0; diff(evtsubT.genx(gz_start:end)) > 100] &...
                    [0; diff(evtsubT.genx(gz_start:end)) < 800] &...
                    evtsubT.genx(gz_start:end) > ABS(abs_idx)+100, 1, 'first')...
                    + (gz_start-1);
                
                
                if isempty(gz_start) 
                    gz_start =  find([0; diff(evtsubT.genx) < -50] &...
                        [0; diff(evtsubT.genx) > -800] & evtsubT.genx < (980-50),1,...
                        'first');
                    gz_return = find([0; diff(evtsubT.genx(gz_start:end)) > 100] &...
                        [0; diff(evtsubT.genx(gz_start:end)) < 800] &...
                        evtsubT.genx(gz_start:end) > ABS(abs_idx)+100, 1, 'first')...
                        + (gz_start-1);
                end
                
            
                if isempty(gz_return)
                    gz_return = find([0; diff(evtsubT.genx(gz_start:end)) > 50] &...
                        [0; diff(evtsubT.genx(gz_start:end)) < 800] & evtsubT.genx(gz_start:end)...
                        > ABS(abs_idx)+50, 1, 'first') + (gz_start-1);
                end
                
                if isempty(gz_start)
                    gz_start = find(evtsubT.ampl_xy>2,1,'first');
                    gz_return = find([0; diff(evtsubT.genx(gz_start:end)) > 50] &...
                        [0; diff(evtsubT.genx(gz_start:end)) < 800] &...
                        evtsubT.genx(gz_start:end) > ABS(abs_idx)+50, 1, 'first')...
                        + (gz_start-1);
                end
               


                if isempty(gz_start) || isempty(gz_return) || isnan(gz_start)...
                        || isnan(gz_return)
                    dud_trial = true;
                end
                
                
                try
                % Calculate Saccade Latency
                    latency = evtsubT.sttime(gz_start)-evtsubT.sttime(1);
                catch
                    latency = NaN;
                end
                

                if ~dud_trial
                % Get minnimum X and Y gaze coordinates
                [min_gazex, Imng] = min(evtsubT.genx(1:end));
                 min_gazey = evtsubT.geny(Imng);
                
                % Get saccade amplitudes
                [max_amplx, Imxa] = min(evtsubT.ampl_x(1:end));
                [max_amplxy, ~] = max(evtsubT.ampl_xy(Imxa));
                
%                 boundtol = 100;
%                 for hh = 40:boundtol
%                 bound4fix = find(evtsubT.genx < ABS(abs_idx)+hh & evtsubT.genx ~= 0);
%                     if ~isempty(bound4fix)
%                         break;
%                     end    
%                 end
                bound4fix = gz_start+1:gz_return-1
                [max_fix, land_idx] = max(evtsubT.evtdur(bound4fix));
                land_idx = bound4fix(land_idx);
                peak_vel = max(evtsubT.pvel(gz_start:gz_return-1));
                num_sacc = sum(strcmp('ENDSACC',...
                    evtsubT.codestring(gz_start:gz_return-1)));
                
                try
                    land_gzx = evtsubT.gavx(land_idx);
                    land_gzy = evtsubT.gavy(land_idx);
                catch
                    land_gzx = NaN;
                    land_gzy = NaN;
                    dud_trial = true;
                end
                
                    
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
                if isempty(land_gzx)
                    land_gzx = NaN;
                end
                if isempty(land_gzy)
                    land_gzy = NaN;
                end
                if isempty(latency)
                    latency = NaN;
                end
                
                try
                    acry_pix = sqrt((land_gzx-ABS(abs_idx)).^2 + (land_gzy-590).^2);
                    acry_deg = acry_pix * degpix;
                catch
                    warning('something wrong with acry_pix and acry_deg')
                    acry_pix = NaN;
                    acry_deg = NaN;
                end
                end    

                % Calculate saccadic accuracy in pixels and degrees
                if acry_deg > 20
                   acry_deg = NaN;
                end
                
                try
                % Define Key Events for Results Table %
                key_evts(1,kk)=min_gazex; key_evts(2,kk)=num_sacc;
                key_evts(3,kk)=max_amplx; key_evts(4,kk)=max_amplxy;
                key_evts(5,kk)=peak_vel; key_evts(6,kk)=land_gzx;
                key_evts(7,kk)=land_gzy; key_evts(8,kk)=max_fix;
                key_evts(9,kk)=acry_deg; key_evts(10,kk)=latency;
                catch
                warning('There was a problem parsing this trial...skipping to next.')
                dud_trial = true;
                end
                
            elseif any(jj == CDidxVals)
                
                if any(jj == CDidxVals(1:2:end-1)) 
                   abs_idx = 3;
                else
                   abs_idx = 4;
                end
                
                inf_trial = find(evtsubT.gstx==100000000 | evtsubT.gstx > 2000);
                evtsubT(inf_trial, :) = [];
                gz_start =  find([0; diff(evtsubT.genx) > 100] & [0; diff(evtsubT.genx)...
                    < 800] & evtsubT.genx > (980+100),1, 'first');
                gz_return = find([0; diff(evtsubT.genx(gz_start:end)) < -100] &...
                    [0; diff(evtsubT.genx(gz_start:end)) > -800] & evtsubT.genx(gz_start:end)...
                    < ABS(abs_idx)-100, 1, 'first') + (gz_start-1);
                
                if isempty(gz_start) 
                    gz_start =  find([0; diff(evtsubT.genx) > 50] & [0; diff(evtsubT.genx)...
                        < 800] & evtsubT.genx > (980+50),1, 'first');
                end
                
                if isempty(gz_return)
                    gz_return = find([0; diff(evtsubT.genx(gz_start:end)) < -50] &...
                        [0; diff(evtsubT.genx(gz_start:end)) > -800] & evtsubT.genx(gz_start:end)...
                        < ABS(abs_idx)-50, 1, 'first');
                end
                
                if isempty(gz_start) 
                    gz_start = find(evtsubT.ampl_xy>2,1,'first');
                    gz_return = find([0; diff(evtsubT.genx(gz_start:end)) < -50] &...
                        [0; diff(evtsubT.genx(gz_start:end)) > -800] &...
                        evtsubT.genx(gz_start:end) < ABS(abs_idx)-50, 1, 'first');
                end
                
                if isempty(gz_start) || isempty(gz_return) ||...
                   isnan(gz_start) || isnan(gz_return)
                   dud_trial = true
                end
                
                try
                % Calculate Saccade Latency
                    latency = evtsubT.sttime(gz_start)-evtsubT.sttime(1);
                catch
                    latency = NaN;
                end
                
                if ~dud_trial
              
                    [max_gazex, Imng] = max(evtsubT.genx);
                    max_gazey = evtsubT.geny(Imng);
                    [max_amplx, Imxa] = max(evtsubT.ampl_x);
                    max_amplxy = max(evtsubT.ampl_xy(Imxa));
                    
                    
%                     boundtol = 200;
%                     for hh = 40:boundtol
%                         bound4fix = find(evtsubT.genx...
%                             > ABS(abs_idx)-hh & evtsubT.genx ~= 0); 
%                         if ~isempty(bound4fix) 
%                             break;
%                         end
%                     end
                    bound4fix = gz_start+1:gz_return-1

                    [max_fix, land_idx] = max(evtsubT.evtdur(bound4fix));
                    land_idx = bound4fix(land_idx);
                    peak_vel = max(evtsubT.pvel(1:gz_return-1));
                    num_sacc = sum(strcmp('ENDSACC',...
                    evtsubT.codestring(gz_start:gz_return-1)));
                    land_gzx = evtsubT.gavx(land_idx);
                    land_gzy = evtsubT.gavy(land_idx);
                end
                
                if isempty(max_fix)
                    max_fix = NaN;
                end
                
                if isempty(max_gazex)
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
                if isempty(land_gzx)
                    land_gzx = NaN;
                end
                if isempty(land_gzy)
                    land_gzy = NaN;
                end
                if isempty(latency)
                    latency = NaN;
                end
                    
                try
                    acry_pix = sqrt((land_gzx-ABS(abs_idx)).^2 + (land_gzy-590).^2);
                    acry_deg = acry_pix * degpix;
                catch
                    warning('something wrong with acry_pix and acry_deg')
                    acry_pix = NaN;
                    acry_deg = NaN;
                end     
                
                if acry_deg > 20
                   acry_deg = NaN;
                end

                try
                    key_evts(1,kk)=max_gazex;
                    key_evts(2,kk)=num_sacc;
                    key_evts(3,kk)=max_amplx;
                    key_evts(4,kk)=max_amplxy;
                    key_evts(5,kk)=peak_vel;
                    key_evts(6,kk)=land_gzx;
                    key_evts(7,kk)=land_gzy;
                    key_evts(8,kk)=max_fix;
                    key_evts(9,kk)=acry_deg;
                    key_evts(10,kk)=latency;
                catch
                    dud_trial = true;
                end               
            end
            
        if dud_trial
            key_evts(1:10,kk) = NaN;
        elseif any(key_evts(1:10,kk)>2200)
            key_evts(1:10,kk) = NaN;
        end
        
        end
        idx1 = ii;
        idx2 = (jj*10)-10+1;
        MATRIX(idx1,idx2:idx2+9) = transpose(nanmean(key_evts,2));
        MATRIX2(ii,jj) = sum(isnan(key_evts(1,:)));
    end
              
end
        

vnames2 = {'a_min_gzx_c1', 'a_nsacc_c1', 'a_max_ampx_c1',...
    'a_max_amp_xy_c1', 'a_pvel_c1' 'a_xfix_c1', 'a_yfix_c1', 'a_maxfix_c1',...
    'a_acrcy_c1', 'a_latncy_c1', 'b_min_gzx_c1', 'b_nsacc_c1', 'b_max_ampx_c1',...
    'b_max_amp_xy_c1', 'b_pvel_c1' 'b_xfix_c1', 'b_yfix_c1', 'b_maxfix_c1',...
    'b_acrcy_c1', 'b_latncy_c1','c_max_gzx_c1', 'c_nsacc_c1', 'c_max_ampx_c1',...
    'c_max_amp_xy_c1', 'c_pvel_c1' 'c_xfix_c1', 'c_yfix_c1', 'c_maxfix_c1',...
    'c_acrcy_c1', 'c_latncy_c1','d_max_gzx_c1', 'd_nsacc_c1', 'd_max_ampx_c1',...
    'd_max_amp_xy_c1', 'd_pvel_c1' 'd_xfix_c1', 'd_yfix_c1', 'd_maxfix_c1',...
    'd_acrcy_c1', 'd_latncy_c1','a_min_gzx_c2', 'a_nsacc_c2', 'a_max_ampx_c2',...
    'a_max_amp_xy_c2', 'a_pvel_c2' 'a_xfix_c2', 'a_yfix_c2', 'a_maxfix_c2',...
    'a_acrcy_c2', 'a_latncy_c2', 'b_min_gzx_c2', 'b_nsacc_c2', 'b_max_ampx_c2',...
    'b_max_amp_xy_c2', 'b_pvel_c2' 'b_xfix_c2', 'b_yfix_c2', 'b_maxfix_c2',...
    'b_acrcy_c2', 'b_latncy_c2','c_max_gzx_c2', 'c_nsacc_c2', 'c_max_ampx_c2',...
    'c_max_amp_xy_c2', 'c_pvel_c2' 'c_xfix_c2', 'c_yfix_c2', 'c_maxfix_c2',...
    'c_acrcy_c2', 'c_latncy_c2','d_max_gzx_c2', 'd_nsacc_c2', 'd_max_ampx_c2',...
    'd_max_amp_xy_c2', 'd_pvel_c2' 'd_xfix_c2', 'd_yfix_c2', 'd_maxfix_c2',...
    'd_acrcy_c2', 'd_latncy_c2','a_min_gzx_c3', 'a_nsacc_c3', 'a_max_ampx_c3',...
    'a_max_amp_xy_c3', 'a_pvel_c3' 'a_xfix_c3', 'a_yfix_c3', 'a_maxfix_c3',...
    'a_acrcy_c3', 'a_latncy_c3', 'b_min_gzx_c3', 'b_nsacc_c3', 'b_max_ampx_c3',...
    'b_max_amp_xy_c3', 'b_pvel_c3' 'b_xfix_c3', 'b_yfix_c3', 'b_maxfix_c3',...
    'b_acrcy_c3', 'b_latncy_c3','c_max_gzx_c3', 'c_nsacc_c3', 'c_max_ampx_c3',...
    'c_max_amp_xy_c3', 'c_pvel_c3' 'c_xfix_c3', 'c_yfix_c3', 'c_maxfix_c3',...
    'c_acrcy_c3', 'c_latncy_c3','d_max_gzx_c3', 'd_nsacc_c3', 'd_max_ampx_c3',...
    'd_max_amp_xy_c3', 'd_pvel_c3' 'd_xfix_c3', 'd_yfix_c3', 'd_maxfix_c3',...
    'd_acrcy_c3', 'd_latncy_c3'};

vtypes2 = cellstr(repmat("double",1,120));
newvars = num2cell(MATRIX);
newvars = cell2table(newvars)
newvars.Properties.VariableNames=vnames2

T4 = [T newvars];
T2.PID = unique(extractBefore(fileIDs, '_'));


ad_acrry_c1 = nanmean([T4.a_acrcy_c1, T4.d_acrcy_c1], 'all')
bc_acrry_c1 = nanmean([T4.b_acrcy_c1, T4.c_acrcy_c1], 'all')
all_acrry_c1 = nanmean([T4.a_acrcy_c1, T4.b_acrcy_c1, T4.c_acrcy_c1,...
    T4.d_acrcy_c1], 2);
all_acrry_c2 = nanmean([T4.a_acrcy_c2, T4.b_acrcy_c2, T4.c_acrcy_c2,...
    T4.d_acrcy_c2], 2);
all_acrry_c3 = nanmean([T4.a_acrcy_c3, T4.b_acrcy_c3, T4.c_acrcy_c3,...
    T4.d_acrcy_c3], 2);
all_pvel_c1 = nanmean([T4.a_pvel_c1, T4.b_pvel_c1, T4.c_pvel_c1,...
    T4.d_pvel_c1], 2);
all_pvel_c2 = nanmean([T4.a_pvel_c2, T4.b_pvel_c2, T4.c_pvel_c2,...
    T4.d_pvel_c2], 2);
all_pvel_c3 = nanmean([T4.a_pvel_c3, T4.b_pvel_c3, T4.c_pvel_c3,...
    T4.d_pvel_c3], 2);
all_latency_c1 = nanmean([T4.a_latncy_c1, T4.b_latncy_c1, T4.c_latncy_c1,...
    T4.d_latncy_c1],2);
all_latency_c2 = nanmean([T4.a_latncy_c2, T4.b_latncy_c2, T4.c_latncy_c2,...
    T4.d_latncy_c2],2);
all_latency_c3 = nanmean([T4.a_latncy_c3, T4.b_latncy_c3, T4.c_latncy_c3,...
    T4.d_latncy_c3],2);
T4 = addvars(T4, all_acrry_c1, all_pvel_c1, all_latency_c1, 'After',44);  
T4 = addvars(T4, all_acrry_c2, all_pvel_c2, all_latency_c2, 'After',87);  
T4 = addvars(T4, all_acrry_c3, all_pvel_c3, all_latency_c3, 'After',130);

