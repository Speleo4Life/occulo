% Closed-Eyes EyeLink Data Analysis
% Raymond MacNeil, 2020-February
% Vision Lab, Department of Psychology, University of British Columbia

% Define file names for analysis
cd '/Users/ray.macneil/Nextcloud/ubc_vision/occulo_imagery/data/Real_Participant_Data/EDFS'
ext = '.edf'
file_names = ls;
file_names = transpose(sort(string(regexp(file_names, ['\w*' ext],...
    'match'))));
subj_count = numel(unique(extractBefore(file_names, '_')));


%Import files batchwise
for ii = 1:length(file_names)
   
   try
       subj = extractBefore(file_names(ii), '_');
       cond = extractBetween(file_names(ii), '_', ext);
       edfm = Edf2Mat(char(file_names(ii)));
       edfall.(subj).(cond) = edfm;
       edfevt.(subj).(cond) = edfm.RawEdf.FEVENT;
   catch
       warning('There was a problem with file: %s', file_names(ii));
   end

end

% Declare out main results table, initially setting up variables for
% calibration info
T = table('Size', [subj_count 4], 'VariableTypes',...
     {'cellstr', 'string', 'string', 'string'}, 'VariableNames',...
     {'pid', 'c1cal', 'c2cal', 'c3cal'});
T.Properties.Description = 'Closed-Eyes EyeLink Data | MacNeil | Vision Lab | 2020'
ssize = numel(fieldnames(edfevt));
field_names = char(fieldnames(edfevt));
pids = unique(extractBefore(file_names, '_'))
T.pid = categorical(pids);

% Extract calibration information
for ii = 1:ssize
    fname = field_names(ii,:);
    cfnames = fieldnames(edfevt.(field_names(ii,:)));
    
    for jj = 1:numel(cfnames)
        idx = find(([edfevt.(fname).(cfnames{jj}).parsedby] == 188),1, 'last');
        if ~isempty(idx)
            T(ii,1+jj) = cellstr(edfevt.(fname).(cfnames{jj})(idx).message); 
        end
    end
end


% Initialize variable names and types for our table that will subset trial events  
new_var_names = varmaker(["C1", "C2", "C3"],["A", "B", "C", "D"], 'false', [],...
    ["ST", "EN"], "IDX");
vartypes = cellstr(repmat("cell",1,numel(new_var_names))); 
evt_indices = table('Size', [subj_count numel(new_var_names)],... 
    'VariableTypes', vartypes, 'VariableNames', new_var_names); 

% Other key information for subsetting trial events
letters = {'A'; 'B'; 'C'; 'D'}; 
seq = {'start'; 'end'};


for ii = 1:ssize
    fname = field_names(ii,:);
    cfnames = fieldnames(edfevt.(field_names(ii,:)));
    for jj = 1:numel(cfnames)
        get = struct2table(edfevt.(fname).(cfnames{jj}));
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

vnames = {'c1a_idx', 'c1b_idx', 'c1c_idx', 'c1d_idx',...
    'c2a_idx', 'c2b_idx', 'c2c_idx', 'c2d_idx',...
    'c3a_idx', 'c3b_idx', 'c3c_idx', 'c3d_idx'};
vtypes = cellstr(repmat("cell", 1, numel(vnames)));
tevt = table('Size', [subj_count, numel(vnames)], 'VariableTypes', vtypes,...
    'VariableNames', vnames);



for ii = 1:ssize
    d = 1;
    for jj = 1:2:width(evt_indices)-1
        holder1 = table2array(evt_indices(ii,jj));
        holder2 = table2array(evt_indices(ii,jj+1));
        loopsz = numel(holder1{:});
        temp = cell(loopsz,1);
        
        for kk = 1:loopsz
            temp{kk}= holder1{:}(kk):holder2{:}(kk);  
        end
        tevt{ii,d} = {temp};
        d=d+1;
    end
end

MATRIX = zeros(subj_count,120);
MATRIX2 = zeros(subj_count,12);

%Declare variables for accuracy calculations
ABidxVals = [1 2 5 6 9 10];
CDidxVals = [3 4 7 8 11 12];
smX_eccentric_offset = 161.40665321197295
lgX_eccentric_offset = 326.4192270718748
degpix = 0.06537277453582728
ABS = [980-lgX_eccentric_offset, 980-smX_eccentric_offset, 980+smX_eccentric_offset, 980+lgX_eccentric_offset]'


%% Main Loop %%

for ii = 1:ssize
    fname = field_names(ii,:)
    cfnames = fieldnames(edfevt.(field_names(ii,:)))
    cond_loop_tags = "";
    
    if numel(cfnames) < 3
       cond_loop_tags = string(repmat(cfnames{1}, numel(letters),1))
       cond_loop_tags = [cond_loop_tags; strrep(cond_loop_tags, '1', '2')];  %#ok<AGROW>
    else
       cond_loop_tags = string(repmat(cfnames{1}, numel(letters),1)); 
       cond_loop_tags = [cond_loop_tags; strrep(cond_loop_tags, '1', '2');...
           strrep(cond_loop_tags, '1', '3')];  %#ok<AGROW>
    end

    for jj = 1:numel(cond_loop_tags)
        dud_trial = false;
        if (jj == 1) || cond_loop_tags{jj}(2) ~= cond_loop_tags{jj-1}(2) 
            get = struct2table(edfevt.(fname).(cond_loop_tags{jj}));
            hget = NaN(height(get),1);
            ampl_x = hget;
            ampl_y = hget;
            ampl_xy = hget;
            
            for mm = 1:height(get)
                if strcmp('ENDSACC', get.codestring(mm)) || strcmp('ENDFIX',...
                        get.codestring(mm))
                    ampl_x(mm) = (get.genx(mm) - get.gstx(mm))...
                        / ((get.eupd_x(mm) + get.supd_x(mm))/2);
                    ampl_y(mm) = (get.geny(mm) - get.gsty(mm))...
                        / ((get.eupd_y(mm) + get.supd_y(mm))/2);
                    ampl_xy(mm) = sqrt(ampl_x(mm).^2 + ampl_y(mm).^2);
                end
                
            end
            
            get = addvars(get, ampl_x, 'After', 'genx');
            get = addvars(get, ampl_y, ampl_xy, 'After', 'eupd_y');
            evtdur = get.entime - get.sttime;
            get = addvars(get, evtdur, 'After', 'entime');
        end
        
        idx = tevt{ii,jj};
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
T2.pid = unique(extractBefore(file_names, '_'));


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

