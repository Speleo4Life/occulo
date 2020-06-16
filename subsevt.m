function tevp = subsevt(edfName, trlType, evtType, ev_ts, ev_dat, xdf, ampf)
% tevp = subsevt(edfName, [,trlType] [,evtType], ev_ts, ev_dat, xdf, ampf)
% edfName (char): name of .edf file, including the extenesion
% trlType (char): 'A', 'B', 'C', 'D', or any comhination {Default 'ABCD'} 
% evtType (str): events to be subsetted - options are: STARTBLINK, 
%                 ENDBLINK, STARTFIX, ENDFIX, STARTSACC, ENDSACC   
%                 {Deafults: ENDFIX, ENDSACC}
% ev_ts   (dble): generated from AnalyzeXDF func
% ev_dat  (dble): generated from AnalyzeXDF func
% xdf   (struct): generated from AnalyzeXDF func
% ampf    (dble): if ENDSACC or ENDFIX are among the events specified for
% subsetting, the value set for ampf will determine the threshold amplitude
% required for the inclusion of an ENDSACC or ENDFIX event in the output
% table. If empty or omitted, defaults to 30 arc minutes (0.5 deg.).
% 
%
% The following script subsets events based on the trial timestamp data 
% stored in the XDF file. At current, event marker timestamps from LSL are
% more accurate for trial-by-event parsing than those supplied directly 
% over the LINK.
% 
% If 'ev_ts', 'ev_dat', and/or 'xdf' are unspecified, or otherwise empty,
% the subsevt() will attempt to load the variables from MATLAB's base
% workspace. Failing that, the user is prompted to run the AnalyzeXDF
% function, which saves relevant (and a few extra) variables into the base
% workspace. This function requires a somewhat modified version of AnalyzeXDF
% so as to assign the needed variables into the workspace. Ideally, run
% AnalyzeXDF prior to calling the subsevt (subset event) function. 
%
%
% Table output also includes x, y, and xy, amplitudes for 'ENDSACC' and 
% 'ENDFIX' events. Sample durations is calculated for all events.  
%
%
% Author: Raymond MacNeil, Vision Lab 
% For 'Closed-Eyes', a.k.a, "The Saccade Study" data.
% Study PIs: RRM HG JTE MC CD JD 
% Research Assistant: CZ
% June, 2019

disp('Processing request...')
validin_chr_trlType = perms(['A', 'B', 'C', 'D']);
validin_str_trlType = string(validin_chr_trlType(1:end,:));

if nargin < 1 || ~exist('edfName', 'var') || isempty(edfName)
    error('You must provide ''edfName''!');
end

if nargin < 2 || ~exist('trlType', 'var') || isempty(trlType)...
        || sum(~contains(validin_str_trlType,...
        string(trlType)) == numel(validin_str_trlType))
    trlType = ['A', 'B', 'C', 'D'];
end

if nargin < 3 || ~exist('evtType', 'var') || isempty(evtType)
    evtType = ["ENDFIX", "ENDSACC"];
end
 
if nargin < 6 || ~exist('ev_ts', 'var') || isempty(ev_ts)...
        || ~exist('ev_dat', 'var') || isempty(ev_dat) || ~exist('xdf', 'var')...
        || isempty(xdf) 
    try
        ev_ts = evalin('base', 'ev_ts');
        ev_dat = evalin('base', 'ev_dat');
        xdf = evalin('base', 'xdf');
        noXDF = false;
    catch
        noXDF = true;
    end
else
   noXDF = false;
end

if nargin < 7 || ~exist('ampf', 'var') || isempty(ampf) 
    ampf = 0.50;
elseif ~isfloat(ampf) || ampf < 0
   ampf = 0.50;
   warning(['Invalid input argument for ''ampf''...reverting to',...
       'default amplitude threshold: %.2f' ampf]);
end

while noXDF
        
runXDF = input(['The "ev_ts", "ev_dat", or "xdf", argument(s) was/were not found.',...
'Do you want to run the "AnalyzeXDF" function (y/n)? '], 's');

switch runXDF
    
    case 'y'
       try
          AnalyzeXDF();
          WaitSecs(1);
          ev_ts = evalin('base', 'ev_ts');
          ev_dat = evalin('base', 'ev_dat');
          xdf = evalin('base', 'xdf');
          noXDF = false;
       catch
           error(['Failed to call "AnalyzeXDF": Ensure that AnalyzeXDF.m ',...
               'and its dependencies have been added to your Matlab file path.']);
       end
       
    case 'n'
        error(['If "ev_ts", "ev_dat", and "xdf" have already been defined, then ensure',...
              'they are loaded into your workspace.']);
    
    otherwise
        disp('Incorrect input. Indicate ''y'' (for yes) or ''n'' (for no)');

end
end


% convert from .Edf to .Mat 
edf = Edf2Mat(edfName);
assignin('base', 'edf', edf);

% modify the ev_dat cell array so that the trial IDs independently reflect
% if they are of the calibration phase or experimental phase
cal = find(contains(ev_dat, 'cal'));
exp = find(contains(ev_dat, 'exp'));

% only condition 1 will contain calibration trials, so need to implement an 
% if statement for the exp string concactenation 
if ~isempty(cal)
calidx_rng = cal(1)+1 : cal(2)-1;
ev_dat(calidx_rng) = strcat('cal_', ev_dat(calidx_rng));
end
expidx_rng = exp(1)+1 : exp(2)-1;
ev_dat(expidx_rng) = strcat('exp_', ev_dat(expidx_rng));    

% get EyeLink timeseries equivelent index positions of ev_ts data
% Identify indices of 'A', 'B', 'C', and 'D' trial timestamps from the
% ev_dat cell array. In a previous version of the task, the code was such 
% that'A', 'B', 'C', and 'D' corresponded with strings char segments '_0_', 
% '_1_', '_2_', and '_3_', respectively. This was later changed for 
% clarity.
% The ev_dat cell array is derived from the .XDF file. 

if sum(contains(ev_dat, '_0_')) ~= 0
tstmp_trls_A = ev_ts(contains(ev_dat, '_0_'))';
tstmp_trls_B = ev_ts(contains(ev_dat, '_1_'))';
tstmp_trls_C = ev_ts(contains(ev_dat, '_2_'))';
tstmp_trls_D = ev_ts(contains(ev_dat, '_3_'))';
elseif sum(contains(ev_dat, '_A_')) ~= 0
tstmp_trls_A = ev_ts(contains(ev_dat, '_A_'))';
tstmp_trls_B = ev_ts(contains(ev_dat, '_B_'))';
tstmp_trls_C = ev_ts(contains(ev_dat, '_C_'))';
tstmp_trls_D = ev_ts(contains(ev_dat, '_D_'))';
end

ev_dat_Aidx = find(ismember(ev_ts, tstmp_trls_A));
ev_dat_Bidx = find(ismember(ev_ts, tstmp_trls_B));
ev_dat_Cidx = find(ismember(ev_ts, tstmp_trls_C));
ev_dat_Didx = find(ismember(ev_ts, tstmp_trls_D));


% determine max trial count for appropriate preallocation and recursion
mx_trial_count = max([length(tstmp_trls_A), length(tstmp_trls_B),...
                      length(tstmp_trls_C), length(tstmp_trls_D)])/2;
% ts_lsl = xdf{1, 1}.time_series(9,:)';  
% tolerance = 0.0000000032;

% determine the index of the EyeLink event stream and then the 
% indices of the respective trial timestamps in the overall
% LSL time stamp array
for i=1:numel(xdf)
 if strcmp({xdf{1,i}.info.name}, {'EyeLink'})
    el_strmidx = i;
 end
end


ts_lsl = xdf{1, el_strmidx}.time_stamps'; % LSL timestamp vector
ts_edf = xdf{1, el_strmidx}.time_series(8,:)'; % EDF timestamp in each LSL sample

% Get closest Eyelink LSL sample indices to event LSL timestamps
[~, lsl_ts_trlA_idx] = arrayfun(@(x) min(abs(x - ts_lsl)), tstmp_trls_A);
[~, lsl_ts_trlB_idx] = arrayfun(@(x) min(abs(x - ts_lsl)), tstmp_trls_B);
[~, lsl_ts_trlC_idx] = arrayfun(@(x) min(abs(x - ts_lsl)), tstmp_trls_C);
[~, lsl_ts_trlD_idx] = arrayfun(@(x) min(abs(x - ts_lsl)), tstmp_trls_D);
% Index into EDF timestamp vector to get event times in EDF file
trlA_ts_evtIdx = ts_edf(lsl_ts_trlA_idx); 
trlB_ts_evtIdx = ts_edf(lsl_ts_trlB_idx); 
trlC_ts_evtIdx = ts_edf(lsl_ts_trlC_idx);
trlD_ts_evtIdx = ts_edf(lsl_ts_trlD_idx); 


% preallocate struct variable for organizing the start and end trials
tstmp_trls(1:mx_trial_count) = struct('Astrt', {NaN}, 'Aend', {NaN}, 'Bstrt',...
    {NaN}, 'Bend', {NaN}, 'Cstrt', {NaN}, 'Cend', {NaN}, 'Dstrt',...
    {NaN}, 'Dend', {NaN});
    
% segregate start and end into a struct and parse by trial (even=st,odd=end)
Astrt = num2cell(trlA_ts_evtIdx(1:2:end));
Aend  = num2cell(trlA_ts_evtIdx(2:2:end));
Bstrt = num2cell(trlB_ts_evtIdx(1:2:end));
Bend  = num2cell(trlB_ts_evtIdx(2:2:end));
Cstrt = num2cell(trlC_ts_evtIdx(1:2:end));
Cend  = num2cell(trlC_ts_evtIdx(2:2:end));
Dstrt = num2cell(trlD_ts_evtIdx(1:2:end));
Dend  = num2cell(trlD_ts_evtIdx(2:2:end));
[tstmp_trls.Astrt] = Astrt{:};
[tstmp_trls.Aend]  = Aend{:};
[tstmp_trls.Bstrt] = Bstrt{:};
[tstmp_trls.Bend]  = Bend{:};
[tstmp_trls.Cstrt] = Cstrt{:};
[tstmp_trls.Cend]  = Cend{:};
[tstmp_trls.Dstrt] = Dstrt{:};
[tstmp_trls.Dend]  = Dend{:};

trlID_track_Astrt = horzcat(Astrt, ev_dat(ev_dat_Aidx(1:2:end))');
trlID_track_Aend = horzcat(Aend, ev_dat(ev_dat_Aidx(2:2:end))');
trlID_track_Bstrt = horzcat(Bstrt, ev_dat(ev_dat_Bidx(1:2:end))');
trlID_track_Bend = horzcat(Bend, ev_dat(ev_dat_Bidx(2:2:end))');
trlID_track_Cstrt = horzcat(Cstrt, ev_dat(ev_dat_Cidx(1:2:end))');
trlID_track_Cend = horzcat(Cend, ev_dat(ev_dat_Cidx(2:2:end))');
trlID_track_Dstrt = horzcat(Dstrt, ev_dat(ev_dat_Didx(1:2:end))');
trlID_track_Dend = horzcat(Dend, ev_dat(ev_dat_Didx(2:2:end))');

% preallocate struct array
idx_trls_rng(1:mx_trial_count) = struct('A', {zeros}, 'A_trl_ID', {''},...
    'B', {zeros}, 'B_trl_ID', {''}, 'C', {zeros}, 'C_trl_ID', {''}, 'D',...
    {zeros}, 'D_trl_ID', {''});

offadj = 100;

i=1;
while i < numel(idx_trls_rng)+1
idx_trls_rng(i).A = tstmp_trls(i).Astrt-offadj : tstmp_trls(i).Aend+offadj;
idx_trls_rng(i).A_trl_ID = trlID_track_Astrt{i,2}(1:end-6);
idx_trls_rng(i).B = tstmp_trls(i).Bstrt-offadj: tstmp_trls(i).Bend+offadj;
idx_trls_rng(i).B_trl_ID = trlID_track_Bstrt{i,2}(1:end-6);
idx_trls_rng(i).C = tstmp_trls(i).Cstrt-offadj : tstmp_trls(i).Cend+offadj;
idx_trls_rng(i).C_trl_ID = trlID_track_Cstrt{i,2}(1:end-6);
idx_trls_rng(i).D = tstmp_trls(i).Dstrt-offadj : tstmp_trls(i).Dend+offadj;
idx_trls_rng(i).D_trl_ID = trlID_track_Dstrt{i,2}(1:end-6);
i = i + 1;
end

lsl_elevt_tsrng = idx_trls_rng;

% concactenate trial sample ranges to ready for subsequent
% event-by-trial parsing of the EDF
lsl_elevt_tscat.A = horzcat(lsl_elevt_tsrng(:).A)';
lsl_elevt_tscat.B = horzcat(lsl_elevt_tsrng(:).B)';
lsl_elevt_tscat.C = horzcat(lsl_elevt_tsrng(:).C)';
lsl_elevt_tscat.D = horzcat(lsl_elevt_tsrng(:).D)';


% 'sttime' is chosen as the field to index the events by trial because the
% 'sttime' is relisted in the END'event' tstmps. In the edf struct array,
% the sttime field values are stored as uint32, but conversion to double is
% required to pass the array to the ismembertol func, which we will use to
% identify the indices


edf_elevt_tstmp = double([edf.RawEdf.FEVENT.sttime]'); 

% get the event indices for each trial type
% If there are occurences of the same element (timestamp values in this 
% instance) at different indices, we need to pass the 'OutputAllIndices'
% argument to the function, otherwise it defaults to the indx of the first
% occurence. However, the output is stored as a cell array, so we must 
% convert to matrix using cell2mat, which conveniently "unbundles" cells
% containing >1 numeric element.

tolerance = 0;
[~, trla_elevt_idx] = ismembertol(lsl_elevt_tscat.A,...
    edf_elevt_tstmp, tolerance, 'ByRows', false, 'OutputAllIndices', true);
trla_el_eidx = nonzeros(cell2mat(trla_elevt_idx)); 
[~, trlb_elevt_idx] = ismembertol(lsl_elevt_tscat.B,...
    edf_elevt_tstmp, tolerance, 'ByRows', false, 'OutputAllIndices', true);
trlb_el_eidx = nonzeros(cell2mat(trlb_elevt_idx)); 
[~, trlc_elevt_idx] = ismembertol(lsl_elevt_tscat.C,...
    edf_elevt_tstmp, tolerance, 'ByRows', false, 'OutputAllIndices', true);
trlc_el_eidx = nonzeros(cell2mat(trlc_elevt_idx)); 
[~, trld_elevt_idx] = ismembertol(lsl_elevt_tscat.D,...
    edf_elevt_tstmp, tolerance, 'ByRows', false, 'OutputAllIndices', true);
trld_el_eidx = nonzeros(cell2mat(trld_elevt_idx)); 


% Create the final subset products 
trlA_el_evts = edf.RawEdf.FEVENT(trla_el_eidx);
trlB_el_evts = edf.RawEdf.FEVENT(trlb_el_eidx);
trlC_el_evts = edf.RawEdf.FEVENT(trlc_el_eidx);
trlD_el_evts = edf.RawEdf.FEVENT(trld_el_eidx);

% preallocate trialID field
[trlA_el_evts(1:end,1).trlID] = char;
[trlB_el_evts(1:end,1).trlID] = char; 
[trlC_el_evts(1:end,1).trlID] = char; 
[trlD_el_evts(1:end,1).trlID] = char; 


i = 1;
while i <= numel(idx_trls_rng) 
idxA = ismember([trlA_el_evts(:).sttime], idx_trls_rng(i).A);
idxB = ismember([trlB_el_evts(:).sttime], idx_trls_rng(i).B);
idxC = ismember([trlC_el_evts(:).sttime], idx_trls_rng(i).C);
idxD = ismember([trlD_el_evts(:).sttime], idx_trls_rng(i).D);

    if sum(idxA) > 0
       [trlA_el_evts(idxA).trlID] = deal(idx_trls_rng(i).A_trl_ID);
    end
    if sum(idxB) > 0
       [trlB_el_evts(idxB).trlID] = deal(idx_trls_rng(i).B_trl_ID);
    end
    if sum(idxC) > 0
       [trlC_el_evts(idxC).trlID] = deal(idx_trls_rng(i).C_trl_ID);
    end
    if sum(idxC) > 0
       [trlD_el_evts(idxD).trlID] = deal(idx_trls_rng(i).D_trl_ID);
    end
i = i+1;   
end


tevp = table;

for i = 1:length(trlType)
    
    if strcmp(trlType(i), 'A')
       user_trlA_evts = trlA_el_evts(ismember({trlA_el_evts(:).codestring},...
       {evtType{:}}));
       tevpA = struct2table(user_trlA_evts);
       trltype = repmat("A", height(tevpA), 1);
       tevpA = addvars(tevpA, trltype, 'After', 'codestring'); 
       tevp = vertcat(tevp, tevpA); %#ok<AGROW>
       
    elseif strcmp(trlType(i), 'B')
  
        user_trlB_evts = trlB_el_evts(ismember({trlB_el_evts(:).codestring},... 
            {evtType{:}}));
        tevpB = struct2table(user_trlB_evts);
        trltype = repmat("B", height(tevpB), 1);
        tevpB = addvars(tevpB, trltype, 'After', 'codestring'); 
        tevp = vertcat(tevp, tevpB); %#ok<AGROW>
        
    elseif strcmp(trlType(i), 'C')
        
        user_trlC_evts = trlC_el_evts(ismember({trlC_el_evts(:).codestring},...
            {evtType{:}}));
        tevpC = struct2table(user_trlC_evts);
        trltype = repmat("C", height(tevpC), 1);
        tevpC = addvars(tevpC, trltype, 'After', 'codestring'); 
        tevp = vertcat(tevp, tevpC); %#ok<AGROW>
        
    elseif strcmp(trlType(i), 'D')

        user_trlD_evts = trlD_el_evts(ismember({trlD_el_evts(:).codestring},...
            {evtType{:}}));
        tevpD = struct2table(user_trlD_evts);
        trltype = repmat("D", height(tevpD), 1);
        tevpD = addvars(tevpD, trltype, 'Before', 'codestring'); 
        tevp = vertcat(tevp, tevpD); %#ok<AGROW>
    else
        error('Input arguments unrecognized.')
    end
end


 tevp = removevars(tevp, {'time', 'type', 'read', 'hstx',... 
            'hsty', 'sta', 'ava' 'henx', 'heny', 'ena', 'havx', 'havy',... 
            'status', 'flags', 'input', 'buttons', 'eye', 'parsedby'});
 evtdur = tevp.entime - tevp.sttime;
 tevp = addvars(tevp, evtdur, 'After', 'entime');
 trialID = string(tevp.trlID);
 evt_type = string(tevp.codestring);
 tevp = addvars(tevp, trialID, evt_type);
 tevp = movevars(tevp, {'trltype', 'trialID', 'evt_type'},'Before', 1);
 tevp = movevars(tevp, 'genx', 'After', 'gstx');
 tevp = movevars(tevp, 'geny', 'After', 'gsty');
 tevp = movevars(tevp, {'gavx', 'gavy'}, 'After', 'evel');

 
 % if 'ENDSACC' or 'ENDFIX' are amongs the event types for subsetting, 
 % create amplitude variables accordingly
 if sum(ismember({'ENDSACC', 'ENDFIX'}, {evtType{:}})) > 0
    
    i = 1;
    htevp = NaN(height(tevp),1);
    ampl_x = htevp; 
    ampl_y = htevp;
    ampl_xy = htevp;
    
        for i = 1:height(tevp)
            if strcmp('ENDSACC', tevp.codestring(i)) || strcmp('ENDFIX', tevp.codestring(i))
            ampl_x(i) = (tevp.genx(i) - tevp.gstx(i))...
                / ((tevp.eupd_x(i) + tevp.supd_x(i))/2); 
            ampl_y(i) = (tevp.geny(i) - tevp.gsty(i))...
                / ((tevp.eupd_y(i) + tevp.supd_y(i))/2); 
            ampl_xy(i) = sqrt(ampl_x(i).^2 + ampl_y(i).^2);
            end
        end

 tevp = addvars(tevp, ampl_x, 'After', 'genx');
 tevp = addvars(tevp, ampl_y, ampl_xy, 'After', 'eupd_y'); 
 ts1 = find(contains(tevp.codestring(:), 'ENDSACC')); 
 ts2 = find(abs(tevp.ampl_x) < ampf);
 rows_removed = tevp(intersect(ts1, ts2), :);
 tevp(intersect(ts1, ts2), :) = [];
 tevp = sortrows(tevp,{'sttime','trltype'},'ascend');
 assignin('base', 'rows_removed', rows_removed);  
 tevp = removevars(tevp, {'trlID', 'codestring'});
 end
 disp('Sucess! Now...pitter patter, let''s get at ''er!')
end
 
  
        