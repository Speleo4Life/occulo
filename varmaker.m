% function VarNamesOut = VarMaker(condLabel, triaLabel, loopTrial,
% numTrialRep, lab3, sortcorL3, lab4)
% condLabel <STR> A string array containing the condition/group labels
% triaLabel <STR> A string array containing the trial labels
% sortcorLT <CHR> A character vector, 'true' OR 'false', to indicate if
%                   sort correction will be need for the trial label.
% loopTrial <CHR> A character vector, 'true' OR 'false', to indicate if
%                   there will numeric repetitions
% numTrialRep <POS INTEGER> Number of trial type variations
% lab3 <STR> A string array representing another component of the trials
% srtcorL3 <CHR> A character vector, 'true' OR 'false', to indicate if
%                sort correction is required for label3. 
% lab4 <STR> A single string acting as an optional universal tag
% Example: VarNamesOut = varmaker(["C1", "C2", "C3"],["A", "B", "C", "D"],... 
%               'true', 2, ["ST", "EN"], "IDX") 
%
% Useage Notes: This was made as a semi-dynamic function though there is
% some specificity built into it. To get the variable label order correct
% in some case will necessitate the use of a different sorting pattern.
% Currently there is a sorting mechanism in place that will preserve the order 
% for the input values provided for lab3, and it should be straightforward to 
% eventually expand this functionality for the other input arguments. 
% 
% Ray R. MacNeil                              
% UBC, Department of Psychology, Vision Lab   
% Email: raymond.macneil@mail.utoronto.ca 
function VarNamesOut = varmaker(condLabel, triaLabel, sortcorLT, loopTrial, numTrialRep, lab3, sortcorL3, lab4)


%% Exception Handling and Input Error

if  nargin < 4 || ~exist('sortcorLT','var') || isempty(sortcorL3) || ~ischar(sortcorL3)
    sortcorLT = 'false';    
end

if  nargin < 5 || ~exist('numTrialRep', 'var') || isempty(numTrialRep)...
        || ~isnumeric(numTrialRep) || numTrialRep ~= round(numTrialRep)...
        || ~isscalar(numTrialRep)
    numTrialRep = 0;
end

if  nargin < 6 || ~exist('lab3', 'var') || isempty(lab3) || ~isstring(lab3)
    lab3 = 1;
end

if  nargin < 7 || ~exist('sortcorL3','var') || isempty(sortcorL3) || ~ischar(sortcorL3)
    sortcorL3 = 'false';    
end
    
if  nargin < 8 || ~exist('lab4', 'var') || isempty(lab4) || ~isstring(lab4)
    lab4 = NaN;
end
    
%% Main Activity

numTLs = numel(triaLabel);

if ~isempty(loopTrial) && strcmp(loopTrial, 'true')
    
    rows2iter = numTrialRep * numel(triaLabel);  
    VarNamesOut = strings(rows2iter, numel(condLabel)); 
    
    for jj = 1:numTLs

        for kk = 1:numTrialRep
            repnum = num2str(kk);
            index = kk + (numTrialRep * (jj-1));
            VarNamesOut(index,:) = condLabel + "_" +...
                triaLabel(jj) + "_" + repnum; 
        end     
    end

else
    rows2iter = numTLs;
    VarNamesOut = strings(rows2iter, numel(condLabel)); 

    for jj = 1:numTLs
            VarNamesOut(jj,:) = condLabel + "_" + triaLabel(jj);
    end
end    
    
if isstring(lab3)
    
    numL3s = numel(lab3)
    VarNamesOut = repmat(VarNamesOut, numL3s, 1);
    
    for ii = 1:numL3s
        idxSTrng = (rows2iter * ii) - (rows2iter - 1);
        idxENrng = rows2iter * ii;
        VarNamesOut(idxSTrng:idxENrng,:) = VarNamesOut(idxSTrng:idxENrng,:) +...
            "_" + lab3(ii);
    end
end

if isstring(lab4)
    VarNamesOut = VarNamesOut + "_" + lab4;
end

VarNamesOut = sortrows(VarNamesOut);


% If requested perform sorting on trial labels so that the input order is
% preserved in the final output.
if strcmp(sortcorLT, 'true')
    % srtInterv = numTLs * numTrialRep * numL3s;
    srtPatern = repelem(triaLabel, numTrialRep*numL3s);  
    %[M,~] = size(VarNamesOut);
    [~, srtIDXfwd] = sort(srtPatern);
    [~, srtIDXrev] = sort(srtIDXfwd);
    %srtVec = 0:1:(srtInterv-1)
    VarNamesOut(:,:) = VarNamesOut(srtIDXrev,:);
   % for ii = 1:srtInterv:M
       % VarNamesOut(srtVec+ii,:) = VarNamesOut(srtIDXrev + (ii-1),:);
   % end
end


% If requested perform sorting on trial labels so that the input order is
% preserved in the final output.
if strcmp(sortcorL3, 'true')
    
    [M,~] = size(VarNamesOut);
    [~, srtIDXfwd] = sort(lab3);
    [~, srtIDXrev] = sort(srtIDXfwd); 
    srtVec = 0:1:numL3s-1

    for ii = 1:numL3s:M
        VarNamesOut(srtVec+ii,:) = VarNamesOut(srtIDXrev + (ii-1),:);
    end
end

VarNamesOut = cellstr(reshape(VarNamesOut, 1, []));


end