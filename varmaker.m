function VarNamesOut = varmaker(condLabel, triaLabel, loopTrial, numTrialRep, lab3, srtcor, lab4) %#ok<INUSD>
% function VarNamesOut = VarMaker(condLabel, triaLabel, loopTrial,
% numTrialRep, lab3, lab4)
% condLabel <STR> A string array containing the condition/group labels
% triaLabel <STR> A string array containing the trial labels
% loopTrial <CHR> A character vector, 'true' OR 'false', to indicate if it
% there will be repeated trials/measures
% numTrialRep <POS INTEGER> Number of trial type variations
% lab3 <STR> A string array representing another component of the trials
% lab4 <STR> A single string acting as an optional universal tag
% Example: VarNamesOut = varmaker(["C1", "C2", "C3"],["A", "B", "C", "D"],... 
%               'true', 3, ["ST", "EN"], "IDX") 

% % % % % % % % % % % % % % % % % % % % % % % %
% Ray R. MacNeil                              %
% UBC, Department of Psychology, Vision Lab   %
% Email: raymond.macneil@mail.utoronto.ca     %
% % % % % % % % % % % % % % % % % % % % % % % %
%
% Useage Notes: This was made as a semi-dynamic function though there is
% some specificity built into. To get the variable label order correct a
% a different sorting pattern will likely be required. You can express the
% sorting pattern as a vector of row indices. See below for more details.

%% Exception Handling and Input Error

if  nargin < 4 || ~exist('numTrialRep', 'var') || isempty(numTrialRep) ...
        || ~isinteger(numTrialRep)
    numTrialRep = 0;
end
% 
if  nargin < 5 || ~exist('lab3', 'var') || isempty(lab3) || ~isstring(lab3)
    lab3 = 1;
end
%  %
if  nargin < 6 || ~exist('srtcor','var') || isempty(srtcor) || ~ischar(srtcor)
    srtcor = 'false';    
end
    
if  nargin < 7 || ~exist('lab4', 'var') || isempty(lab4) || ~isstring(lab4)
    lab4 = NaN;
end


    
 

%% Main Activity

if ~isempty(loopTrial) && strcmp(loopTrial, 'true')
    
    rows2iter = numTrialRep * numel(triaLabel);  
    VarNamesOut = strings(rows2iter, numel(condLabel)); 
    
    for jj = 1:numel(triaLabel)

        for kk = 1:numTrialRep
            repnum = num2str(kk);
            index = kk + (numTrialRep * (jj-1));
            VarNamesOut(index,:) = condLabel + "_" +...
                triaLabel(jj) + "_" + repnum; 
        end     
    end

else
    rows2iter = numel(triaLabel);
    VarNamesOut = strings(rows2iter, numel(condLabel)); 

    for jj = 1:numel(triaLabel)
            VarNamesOut(jj,:) = condLabel + "_" + triaLabel(jj);
    end
end    
    
if isstring(lab3)
    VarNamesOut = repmat(VarNamesOut, numel(lab3), 1);
    
    for ii = 1:numel(lab3)
        idxSTrng = (rows2iter * ii) - (rows2iter - 1);
        idxENrng = rows2iter * ii;
        VarNamesOut(idxSTrng:idxENrng,:) = VarNamesOut(idxSTrng:idxENrng,:) +...
            "_" + lab3(ii);
    end
end

if isstring(lab4)
    VarNamesOut = VarNamesOut + "_" + lab4;
end

% You May Need to Edit this to fit your specific sorting requirements
% determine the length of your sort vector containing the ii+x index with 
% x >= 0 that represents the swapping pattern 
VarNamesOut = sortrows(VarNamesOut);

if strcmp(srtcor, 'true')
    
    
    [~, srtIDXfwd] = sort(lab3);
    [~, srtIDXrev] = sort(srtIDXfwd);
    srtVec = 0:1:numel(lab3)

    for ii = 1:numel(lab3):length(VarNamesOut)
        VarNamesOut(srtVec+ii,:) = VarNamesOut(srtIDXrev + (ii-1),:);
    end
end

VarNamesOut = cellstr(reshape(VarNamesOut, 1, []));


end


