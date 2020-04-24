% condLabel = ["C1", "C2", "C3"];
% triaLabel = ["A", "B", "C", "D"];
% loopTrial = 'true';
% numTrialRep = 2;
% lab3 = ["ST", "EN"];
% srtcor = 'true';
% lab4 = "IDX";

function VarNamesOut = varmaker(condLabel, triaLabel, loopTrial, numTrialRep, lab3, srtcorL3, lab4) 
% function VarNamesOut = VarMaker(condLabel, triaLabel, loopTrial,
% numTrialRep, lab3, lab4)
% condLabel <STR> A string array containing the condition/group labels
% triaLabel <STR> A string array containing the trial labels
% loopTrial <CHR> A character vector, 'true' OR 'false', to indicate if
% there will be repeated trials/measures
% numTrialRep <POS INTEGER> Number of trial type variations
% lab3 <STR> A string array representing another component of the trials
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

%% Exception Handling and Input Error

if  nargin < 4 || ~exist('numTrialRep', 'var') || isempty(numTrialRep)...
        || ~isnumeric(numTrialRep) || numTrialRep ~= round(numTrialRep)...
        || ~isscalar(numTrialRep)
    numTrialRep = 0;
end

if  nargin < 5 || ~exist('lab3', 'var') || isempty(lab3) || ~isstring(lab3)
    lab3 = 1;
end

if  nargin < 6 || ~exist('srtcor','var') || isempty(srtcorL3) || ~ischar(srtcorL3)
    srtcorL3 = 'false';    
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
    
    L3numel = numel(lab3)
    VarNamesOut = repmat(VarNamesOut, L3numel, 1);
    
    for ii = 1:L3numel
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

if strcmp(srtcorL3, 'true')
    

    [M,~] = size(VarNamesOut);
    [~, srtIDXfwd] = sort(lab3);
    [~, srtIDXrev] = sort(srtIDXfwd);
    srtVec = 0:1:L3numel-1

    for ii = 1:L3numel:M
        VarNamesOut(srtVec+ii,:) = VarNamesOut(srtIDXrev + (ii-1),:);
    end
end

VarNamesOut = cellstr(reshape(VarNamesOut, 1, []));


end