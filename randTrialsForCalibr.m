% function rand_trials = randTrialsForCalibr()
%
% INPUT:
% participants = cell array of char vectors representing participant ids
%    ex. participants = {'P019','P020','P021','P022','P023','P024','P025'};
% printCSV = if true, then print a formatted CSV file of the rand trials
%
%
% OUTPUT:
% rand_trials = matrix of rand trial indices
%
%
% Generates CSV with random trial list needed for calibration (P19-25)
% by Jamie Dunkle, June 2020, UBC Vision Lab

function rand_trials = randTrialsForCalibr(participants, printCSV)

    if nargin < 1 || ~exist('participants', 'var') || isempty(participants) ||...
            ~iscellstr(participants)
        warning(['participants: no argument or an invalid argument. Defaulting'...
            ' to ''P19 - P25''.']);
        participants = {'P019','P020','P021','P022','P023','P024','P025'};
    end 
    
    if nargin < 2 || ~exist('printCSV', 'var') || isempty(printCSV) ||...
            ~islogical(printCSV)
        warning(['printCSV: no argument or an invalid argument. Defaulting'...
            ' to ''false''.']);
        printCSV = false;
    end 

    num_parts = length(participants);
    num_trials_per_letter = 10;
    
    % get all possible permutations of total num of trials
    M=perms(1:num_trials_per_letter);  %tabulate once only
    len_M = factorial(num_trials_per_letter); % number of permutations in M
    
    % for each participant, select random row of the permutation matrix
    J=M(randi(len_M, num_parts, 1), :);

    % take only the first three numbers in each permutation row,
    % this gives us three random trial indices per each participant
    rand_trials = J(:,1:3);

    if printCSV
        out_csv = '';
        for i=1:num_parts
            row_str = sprintf('%d,' , rand_trials(i,:));
            out_csv = [out_csv participants{i} ',' row_str(1:end-1) '\n'];
        end
        
        fileID = fopen('rand_calib_trials.csv','w+');
        fprintf(fileID, 'Participant,Rand1,Rand2,Rand3\n');
        fprintf(fileID, out_csv);
        fclose(fileID);
    end
end