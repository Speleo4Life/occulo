function plotbars(data, FolderName)
% Create Plots folder
if ~exist(FolderName, 'dir');
   mkdir(FolderName);
end
    
%%% With outliers 

% Magnitude 
figure(1);
y = data.mags.y;
y_error = data.mags.error;

bar(y,'grouped')
hold on
%errorbar(y,y_error)

% Find the number of groups and the number of bars in each group
ngroups = size(y, 1);
nbars = size(y_error, 2);

% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), y_error(:,i), 'k', 'linestyle', 'none');
end

set(gca,'xticklabel',{'C01','C02','C03','C04'});
yline(-22,'-.k', 'A : -22 deg.');
yline(-11,'-.k', 'B : -11 deg.');
yline(22,'-.k', 'D : 22 deg.');
yline(11,'-.k', 'C : 11 deg.');

title('Variation of sccade magnitudes in each condition : with outliers');
xlabel('Experiment condition'); 
ylabel('Saccade magnitude in degrees'); 

hold off 

legend({'A - EOG','A - EL','B - EOG','B - EL','C - EOG','C - EL','D - EOG','D - EL'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accuracy 
figure(2);

y = data.acc.y;


bar(y,'grouped')
hold on


set(gca,'xticklabel',{'C01','C02','C03','C04'});
yline(10,'-.k', '10');
yline(15,'-.k', '15');
yline(18,'-.k', '18');
yline(19,'-.k', '19');

title('Variation of average sccadic error in each condition : with outliers');
xlabel('Experiment condition'); 
ylabel('Average error'); 

hold off 
legend({'A - EOG','A - EL','B - EOG','B - EL','C - EOG','C - EL','D - EOG','D - EL'});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Peak Velocity
figure(3);

y = data.pvel.y;
y_error = data.pvel.error;

bar(y,'grouped')
hold on
%errorbar(y,y_error)

% Find the number of groups and the number of bars in each group
ngroups = size(y, 1);
nbars = size(y_error, 2);

% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), y_error(:,i), 'k', 'linestyle', 'none');
end

set(gca,'xticklabel',{'C01','C02','C03','C04'});

title('Variation of sccade peak velocities in each condition : with outliers');
xlabel('Experiment condition'); 
ylabel('Saccade peak velosity in degrees per seconds'); 

hold off 

legend({'A - EOG','A - EL','B - EOG','B - EL','C - EOG','C - EL','D - EOG','D - EL'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Latency
figure(4);

y = data.lat.y;
y_error = data.lat.error;

bar(y,'grouped')
hold on
%errorbar(y,y_error)

% Find the number of groups and the number of bars in each group
ngroups = size(y, 1);
nbars = size(y_error, 2);

% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), y_error(:,i), 'k', 'linestyle', 'none');
end

set(gca,'xticklabel',{'C01','C02','C03','C04'});

title('Variation of sccade latency in each condition : with outliers');
xlabel('Experiment condition'); 
ylabel('Saccade latency in milliseconds'); 

hold off 

legend({'A - EOG','A - EL','B - EOG','B - EL','C - EOG','C - EL','D - EOG','D - EL'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Without outliers 

% Magnitude 
figure(5);
y = data.mags_out.y;
y_error = data.mags_out.error;

bar(y,'grouped')
hold on
%errorbar(y,y_error)

% Find the number of groups and the number of bars in each group
ngroups = size(y, 1);
nbars = size(y_error, 2);

% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), y_error(:,i), 'k', 'linestyle', 'none');
end

set(gca,'xticklabel',{'C01','C02','C03','C04'});
yline(-22,'-.k', 'A : -22 deg.');
yline(-11,'-.k', 'B : -11 deg.');
yline(22,'-.k', 'D : 22 deg.');
yline(11,'-.k', 'C : 11 deg.');

title('Variation of sccade magnitudes in each condition : without outliers');
xlabel('Experiment condition'); 
ylabel('Saccade magnitude in degrees'); 

hold off 

legend({'A - EOG','A - EL','B - EOG','B - EL','C - EOG','C - EL','D - EOG','D - EL'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accuracy 
figure(6);

y = data.acc_out.y;

bar(y,'grouped')
hold on


set(gca,'xticklabel',{'C01','C02','C03','C04'});
yline(10,'-.k', '10');
yline(15,'-.k', '15');
yline(18,'-.k', '18');
yline(19,'-.k', '19');

title('Variation of average sccadic error in each condition : without outliers');
xlabel('Experiment condition'); 
ylabel('Average error'); 

hold off 
legend({'A - EOG','A - EL','B - EOG','B - EL','C - EOG','C - EL','D - EOG','D - EL'});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Peak Velocity
figure(7);

y = data.pvel_out.y;
y_error = data.pvel_out.error;

bar(y,'grouped')
hold on
%errorbar(y,y_error)

% Find the number of groups and the number of bars in each group
ngroups = size(y, 1);
nbars = size(y_error, 2);

% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), y_error(:,i), 'k', 'linestyle', 'none');
end

set(gca,'xticklabel',{'C01','C02','C03','C04'});

title('Variation of sccade peak velocities in each condition : without outliers');
xlabel('Experiment condition'); 
ylabel('Saccade peak velosity in degrees per seconds'); 

hold off 

legend({'A - EOG','A - EL','B - EOG','B - EL','C - EOG','C - EL','D - EOG','D - EL'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Latency
figure(8);

y = data.lat_out.y;
y_error = data.lat_out.error;


bar(y,'grouped')
hold on
%errorbar(y,y_error)

% Find the number of groups and the number of bars in each group
ngroups = size(y, 1);
nbars = size(y_error, 2);

% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), y_error(:,i), 'k', 'linestyle', 'none');
end

set(gca,'xticklabel',{'C01','C02','C03','C04'});

title('Variation of sccade latency in each condition : without outliers');
xlabel('Experiment condition'); 
ylabel('Saccade latency in milliseconds'); 

hold off 

legend({'A - EOG','A - EL','B - EOG','B - EL','C - EOG','C - EL','D - EOG','D - EL'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Improvement study : difference between EOG and EL plots 
figure(9);

y = data.compare.y;
y_error = data.compare.error;


bar(y,'grouped')
hold on
%errorbar(y,y_error)

% Find the number of groups and the number of bars in each group
ngroups = size(y, 1);
nbars = size(y_error, 2);

% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), y_error(:,i), 'k', 'linestyle', 'none');
end

set(gca,'xticklabel',{'With outlers','Without outlers'});

title('Improvement study : difference between EOG and EL plots');
xlabel('Condition'); 
ylabel('Differece between EOG and EL'); 

hold off 
legend({'Magnitude','Accuracy'});

figure(10);

y = data.diffs.y;
y_error = data.diffs.error;


bar(y,'grouped')
hold on
%errorbar(y,y_error)

% Find the number of groups and the number of bars in each group
ngroups = size(y, 1);
nbars = size(y_error, 2);

% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), y_error(:,i), 'k', 'linestyle', 'none');
end

set(gca,'xticklabel',{'C01','C02','C03','C04'});

title('Variation of difference between Elink and Eog');
xlabel('Experiment condition'); 
ylabel('Difference in angle'); 

hold off 

legend({'A','A - w/o outliers','B','B - w/o outliers','C' ,'C - w/o outliers'});


FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  savefig(fullfile(FolderName, [FigName '.fig']));
end



