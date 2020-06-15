empty = zeros(4,4);
data_struct = struct("mags", struct("y", empty, "error", empty), ...
    "mags_out", struct("y", empty, "error", empty), ...,
    "acc", struct("y", empty, "error", empty), ...,
    "acc_out", struct("y", empty, "error", empty), ...,
    "pvel", struct("y", empty, "error", empty), ...,
    "pvel_out", struct("y", empty, "error", empty), ...,
    "lat", struct("y", empty, "error", empty), ...,
    "lat_out", struct("y", empty, "error", empty), ...,
    "diffs", struct("y", empty, "error", empty), ...,
    "compare", struct());

% EOG Eyelink for A B C D, one row for each condition
% This is kind of ugly but can optimise later
for i = 1:4 % Loop through conditions
    for j = 1:4 % Loop through points
        data_struct.mags.y(i, 2*(j-1)+1) = double(averages_table(4 + num_participants + (i-1)*(num_participants+4), 3+4*(j-1)));
        data_struct.mags.y(i, 2*(j-1)+2) = double(averages_table(4 + num_participants + (i-1)*(num_participants+4), 20+4*(j-1)));
        data_struct.mags.error(i, 2*(j-1)+1) = double(averages_table(5 + num_participants + (i-1)*(num_participants+4), 3+4*(j-1)));
        data_struct.mags.error(i, 2*(j-1)+2) = double(averages_table(5 + num_participants + (i-1)*(num_participants+4), 20+4*(j-1)));
        data_struct.mags_out.y(i, 2*(j-1)+1) = double(averages_table_out(4 + num_participants + (i-1)*(num_participants+4), 3+4*(j-1)));
        data_struct.mags_out.y(i, 2*(j-1)+2) = double(averages_table_out(4 + num_participants + (i-1)*(num_participants+4), 20+4*(j-1)));
        data_struct.mags_out.error(i, 2*(j-1)+1) = double(averages_table_out(5 + num_participants + (i-1)*(num_participants+4), 3+4*(j-1)));
        data_struct.mags_out.error(i, 2*(j-1)+2) = double(averages_table_out(5 + num_participants + (i-1)*(num_participants+4), 20+4*(j-1)));
        
        data_struct.acc.y(i, 2*(j-1)+1) = double(averages_table(4 + num_participants + (i-1)*(num_participants+4), 4+4*(j-1)));
        data_struct.acc.y(i, 2*(j-1)+2) = double(averages_table(4 + num_participants + (i-1)*(num_participants+4), 21+4*(j-1)));
        data_struct.acc.error(i, 2*(j-1)+1) = double(averages_table(5 + num_participants + (i-1)*(num_participants+4), 4+4*(j-1)));
        data_struct.acc.error(i, 2*(j-1)+2) = double(averages_table(5 + num_participants + (i-1)*(num_participants+4), 21+4*(j-1)));
        data_struct.acc_out.y(i, 2*(j-1)+1) = double(averages_table_out(4 + num_participants + (i-1)*(num_participants+4), 4+4*(j-1)));
        data_struct.acc_out.y(i, 2*(j-1)+2) = double(averages_table_out(4 + num_participants + (i-1)*(num_participants+4), 21+4*(j-1)));
        data_struct.acc_out.error(i, 2*(j-1)+1) = double(averages_table_out(5 + num_participants + (i-1)*(num_participants+4), 4+4*(j-1)));
        data_struct.acc_out.error(i, 2*(j-1)+2) = double(averages_table_out(5 + num_participants + (i-1)*(num_participants+4), 21+4*(j-1)));
        
        data_struct.pvel.y(i, 2*(j-1)+1) = double(averages_table(4 + num_participants + (i-1)*(num_participants+4), 5+4*(j-1)));
        data_struct.pvel.y(i, 2*(j-1)+2) = double(averages_table(4 + num_participants + (i-1)*(num_participants+4), 22+4*(j-1)));
        data_struct.pvel.error(i, 2*(j-1)+1) = double(averages_table(5 + num_participants + (i-1)*(num_participants+4), 5+4*(j-1)));
        data_struct.pvel.error(i, 2*(j-1)+2) = double(averages_table(5 + num_participants + (i-1)*(num_participants+4), 22+4*(j-1)));
        data_struct.pvel_out.y(i, 2*(j-1)+1) = double(averages_table_out(4 + num_participants + (i-1)*(num_participants+4), 5+4*(j-1)));
        data_struct.pvel_out.y(i, 2*(j-1)+2) = double(averages_table_out(4 + num_participants + (i-1)*(num_participants+4), 22+4*(j-1)));
        data_struct.pvel_out.error(i, 2*(j-1)+1) = double(averages_table_out(5 + num_participants + (i-1)*(num_participants+4), 5+4*(j-1)));
        data_struct.pvel_out.error(i, 2*(j-1)+2) = double(averages_table_out(5 + num_participants + (i-1)*(num_participants+4), 22+4*(j-1)));
        
        data_struct.lat.y(i, 2*(j-1)+1) = double(averages_table(4 + num_participants + (i-1)*(num_participants+4), 6+4*(j-1)));
        data_struct.lat.y(i, 2*(j-1)+2) = double(averages_table(4 + num_participants + (i-1)*(num_participants+4), 23+4*(j-1)));
        data_struct.lat.error(i, 2*(j-1)+1) = double(averages_table(5 + num_participants + (i-1)*(num_participants+4), 6+4*(j-1)));
        data_struct.lat.error(i, 2*(j-1)+2) = double(averages_table(5 + num_participants + (i-1)*(num_participants+4), 23+4*(j-1)));
        data_struct.lat_out.y(i, 2*(j-1)+1) = double(averages_table_out(4 + num_participants + (i-1)*(num_participants+4), 6+4*(j-1)));
        data_struct.lat_out.y(i, 2*(j-1)+2) = double(averages_table_out(4 + num_participants + (i-1)*(num_participants+4), 23+4*(j-1)));
        data_struct.lat_out.error(i, 2*(j-1)+1) = double(averages_table_out(5 + num_participants + (i-1)*(num_participants+4), 6+4*(j-1)));
        data_struct.lat_out.error(i, 2*(j-1)+2) = double(averages_table_out(5 + num_participants + (i-1)*(num_participants+4), 23+4*(j-1)));
        
        data_struct.diffs.y(i, 2*(j-1)+1) = double(averages_table(4 + num_participants + (i-1)*(num_participants+4), 36+j));
        data_struct.diffs.y(i, 2*(j-1)+2) = double(averages_table_out(4 + num_participants + (i-1)*(num_participants+4), 36+j));
        data_struct.diffs.error(i, 2*(j-1)+1) = double(averages_table(5 + num_participants + (i-1)*(num_participants+4), 36 + j));
        data_struct.diffs.error(i, 2*(j-1)+2) = double(averages_table_out(5 + num_participants + (i-1)*(num_participants+4), 36 + j));
    end
end


diff_mag = abs(data_struct.mags.y(:, 2:2:end) - data_struct.mags.y(:, 1:2:end));
diff_acc = abs(data_struct.acc.y(:, 2:2:end) - data_struct.acc.y(:, 1:2:end));
diff_mag_out = abs(data_struct.mags_out.y(:, 2:2:end) - data_struct.mags_out.y(:, 1:2:end));
diff_acc_out = abs(data_struct.acc_out.y(:, 2:2:end) - data_struct.acc_out.y(:, 1:2:end));

data_struct.compare.y = [nanmean(diff_mag, 'all'), nanmean(diff_acc, 'all');
    nanmean(diff_mag_out, 'all'), nanmean(diff_acc_out, 'all')];
data_struct.compare.error = [nanstd(diff_mag, 0, 'all'), nanstd(diff_acc, 0, 'all');
    nanstd(diff_mag_out, 0, 'all'), nanstd(diff_acc_out, 0, 'all')];

plotbars(data_struct, plot_folder_name);

