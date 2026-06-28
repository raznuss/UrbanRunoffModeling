
clc
clear
close all

directory = 'D:\Development\RESEARCH\Raanana\data\rain_radar\Row_radar_MAT_files';
% directory = 'C:\Users\raznu\Development\RESEARCH\Raanana\data\rain_radar\Row_radar_correct_MAT_files';

%% Plots the differences bewtween time slices
matFiles = dir(fullfile(directory, '*.mat'));
numFiles = length(matFiles);
% The number of rows and columns for subplots
numRows = ceil(sqrt(numFiles));
numCols = ceil(numFiles / numRows);
% Create the big figure
bigFigure = figure;
% Iterate through each MAT file
for i = 1:numFiles
    load(fullfile(directory, matFiles(i).name));
    timeDiff = minute(diff(event_file.time));
    
    % Find the indices where time difference is bigger than 5
    missing_time = find(timeDiff > 5);
    
    % Retrieve the corresponding time values
    eventTimes = event_file.time(missing_time);
    
%     Display the time values
    disp((eventTimes));

% 
    subplot(numRows, numCols, i);
    plot(event_file.time(2:end), timeDiff);
    [~, fileName, ~] = fileparts(matFiles(i).name);
    title(fileName);
%     
    xlabel('Time');
    ylabel('Time Difference (minutes)');
end

sgtitle('Time Differences for MAT Files');


%% Create the missing rain data by doing linear interpolation

% close all
for i = 1:numFiles
    load(fullfile(directory, matFiles(i).name));
    time_step = 5; % in minutes
    rain_field = event_file.rain; % Assuming 'rain' is the variable name in the MAT file
    time_field = event_file.time; % Assuming 'time' is the variable name in the MAT file

    t1 = datetime(datestr(time_field(1)));
    t2 = datetime(datestr(time_field(end)));

    full_time_vector = t1:minutes(5):t2;
    full_time_vector_datenum = datenum(full_time_vector);
    interpolated_rain_field = zeros(112, 39, length(full_time_vector_datenum));

    for time = 1:length(time_field)
        if strcmp(matFiles(i).name, '20191227.mat'); % This file doesn't need correction;
            continue; % Skip this iteration
        end
        display(datestr(time_field(1)));
        display(matFiles(i).name);

        index = find(abs((full_time_vector_datenum - time_field(time))) < datenum(minutes(1)));
        
        % this helps to avoid empty index cause by intervals less the 5 min
        if ~isempty(index) % Check if index is not empty
            for j = 1:length(index)
                interpolated_rain_field(:,:,index(j)) = rain_field(:,:,time);
            end
        end
%         disp(index)
%         interpolated_rain_field(:,:,index) = rain_field(:,:,time);
    end

    % Required number of time steps
    % Create a new rain field copy
    time_diff = minute(diff(event_file.time)); % Compute the time differences
    missing_indices = find(time_diff > 5); % Convert 5 minutes to seconds

    % Perform linear interpolation to fill the missing rain values
    for j = 1:length(missing_indices)
        index_temp = missing_indices(j);
        start_time = time_field(index_temp);
        end_time = time_field(index_temp + 1);
        missing_time = minute(end_time - start_time);

        missing_time_steps = (missing_time / time_step) - 1;
        missing_steps_vector = ones(112, 39, ceil(missing_time_steps));
        value_interval = (rain_field(:,:,index_temp + 1) - rain_field(:,:,index_temp)) / (missing_time_steps + 1);

        a = [1:missing_time_steps]';
        B = bsxfun(@times, missing_steps_vector, reshape(a, 1, 1, [])); % X
        time_steps_to_add = rain_field(:,:,index_temp) + (B .* value_interval);

        index_in_interpolate_field = find(abs((full_time_vector_datenum - start_time)) < minutes(1));
        interpolated_rain_field(:,:,index_in_interpolate_field + 1 : (index_in_interpolate_field + missing_time_steps)) = time_steps_to_add;

    end

    % Save the file with corrected rain values
    event_file = event_file;
    event_file.rain = interpolated_rain_field;
    event_file.time = full_time_vector_datenum;
    [~, filename, ~] = fileparts(matFiles(i).name);
    save(fullfile('D:\Development\RESEARCH\Raanana\data\rain_radar\Row_radar_correct_MAT_files_SA_storms', [filename]), 'event_file');

end


% The interpolated_rain_field now contains the rain field with the missing values filled using interpolation



