clc
clear

% Source and destination directories
%%% use for future data:
% sourceDir = 'S:\hydrolab\home\moshe\Research\wrf\PGW\functions\savedMat\singleEvents\V3\PGW_All8';
% destDir = 'S:\hydrolab\home\Raz\WRF\future';
% filePattern = fullfile(sourceDir, 'WRFrain_PGW_V3_HPE*.mat');

%%% use for historical data:
sourceDir = 'S:\hydrolab\home\moshe\Research\wrf\PGW\functions\savedMat\singleEvents\V3\Control';
destDir = 'S:\hydrolab\home\Raz\WRF\historical';
filePattern = fullfile(sourceDir, 'WRFrain_noPGW_HPE*.mat');

% List all .mat files in the source directory
fileList = dir(filePattern);

disp(['Found ', num2str(numel(fileList)), ' files to process.']);

% Iterate over each file
for fileIdx = 1:numel(fileList)
    % Construct the full file path
    filePath = fullfile(sourceDir, fileList(fileIdx).name);
    
    % Load the original file
    originalFile = load(filePath);
    
    % Extract the times from the cell array
    timesCell = originalFile.times;    %% for future use timesPGW
    
    % Preallocate a numeric array for datenum values
    timesList = zeros(size(timesCell));
    
    % Convert each datetime string to datenum
    for i = 1:numel(timesCell)
        timesList(i) = datenum(timesCell{i});
    end
    
    % Remove timesPGW field from the original structure
    originalFile = rmfield(originalFile, 'times');  %% for future use timesPGW
    
    % Add timesList as a new field to the original structure
    originalFile.timesList = timesList;
    
    % Construct the new file path in the destination directory
    [~, fileName, fileExt] = fileparts(filePath);
    newFilePath = fullfile(destDir, [fileName, fileExt]);
    
    % Save the modified structure to a new .mat file
    save(newFilePath, '-struct', 'originalFile', '-v7.3');
    
    disp(['Processed: ', fileList(fileIdx).name, ', Saved to: ', fileName, fileExt]);
end

disp('All files processed and saved.');
