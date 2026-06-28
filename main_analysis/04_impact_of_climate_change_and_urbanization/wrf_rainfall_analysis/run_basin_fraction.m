% Run basin and domain wet-area fraction calculations
clc ; clear; 
% ----- paths -----
shpPath       = 'D:\Development\RESEARCH\Raanana\gis\GIS\watershead_for_radar\Raanana_wsp_merged_wgs84.shp';
eventsDirHist = '\\vscifs\hydrolab1\hydrolab\home\moshe\Research\wrf\PGW\functions\savedMat\singleEvents\V3\Control\';
eventsDirPGW  = '\\vscifs\hydrolab1\hydrolab\home\moshe\Research\wrf\PGW\functions\savedMat\singleEvents\V3\PGW_All8\';

% ----- thresholds -----
% thresholds = [0.1 5 80];   % mm/hr, adjust as needed
thresholds = [0.1 1 2 5 7.5 10 20 50 60 80 100 120 150];

% ----- options -----
save_dir_basin  = fullfile(pwd, 'basin_wet_fraction');
save_dir_domain = fullfile(pwd, 'domain_wet_fraction');
events_idx      = [];    % [] -> all events
verbose         = true;
show_waitbar    = true;

% % ----- basin wet fraction -----
% fprintf('=== Computing BASIN wet fraction ===\n');
% results_basin = compute_basin_wet_fraction( ...
%     thresholds, shpPath, eventsDirHist, eventsDirPGW, ...
%     save_dir_basin, events_idx, verbose, show_waitbar);

% ----- domain wet fraction -----
fprintf('\n=== Computing DOMAIN wet fraction ===\n');
results_domain = compute_domain_wet_fraction( ...
    thresholds, eventsDirHist, eventsDirPGW, ...
    save_dir_domain, events_idx, verbose, show_waitbar);

fprintf('\nAll done.\n');
