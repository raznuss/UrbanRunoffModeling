function results = wrf_basin_analysis(verbose, show_waitbar, save_dir, basin_shp, events_idx, samples_per_side)
% Basin focused rainfall analysis over a single basin polygon in WGS84.
% Compares HIST vs PGW using independent maxima inside the basin
% for every box size including a single pixel.
% CSV column order is unchanged. Labels in block_size are:
%   'Raanana Basin'
%   'Raanana Basin X pixel'  where X is the box side length in pixels

if nargin < 1, verbose = true; end
if nargin < 2, show_waitbar = true; end
if nargin < 3 || isempty(save_dir), save_dir = fullfile(pwd,'wrf_basin_analysis'); end
if nargin < 4 || isempty(basin_shp)
    basin_shp = 'D:\Development\RESEARCH\Raanana\gis\GIS\watershead_for_radar\Raanana_wsp_merged_wgs84.shp';
end
if nargin < 5, events_idx = [ ]; end
if nargin < 6 || isempty(samples_per_side), samples_per_side = 7; end
if ~isfolder(save_dir), mkdir(save_dir); end

% Configuration
eventsDirHist = '\\vscifs\hydrolab1\hydrolab\home\moshe\Research\wrf\PGW\functions\savedMat\singleEvents\V3\Control\';
eventsDirPGW  = '\\vscifs\hydrolab1\hydrolab\home\moshe\Research\wrf\PGW\functions\savedMat\singleEvents\V3\PGW_All8\';
durations_min = [10 20 30 60 120 240 360 720 1440];

box_sizes      = [1 2 3];       % compute for these box sizes, including single pixel
overlap_thresh = 0.0;           % a grid node is in basin mask if W_basin > overlap_thresh
basin_label    = "Raanana Basin";

% Discover events
files_hist_all = basin_list_mat_files(eventsDirHist);
files_pgw_all  = basin_list_mat_files(eventsDirPGW);
nEventsAll = min(numel(files_hist_all), numel(files_pgw_all));

if isempty(events_idx)
    events_idx = 1:nEventsAll;
else
    events_idx = unique(events_idx(:)');
    if any(events_idx < 1 | events_idx > nEventsAll)
        error('events_idx must be within 1..%d', nEventsAll);
    end
end

files_hist = files_hist_all(events_idx);
files_pgw  = files_pgw_all(events_idx);

nEvents = numel(events_idx);
nDur    = numel(durations_min);

% Read basin geometry once
try
    GT = readgeotable(basin_shp);
    basinShape = GT.Shape(1);
    if any(strcmpi(GT.Properties.VariableNames,'name'))
        basin_name_meta = string(GT.name(1));
    elseif any(strcmpi(GT.Properties.VariableNames,'Name'))
        basin_name_meta = string(GT.Name(1));
    else
        basin_name_meta = "BASIN";
    end
catch ME
    error('Failed to read basin shapefile: %s', ME.message);
end

% Outputs for classic basin averages
I_hist = nan(nEvents, nDur);   % basin average intensity
I_pgw  = nan(nEvents, nDur);
T_hist = nan(nEvents, 1);      % basin average total
T_pgw  = nan(nEvents, 1);

hist_kept = cell(0,1);
pgw_kept  = cell(0,1);

W_basin    = [];   % area fraction weights
mask_basin = [];   % logical mask of nodes that belong to the basin

rows = {};         % accumulates CSV rows

if show_waitbar, h = waitbar(0,'Processing events...'); end

% Event loop
for j = 1:nEvents
    src_idx = events_idx(j);
    if verbose
        fprintf('Event %02d/%02d (source #%d)\n', j, nEvents, src_idx);
        fprintf('  Hist file: %s\n', files_hist{j});
        fprintf('  PGW  file: %s\n', files_pgw{j});
    end
    t0 = tic;

    % Load rain cubes and time vectors
    try
        [R_hist, t_hist, lat, lon] = basin_load_event_with_latlon(fullfile(eventsDirHist, files_hist{j}));
        [R_pgw , t_pgw , ~  ,  ~ ] = basin_load_event_with_latlon(fullfile(eventsDirPGW , files_pgw{j}));
    catch ME
        warning('Skipping event %d (src %d): %s', j, src_idx, ME.message);
        if show_waitbar, waitbar(j/nEvents, h, sprintf('Event %d skipped', j)); end
        continue
    end

    % Drop the first time step for parity with previous workflows
    if size(R_hist,3) > 0, R_hist(:,:,1) = []; t_hist(1) = []; end
    if size(R_pgw ,3) > 0, R_pgw (:,:,1) = []; t_pgw(1)  = []; end

    % Align by overlapping timestamps
    [R_hist_a, R_pgw_a, t_aligned] = basin_align_time_series_3d(R_hist, t_hist, R_pgw, t_pgw);
    if isempty(t_aligned)
        warning('Event %d (src %d): no overlapping timestamps', j, src_idx);
        if show_waitbar, waitbar(j/nEvents, h, sprintf('Event %d no overlap', j)); end
        continue
    end

    % Build basin weights and mask once using the first processed event
    if isempty(W_basin)
        if isempty(lat) || isempty(lon)
            error('lat or lon grid missing in MAT files.');
        end
        [W_basin, ~] = basin_build_weights(basinShape, lat, lon, samples_per_side, verbose);
        mask_basin = W_basin > overlap_thresh;
        try
            basin_sanity_plot(basinShape, lat, lon, save_dir);
        catch MEp
            warning('Sanity plot failed: %s', MEp.message);
        end
    end

    % Per pixel event totals [mm]
    dt_min = basin_infer_dt_minutes(t_aligned);
    dt_hr  = dt_min/60;
    H_hist = sum(R_hist_a, 3) * dt_hr;   % [ny x nx]
    H_pgw  = sum(R_pgw_a , 3) * dt_hr;

    % Basin average totals via area weighted mean
    T_hist(j) = sum(W_basin .* H_hist, 'all');
    T_pgw (j) = sum(W_basin .* H_pgw , 'all');

    % CSV row: basin average total
    rows(end+1,:) = { ...
        j, 'total', 0, char(basin_label), false, ...
        files_hist{j}, files_pgw{j}, ...
        T_hist(j), T_pgw(j), T_pgw(j)-T_hist(j), T_pgw(j)./T_hist(j), 'mm'};

    % Independent maxima for totals for each box size including single pixel
    for bx = box_sizes
        [hist_box_max_total, pgw_box_max_total] = max_box_total_both(H_hist, H_pgw, mask_basin, bx);
        rows(end+1,:) = { ...
            j, 'total', 0, sprintf('%s %d pixel', basin_label, bx), false, ...
            files_hist{j}, files_pgw{j}, ...
            hist_box_max_total, pgw_box_max_total, ...
            pgw_box_max_total - hist_box_max_total, pgw_box_max_total ./ hist_box_max_total, 'mm'};
    end

    % Intensities for each duration D
    for d = 1:nDur
        D = durations_min(d);
        w = max(1, round(D / dt_min));
        Kt = ones(1,1,w);

        % Temporal means over full windows only
        M_hist = convn(R_hist_a, Kt, 'valid') / w;   % [ny x nx x Tvalid]
        M_pgw  = convn(R_pgw_a , Kt, 'valid') / w;
        Tvalid = size(M_hist,3);

        % Basin average intensity: area weighted mean at each time then max over time
        if Tvalid == 0
            I_hist(j,d) = NaN; I_pgw(j,d) = NaN;
        else
            zH = -Inf; zP = -Inf;
            for tt = 1:Tvalid
                zH = max(zH, sum(W_basin .* M_hist(:,:,tt), 'all'));
                zP = max(zP, sum(W_basin .* M_pgw (:,:,tt), 'all'));
            end
            I_hist(j,d) = zH; I_pgw(j,d) = zP;
        end
        rows(end+1,:) = { ...
            j, 'intensity', D, char(basin_label), false, ...
            files_hist{j}, files_pgw{j}, ...
            I_hist(j,d), I_pgw(j,d), I_pgw(j,d) - I_hist(j,d), I_pgw(j,d) ./ I_hist(j,d), 'mm_per_hr'};

        % Independent maxima for intensity for each box size including single pixel
        for bx = box_sizes
            if Tvalid == 0
                hist_box_max_int = NaN; pgw_box_max_int = NaN;
            else
                [hist_box_max_int, pgw_box_max_int] = max_box_intensity_both(M_hist, M_pgw, mask_basin, bx);
            end
            rows(end+1,:) = { ...
                j, 'intensity', D, sprintf('%s %d pixel', basin_label, bx), false, ...
                files_hist{j}, files_pgw{j}, ...
                hist_box_max_int, pgw_box_max_int, ...
                pgw_box_max_int - hist_box_max_int, pgw_box_max_int ./ hist_box_max_int, 'mm_per_hr'};
        end
    end

    hist_kept{end+1,1} = files_hist{j}; %#ok<AGROW>
    pgw_kept {end+1,1} = files_pgw{j};  %#ok<AGROW>

    if verbose, fprintf('  Done event %02d in %.1fs\n', j, toc(t0)); end
    if show_waitbar, waitbar(j/nEvents, h, sprintf('Event %d of %d', j, nEvents)); end
end

if show_waitbar, close(h); end

% Keep only processed events for classic outputs
nKept    = numel(hist_kept);
I_hist   = I_hist(1:nKept,:);
I_pgw    = I_pgw (1:nKept,:);
T_hist   = T_hist(1:nKept);
T_pgw    = T_pgw (1:nKept);

% Pack and save
results.durations_min    = durations_min;
results.basin_name       = basin_name_meta;
results.W_basin          = W_basin;
results.mask_basin       = mask_basin;
results.intensity_hist   = I_hist;
results.intensity_pgw    = I_pgw;
results.total_hist_mm    = T_hist;
results.total_pgw_mm     = T_pgw;
results.event_hist_files = hist_kept;
results.event_pgw_files  = pgw_kept;
results.box_sizes        = box_sizes;
results.overlap_thresh   = overlap_thresh;

mat_path = fullfile(save_dir, 'basin_stats_results.mat');
save(mat_path, 'results');

csv_path = basin_write_results_csv(save_dir, rows);

if verbose
    fprintf('Saved MAT: %s\n', mat_path);
    fprintf('Saved CSV: %s\n', csv_path);
end
end

% ===== Helpers =====

function files = basin_list_mat_files(dir_path)
% Return sorted list of .mat files in a folder
d = dir(fullfile(dir_path, '*.mat'));
d = d(~[d.isdir]);
files = sort({d.name});
end

function [R, tvec, lat, lon] = basin_load_event_with_latlon(fullpath_mat)
% Load rain cube and time vector and grid nodes if present
S = load(fullpath_mat);
lat = []; lon = [];
% rain
if isfield(S, 'rainRate')
    R = S.rainRate;
elseif isfield(S, 'rainRatePGW')
    R = S.rainRatePGW;
else
    error('No rain rate variable in %s', fullpath_mat);
end
% time
if isfield(S, 'times')
    tvec = basin_normalize_times(S.times);
elseif isfield(S, 'timesPGW')
    tvec = basin_normalize_times(S.timesPGW);
elseif isfield(S, 'time')
    tvec = basin_normalize_times(S.time);
else
    error('No time variable in %s', fullpath_mat);
end
% grid
if isfield(S, 'lat'), lat = double(S.lat); end
if isfield(S, 'lon'), lon = double(S.lon); end
end

function tvec = basin_normalize_times(times_in)
% Normalize time input into a column datetime vector
if iscell(times_in)
    if ~isempty(times_in) && isa(times_in{1}, 'datetime')
        tvec = vertcat(times_in{:});
    else
        tvec = datetime(cell2mat(times_in(:)), 'ConvertFrom', 'datenum');
    end
elseif isa(times_in, 'datetime')
    tvec = times_in(:);
elseif isnumeric(times_in)
    tvec = datetime(times_in(:), 'ConvertFrom', 'datenum');
else
    error('Unsupported time data type');
end
end

function [R1a, R2a, t_common] = basin_align_time_series_3d(R1, t1, R2, t2)
% Intersect timestamps and align the two rain cubes accordingly
[t_common, ia, ib] = intersect(t1, t2, 'stable');
R1a = R1(:,:,ia);
R2a = R2(:,:,ib);
end

function dt_min = basin_infer_dt_minutes(tvec)
% Infer time step in minutes from a datetime vector
if numel(tvec) < 2
    error('Cannot infer time step from fewer than two timestamps');
end
dt = median(diff(tvec));
dt_min = minutes(dt);
if dt_min <= 0
    error('Non positive time step inferred');
end
end

function [W_norm, info] = basin_build_weights(basinShape, lat, lon, s, verbose)
% Build area fraction weights W on the model grid using uniform subcell sampling
if s < 3, s = 3; end
[ny, nx] = size(lat);
if ~isequal(size(lon), [ny nx]), error('lat lon size mismatch'); end

W = zeros(ny, nx);
u = (0.5:s-0.5)/s; v = u; [U, V] = meshgrid(u, v);

% Limit computations to the basin bounding box for speed
try
    [latlim, lonlim] = geobounds(basinShape);
    latmin = latlim(1); latmax = latlim(2);
    lonmin = lonlim(1); lonmax = lonlim(2);
catch
    latmin = -Inf; latmax = Inf; lonmin = -Inf; lonmax = Inf;
end

for i = 1:(ny-1)
    for j = 1:(nx-1)
        % fast reject using cell bounding box
        lat_cell = [lat(i,j) lat(i,j+1) lat(i+1,j) lat(i+1,j+1)];
        lon_cell = [lon(i,j) lon(i,j+1) lon(i+1,j) lon(i+1,j+1)];
        if max(lat_cell) < latmin || min(lat_cell) > latmax || ...
           max(lon_cell) < lonmin || min(lon_cell) > lonmax
            continue
        end

        % four corners
        lat00 = lat(i  , j  ); lon00 = lon(i  , j  );
        lat01 = lat(i  , j+1); lon01 = lon(i  , j+1);
        lat10 = lat(i+1, j  ); lon10 = lon(i+1, j  );
        lat11 = lat(i+1, j+1); lon11 = lon(i+1, j+1);

        % bilinear interpolation inside the cell
        lat_uv = (1-U).*(1-V).*lat00 + U.*(1-V).*lat01 + (1-U).*V.*lat10 + U.*V.*lat11;
        lon_uv = (1-U).*(1-V).*lon00 + U.*(1-V).*lon01 + (1-U).*V.*lon10 + U.*V.*lon11;

        % fraction of sub points inside the basin
        inside = isinterior(basinShape, geopointshape(lat_uv(:), lon_uv(:)));
        alpha  = mean(inside);

        % distribute the fraction to the four corner nodes
        W(i  , j  ) = W(i  , j  ) + 0.25*alpha;
        W(i  , j+1) = W(i  , j+1) + 0.25*alpha;
        W(i+1, j  ) = W(i+1, j  ) + 0.25*alpha;
        W(i+1, j+1) = W(i+1, j+1) + 0.25*alpha;
    end
end

Wsum = sum(W, 'all');
if Wsum <= 0
    error('Basin polygon does not overlap the model grid.');
end
W_norm = W / Wsum;

info.s    = s;
info.sumW = sum(W_norm,'all');
info.nnz  = nnz(W_norm>0);
if verbose
    fprintf('  Basin weights built: %d grid points with positive weight. Sum(W)=%.6f\n', info.nnz, info.sumW);
end
end

function basin_sanity_plot(basinShape, lat, lon, save_dir)
% Quick sanity plot of basin outline versus model domain box
f = figure('Color','w','Units','normalized','Position',[0.15 0.15 0.6 0.6]);
gx = geoaxes(f); hold(gx,'on');
try
    L = shaperead('landareas','UseGeoCoords',true);
    for k = 1:numel(L)
        geoplot(gx, L(k).Lat, L(k).Lon, 'Color',[0.7 0.7 0.7], 'LineWidth', 0.4);
    end
catch
end
latmin = min(lat(:)); latmax = max(lat(:));
lonmin = min(lon(:)); lonmax = max(lon(:));
geoplot(gx, [latmin latmin latmax latmax latmin], [lonmin lonmax lonmax lonmin lonmin], 'k--', 'LineWidth', 1.2);
geoplot(gx, basinShape, 'EdgeColor', [0.85 0.2 0.2], 'FaceColor','none', 'LineWidth', 1.8);
geolimits(gx, [latmin latmax], [lonmin lonmax]);
geobasemap(gx, 'streets-light');
title(gx, 'Basin vs WRF domain');
saveas(f, fullfile(save_dir,'basin_sanity_plot.png'));
close(f);
end

function [hist_max, pgw_max] = max_box_total_both(H_hist, H_pgw, mask_basin, bx)
% Independent box mean maxima for totals for HIST and for PGW
% Consider only boxes that touch the basin mask
K = ones(bx)/(bx*bx);
S_hist = conv2(H_hist, K, 'valid');   % map of box means
S_pgw  = conv2(H_pgw , K, 'valid');

C = conv2(double(mask_basin), ones(bx), 'valid') > 0;  % candidate boxes that touch the basin
S_hist(~C) = -Inf;
S_pgw (~C) = -Inf;

hist_max = max(S_hist(:));
pgw_max  = max(S_pgw (:));
if ~isfinite(hist_max), hist_max = NaN; end
if ~isfinite(pgw_max ), pgw_max  = NaN; end
end

function [hist_max, pgw_max] = max_box_intensity_both(M_hist, M_pgw, mask_basin, bx)
% Independent box mean maxima over time for intensities for HIST and for PGW
% Consider only boxes that touch the basin mask
ny = size(M_hist,1); nx = size(M_hist,2); Tvalid = size(M_hist,3);
Kb = ones(bx)/(bx*bx);
C  = conv2(double(mask_basin), ones(bx), 'valid') > 0;

Qhist = -Inf(ny-bx+1, nx-bx+1);  % max over time of box means for HIST
Qpgw  = -Inf(ny-bx+1, nx-bx+1);  % max over time of box means for PGW

for tt = 1:Tvalid
    Sh = conv2(M_hist(:,:,tt), Kb, 'valid');
    Sp = conv2(M_pgw (:,:,tt), Kb, 'valid');
    Sh(~C) = -Inf;
    Sp(~C) = -Inf;
    Qhist = max(Qhist, Sh);
    Qpgw  = max(Qpgw , Sp);
end

hist_max = max(Qhist(:));
pgw_max  = max(Qpgw (:));
if ~isfinite(hist_max), hist_max = NaN; end
if ~isfinite(pgw_max ), pgw_max  = NaN; end
end

function csv_path = basin_write_results_csv(save_dir, rows)
% Write rows to a single CSV with fixed column order
T = cell2table(rows, 'VariableNames', { ...
    'event_id','measure_type','duration_min','block_size','is_full_domain', ...
    'hist_file','pgw_file','hist_value','pgw_value','diff','ratio','units'});
csv_path = fullfile(save_dir, 'basin_stats_events.csv');
writetable(T, csv_path);
end
