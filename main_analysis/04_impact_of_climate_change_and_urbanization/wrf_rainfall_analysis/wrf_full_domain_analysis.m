    function results = compute_unconditional_blockmax_fast(verbose, show_waitbar, save_dir, events_idx)
% Regional block maximum of UNCONDITIONAL rain rate over multiple durations,
% plus UNCONDITIONAL totals per spatial block size (all in ONE output CSV).
%
% Assumes rainRate cubes have finite non-negative values.
%
% Intensities (per duration D, per block size b):
%   time mean over full windows (convn 'valid') -> spatial block-mean -> spatial max per time -> max over time.
%   If an event is shorter than D (no full windows), intensity values are left as NaN and STILL written to CSV.
%
% Totals (per block size b):
%   per-pixel time accumulation [mm] over the event -> spatial block-mean -> spatial max.
%
% inputs:
%   verbose      (logical) print progress (default true)
%   show_waitbar (logical) show waitbar (default true)
%   save_dir     (char/str) output folder (default <pwd>/wrf_rainfall_analysis)
%   events_idx   (vector) 1-based indices of events to run; [] or omitted => run all
%
% results = compute_unconditional_blockmax_fast(verbose, show_waitbar, save_dir, events_idx)

if nargin < 1, verbose = true; end
if nargin < 2, show_waitbar = true; end
if nargin < 3 || isempty(save_dir), save_dir = fullfile(pwd,'wrf_rainfall_analysis'); end
if nargin < 4, events_idx = []; end
if ~isfolder(save_dir), mkdir(save_dir); end

% -----------------------------
% Settings
% -----------------------------
eventsDirHist = '\\vscifs\hydrolab1\hydrolab\home\moshe\Research\wrf\PGW\functions\savedMat\singleEvents\V3\Control\';
eventsDirPGW  = '\\vscifs\hydrolab1\hydrolab\home\moshe\Research\wrf\PGW\functions\savedMat\singleEvents\V3\PGW_All8\';

% Temporal windows [minutes]  (FULL WINDOWS ONLY)
durations_min = [10 20 30 60 120 240 360 720 1440];

% Spatial block sizes [pixels] (square neighborhoods)
block_sizes = [1 3 6 10 25 50 100 150 200 250 350];

% Include a full domain block (whole grid at once; last column with block_size=0)
include_full_domain_block = true;

% Spatial block mode: true -> 'valid' (only full neighborhoods), false -> 'same'
spatial_mode_valid = true;

% -----------------------------
% Discover events (+ optional subset)
% -----------------------------
files_hist_all = list_mat_files(eventsDirHist);
files_pgw_all  = list_mat_files(eventsDirPGW);

nEventsAll = min(numel(files_hist_all), numel(files_pgw_all));
if isempty(events_idx)
    events_idx = 1:nEventsAll;                 % run all
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
nB0     = numel(block_sizes);
nB      = nB0 + double(include_full_domain_block);  % last column is FULL domain if enabled

% Preallocate outputs
RM_hist  = nan(nEvents, nDur, nB);     % intensities [events x dur x block] mm/hr
RM_pgw   = nan(nEvents, nDur, nB);
TB_hist  = nan(nEvents,        nB);    % totals per block [events x block] mm
TB_pgw   = nan(nEvents,        nB);

% Keep file names
hist_kept = cell(0,1);
pgw_kept  = cell(0,1);

if show_waitbar, h = waitbar(0,'Processing events...'); end

% -----------------------------
% Loop over events
% -----------------------------
for j = 1:nEvents
    src_idx = events_idx(j);
    if verbose
        fprintf('Event %02d/%02d (source #%d)\n', j, nEvents, src_idx);
        fprintf('  Hist file: %s\n', files_hist{j});
        fprintf('  PGW  file: %s\n', files_pgw{j});
    end
    t_event = tic;

    % Load rain and time
    try
        [R_hist, t_hist] = load_event_any(fullfile(eventsDirHist, files_hist{j}));
        [R_pgw , t_pgw ] = load_event_any(fullfile(eventsDirPGW , files_pgw{j}));
    catch ME
        warning('Skipping event %d (source %d): %s', j, src_idx, ME.message);
        if show_waitbar, waitbar(j/nEvents, h, sprintf('Event %d skipped', j)); end
        continue
    end

    % Drop first step for parity with earlier workflow
    if size(R_hist,3) > 0, R_hist(:,:,1) = []; t_hist(1) = []; end
    if size(R_pgw ,3) > 0, R_pgw (:,:,1) = []; t_pgw(1)  = []; end

    % Align by common timestamps
    [R_hist_a, R_pgw_a, t_aligned] = align_time_series_3d(R_hist, t_hist, R_pgw, t_pgw);
    if isempty(t_aligned)
        warning('Skipping event %d (source %d): no overlapping timestamps', j, src_idx);
        if show_waitbar, waitbar(j/nEvents, h, sprintf('Event %d no overlap', j)); end
        continue
    end

    % Event total depth per pixel [mm]
    dt_min = infer_dt_minutes(t_aligned);
    dt_hr  = dt_min/60;
    H_hist = sum(R_hist_a, 3) * dt_hr;   % [ny x nx] mm
    H_pgw  = sum(R_pgw_a , 3) * dt_hr;

    % Totals per spatial block (block-mean of H, then spatial max)
    [ny, nx, ~] = size(R_hist_a);
    for b = 1:nB0
        bs = block_sizes(b);
        if spatial_mode_valid && (ny < bs || nx < bs)
            TB_hist(j,b) = NaN;
            TB_pgw (j,b) = NaN;
        else
            TB_hist(j,b) = block_spatial_mean_max(H_hist, bs, spatial_mode_valid);
            TB_pgw (j,b) = block_spatial_mean_max(H_pgw , bs, spatial_mode_valid);
        end
    end
    if include_full_domain_block
        bIdx = nB0 + 1;
        TB_hist(j,bIdx) = mean(H_hist, 'all');   % FULL domain mean total [mm]
        TB_pgw (j,bIdx) = mean(H_pgw , 'all');
    end

    % Intensities per duration with FULL WINDOWS ONLY (convn 'valid')
    for d = 1:nDur
        w = max(1, round(durations_min(d) / dt_min));  % samples per window

        K = ones(1,1,w);                               % time box filter
        M_hist = convn(R_hist_a, K, 'valid') / w;      % [ny x nx x (T-w+1)], mm/hr
        M_pgw  = convn(R_pgw_a , K, 'valid') / w;

        % If no full windows (event shorter than duration): leave NaN and continue,
        % but DO NOT crash (rows will still be written as NaN in the CSV).
        Tvalid = size(M_hist,3);
        if Tvalid == 0
            if verbose
                fprintf('    Duration %d min: no full windows (event shorter than D) -> writing NaN\n', durations_min(d));
            end
            % Values in RM_* remain NaN as preallocated
            continue
        end

        % Spatial sliding block-mean then spatial max-over-time
        for b = 1:nB0
            bs = block_sizes(b);
            if spatial_mode_valid && (ny < bs || nx < bs)
                RM_hist(j,d,b) = NaN;
                RM_pgw (j,d,b) = NaN;
            else
                RM_hist(j,d,b) = block_regional_max_over_time_fast(M_hist, bs, spatial_mode_valid);
                RM_pgw (j,d,b) = block_regional_max_over_time_fast(M_pgw , bs, spatial_mode_valid);
            end
        end

        % FULL domain intensity: mean over field per time, then max over time
        if include_full_domain_block
            bIdx = nB0 + 1;
            Z_hist = squeeze(mean(mean(M_hist,1),2));   % [Tvalid x 1]
            Z_pgw  = squeeze(mean(mean(M_pgw ,1),2));
            RM_hist(j,d,bIdx) = max(Z_hist);
            RM_pgw (j,d,bIdx) = max(Z_pgw );
        end
    end

    hist_kept{end+1,1} = files_hist{j}; %#ok<AGROW>
    pgw_kept {end+1,1} = files_pgw{j};  %#ok<AGROW>

    if verbose
        fprintf('  Done event %02d in %.1fs\n', j, toc(t_event));
    end
    if show_waitbar, waitbar(j/nEvents, h, sprintf('Event %d of %d', j, nEvents)); end
end

if show_waitbar, close(h); end

% Keep only processed events
nKept   = numel(hist_kept);
RM_hist = RM_hist(1:nKept,:,:);
RM_pgw  = RM_pgw (1:nKept,:,:);
TB_hist = TB_hist(1:nKept,:);
TB_pgw  = TB_pgw (1:nKept,:);

% Labels (append FULL-domain as block_size==0)
block_sizes_out = block_sizes(:);
block_is_full   = false(nB,1);
if include_full_domain_block
    block_sizes_out = [block_sizes_out; 0];
    block_is_full(end) = true;
end

% -----------------------------
% Pack results (MAT for programmatic use)
% -----------------------------
results.durations_min     = durations_min;
results.block_sizes       = block_sizes_out;   % last (==0) means FULL domain
results.block_is_full     = block_is_full;

results.RM_hist           = RM_hist;           % [events x dur x block], mm/hr
results.RM_pgw            = RM_pgw;            % [events x dur x block], mm/hr
results.TOTblock_hist_mm  = TB_hist;           % [events x block], mm
results.TOTblock_pgw_mm   = TB_pgw;            % [events x block], mm

results.event_hist_files  = hist_kept;
results.event_pgw_files   = pgw_kept;

% -----------------------------
% Save outputs (ONE CSV + MAT)
% -----------------------------
mat_path = fullfile(save_dir, 'block_stats_results.mat');
save(mat_path, 'results');

csv_path = write_all_results_csv(save_dir, durations_min, block_sizes_out, block_is_full, ...
    RM_hist, RM_pgw, TB_hist, TB_pgw, hist_kept, pgw_kept);

if verbose
    fprintf('Saved MAT:  %s\n', mat_path);
    fprintf('Saved CSV:  %s\n', csv_path);
end
end

% ============================================================
% Helpers
% ============================================================

function files = list_mat_files(dir_path)
d = dir(dir_path);
d = d(~[d.isdir]);
files = {d.name};
is_mat = endsWith(files, '.mat', 'IgnoreCase', true);
files = files(is_mat);
files = sort(files);
end

function [R, tvec] = load_event_any(fullpath_mat)
S = load(fullpath_mat);
if isfield(S, 'rainRate')
    R = S.rainRate;
elseif isfield(S, 'rainRatePGW')
    R = S.rainRatePGW;
else
    error('No rain rate variable in %s', fullpath_mat);
end

if isfield(S, 'times')
    tvec = normalize_times(S.times);
elseif isfield(S, 'timesPGW')
    tvec = normalize_times(S.timesPGW);
elseif isfield(S, 'time')
    tvec = normalize_times(S.time);
else
    error('No time variable in %s', fullpath_mat);
end
end

function tvec = normalize_times(times_in)
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
    error('Unsupported times data type');
end
end

function [R1a, R2a, t_common] = align_time_series_3d(R1, t1, R2, t2)
[t_common, ia, ib] = intersect(t1, t2, 'stable');
R1a = R1(:,:,ia);
R2a = R2(:,:,ib);
end

function rm = block_regional_max_over_time_fast(M, block_size, spatial_mode_valid)
% Spatial block-mean maps per time, then spatial max of each, then max over time
if isempty(M) || size(M,3) == 0
    rm = NaN; return
end
K = ones(block_size, block_size) / (block_size^2);
rm = -Inf;
if spatial_mode_valid
    for t = 1:size(M,3)
        S = conv2(M(:,:,t), K, 'valid');
        smax = max(S, [], 'all');
        if smax > rm, rm = smax; end
    end
else
    for t = 1:size(M,3)
        S = conv2(M(:,:,t), K, 'same');
        smax = max(S, [], 'all');
        if smax > rm, rm = smax; end
    end
end
if isinf(rm), rm = NaN; end
end

function vmax = block_spatial_mean_max(H, block_size, spatial_mode_valid)
% Spatial block mean over a 2D map H, then the spatial maximum
K = ones(block_size, block_size) / (block_size^2);
if spatial_mode_valid
    S = conv2(H, K, 'valid');
else
    S = conv2(H, K, 'same');
end
vmax = max(S, [], 'all');
end

function dt_min = infer_dt_minutes(tvec)
if numel(tvec) < 2
    error('Cannot infer time step from fewer than two timestamps');
end
dt = median(diff(tvec));
dt_min = minutes(dt);
if dt_min <= 0
    error('Non positive time step inferred');
end
end

% -----------------------------
% ONE CSV writer (intensities + totals)
% -----------------------------
function csv_path = write_all_results_csv(save_dir, durations_min, block_sizes_out, block_is_full, ...
    RM_hist, RM_pgw, TB_hist, TB_pgw, hist_kept, pgw_kept)

rows = {};
[ne, nd, nb] = size(RM_hist);

% Intensities rows
for e = 1:ne
    for d = 1:nd
        for b = 1:nb
            rows(end+1,:) = { ...
                e, 'intensity', durations_min(d), block_sizes_out(b), logical(block_is_full(b)), ...
                hist_kept{e}, pgw_kept{e}, ...
                RM_hist(e,d,b), RM_pgw(e,d,b), ...
                RM_pgw(e,d,b) - RM_hist(e,d,b), RM_pgw(e,d,b) ./ RM_hist(e,d,b), ...
                'mm_per_hr'}; %#ok<AGROW>
        end
    end
end

% Totals rows (duration_min=0)
for e = 1:ne
    for b = 1:nb
        rows(end+1,:) = { ...
            e, 'total', 0, block_sizes_out(b), logical(block_is_full(b)), ...
            hist_kept{e}, pgw_kept{e}, ...
            TB_hist(e,b), TB_pgw(e,b), ...
            TB_pgw(e,b) - TB_hist(e,b), TB_pgw(e,b) ./ TB_hist(e,b), ...
            'mm'}; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', { ...
    'event_id','measure_type','duration_min','block_size','is_full_domain', ...
    'hist_file','pgw_file', ...
    'hist_value','pgw_value','diff','ratio','units'});

csv_path = fullfile(save_dir, 'block_stats_events.csv');
writetable(T, csv_path);
end
