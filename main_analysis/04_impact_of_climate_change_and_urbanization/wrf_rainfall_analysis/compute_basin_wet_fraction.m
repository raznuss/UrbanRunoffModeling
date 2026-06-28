function results = compute_basin_wet_fraction(thresholds_mm_per_hr, shpPath, eventsDirHist, eventsDirPGW, save_dir, events_idx, verbose, show_waitbar)
% Compute mean wet-area fraction over a basin for multiple events, thresholds and scenarios.
%
% For each event and each rain-rate threshold, we:
%   1) Compute, at each time step, the fraction of basin pixels where rainRate > threshold.
%   2) Take the time-mean of that fraction over the full event duration.
%
% Outputs are saved both as a MAT struct and as a CSV in long format.
%
% Inputs:
%   thresholds_mm_per_hr : vector of thresholds in mm/hr (e.g. [0.1 0.5 1 2 5])
%   shpPath              : path to basin shapefile (WGS84)
%   eventsDirHist        : folder with historical event MAT files
%   eventsDirPGW         : folder with PGW event MAT files
%   save_dir             : output folder (optional, default: <pwd>/basin_wet_fraction)
%   events_idx           : which events to process (optional, default: all)
%   verbose              : logical, print progress and sanity checks (optional, default: true)
%   show_waitbar         : logical, show waitbar (optional, default: true)
%
% Event MAT files are expected to contain:
%   - rainRate or rainRatePGW  [ny x nx x T], units mm/hr
%   - lat, lon                 [ny x nx], WGS84
%   - times / timesPGW / time  (datetime or datenum or cell)
%
% Result fields:
%   results.thresholds_mm_per_hr
%   results.mean_wet_hist   [nEvents x nThr]
%   results.mean_wet_pgw    [nEvents x nThr]
%   results.event_hist_files
%   results.event_pgw_files
%   results.mask_basin      [ny x nx] logical
%   results.lat, results.lon

    if nargin < 1 || isempty(thresholds_mm_per_hr)
        thresholds_mm_per_hr = 0.1;  % default single threshold
    end
    if nargin < 2
        error('You must provide shpPath.');
    end
    if nargin < 3
        error('You must provide eventsDirHist.');
    end
    if nargin < 4
        error('You must provide eventsDirPGW.');
    end
    if nargin < 5 || isempty(save_dir)
        save_dir = fullfile(pwd, 'basin_wet_fraction');
    end
    if nargin < 6, events_idx   = [];    end
    if nargin < 7, verbose      = true;  end
    if nargin < 8, show_waitbar = true;  end

    if ~isfolder(save_dir)
        mkdir(save_dir);
    end

    thresholds_mm_per_hr = thresholds_mm_per_hr(:)';  % row vector
    nThr = numel(thresholds_mm_per_hr);

    % ----------------------------------------
    % Discover events
    % ----------------------------------------
    files_hist_all = list_mat_files(eventsDirHist);
    files_pgw_all  = list_mat_files(eventsDirPGW);

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
    if nEvents == 0
        error('No events found to process.');
    end

    % ----------------------------------------
    % Load basin polygon and grid from first historical event
    % ----------------------------------------
    % Basin polygon
    GT = readgeotable(shpPath);
    basinShape = GT.Shape(1);  % assume first polygon is the basin

    % Grid lat/lon
    S0 = load(fullfile(eventsDirHist, files_hist{1}));
    if ~isfield(S0, 'lat') || ~isfield(S0, 'lon')
        error('First event file must contain lat and lon fields.');
    end
    lat = double(S0.lat);
    lon = double(S0.lon);
    [ny, nx] = size(lat);

    % Build basin mask on the grid (pixel centers inside basin)
    if verbose
        fprintf('Building basin mask from shapefile and lat/lon grid...\n');
    end
    P = geopointshape(lat(:), lon(:));
    inside_vec = isinterior(basinShape, P);
    mask_basin = reshape(inside_vec, ny, nx);
    nBasinPix  = sum(mask_basin(:));
    if nBasinPix == 0
        error('Basin mask is empty – check shapefile and grid alignment.');
    end

    if verbose
        fprintf('  Basin mask has %d pixels inside basin.\n', nBasinPix);
    end

    % ----------------------------------------
    % Preallocate outputs
    % ----------------------------------------
    mean_wet_hist = nan(nEvents, nThr);   % [event x threshold], fraction 0–1
    mean_wet_pgw  = nan(nEvents, nThr);

    hist_kept = cell(0,1);
    pgw_kept  = cell(0,1);

    if show_waitbar
        h = waitbar(0, 'Processing events...');
    end

    % ----------------------------------------
    % Loop over events
    % ----------------------------------------
    for j = 1:nEvents
        src_idx = events_idx(j);

        if verbose
            fprintf('Event %02d/%02d (source #%d)\n', j, nEvents, src_idx);
            fprintf('  Hist file: %s\n', files_hist{j});
            fprintf('  PGW  file: %s\n', files_pgw{j});
        end

        t_event = tic;

        % Load historical and PGW cubes
        try
            [R_hist, ~] = load_event_any(fullfile(eventsDirHist, files_hist{j}));
            [R_pgw , ~] = load_event_any(fullfile(eventsDirPGW , files_pgw{j}));
        catch ME
            warning('Skipping event %d (source %d): %s', j, src_idx, ME.message);
            if show_waitbar
                waitbar(j/nEvents, h, sprintf('Event %d skipped', j));
            end
            continue
        end

        % Basic size check
        if size(R_hist,1) ~= ny || size(R_hist,2) ~= nx
            warning('Event %d historical grid size mismatch, skipping.', src_idx);
            continue;
        end
        if size(R_pgw,1) ~= ny || size(R_pgw,2) ~= nx
            warning('Event %d PGW grid size mismatch, skipping.', src_idx);
            continue;
        end

        % Optional: drop first time step (parity with other workflows)
        if size(R_hist,3) > 0
            R_hist(:,:,1) = [];
        end
        if size(R_pgw,3) > 0
            R_pgw(:,:,1) = [];
        end

        % Compute mean wet fraction for each threshold
        for k = 1:nThr
            thr = thresholds_mm_per_hr(k);
            mean_wet_hist(j,k) = compute_mean_wet_fraction(R_hist, mask_basin, thr);
            mean_wet_pgw (j,k) = compute_mean_wet_fraction(R_pgw , mask_basin, thr);
        end

        % Sanity check print: mean wet fractions per threshold
        if verbose
            fprintf('  Sanity check: mean wet-area fraction (historical vs PGW)\n');
            for k = 1:nThr
                thr  = thresholds_mm_per_hr(k);
                vh   = mean_wet_hist(j,k);
                vp   = mean_wet_pgw(j,k);

                if isfinite(vh) && vh > 0
                    ratio = vp / vh;
                else
                    ratio = NaN;
                end

                fprintf('    thr = %.3f mm/hr: hist = %.3f, pgw = %.3f, ratio = %.2f\n', ...
                    thr, vh, vp, ratio);
            end
        end

        hist_kept{end+1,1} = files_hist{j}; %#ok<AGROW>
        pgw_kept {end+1,1} = files_pgw{j};  %#ok<AGROW>

        if verbose
            fprintf('  Done event %02d in %.1fs\n', j, toc(t_event));
        end
        if show_waitbar
            waitbar(j/nEvents, h, sprintf('Event %d of %d', j, nEvents));
        end
    end

    if show_waitbar
        close(h);
    end

    % Keep only processed events (in case some were skipped)
    nKept = numel(hist_kept);
    mean_wet_hist = mean_wet_hist(1:nKept, :);
    mean_wet_pgw  = mean_wet_pgw (1:nKept, :);

    % ----------------------------------------
    % Pack results struct
    % ----------------------------------------
    results.thresholds_mm_per_hr = thresholds_mm_per_hr;
    results.mean_wet_hist        = mean_wet_hist;   % [event x threshold], fraction 0–1
    results.mean_wet_pgw         = mean_wet_pgw;
    results.event_hist_files     = hist_kept;
    results.event_pgw_files      = pgw_kept;
    results.mask_basin           = mask_basin;
    results.lat                  = lat;
    results.lon                  = lon;

    if ~isempty(save_dir)
        mat_path = fullfile(save_dir, 'basin_wet_fraction_results.mat');
        save(mat_path, 'results');

        csv_path = write_wet_fraction_csv(save_dir, thresholds_mm_per_hr, ...
            mean_wet_hist, mean_wet_pgw, hist_kept, pgw_kept);

        if verbose
            fprintf('Saved MAT: %s\n', mat_path);
            fprintf('Saved CSV: %s\n', csv_path);
        end
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
        tvec = [];
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

function mfrac = compute_mean_wet_fraction(R, mask_basin, thr_mm_per_hr)
% R: [ny x nx x T], mm/hr
% mask_basin: [ny x nx] logical
% thr_mm_per_hr: scalar threshold
%
% Returns:
%   mfrac – time-mean fraction of basin pixels where rainRate > thr (0–1)

    [~, ~, T] = size(R);
    if T == 0
        mfrac = NaN;
        return;
    end

    % Broadcast mask over time
    mask3 = repmat(mask_basin, 1, 1, T);

    wet = (R > thr_mm_per_hr) & mask3;          % basin pixels above threshold
    basinPixCount = sum(mask_basin(:));         % constant

    wet_counts = squeeze(sum(sum(wet, 1), 2));  % [T x 1]
    frac_t     = double(wet_counts) / double(basinPixCount);

    mfrac = mean(frac_t);
end

function csv_path = write_wet_fraction_csv(save_dir, thresholds_mm_per_hr, ...
    mean_wet_hist, mean_wet_pgw, hist_kept, pgw_kept)

    rows = {};
    [nEvents, nThr] = size(mean_wet_hist);

    for e = 1:nEvents
        for k = 1:nThr
            thr = thresholds_mm_per_hr(k);

            % Historical row
            rows(end+1,:) = { ...
                e, 'historical', thr, ...
                hist_kept{e}, ...
                mean_wet_hist(e,k), ...
                'fraction'}; %#ok<AGROW>

            % PGW row
            rows(end+1,:) = { ...
                e, 'future', thr, ...
                pgw_kept{e}, ...
                mean_wet_pgw(e,k), ...
                'fraction'}; %#ok<AGROW>
        end
    end

    T = cell2table(rows, 'VariableNames', { ...
        'event_id', 'scenario', 'threshold_mm_per_hr', ...
        'event_file', 'mean_wet_fraction', 'units'});

    csv_path = fullfile(save_dir, 'basin_wet_fraction_events.csv');
    writetable(T, csv_path);
end
