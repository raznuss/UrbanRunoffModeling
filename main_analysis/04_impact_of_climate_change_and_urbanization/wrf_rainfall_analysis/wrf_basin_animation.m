function wrf_basin_analysis()
% Plot an animation of rain intensity for event #1 (historical),
% frames t=1..1000, overlay a WGS84 basin polygon, and draw coastlines LAST.

%% ---------------- User paths ----------------
eventsDirHist = '\\vscifs\hydrolab1\hydrolab\home\moshe\Research\wrf\PGW\functions\savedMat\singleEvents\V3\Control\';
basin_shp     = 'D:\Development\RESEARCH\Raanana\gis\GIS\watershead_for_radar\Raanana_wsp_merged_wgs84.shp';

%% --------------- Load event #1 ---------------
files = list_mat_files(eventsDirHist);
assert(~isempty(files), 'No .mat files found in %s', eventsDirHist);
matfile = fullfile(eventsDirHist, files{1});
S = load(matfile);

% Detect rain + time
if isfield(S, 'rainRate')
    R = S.rainRate;
elseif isfield(S, 'rainRatePGW')
    R = S.rainRatePGW;
else
    error('No rainRate / rainRatePGW variable in %s', matfile);
end

if isfield(S, 'times')
    tvec = normalize_times(S.times);
elseif isfield(S, 'time')
    tvec = normalize_times(S.time);
elseif isfield(S, 'timesPGW')
    tvec = normalize_times(S.timesPGW);
else
    tvec = [];
end

% Lat/Lon (geographic, degrees)
if isfield(S, 'lat') && isfield(S, 'lon')
    lat = double(S.lat);
    lon = double(S.lon);
else
    error('lat/lon not found in %s. Please include lat/lon in the event MAT.', matfile);
end

% Animation range
T = size(R,3);
tEnd = min(144, T);
R = double(R);  % for plotting/math

% Color limits consistent along the animation
vmin = 0;
vmax = max(R(:,:,1:tEnd), [], 'all');  % max over frames 1..tEnd

%% --------------- Make map axes ---------------
latlim = [min(lat(:)) max(lat(:))];
lonlim = [min(lon(:)) max(lon(:))];

f = figure('Color','w','Units','normalized','Position',[0.08 0.08 0.68 0.78]);
ax = worldmap(latlim, lonlim);  % geographic axes
setm(ax, 'ParallelLabel','on','MeridianLabel','on');
hold on

% Light land layer (optional, מתחת לכול)
try
    L = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(ax, L, 'FaceColor', [0.96 0.96 0.96], 'EdgeColor', [0.75 0.75 0.75]);
catch
    warning('Could not load landareas. Proceeding without land fill.');
end

% First frame (so שהסדרי השכבות ייקבעו: רסטר -> אגן -> קווי חוף)
Z1 = R(:,:,1);
hImg = geoshow(ax, lat, lon, Z1, 'DisplayType','texturemap');  % raster bottom
colormap(parula); caxis([vmin, vmax]);
cb = colorbar; ylabel(cb, 'Rain rate (mm/hr)');

% Basin overlay (middle)
try
    Sbasin = shaperead(basin_shp, 'UseGeoCoords', true);
    geoshow(ax, Sbasin, 'FaceColor','none', 'EdgeColor',[0.85 0.2 0.2], 'LineWidth',1.8);
catch ME
    warning('Shapefile read/plot failed: %s', ME.message);
end

% Coastlines LAST (top-most)
try
    C = load('coastlines'); % coastlat, coastlon
    geoshow(ax, C.coastlat, C.coastlon, 'Color', [0 0 0], 'LineWidth', 0.75);
catch
    warning('Could not load coastlines.');
end

gridm on; mlabel on; plabel on;

% Title
if ~isempty(tvec)
    title(sprintf('Event #1 (historical) — frame 1/%d (%s)', tEnd, string(tvec(1))));
else
    title(sprintf('Event #1 (historical) — frame 1/%d', tEnd));
end

%% --------------- Animate frames ---------------
delay_s = 0.12;   % זמן בין פריימים
for t = 1:tEnd
    set(hImg, 'CData', R(:,:,t));  % update raster only; layers stay ordered

    if ~isempty(tvec)
        title(sprintf('Event #1 (historical) — frame %d/%d (%s)', t, tEnd, string(tvec(t))));
    else
        title(sprintf('Event #1 (historical) — frame %d/%d', t, tEnd));
    end

    drawnow;
    pause(delay_s);
end

fprintf('Animation complete: %s, frames 1..%d\n', matfile, tEnd);
end

%% ----------------- Helpers -----------------
function files = list_mat_files(dir_path)
d = dir(fullfile(dir_path, '*.mat'));
files = {d.name};
files = sort(files);
end

function tvec = normalize_times(times_in)
if iscell(times_in)
    if ~isempty(times_in) && isa(times_in{1}, 'datetime')
        tvec = vertcat(times_in{:});
    else
        tvec = datetime(cell2mat(times_in(:)), 'ConvertFrom','datenum');
    end
elseif isa(times_in, 'datetime')
    tvec = times_in(:);
elseif isnumeric(times_in)
    tvec = datetime(times_in(:), 'ConvertFrom','datenum');
else
    tvec = [];
end
end




