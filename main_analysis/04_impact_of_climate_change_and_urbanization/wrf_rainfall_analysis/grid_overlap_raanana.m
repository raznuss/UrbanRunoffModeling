%% ===== קלטים =====
shpPath = 'D:\Development\RESEARCH\Raanana\gis\GIS\watershead_for_radar\Raanana_wsp_merged_wgs84.shp';
matPath = '\\vscifs\hydrolab1\hydrolab\home\moshe\Research\wrf\PGW\functions\savedMat\singleEvents\V3\Control\WRFrain_noPGW_HPE01.mat';

samplesPerSide = 7;      % דיוק דגימה בתוך תא
overlapThresh  = 0.0;    % סף חפיפה של תא עם האגן כדי להיחשב רלוונטי
boxSize        = 3;      % גודל בוקס בבסיס תאים של הרשת

%% ===== טעינת רשת ושייפ =====
S = load(matPath,'lat','lon');
assert(isfield(S,'lat') && isfield(S,'lon'), 'בקובץ MAT חייבות להיות המטריצות lat ו lon');

lat = double(S.lat);
lon = double(S.lon);
[ny,nx] = size(lat);
assert(isequal(size(lon),[ny,nx]), 'lat ו lon חייבות להיות באותו גודל');

GT = readgeotable(shpPath);
basinShape = GT.Shape(1);      % geopolyshape של האגן

%% ===== חישוב חפיפה לכל תא רשת מול האגן =====
s = max(3, samplesPerSide);
u = (0.5:s-0.5)/s;  v = u;  [U,V] = meshgrid(u,v);

alpha = zeros(ny-1, nx-1);     % חפיפה של תאי הרשת

try
    [latlim, lonlim] = geobounds(basinShape);
    latminBB = latlim(1); latmaxBB = latlim(2);
    lonminBB = lonlim(1); lonmaxBB = lonlim(2);
catch
    latminBB = -Inf; latmaxBB = Inf; lonminBB = -Inf; lonmaxBB = Inf;
end

for i = 1:(ny-1)
    for j = 1:(nx-1)
        % דילוג מהיר אם התא מחוץ לתיבה המקיפה
        lat_cell = [lat(i,j) lat(i,j+1) lat(i+1,j) lat(i+1,j+1)];
        lon_cell = [lon(i,j) lon(i,j+1) lon(i+1,j) lon(i+1,j+1)];
        if max(lat_cell) < latminBB || min(lat_cell) > latmaxBB || ...
           max(lon_cell) < lonminBB || min(lon_cell) > lonmaxBB
            continue
        end

        % ארבע פינות התא
        lat00 = lat(i  , j  ); lon00 = lon(i  , j  );
        lat01 = lat(i  , j+1); lon01 = lon(i  , j+1);
        lat10 = lat(i+1, j  ); lon10 = lon(i+1, j  );
        lat11 = lat(i+1, j+1); lon11 = lon(i+1, j+1);

        % דגימות משנה ביליניאריות בתוך התא
        lat_uv = (1-U).*(1-V).*lat00 + U.*(1-V).*lat01 + (1-U).*V.*lat10 + U.*V.*lat11;
        lon_uv = (1-U).*(1-V).*lon00 + U.*(1-V).*lon01 + (1-U).*V.*lon10 + U.*V.*lon11;

        inside = isinterior(basinShape, geopointshape(lat_uv(:), lon_uv(:)));
        alpha(i,j) = mean(inside);
    end
end

maskCells = alpha > overlapThresh;    % תאי רשת שנוגעים באגן

%% ===== מציאת כל הבוקסים האפשריים בגודל 3 על 3 שנוגעים באגן =====
% מסכת מועמדים: בוקס נחשב אם לפחות תא אחד בתוך האגן
C = conv2(double(maskCells), ones(boxSize), 'valid') > 0;   % גודל C הוא [(ny-1)-(boxSize)+1] על [(nx-1)-(boxSize)+1]
[iBox, jBox] = find(C);
nBoxes = numel(iBox);

%% ===== בניית מצולעים לציור התאים הרלוונטיים =====
latList = []; lonList = [];
for i = 1:(ny-1)
    for j = 1:(nx-1)
        if ~maskCells(i,j), continue; end
        latpoly = [lat(i,j) lat(i,j+1) lat(i+1,j+1) lat(i+1,j) lat(i,j) NaN];
        lonpoly = [lon(i,j) lon(i,j+1) lon(i+1,j+1) lon(i+1,j) lon(i,j) NaN];
        latList = [latList latpoly]; %#ok<AGROW>
        lonList = [lonList lonpoly]; %#ok<AGROW>
    end
end

%% ===== ציור מפה =====
f = figure('Color','w','Units','normalized','Position',[0.15 0.15 0.7 0.7]);
gx = geoaxes(f); hold(gx,'on');

% תאי רשת רלוונטיים
hCells = geoplot(gx, latList, lonList, 'LineWidth', 0.7);

% האגן עצמו — שימוש ב EdgeColor ולא ב Color כדי למנוע את השגיאה
hBasin = geoplot(gx, basinShape, 'EdgeColor', [0.85 0.2 0.2], 'FaceColor', 'none', 'LineWidth', 1.8);

% ציור כל הבוקסים האפשריים 3 על 3 הצמודים לרשת
% כל בוקס מצויר כמלבן לפי ארבעת צמתי הפינות שלו
if nBoxes > 0
    cmap = lines(max(nBoxes, 7));   % מפה של צבעים לחזרה מחזורית
    for k = 1:nBoxes
        ii = iBox(k); jj = jBox(k);          % אינדקס פינת שמאל עליון במרחב של תאי רשת
        % פינות המלבן של בוקס 3 על 3 תאים משתמשות בצמתים [ii .. ii+3] ו [jj .. jj+3]
        latBox = [lat(ii    , jj    ), lat(ii    , jj+3), lat(ii+3, jj+3), lat(ii+3, jj    ), lat(ii    , jj    )];
        lonBox = [lon(ii    , jj    ), lon(ii    , jj+3), lon(ii+3, jj+3), lon(ii+3, jj    ), lon(ii    , jj    )];
        col = cmap(mod(k-1,size(cmap,1))+1, :);
        geoplot(gx, latBox, lonBox, 'LineWidth', 1.2, 'Color', col);
        % אם תרצה גם מילוי שקוף בתוך הבוקס אפשר להשתמש ב geopolyshape:
        % shpBox = geopolyshape(latBox, lonBox);
        % geoplot(gx, shpBox, 'FaceColor', col, 'FaceAlpha', 0.08, 'EdgeColor', col, 'LineWidth', 1.0);
    end
end

% תחומי תצוגה ומפה בסיסית
latmin = min(lat(:)); latmax = max(lat(:));
lonmin = min(lon(:)); lonmax = max(lon(:));
geolimits(gx, [latmin latmax], [lonmin lonmax]);
try, geobasemap(gx,'streets-light'); catch, end

title(gx, sprintf('Cells touching Raanana Basin and all %dx%d candidate boxes', boxSize, boxSize));
legend(gx, [hCells hBasin], {'Relevant grid cells','Basin outline'}, 'Location','best');
