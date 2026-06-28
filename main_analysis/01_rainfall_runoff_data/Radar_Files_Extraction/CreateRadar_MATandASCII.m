%% This script convert RADAR .mat file to ASCII after storms where boudnded using "BoundStormsInRadarData" algorithm has been used.
%% The code reads the .mat structure and create ASCII files day by day.

% The field peak_time is the peak time date and the field time is the day
% with rainfall the contribute to this peak. Therefore, events can be
% distingushied by sharing the same !!peak_time!!

clear 
clc
path_out = 'D:\Development\RESEARCH\Raanana\data\rain_radar\Row_radar_ASC_files\';
path_radar = 'S:\hydrolab\ShareData\radarDatabase\IMS\IMS_archive_cleaned_Adj_NEW\daily_SR';
time_step = 5 ; % radar time step in minutes;

%%Load one radar data just for grd object.
load('S:\hydrolab\ShareData\radarDatabase\IMS\IMS_archive_cleaned_Adj_NEW\daily_SR\2007-2008\20071029.mat')
raanana_watershed = shaperead('D:\Development\RESEARCH\Raanana\gis\GIS\watershead_for_radar\Raanana_wsp.shp');

x = grd.xllcorner:grd.cellsize:(grd.ncols-1)*grd.cellsize+grd.xllcorner;
y = flip(grd.yllcorner:grd.cellsize:(grd.nrows-1)*grd.cellsize+grd.yllcorner);

[X,Y] = meshgrid(x,y);
raanana_box = raanana_watershed(end).BoundingBox;

xline = [raanana_box(1,1),raanana_box(1,1),raanana_box(2,1),raanana_box(2,1),raanana_box(1,1)]';
yline = [raanana_box(1,2),raanana_box(2,2),raanana_box(1,2),raanana_box(2,2),raanana_box(1,2)]';
temp(:,1)=xline;
temp(:,2)=yline;
temp=sortrows(temp,1);
raanana_box_area = polyarea(temp(:,1),temp(:,2));


[~,xindexmin] = min(abs(x-raanana_box(1,1)));
[~,xindexmax] = min(abs(x-raanana_box(2,1)));

[~,yindexmin] = min(abs(y-raanana_box(1,2)));
[~,yindexmax] = min(abs(y-raanana_box(2,2)));



load('D:\MY_CODES\Raanana_urban_floods\Radar_Files_Extraction\event_structure_bigger_6.5to13dis.mat')

for i = 1:size(event_struc,2)


   
    event_start_date = event_struc(i).start;
    event_end_date = event_struc(i).end;
    event_peak_date = event_struc(i).peak_time;
    event_days = [datetime(datestr(event_struc(i).start)),datetime(datestr(event_struc(i).end))];

    %% Make sure you take the correct day (days in the radar are 8 am to 20)
    if hour(event_start_date) < 8 
        event_days(1) = event_days(1)-1;
    end
    if hour(event_end_date) <8
        event_days(end) = event_days(end)-1;
    end
    
    event_days = [dateshift(datetime(datestr(event_days(1))),'start','day'):1:dateshift(datetime(datestr(event_days(end))),'start','day')];
    event_days_datevec = datevec(event_days);
    dir_path = [sprintf('%02d',event_days_datevec(1,1)), sprintf('%02d', event_days_datevec(1,2)),...
                sprintf('%02d',event_days_datevec(1,3)), '-',  sprintf('%02d',event_days_datevec(end,1)),...
                 sprintf('%02d',event_days_datevec(end,2)),  sprintf('%02d',event_days_datevec(end,3))];
    mkdir([path_out dir_path])
    
    
    
    event_start_flag = 1;
    radar_rain_out= 0;
    radar_rain_in= 0;
    mat_clipped_all = [];
    event_file = struct();
    cnt=1;
      
    for j = 1:length(event_days)
        %% directory is Sep first year until June last year. i.e diretory 2010-2011 is sep 2010 - june 2011
        date_datevec = event_days_datevec(j,:);
        formatOut = 'yyymmdd';
        year =  date_datevec(1,1);
        month =  date_datevec(1,2);
        day =  date_datevec(1,3);

        if month < 9 % decide which hydrological year the data should be in
            path_day = [path_radar '\' num2str(year-1) '-' num2str(year) '\' num2str(year) sprintf('%02d', month) sprintf('%02d',day) '.mat'];
            path_day_short = [num2str(year-1) '-' num2str(year) '\' num2str(year) sprintf('%02d', month) sprintf('%02d',day) '.mat'];
        else
            path_day = [path_radar '\' num2str(year) '-' num2str(year+1) '\' num2str(year) sprintf('%02d', month) sprintf('%02d',day) '.mat'];
            path_day_short =['\' num2str(year) '-' num2str(year+1) '\' num2str(year) sprintf('%02d', month) sprintf('%02d',day) '.mat'];
        end
        
        load(path_day)
%% now loops on times
        if event_start_flag==1
            [diffrence,start_index] = min(abs(time-event_start_date));    %% Take start of the rain  as after 5 precent already precipitate.
            if diffrence>datenum(minutes(10)) % not in this day
                continue
            end
            event_start_flag=0;
            time_out=time(start_index:end);
            rainrate_mmh_out = rainrate_mmh(:,:,start_index:end);
            
            % mat_clipped_all = rainrate_mmh_out(yindexmax-10:yindexmin+4,xindexmin-2:xindexmax+2,:); 
            mat_clipped_all = rainrate_mmh_out(yindexmax-50:yindexmin+50, xindexmin-15:xindexmax+15, :);

            time_clipped_all = time_out;
        elseif j == length(event_days) % means this is the last day of the event
            [diffrence,end_index] = min(abs(time-event_end_date));
            if diffrence>datenum(minutes(10))
                continue
            end
            time_out=time(1:end_index);
            rainrate_mmh_out = rainrate_mmh(:,:,1:end_index);
            
            mat_clipped_all = cat(3,mat_clipped_all,rainrate_mmh_out(yindexmax-50:yindexmin+50, xindexmin-15:xindexmax+15, :)); % take time series only until end index
            time_clipped_all = cat(1,time_clipped_all, time_out);
        else
            time_out=time;
            rainrate_mmh_out = rainrate_mmh;
            
            mat_clipped_all = cat(3,mat_clipped_all,rainrate_mmh_out(yindexmax-50:yindexmin+50, xindexmin-15:xindexmax+15, :)); %% take all time series
            time_clipped_all = cat(1,time_clipped_all, time_out);
        end            
    end
    
    convert_to_depth_factor = 60/(time_step);

    %%
    data_out =  mat_clipped_all;
    time_out =   time_clipped_all;

    event_struc(i).radar_rain =  data_out;
    event_struc(i).radar_times = time_out;

    kvec=[];
    cellsize = grd.cellsize;
    nanvalue = -9999;
    xll = x(xindexmin-15);
    yll = y(yindexmin+50);
    
    event_file.rain = data_out;
    event_file.time = time_clipped_all;
    event_file.cellsize = 500;
    event_file.xll = xll;
    event_file.yll = yll;
    save(['D:\Development\RESEARCH\Raanana\data\rain_radar\SA_rain_radar_transformation\Row_radar_MAT_files\' num2str(yyyymmdd(datetime(datestr(event_peak_date)))) '.mat'],'event_file')
    
    time_vec_datevec = datevec(time_out);
    for k = 1:length(time_out)
        mat = data_out(:,:,k);
        mat = mat./convert_to_depth_factor;
        mat(isnan(mat))=0;
        year =  time_vec_datevec(k,1);
        month = time_vec_datevec(k,2);
        day =  time_vec_datevec(k,3);
        hours = time_vec_datevec(k,4);
        minu = time_vec_datevec(k,5);
        file_full_path = [path_out, dir_path,'\',sprintf('%02d', year),...
            sprintf('%02d', month),sprintf('%02d', day),'_',sprintf('%02d', hours), sprintf('%02d', minu)];
        Path_out = file_full_path;
        cnt=cnt+1;
            hold on;scatter([xindexmin,xindexmin,xindexmax,xindexmax],[yindexmin,yindexmax,yindexmin,yindexmax])
        CreateASCII(mat,Path_out,nanvalue,xll,yll,cellsize)
    end   
end

