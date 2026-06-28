
%% load structure with all the dates
% The field peak_time is the peak time date and the field time is the day
% with rainfall the contribute to this peak. Therefore, events can be
% distingushied by sharing the same !!peak_time!!
clear 
clc

% load('D:\MY_CODES\Raanana_urban_floods\Radar_Files_Extraction\Gauges_Radar_raanana_byevent_bigger_than6.5cms.mat')
load('D:\MY_CODES\Raanana_urban_floods\Radar_Files_Extraction\Gauges_Radar_raanana_byevent_bigger_than6.5to13dis.mat')
% load('D:\MY_CODES\Raanana_urban_floods\Radar_Files_Extraction\high_dischrage_and_dates_raanana.csv')

%%Load one radar data just for grd object.
load('S:\hydrolab\ShareData\radarDatabase\IMS\IMS_archive_cleaned_Adj_NEW\daily_SR\2007-2008\20071029.mat')
path_out = 'D:\Development\RESEARCH\Raanana\data\rain_radar\Row_radar_ASC_files\';
path_radar = 'S:\hydrolab\ShareData\radarDatabase\IMS\IMS_archive_cleaned_Adj_NEW\daily_SR';
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

event_dates = unique([data_gague_radar.peak_time]);
convert_to_depth_factor = 12; % (5min minutes in 1 hour)

event_struc = struct();
for i = 1:length(event_dates)
% for i = 10
    i
    event_struc(i).peak_time = event_dates(i); %% this is the real date
    event_struc(i).mm_ws = []; %% milimiter over all the watershed for each time step
    event_struc(i).time_ws = []; %% time vector of the milimeter over the watersehd
    event_struc(i).split_flag = 0;  % if it will be one, I split the storm within a day

    rain_amout_flag = 1;
    day_temp = event_dates(i);
    date_datevec = datevec(day_temp);
    if date_datevec(1,4)<8 % hour < 8am  % if true the peak day is the day before in the radar file  because each day is from 8am to 8am;
            event_dates(i) = event_dates(i)-1;
            day_temp   = day_temp -1;
            date_datevec = datevec(day_temp);
    end
    
    
    peak_day_flag=1;
    minutes_missing=0;
    while rain_amout_flag == 1 %% loop backward from peak time;
        formatOut = 'yyymmdd';
        year =  date_datevec(1,1);
        month =  date_datevec(1,2);
        day =  date_datevec(1,3);
        hour = date_datevec(1,4);
        
        %% Decide which hydrological year the data should be in
        if month < 9 
            path_day = [path_radar '\' num2str(year-1) '-' num2str(year) '\' num2str(year) sprintf('%02d', month) sprintf('%02d',day) '.mat'];
            path_day_short = [num2str(year-1) '-' num2str(year) '\' num2str(year) sprintf('%02d', month) sprintf('%02d',day) '.mat'];
        else
            path_day = [path_radar '\' num2str(year) '-' num2str(year+1) '\' num2str(year) sprintf('%02d', month) sprintf('%02d',day) '.mat'];
            path_day_short =['\' num2str(year) '-' num2str(year+1) '\' num2str(year) sprintf('%02d', month) sprintf('%02d',day) '.mat'];
        end
        
        if exist(path_day)==0 %not exist    
            rain_amout_flag=0;
            event_struc(i).in_radar=0;
            continue
        end

        %% Load the data
        event_struc(i).in_radar=1;
        load(path_day)
        time_temp = time(end);
        %%
        
        %% Count the number of missing minutes.
        if length(time)<288 %% 280 minutes is full data length (12 *24) 
            minutes_missing = minutes_missing+[288-length(time)];
        end
    
        [cm_over_raanana] = squeeze(nansum(nansum(rainrate_mmh(yindexmax:yindexmin,xindexmin:xindexmax,:)./convert_to_depth_factor.*grd.cellsize^2./1000))); % m^3 [cubic meter] of rain;
        mm_over_raanana = (cm_over_raanana/raanana_box_area)*1000; % for every 5 min.
        acc_rain = cumsum(mm_over_raanana)./sum(mm_over_raanana); 
        
        
        %% DAILY CONDITION 
        %    5 milimiter thershold.  Check what is the portion of the watersehd with > 5mm in a day
        above_threshold_precent = sum(sum(daily_rain_mm(yindexmax:yindexmin,xindexmin:xindexmax) > 5)) / (size(daily_rain_mm(yindexmax:yindexmin,xindexmin:xindexmax) ,1)*size(daily_rain_mm(yindexmax:yindexmin,xindexmin:xindexmax) ,2));
        
        
        
        %% Going backward from the peak day.
        %   Condition #1 -  break out if :
        %                     <10% of the watershed get rain above >5 mm in total in this day.
        %              OR,
        %                     10 hours from the end of the day backwards the avergae precipitaton over the watershed is less than 5 mm
        %              OR,
        %                      If it is the peakday, go backwards no matter
        %                      what.
        if (((above_threshold_precent <0.1) || (...
                                    sum(mm_over_raanana(end-120:end))<5)) && (peak_day_flag~=1)) 
                peak_day_flag=0;
               break 
       
       %    Condition #2
       %        IF within 10 hours interval, the cummulutive rainfall over the
       %          rain in watershed is  less than 10 mm, bound this storm with the last       date where is was  less than 10 mm. 
       %        OR, 
       %             
      
        elseif min(conv(mm_over_raanana,ones(1,120),'same'))<5 && (peak_day_flag~=1)
            start_index = find(conv(mm_over_raanana,ones(1,120),'same')<10,1,'first');
            event_struc(i).start = time(start_index);
            event_struc(i).mm_ws = [event_struc(i).mm_ws; mm_over_raanana(start_index:end)]; 
            event_struc(i).time_ws = [event_struc(i).time_ws; time(start_index:end)]; 
            event_struc(i).split_flag = 1; 
            peak_day_flag = 0;
            break  
        else %% This day with alot of rain/peak day. therefore, save the start of the day as when 10% precent already precipitate and go check the previous day
            day_temp = day_temp-1;
            date_datevec = datevec(day_temp);
            start_index = find(acc_rain<0.01, 1, 'last' );        %% Take start of the rain  as after 10 precent already precipitate.
            if isempty(start_index) 
                start_index=1;
            end
            event_struc(i).start = time(start_index);
            if event_struc(i).start > event_struc(i).peak_time
                start_index=1;
                event_struc(i).start = time(start_index);
            end         
            event_struc(i).mm_ws = [event_struc(i).mm_ws; mm_over_raanana(start_index:end)]; 
            event_struc(i).time_ws = [event_struc(i).time_ws; time(start_index:end)];
            peak_day_flag = 0;
        end
    end
        
    rain_amout_flag = 1;
    peak_day_flag=1;
    day_temp = event_dates(i); % take 
    date_datevec = datevec(day_temp);
    [~ ,peak_index] = min(abs(time-event_dates(i)));
    
    %% Going forward from the peak day.
    while rain_amout_flag == 1 %% loop backward from peak time;
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
        
        try
        load(path_day)
        catch
            rain_amout_flag=0;
            continue
        end
%         
        if length(time)<280
            minutes_missing = minutes_missing+[288-length(time)];
        end
        
        time_temp = time(1);
        [cm_over_raanana] = squeeze(nansum(nansum(rainrate_mmh(yindexmax:yindexmin,xindexmin:xindexmax,:)./convert_to_depth_factor.*grd.cellsize^2./1000))); % m^3 of rain;
        mm_over_raanana = (cm_over_raanana/raanana_box_area)*1000; % for every 5 min.
        acc_rain = cumsum(mm_over_raanana)./sum(mm_over_raanana);
        above_threshold_precent = sum(sum(daily_rain_mm(yindexmax:yindexmin,xindexmin:xindexmax) > 5)) / (size(daily_rain_mm(yindexmax:yindexmin,xindexmin:xindexmax) ,1)*size(daily_rain_mm(yindexmax:yindexmin,xindexmin:xindexmax) ,2));
       
             
        %% If end day is the peak time day, take 5 hours after the peak and check if > average 5 mm has been precipitate. 
        %% If yes and if there is enough time until the end of the day (the peak time is not at the end of the day), take this time + 5 hours as the end of the flow.
        if peak_day_flag %% 
            peak_day_flag=0;
            if sum(mm_over_raanana(peak_index:end)) <5 && (length(time)-peak_index)>60
                event_struc(i).end = event_struc(i).peak_time + datenum(hours(5));
                event_struc(i).duration_h = (event_struc(i).end - event_struc(i).start)*24;
                break
            else                %% Peak Dischrge day but there is still much rain in this day so go to the next dayt.             
                day_temp = day_temp+1;
                date_datevec = datevec(day_temp);
                continue
            end
        end
        if (above_threshold_precent <0.1 ||...
                                    sum(mm_over_raanana(1:72))<5) %% Not the peak dischrage day so check if on on average, more than 5 mm has been precipitate in 6 hours.
                                
            event_struc(i).end = time(1);
            event_struc(i).duration_h = (event_struc(i).end - event_struc(i).start)*24;
            break 
        else % This day with alot of rainfall, therefore save the time where 0.9 of the rainfall has already precipitate and go check next day
            day_temp = day_temp+1;
            date_datevec = datevec(day_temp);
            end_index = find(acc_rain<0.95, 1, 'last' );        %% Take end of the rain  as after 95 precent already precipitate.
            event_struc(i).end = time(end_index);
            event_struc(i).duration_h = (event_struc(i).end - event_struc(i).start)*24;
        end
    end
    
    %% Just for sorting the time vector
%     temp_a(:,1) = event_struc(i).time_ws;
%     temp_a(:,2) = event_struc(i).mm_ws;
% %     temp_a = sortrows(temp_a,1);
%     event_struc(i).time_ws = temp_a(:,1);
%     event_struc(i).mm_ws =  temp_a(:,2);
    temp_a=[];
    event_struc(i).missing_min = minutes_missing;
    
end

save('D:\MY_CODES\Raanana_urban_floods\Radar_Files_Extraction\event_structure_bigger_6.5to13dis.mat','event_struc');
