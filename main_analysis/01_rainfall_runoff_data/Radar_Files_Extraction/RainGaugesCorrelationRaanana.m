clear
clc
flow_data = readtable('D:\MY_CODES\Raanana_urban_floods\Radar_Files_Extraction\high_dischrage_and_dates_raanana_6.5to13dis.csv');
load('D:\MY_CODES\Raanana_urban_floods\Radar_Files_Extraction\gaugeData_updated_042023.mat')

%% Clip gauges within RM area
xcordmin = 185000;
xcordmax = 197000;

ycordmin = 668000;
ycordmax = 685000;

raanana_watershed = shaperead('D:\Development\RESEARCH\Raanana\gis\GIS\watershead_for_radar\Raanana_wsp.shp');

indexim = find([gaugeData.x]> xcordmin & [gaugeData.x] < xcordmax & [gaugeData.y]> ycordmin & [gaugeData.y] < ycordmax);
gauges_names = {gaugeData(indexim).name};
f1 = figure;scatter([gaugeData(indexim).x],[gaugeData(indexim).y]); hold on; text([gaugeData(indexim).x],[gaugeData(indexim).y],gauges_names)
hold on;

% print(f1, '\Radar_Gauges_Calibration\radar_domain.jpg', '-djpeg', '-r600')
plot(raanana_watershed(end).X,raanana_watershed(end).Y);hold on
plot([raanana_watershed(end).BoundingBox(:,1)],[raanana_watershed(end).BoundingBox(:,2)]);hold on


dates_above_thrshold = x2mdate(flow_data{:,1});
%% 
%%% Take all days after and before the peak discrharge with rain above threshold 
all_days=[];
cnt=1;
% raz = [];
for j =1:length(dates_above_thrshold)
    
    date_matnum = dates_above_thrshold(j); %peak discharge date
    
    for i = 1:length(indexim) % indexim of gauges in RM, Loope on gauge
%         raz = [raz cnt];
        gauge_index = indexim(i);
        gauge_xcord = gaugeData(gauge_index).x;
        gauge_ycord = gaugeData(gauge_index).y;

     
        data_gague_radar(cnt).id = gaugeData(gauge_index).id;
        data_gague_radar(cnt).name = gaugeData(gauge_index).name;
        data_gague_radar(cnt).x = gaugeData(gauge_index).x;
        data_gague_radar(cnt).y = gaugeData(gauge_index).y;
        data_gague_radar(cnt).time = date_matnum;

        %% get the rain from the station
        gauge_time_vec = [gaugeData(gauge_index).time];
        gauge_rain_vec = [gaugeData(gauge_index).vals];
        date_matnum_temp = date_matnum;
        time_index = find(( date_matnum_temp-gauge_time_vec)<1 & ( date_matnum_temp-gauge_time_vec)>=0 );
        if isempty(time_index)
             date_matnum_temp=  date_matnum_temp-1;
            time_index = find(( date_matnum_temp-gauge_time_vec)<1 & ( date_matnum_temp-gauge_time_vec)>=0 );
        end
        if isempty(time_index)
            disp('out')
            data_gague_radar(cnt).in_gauge = NaN;
            continue
        end
        data_gague_radar(cnt).in_gauge = 1;
        
        rain_flag=1;
        % find days with rain before the peak discharge
        while rain_flag==1

            disp('in')
            time_index = find((date_matnum_temp-gauge_time_vec)<1 & (date_matnum_temp-gauge_time_vec)>=0 );
            if isempty(time_index)
                break
            end
            if gauge_rain_vec(time_index) > 5 %mm
                data_gague_radar(cnt).gauge_rain =  gauge_rain_vec(time_index);
                data_gague_radar(cnt).id = gaugeData(gauge_index).id;
                data_gague_radar(cnt).name = gaugeData(gauge_index).name;
                data_gague_radar(cnt).x = gaugeData(gauge_index).x;
                data_gague_radar(cnt).y = gaugeData(gauge_index).y;
                data_gague_radar(cnt).time = date_matnum_temp;
                data_gague_radar(cnt).peak_time = date_matnum;
                data_gague_radar(cnt).in_gauge = 1;
                try
                    data_gague_radar(cnt).gauge_flag = gaugeData(gauge_index).flag(time_index);
                catch
                    data_gague_radar(cnt).gauge_flag = 8;
                end
                cnt = cnt+1;

            end
            date_matnum_temp= date_matnum_temp-1;

        end
        % find days with rain after the peak discharge
        date_matnum_temp = date_matnum;
        while rain_flag==1

            time_index = find((date_matnum_temp-gauge_time_vec)<1 & (date_matnum_temp-gauge_time_vec)>=0 );
            if isempty(time_index)
                break
            end
            if gauge_rain_vec(time_index) > 5 %mm
                data_gague_radar(cnt).gauge_rain =  gauge_rain_vec(time_index);
                data_gague_radar(cnt).id = gaugeData(gauge_index).id;
                data_gague_radar(cnt).name = gaugeData(gauge_index).name;
                data_gague_radar(cnt).x = gaugeData(gauge_index).x;
                data_gague_radar(cnt).y = gaugeData(gauge_index).y;
                data_gague_radar(cnt).in_gauge = 1;
                data_gague_radar(cnt).time = date_matnum_temp;
                data_gague_radar(cnt).peak_time = date_matnum;

                try
                    data_gague_radar(cnt).gauge_flag = gaugeData(gauge_index).flag(time_index);
                catch
                    data_gague_radar(cnt).gauge_flag = 8;
                end
                cnt = cnt+1;

            end
            date_matnum_temp= date_matnum_temp+1;


        end
    end
end

%% Loop over all dates for radar data.
cnt=1;
firstflag=1; % for coordinate;
times = unique([data_gague_radar.time]);
[data_gague_radar(:).in_radar]= deal(0);
for j = 1:length([times]) %% Loope on date
    j
    date_matnum = times(j); %
    date_datevec = datevec(date_matnum);
    indexim = find([data_gague_radar.time]==date_matnum ); %% all gauges that have data on this day.
    path_radar = 'S:\hydrolab\ShareData\radarDatabase\IMS\IMS_archive_cleaned_Adj_NEW\daily_SR';
    %% directory is Sep first year until June last year. i.e diretory 2010-2011 is sep 2010 - june 2011
    formatOut = 'yyymmdd';
    year =  date_datevec(1,1);
    month =  date_datevec(1,2);
    day =  date_datevec(1,3);

    if month < 9 % decide which hydrological year the data should be in
        path_day = [path_radar '\' num2str(year-1) '-' num2str(year) '\' num2str(year) sprintf('%02d', month) sprintf('%02d',day) '.mat'];
    else
        path_day = [path_radar '\' num2str(year) '-' num2str(year+1) '\' num2str(year) sprintf('%02d', month) sprintf('%02d',day) '.mat'];
    end

    if exist(path_day)==0 %not exist    
        continue
    end

    try
        load(path_day)
        [data_gague_radar(indexim).in_radar]=deal(1);
    catch
        continue
    end

    for i= 1:length(indexim)
        gauge_index = indexim(i);
        gauge_xcord = data_gague_radar(gauge_index).x;
        gauge_ycord = data_gague_radar(gauge_index).y;

        if firstflag==1
            xcords = linspace(grd.xllcorner,grd.xllcorner+(grd.cellsize*grd.ncols),grd.ncols);
            ycords = fliplr(linspace(grd.yllcorner,grd.yllcorner+(grd.cellsize*grd.nrows),grd.nrows));
            [X,Y] = meshgrid(xcords,ycords);
            firstflag=0;
        end

        distance_diff = ([X-gauge_xcord].^2+[Y-gauge_ycord].^2).^0.5;
        [~, ind] = min(distance_diff(:));
        radar_rain = daily_rain_mm(ind);
        data_gague_radar(gauge_index).radar_rain =  radar_rain;
    end
end

% Cut low quality gauge data
% indexim = find([data_gague_radar.gauge_flag]>=1);
% data_gague_radar(indexim)=[];

indexim =  find([data_gague_radar.in_radar]==0);
data_gague_radar(indexim)=[];
% dates_above_thrshold = DaliyaReport.MaxDis_date(moderate_indexim);

%% Unified the data per events (days before and after peak flow)
data_gague_radar([data_gague_radar(:).radar_rain]==0)=[];
gauges = unique([data_gague_radar.id]);

data_gague_radar_unified= struct();
cnt=1
for i = 1:length(dates_above_thrshold) 
    for j = 1:length(gauges)
        gauge_index = gauges(j);
        date_in = dates_above_thrshold(i);
        indexim = find(abs(date_in-[data_gague_radar.time])<=1 & [data_gague_radar.id] == gauge_index );
        
        if isempty(indexim) %| length(indexim)<3
            disp(['date' num2str(i) ' gauge ' num2str(j)])
            continue
        end
        
        radar_rain = sum([data_gague_radar(indexim).radar_rain]);
        gauge_rain = sum([data_gague_radar(indexim).gauge_rain]);

        data_gague_radar_unified(cnt).id = gauge_index;
        data_gague_radar_unified(cnt).name = data_gague_radar(indexim(1)).name;
        data_gague_radar_unified(cnt).x = data_gague_radar(indexim(1)).x;
        data_gague_radar_unified(cnt).y = data_gague_radar(indexim(1)).y;
        data_gague_radar_unified(cnt).time = date_in;
        data_gague_radar_unified(cnt).radar_rain = radar_rain;
        data_gague_radar_unified(cnt).gauge_rain = gauge_rain;
        cnt=cnt+1;

    end
end 

%% Plot the full data - WITHOUT BIAS CORRECTION  !data_gague_radar_unified!

[FSE,BIAS,CORR,RMSE,NSD,NASH]=statCalc([data_gague_radar_unified.gauge_rain],[data_gague_radar_unified.radar_rain]);
mdl = fitlm([data_gague_radar_unified.radar_rain],[data_gague_radar_unified.gauge_rain,]);
r2 = mdl.Rsquared.Ordinary;
rmsd = rms([data_gague_radar_unified.radar_rain]-[data_gague_radar_unified.gauge_rain,]);

% without statcalc plot from here[data_gague_radar.radar_rain,]
close all

f1 = figure;plot1= scatter([data_gague_radar_unified.gauge_rain],[data_gague_radar_unified.radar_rain],20,[0.5,0.4,0.8],'filled');hold on
% with bias correction

plot2 = plot([0,120],[0,120],'k','LineWidth',2);
xlabel('Gauge');ylabel('Radar')
plot2.Color(4) = 0.4;
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.5         );
legend([plot1,plot2],{'Daily amount [mm]','1:1'},'Location','northwest','FontSize',10);
ax = gca;
ax.FontSize = 20; 
title('radar - gauge correlation per event','Interpreter','none','FontSize',14)
xlim([0,120])
ylim([0,120])
num_of_events = length(unique([data_gague_radar_unified.time]));
text(0.7,0.2,[{['R^{2} = ' num2str(round(r2,2))]},{['RMSD = ' num2str(round(rmsd,1))]},{['bias = ' num2str(round(BIAS,3))]}],...
            'FontSize',14,'Color',[0.5 0.5 0.5],'Unit','normalized')
print(f1, 'C:\Users\raznu\Development\RESEARCH\Raanana\figures\Radar_Gauges_Calibration\radar_gauga_perEvent_with_bias.jpg', '-djpeg', '-r600')


%% Plot the full data - WITHOUT BIAS CORRECTION 

[FSE,BIAS,CORR,RMSE,NSD,NASH]=statCalc([data_gague_radar.gauge_rain],[data_gague_radar.radar_rain]);
mdl = fitlm([data_gague_radar.radar_rain],[data_gague_radar.gauge_rain,]);
r2 = mdl.Rsquared.Ordinary;
rmsd = rms([data_gague_radar.radar_rain]-[data_gague_radar.gauge_rain,]);

% without statcalc plot from here[data_gague_radar.radar_rain,]
% close all

f1 = figure;plot1= scatter([data_gague_radar.gauge_rain],[data_gague_radar.radar_rain],20,[0.5,0.4,0.8],'filled');hold on
% with bias correction

plot2 = plot([0,120],[0,120],'k','LineWidth',2);
xlabel('Gauge');ylabel('Radar')
plot2.Color(4) = 0.4;
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.5         );
legend([plot1,plot2],{'Daily amount [mm]','1:1'},'Location','northwest','FontSize',10);
ax = gca;
ax.FontSize = 20; 
title('radar - gauge correlation per day','Interpreter','none','FontSize',14)
xlim([0,120])
ylim([0,120])
num_of_events = length(unique([data_gague_radar.time]));
text(0.7,0.2,[{['R^{2} = ' num2str(round(r2,2))]},{['RMSD = ' num2str(round(rmsd,1))]},{['bias = ' num2str(round(BIAS,3))]}],...
            'FontSize',14,'Color',[0.5 0.5 0.5],'Unit','normalized')
print(f1, 'C:\Users\raznu\Development\RESEARCH\Raanana\figures\Radar_Gauges_Calibration\radar_gauga_perDay_with_bias.jpg', '-djpeg', '-r600')




%% Plot the full data - WITH BIAS CORRECTION
radar_correction = [data_gague_radar.radar_rain]./BIAS

[FSE,BIAS,CORR,RMSE,NSD,NASH]=statCalc([data_gague_radar.gauge_rain],radar_correction);
mdl = fitlm(radar_correction,[data_gague_radar.gauge_rain,]);
r2 = mdl.Rsquared.Ordinary;
rmsd = rms(radar_correction-[data_gague_radar.gauge_rain,]);
% close all

f1 = figure;plot1= scatter([data_gague_radar.gauge_rain],radar_correction,20,[0.5,0.4,0.8],'filled');hold on
plot2 = plot([0,120],[0,120],'k','LineWidth',2);
xlabel('Gauge');ylabel('Radar')
plot2.Color(4) = 0.4;
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      ,...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.5         );
legend([plot1,plot2],{'Daily amount [mm]','1:1'},'Location','northwest','FontSize',10);

ax = gca;
ax.FontSize = 20; 
title('radar - gauge correlation per day; bias correction','Interpreter','none','FontSize',14)                                                                             
xlim([0,120])
ylim([0,120])
num_of_events = length(unique([data_gague_radar.time]));
text(0.7,0.2,[{['R^{2} = ' num2str(round(r2,2))]},{['RMSD = ' num2str(round(rmsd,1))]},{['bias = ' num2str(round(BIAS,3))]}],...
            'FontSize',14,'Color',[0.5 0.5 0.5],'Unit','normalized')
print(f1, 'C:\Users\raznu\Development\RESEARCH\Raanana\figures\Radar_Gauges_Calibration\radar_gauga_perDay_with_bias_after_bias_correction.jpg', '-djpeg', '-r600')
%% Event radar-gauge multi-plot 

events = unique([data_gague_radar.peak_time]);
event_struct = struct();
cnt_fig=1;
f2 = figure('units','normalized','outerposition',[0.05 0.02 0.95 0.99]);
subplot(3,7,1);hold on;

for i = 1:length(events)
    event = events(i);
    indexim = find([data_gague_radar.peak_time] ==  event);

    event_struct(i).date=event;
    event_struct(i).gauges = [{data_gague_radar(indexim).name}];
    event_struct(i).gaugues_rain = [data_gague_radar(indexim).gauge_rain];
    event_struct(i).radar_rain = [data_gague_radar(indexim).radar_rain];

    [FSE,BIAS,CORR,RMSE,NSD,NASH]=statCalc([data_gague_radar(indexim).gauge_rain],[data_gague_radar(indexim).radar_rain]);
    
    event_struct(i).bias =BIAS;
    mdl = fitlm([data_gague_radar(indexim).radar_rain],[data_gague_radar(indexim).gauge_rain,]);
    r2 = mdl.Rsquared.Ordinary;
    event_struct(i).rsquare = r2;
    event_struct(i).rmsd = RMSE;
    subplot(6,6,cnt_fig);hold on;     
    plot1= scatter([data_gague_radar(indexim).gauge_rain,],[data_gague_radar(indexim).radar_rain],20,[0.5,0.4,0.8],'filled');hold on
    plot2 = plot([0,150],[0,150],'k','LineWidth',2);
    xlabel('Gauge');ylabel('Radar')
    plot2.Color(4) = 0.4;
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      ,...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'LineWidth'   , 1.5         );
%     legend([plot1,plot2],{'Daily amount [mm]','1:1'},'Location','northwest')
    ax = gca;
%     ax.FontSize = 10; 
%     title('radar - gauge correlation per day; bias correction','Interpreter','none','FontSize',14) 
    xlim([0,max(max([data_gague_radar(indexim).radar_rain],[data_gague_radar(indexim).gauge_rain,]))+2])
    ylim([0,max(max([data_gague_radar(indexim).radar_rain],[data_gague_radar(indexim).gauge_rain,]))+2])
    text(0.69,0.28,[{['R^{2} = ' num2str(round(r2,2))]},{['RMSD = ' num2str(round(RMSE,1))]},{['bias = ' num2str(round(BIAS,2))]}],...
                'FontSize',8,'Color',[0.5 0.5 0.5],'Units','normalized')
    title(datestr(event));
    cnt_fig=cnt_fig+1;

end

print(f2, 'C:\Users\raznu\Development\RESEARCH\Raanana\figures\Radar_Gauges_Calibration\Multiplot_radar_gauga_perEvent_with_bias.jpg', '-djpeg', '-r600')

save('D:\MY_CODES\Raanana_urban_floods\Radar_Files_Extraction\Gauges_Radar_raanana_byevent_bigger_than6.5to13dis.mat','data_gague_radar')
tt = datevec([data_gague_radar.time]);

%% 

%% Plot the full data - WITHOUT BIAS CORRECTION  !data_gague_radar_nozero!
%%   
% data_gague_radar_unified_origi =data_gague_radar_unified;
% data_gague_radar = data_gague_radar_unified;

% load('D:\Yuval\Phd\BARD\WatershedScale\data\Gauges_Radar_RM_byevent.mat')

% data_gague_radar_origi = data_gague_radar;

% Logical indexing to filter rows

% %%  filter the out no rain events:
% % Convert the structure to a table
% data_gague_radar_nozero = struct2table(data_gague_radar_unified);
% % Filter out rows where 'gauge_rain' is 0
% data_gague_radar_nozero = data_gague_radar_nozero(data_gague_radar_nozero.gauge_rain ~= 0, :);

% [FSE,BIAS,CORR,RMSE,NSD,NASH]=statCalc([data_gague_radar_nozero.gauge_rain],[data_gague_radar_nozero.radar_rain]);
% mdl = fitlm([data_gague_radar_nozero.radar_rain],[data_gague_radar_nozero.gauge_rain,]);
% r2 = mdl.Rsquared.Ordinary;
% rmsd = rms([data_gague_radar_nozero.radar_rain]-[data_gague_radar_nozero.gauge_rain,]);
% 
% % without statcalc plot from here[data_gague_radar.radar_rain,]
% close all
% 
% f1 = figure;plot1= scatter([data_gague_radar_nozero.gauge_rain],[data_gague_radar_nozero.radar_rain],20,[0.5,0.4,0.8],'filled');hold on
% % with bias correction

% plot2 = plot([0,120],[0,120],'k','LineWidth',2);
% xlabel('Gauge');ylabel('Radar')
% plot2.Color(4) = 0.4;
% set(gca, ...
%   'Box'         , 'on'     , ...
%   'TickDir'     , 'in'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'XMinorTick'  , 'off'      , ...
%   'YMinorTick'  , 'off'      , ...
%   'YGrid'       , 'off'      , ...
%   'XGrid'       , 'off'      ,...
%   'XColor'      , [.3 .3 .3], ...
%   'YColor'      , [.3 .3 .3], ...
%   'LineWidth'   , 1.5         );
% legend([plot1,plot2],{'Daily amount [mm]','1:1'},'Location','northwest','FontSize',10);
% ax = gca;
% ax.FontSize = 20; 
% title('radar - gauge correlation','Interpreter','none','FontSize',14)
% xlim([0,120])
% ylim([0,120])
% num_of_events = length(unique([data_gague_radar.time]));
% text(0.7,0.2,[{['R^{2} = ' num2str(round(r2,2))]},{['RMSD = ' num2str(round(rmsd,1))]},{['bias = ' num2str(round(BIAS,3))]}],...
%             'FontSize',14,'Color',[0.5 0.5 0.5],'Unit','normalized')
% print(f1, 'C:\Users\raznu\Development\RESEARCH\Raanana\Radar_IMS\radar_gauga_with_bias.jpg', '-djpeg', '-r600')






%% load gauge data from IMS to local 
% raz = [];

% Insert manually the gague data from IMS xlxs files

%% 
% raz(:,1) = x2mdate(raz(:,1));
% station_id = 134000;
% 
% 
% sation_index = find([gaugeData.id] ==  station_id);
% old_rain_data  = gaugeData(sation_index).vals;
% old_time_data  = gaugeData(sation_index).time;
% 
% new_rain_data = [old_rain_data ;raz(:,2)];
% new_time_data = [old_time_data ;raz(:,1)];
% 
% gaugeData(sation_index).vals = new_rain_data;
% gaugeData(sation_index).time = new_time_data;
% 