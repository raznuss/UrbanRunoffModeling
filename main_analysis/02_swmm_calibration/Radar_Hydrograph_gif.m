
clear

% CREATE GIF ANIMATION OF RADAR images and hydrograph over watershed

data = readtable("D:\Development\RESEARCH\Raanana\20191212_storm_data.csv");
data.Time = datenum(data.Time);
data = table2struct(data);

% Save the data to a .mat file
save("D:\Development\RESEARCH\Raanana\20191212_storm_data.mat", 'data');


%%
load("D:\Development\RESEARCH\Raanana\data\rain_radar\Row_radar_correct_MAT_files\20191213.mat")
load("D:\Development\RESEARCH\Raanana\20191212_storm_data.mat")
% data = data'
close all
timestep = 5;


latmin = event_file.yll - 9000;
latmax = latmin+1300*(size(event_file.rain,1)-1);

lonmin = event_file.xll-9000;
lonmax = lonmin+1300*(size(event_file.rain,2)-1);


rastersize= [size(event_file.rain,1),size(event_file.rain,2)];
R = maprefpostings([lonmin,lonmax],[latmin,latmax],rastersize);
R.ProjectedCRS= projcrs(2039);
W = worldFileMatrix(R);
R = maprasterref(W,rastersize,"cells");
S = shaperead("D:\Development\RESEARCH\Raanana\gis\GIS\watershead_for_radar\Raanana_wsp.shp");



% flowtimes = data(:,1);
% flow = data(:,2);
% raingauge = data(:,4);

flowtimes = [data.Time];
flow = [data.SWMMRunoff_CMS_];
raingauge = [data.rainfall_mm_h_];
rain = event_file.rain; % Spatail radar data

radar_times  = event_file.time;
t1 = radar_times(1);
t2 = radar_times(end);
times = datetime(datevec(t1)):minutes(timestep):datetime(datevec(t2));
times = datenum(times);

conve_flow = conv(flow,[1,1,1,1,1,1,1,1,1,1],'same')./10;

virtual_gauge = [];
new_flow_vec = [];
cnt=0;
for k=1:length(times)
    [val,indx] = min(abs(flowtimes-times(k)));
    if val <=0.005
        new_flow_vec(k)=flow(indx);
        virtual_gauge(k)=raingauge(indx);
        cnt=cnt+1;
    elseif k >1
        new_flow_vec(k)=nan;
        virtual_gauge(k)=nan;
        cnt=cnt+1;
    end
end

    
event_date = radar_times(1);
event_date_str = num2str(yyyymmdd(datetime(datestr(event_date))));
cnt=1;
h = figure('Renderer', 'painters', 'Position', [10 20 1000 800]);
axis tight manual % this ensures that getframe() returns a consistent size
set(gcf,'color','w');

outletcolor = 'red';
gaugecolor = 'blue';
rain_threshold_to_plot = 5; 
filename = ['D:\Development\RESEARCH\Raanana\figures\gif\' event_date_str  '.gif']; % SET YOU GIFF PATH
plot_jumps = 1; %plot in jumpts of two time steps
for j = 1:plot_jumps:length(new_flow_vec)%:1:size(rain,3)
    subplot(2,2,[1,2]);
    if j > size(rain,3)
        KK.Z(:)= 0;
    else
        KK.Z = rain(:,:,j);
    end
   
   
    rain = event_file.rain(:,:,j);
    rain(rain<rain_threshold_to_plot) = NaN;
    hh(2) =  mapshow(rain,R,'DisplayType','surface');hold on
    xlim([lonmin,lonmax])
    ylim([latmin,latmax])
    mapshow(S,'FaceColor','none');

    xticklabels('')
    yticklabels('')
    clim([rain_threshold_to_plot,40])
    ff = colorbar;
    ylabel(ff, 'Rain intensity [mm h^{-1}]')
    title(datestr(times(j)))
    ax=gca;
    ax.FontSize=20;
    subplot(2,2,[3]);
    plot(new_flow_vec(1:j),'k','LineStyle','-');hold on
    xlim([0,round((times(end)-times(1))*24*(60/timestep),0)])
    ylim([0,max(new_flow_vec+1)]);
    ylabel('Discharge [m^{3} s^{-1}]','FontSize' ,16);
    xlabel('Time [min]','FontSize' ,16);
    title('Discharge at outlet','FontSize',16,'Color',outletcolor);
    ax=gca;
    ax.FontSize=20;
    subplot(2,2,4);
    plot(virtual_gauge(1:j),'k','LineStyle','-');hold on
    xlim([0,round((times(end)-times(1))*24*(60/timestep),0)])
    ylim([0,max(virtual_gauge)+1]);
    ylabel('Rain intensity [mm h^{-1}]','FontSize', 16);
    xlabel('Time [min]','FontSize' ,16);
    title('Virtual rain gauge', 'FontSize',16,'Color', gaugecolor  )
    ax=gca;
    ax.FontSize=20;
    if j ==1
        labels = xticklabels;
        for l = 1:size(labels,1)
            ll = str2num(labels{l});
            labels{l} = num2str(ll*timestep*plot_jumps);
        end
    end
    subplot(2,2,4)
    xticklabels(labels);

    subplot(2,2,3)
    xticklabels(labels);
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
  
    if cnt == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.01); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.01); 
    end 
    cnt=cnt+1;
end
close all
