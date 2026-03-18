
clear all
LoadCommonSettings
%%
% ########## altimeter data to choose ##########
L2_dir = '/home/fehelers/ownCloud/sarsim4rofi/Data/altimetry/ROFI_Level2';
%mission = 'S3B';
mission = 'S3A'
orbit = '370';
%orbit = '279'
freq = 'f_120Hz';
product = 'UF';

latT.('S3A370') = 52.25
latT.('S3B370') = 52

L2_files   = dir(fullfile(L2_dir,[mission,'*_SR_1_SRA_A*',orbit,'*']));

%% ################# load simulation data

if strcmp(mission,'S3A')&strcmp(orbit,'370')
    FName = '/home/fehelers/ownCloud/sarsim4rofi/Data/simulations/S3A370/DFM_OUTPUT_DCSM-FM_S3A370_postprocessed.nc'
    cycles = 38:1:49; % timing/cycles to process
    lat_plot = [52 53];
elseif strcmp(mission,'S3B')&strcmp(orbit,'370')
    FName = '/home/fehelers/ownCloud/sarsim4rofi/Data/simulations/S3B370/DFM_OUTPUT_DCSM-FM_S3B370_postprocessed.nc'
    cycles = 19:1:33; % timing/cycles to process
    lat_plot = [51.8 53];
end

sim=[];

simtimes = ncread(FName,'time'); % hours since 2018-12-01T10:00:00
simtimes = datetime(2018,12,01,simtimes+10,0,0);

sim.LAT = ncread(FName,'lat')
sim.LON = ncread(FName,'lon')
sim.TIME = simtimes;%ncread(FName,'time')
sim.SSHi = ncread(FName,'mesh2d_s1')
sim.SAL = ncread(FName,'mesh2d_sa1')
sim.SAL = squeeze(sim.SAL(20,:,:,:)) % consider surface salinity only
%sim.RHO = ncread(FName,'mesh2d_rho') % omit rho for now
sim.H_STERIC = ncread(FName,'hsteric')


%% ############ load and plot altimetry data and colocate the simulation to it

N = 2500; % approximate number of values in altimeter dataset at this posting rate
M = numel(cycles)

ds = [] % struct with data arrays

varnames = {'LAT','LON','SSHi_alt','SWH_alt','Geoid','SSHi_sim','SAL_sim','H_steric_sim'}

for name = varnames
    ds.(name{1}) = NaN*zeros(M,N);
end
ds.('StartTime') = {};


for k = 1:numel(cycles)
    
    cycle = num2str(cycles(k));
    altdata = load_CS2_altdata(L2_dir, mission, orbit, cycle, freq, product);
    N_alt = numel(altdata.LAT)

    % ################# interpolate sim data to altdata altimetry track

    [val,time_ind] = min(abs(datenum(sim.TIME)-datenum(altdata.StartTime)));
    if 24*60*val < 10
        % interpolate the 
    %     altdata.SSHi_sim        = griddata(sim.LAT,sim.LON,sim.SSHi(time_ind,:,:),altdata.LAT,altdata.LON,'linear');
    %     altdata.H_STERIC_sim    = griddata(sim.LAT,sim.LON,sim.H_STERIC(time_ind,:,:),altdata.LAT,altdata.LON,'linear');
    %     altdata.SAL_sim         = griddata(sim.LAT,sim.LON,sim.SAL(time_ind,:,:),altdata.LAT,altdata.LON,'linear');
    %     
        [lat,lon] = ndgrid(sim.LAT,sim.LON);

        F = griddedInterpolant(lat,lon,squeeze(sim.SSHi(time_ind,:,:)));
        altdata.SSHi_sim = F(altdata.LAT,altdata.LON);

        F = griddedInterpolant(lat,lon,squeeze(sim.H_STERIC(time_ind,:,:)));
        altdata.H_STERIC_sim = F(altdata.LAT,altdata.LON);

        F = griddedInterpolant(lat,lon,squeeze(sim.SAL(time_ind,:,:)));
        altdata.SAL_sim = F(altdata.LAT,altdata.LON);
    else
        warning('no simulation data field with less than 10 min difference found')
    end

    
    % ######## get a Gaussian filtering kernel "simulating" the altimeter
    % resolution in the lat-lon image:
    % 1) ###################### get approximate cross-track direction within the lat-lon-space in
    % lat-lon coords
    clat = median((altdata.LAT)); % center of track Latitude
    clon = median((altdata.LON)); % center of track Longitude
    dlat = median(diff(altdata.LAT));
    dlon = median(diff(altdata.LON));
    
    wgs84 = wgs84Ellipsoid;
    [atrackE,atrackN,trackZ] = geodetic2enu(clat+dlat,clon+dlon,0,clat,clon,0,wgs84);
%     figure;
%     plot([0 atrackE],[0 atrackN])
    
    ctrackE = atrackN;
    ctrackN = -atrackE;
    
    scale = 1500; % ~3 km footprint that determines leading edge (may be tweaked)
    ctrackNn = scale*ctrackN./sqrt(ctrackN^2+ctrackE^2);
    ctrackEn = scale*ctrackE./sqrt(ctrackN^2+ctrackE^2);
%     atrackNn = scale*atrackN./sqrt(atrackN^2+atrackE^2);
%     atrackEn = scale*atrackE./sqrt(atrackN^2+atrackE^2);

    [clat_1,clon_1,~] = enu2geodetic(ctrackEn,ctrackNn,0,clat,clon,0,wgs84);
    [clat_2,clon_2,~] = enu2geodetic(-ctrackEn,-ctrackNn,0,clat,clon,0,wgs84);
    [alat_1,alon_1,~] = enu2geodetic(atrackE/2,atrackN/2,0,clat,clon,0,wgs84);
    [alat_2,alon_2,~] = enu2geodetic(-atrackE/2,-atrackN/2,0,clat,clon,0,wgs84);

    % 2) #################### define a filter kernel according to the
    % lat-lon grid that resembles the resolution of the SAR data
    
    r_lat = mean(diff(sim.LAT));
    r_lon = mean(diff(sim.LON));
    
    c_vec = [((clat_1-clat)/r_lat) ((clon_1-clon)/r_lon)]';
    a_vec = [((alat_1-clat)/r_lat) ((alon_1-clon)/r_lon)]';
    sigma_c = sqrt( c_vec(1).^2  +  c_vec(2).^2 ); % define sigma as the length of the 1.5 km cross-track vector in pixel units
    sigma_a = sqrt( a_vec(1).^2  +  a_vec(2).^2 ); % define sigma as the length of the 330/2 m along-track vector in pixel units

    A = [c_vec  a_vec]; % transformation matrix from cross_track-along_track to lat-lon
    S = diag([1 1]); % matrix of unit variances?
    
    Sfilter = A*S*(A');
    
    xl = 5;
    yl = 20;
    xs = -xl:xl;
    ys = -yl:yl;
    
    K = zeros(size(xs,2),size(ys,2));
    
    for i = xs
        for j = ys
            K(i+xl+1,j+yl+1) = exp(-0.5*[i j]*inv(Sfilter)*[i j]');
        end
    end
    
    K = K./sum(K(:));
    % ##################### filter found
    
    
    % now filter the output fields prior to interpolation onto the
    % satellite groudn track
    if 24*60*val < 10
        % interpolate the 
    %     altdata.SSHi_sim        = griddata(sim.LAT,sim.LON,sim.SSHi(time_ind,:,:),altdata.LAT,altdata.LON,'linear');
    %     altdata.H_STERIC_sim    = griddata(sim.LAT,sim.LON,sim.H_STERIC(time_ind,:,:),altdata.LAT,altdata.LON,'linear');
    %     altdata.SAL_sim         = griddata(sim.LAT,sim.LON,sim.SAL(time_ind,:,:),altdata.LAT,altdata.LON,'linear');
    %     
        [lat,lon] = ndgrid(sim.LAT,sim.LON);
        
        field = squeeze(sim.SSHi(time_ind,:,:));
        field = conv2(field,K,'same');
        F = griddedInterpolant(lat,lon,field);
        altdata.SSHi_smooth = F(altdata.LAT,altdata.LON);

        field = squeeze(sim.H_STERIC(time_ind,:,:));
        field = conv2(field,K,'same');
        F = griddedInterpolant(lat,lon,field);
        altdata.H_STERIC_smooth = F(altdata.LAT,altdata.LON);

        field = squeeze(sim.SAL(time_ind,:,:));
        field = conv2(field,K,'same');
        F = griddedInterpolant(lat,lon,field);
        altdata.SAL_smooth = F(altdata.LAT,altdata.LON);
    else
        warning('no simulation data field with less than 10 min difference found')
    end

    
    
    % make a plot of the intercomparison, including along-track / cross-track resolution and filter:
    fig = figure('units','inch','position',[0,0,18,5]);
    set(gcf,'color','w');
    
    subplot(1,3,1)
    % h_steric and satellite pass
    imagesc(sim.LON,sim.LAT,squeeze(sim.H_STERIC(time_ind,:,:)));hold on
    plot(altdata.LON,altdata.LAT,'k--');hold on
    
    % footprint dimensions
    plot([clon_1 clon_2],[clat_1 clat_2],'r-');hold on
    plot([alon_1 alon_2],[alat_1 alat_2],'g-');hold on
    
    % filter
    imagesc(r_lon*ys+clon-0.2,r_lat*xs+clat,K./max(K(:))*0.1);
    
    legend('satellite ground track','cross-track distance of ~3 km','along-track resolution');
    
    set(gca,'YDir','normal')
    h = colorbar();
    ylabel(h,'steric water level')
    ylabel('Latitude (deg N)')
    xlabel('Longitude (deg E)')
    ylim(lat_plot)
    xlim([3.5 4.7])
    caxis([-0.1,0.1])
    
    
    subplot(1,3,2)
    plot(altdata.SSHi,altdata.LAT);hold on
    plot(altdata.SSHi_sim,altdata.LAT);hold on
    xlim([-1.5 1.5])
    ylim(lat_plot)
    xlabel('instantaneous sea surface height (m)')
    grid on
    legend('altimeter','model')

    title(['cycle: ' cycle])
    
    subplot(1,3,3)
    field = altdata.SSHi(:) - altdata.SSHi_sim(:);
    plot(field,altdata.LAT,'Color','#a6a6a6');hold on
    plot(movmedian(field,120,'omitnan'),altdata.LAT,'k','LineWidth',2);hold on

%     corr = altdata.COR.sea_state_bias_01_ku(altdata.LAT) %+ altdata.COR.inv_bar_cor_01(altdata.LAT)
%     field = altdata.SSHi(:) - corr(:) - altdata.SSHi_sim(:);
%     plot(field,altdata.LAT,'r-');hold on
%     plot(movmedian(field,20),altdata.LAT,'r','LineWidth',2);hold on

%     plot(altdata.SSHi(:) - altdata.SSHi_smooth(:),altdata.LAT,'r-');hold on
%     plot(movmedian(altdata.SSHi(:) - altdata.SSHi_smooth(:),20),altdata.LAT,'r','LineWidth',2);hold on

    plot(altdata.H_STERIC_sim,altdata.LAT,'LineWidth',2);hold on
    ylabel('difference (m)')

    grid on
    legend('','difference of SSHi (alt) - SSHi (sim) with ssb','steric water level','Location','south')
    xlim([-0.2 0.2])
    ylim(lat_plot)
    
    saveas(gcf,['/home/fehelers/ownCloud/sarsim4rofi/Data/comparison_sim_alt/' mission orbit '_' cycle '_' datestr(sim.TIME(time_ind)) '.png'])
    close all
    % ############# save the altdata into ds fields
    %varnames = {'LAT','LON','SSHi_alt','SWH_alt','Geoid','SSHi_sim','SAL_sim','H_steric_sim'}
    
    ds.LAT(k,1:N_alt) = altdata.LAT;
    ds.LON(k,1:N_alt) = altdata.LON;
    ds.SSHi_alt(k,1:N_alt) = altdata.SSHi;
    ds.SWH_alt(k,1:N_alt) = altdata.SWH;
    ds.Geoid(k,1:N_alt) = altdata.Geoid;
    ds.SSHi_sim(k,1:N_alt) = altdata.SSHi_sim;
    ds.SAL_sim(k,1:N_alt) = altdata.SAL_sim;
    ds.H_steric_sim(k,1:N_alt) = altdata.H_STERIC_sim;
    ds.StartTime{k} = altdata.StartTime;
end

%% play a bit with the dataset ds
rgb = [ ...
    94    79   162
    50   136   189
   102   194   165
   171   221   164
   230   245   152
   255   255   191
   254   224   139
   253   174    97
   244   109    67
   213    62    79
   158     1    66  ] / 255;

rgb = interp1(1:11,rgb,1:0.05:11);

mask = (ds.LAT>latT.([mission orbit]))
for name = varnames
    ds.(name{1})(~mask)=NaN
end
% flag out altimeter measurements between -5m and +5 m (much more than tidal range!)
ds.SSHi_alt(((ds.SSHi_alt>5)|(ds.SSHi_alt<-5))) = NaN;

% difference of altimeter and proper model
ds.SSHi_diff = ds.SSHi_alt - ds.SSHi_sim;
% difference of altimeter and model - ~halosteric waterlevel contribution
% (mean removed from steric water level before)

ds.H_halosteric_sim = ds.H_steric_sim - mean(ds.H_steric_sim,2,'omitnan'); % as an approximation, the termosteric contribution is uniform and hence approximated as the mean of the total steric water level
ds.SSHi_diff_mod1 = ds.SSHi_alt - (ds.SSHi_sim - ds.H_halosteric_sim);
ds.SSHi_diff_mod2 = ds.SSHi_alt - (ds.SSHi_sim - ds.H_steric_sim);


x = 1:numel(ds.LAT(1,:));
mask_x = ~isnan(ds.LAT(1,:));
lat_interp = interp1(x(mask_x),ds.LAT(1,mask_x),x,'Linear','extrap');

fig = figure('units','inch','position',[0,0,18,12]);
set(gcf,'color','w');

subplot(3,3,1)
imagesc(lat_interp,cycles,ds.SSHi_sim)
caxis([-1.5 1.5])
title('model')
h = colorbar()
ylabel(h,'SSH (m)')
ylabel('#cycle')
xlabel('Latitude (deg N)')

subplot(3,3,2)
imagesc(lat_interp,cycles,ds.SSHi_alt)
caxis([-1.5 1.5])
title('altimeter')
h = colorbar()
ylabel(h,'SSH (m)')
ylabel('#cycle')
xlabel('Latitude (deg N)')

ax1 = subplot(3,3,4)
field = movmedian(ds.SSHi_diff,60,2,'omitnan')
imagesc(lat_interp,cycles,field(:,1:60:end))
% title({['altimeter - model'],['median bias: ' num2str(median(ds.SSHi_diff(:),'omitnan')*100,2) 'cm'],['1.4826\timesMAD (res=340 m): ' num2str(1.4826*mad(100*ds.SSHi_diff(:),1),2) 'cm'],['1.4826\timesMAD (res=7 km ): ' num2str(1.4826*mad(100*field,1),2) 'cm']})

title({ ['altimeter - model'], ...
        ['median bias: ' num2str(median(ds.SSHi_diff(:),'omitnan')*100,2) 'cm'], ...
        ['1.4826\timesMAD: ' num2str(1.4826*mad(100*ds.SSHi_diff(:),1),2) ' cm (on 340 m), ' num2str(1.4826*mad(100*field(:),1),2) ' cm (on 3.4 km)'], ...
        ['std: ' num2str(std(100*field(:),'omitnan'),2) ' cm (on 3.4 km)'], ...
        })
caxis([-0.1 0.1])
h = colorbar()
ylabel(h,'\DeltaSSH (m)')
colormap(rgb)
ylabel('#cycle')
xlabel('Latitude (deg N)')

ax1 = subplot(3,3,5)
field = movmedian(ds.SSHi_diff,60,2,'omitnan')-median(ds.SSHi_diff(:),'omitnan');
imagesc(lat_interp,cycles,field(:,1:60:end));hold on
title({ ['altimeter - model (median removed)'], ...
        ['1.4826\timesMAD: ' num2str(1.4826*mad(100*ds.SSHi_diff(:),1),2) ' cm (on 340 m), ' num2str(1.4826*mad(100*field(:),1),2) ' cm (on 3.4 km)'], ...
        ['std: ' num2str(std(100*field(:),'omitnan'),2) ' cm (on 3.4 km)'], ...
        })
drawnow;
hold on;
caxis([-0.1 0.1])
h = colorbar()
ylabel(h,'\DeltaSSH (m)')
colormap(rgb)
ylabel('#cycle')
xlabel('Latitude (deg N)')

ax1 = subplot(3,3,3)
imagesc(lat_interp,cycles,ds.H_steric_sim);hold on
title({ ['model: steric water level']})
drawnow;
hold on;
caxis([-0.1 0.1])
h = colorbar()
ylabel(h,'model: steric water level (m)')
colormap(rgb)
ylabel('#cycle')
xlabel('Latitude (deg N)')

ax1 = subplot(3,3,6)
field = movmedian(ds.SWH_alt,60,2,'omitnan');
imagesc(lat_interp,cycles,field(:,1:60:end));hold on
title({ ['altimeter: SWH']})
drawnow;
hold on;
caxis([-1 6])
h = colorbar()
ylabel(h,'significant wave height (m)')
colormap(rgb)
ylabel('#cycle')
xlabel('Latitude (deg N)')

ax1 = subplot(3,3,8)
field = movmedian(ds.SSHi_diff_mod1,60,2,'omitnan')-median(ds.SSHi_diff_mod1(:),'omitnan');
imagesc(lat_interp,cycles,field(:,1:60:end));hold on
title({ ['altimeter - model - H_{halosteric} (median removed)'], ...
        ['1.4826\timesMAD: ' num2str(1.4826*mad(100*ds.SSHi_diff_mod1(:),1),2) ' cm (on 340 m), ' num2str(1.4826*mad(100*field(:),1),2) ' cm (on 3.4 km)'], ...
        ['std: ' num2str(std(100*field(:),'omitnan'),2) ' cm (on 3.4 km)'], ...
        })
drawnow;
hold on;
caxis([-0.1 0.1])
h = colorbar()
ylabel(h,'\DeltaSSH (m)')
colormap(rgb)
ylabel('#cycle')
xlabel('Latitude (deg N)')

ax1 = subplot(3,3,7)
field = movmedian(ds.SSHi_diff_mod2,60,2,'omitnan')-median(ds.SSHi_diff_mod2(:),'omitnan');
imagesc(lat_interp,cycles,field(:,1:60:end));hold on
title({ ['altimeter - model - H_{steric} (median removed)'], ...
        ['1.4826\timesMAD: ' num2str(1.4826*mad(100*ds.SSHi_diff_mod2(:),1),2) ' cm (on 340 m), ' num2str(1.4826*mad(100*field(:),1),2) ' cm (on 3.4 km)'], ...
        ['std: ' num2str(std(100*field(:),'omitnan'),2) ' cm (on 3.4 km)'], ...
        })
drawnow;
hold on;
caxis([-0.1 0.1])
h = colorbar()
ylabel(h,'\DeltaSSH (m)')
colormap(rgb)
ylabel('#cycle')
xlabel('Latitude (deg N)')


saveas(gcf,['/home/fehelers/ownCloud/sarsim4rofi/Data/comparison_sim_alt/statistics_' mission orbit '.png'])

%% plot differences of altimeter and model 
figure;
ax1 = subplot(1,2,1)
% differences altimetry/model with median bias removed
field = movmedian(ds.SSHi_diff,60,2,'omitnan')-median(ds.SSHi_diff(:),'omitnan');
overpass_median = (median(field,2,'omitnan'));
%plot(field'-overpass_median','k')
plot(lat_interp,100*field','Color',[0, 0, 0, 0.2]); hold on
plot(lat_interp(1:60:end),mean(100*field(:,1:60:end),1),'k','Linewidth',3); hold on
ylim([-20 20])
grid on
title('SSH_{altimeter} - SSH_{model}')
ylabel('residual (cm)')
xlabel('Latitude')

ax2 = subplot(1,2,2)
% differences altimetry - (model-steric height) with median bias removed
field = movmedian(ds.SSHi_diff_mod2,60,2,'omitnan')-median(ds.SSHi_diff_mod2(:),'omitnan');
overpass_median = (median(field,2,'omitnan'));
%plot(field'-overpass_median','k')
plot(lat_interp,100*field','Color',[0, 0.5, 0.5, 0.2]);hold on
plot(lat_interp(1:60:end),mean(100*field(:,1:60:end),1),'Color',[0, 0.5, 0.5, 1],'Linewidth',3); hold on
ylim([-20 20])
grid on
title('SSH_{altimeter} - (SSH_{model}-H_{steric})')
xlabel('Latitude')

linkaxes([ax1 ax2],'xy')

saveas(gcf,['/home/fehelers/ownCloud/sarsim4rofi/Data/comparison_sim_alt/slope_' mission orbit '.png'])


%% play around with the difference data and the SWH
y = movmedian(ds.SSHi_diff,60,2,'omitnan');
x = movmedian(ds.SWH_alt,60,2,'omitnan');
scatter(x(:,1:60:end),y(:,1:60:end),'k.')


%% useful function definitions
function altdata = load_CS2_altdata(L2_dir,mission, orbit, cycle, freq, product)
    LoadCommonSettings
    %load L2 data at requested retracking frequency
    L2_file = dir(fullfile(L2_dir,[mission,'*_SR_1_SRA_A*',cycle,'_',orbit,'___*']));
    if ~isempty(L2_file)
        CS2 = load(fullfile(L2_file.folder,L2_file.name));

        altdata = CS2.CS2.(freq).(product);
        altdata.StartTime = CS2.CS2.StartTime;

        % load geoid as interpolant and save into the altdata struct
        load(fullfile(PathDATA,'ROFI','Geoid','NLGEO2018_Netherlands.mat'));
        Geoid = griddedInterpolant(SYN.LATgrd,SYN.LONgrd,SYN.NLGEO2018);
        clear SYN

        altdata.COR = CS2.CS2.COR;
        
        altdata.Geoid = Geoid(altdata.LAT,altdata.LON);
        altdata.SSHi = altdata.SSHi(:) - altdata.Geoid(:);
        
        % apply additional ssb correction as in L2 data file:
        corr = altdata.COR.sea_state_bias_01_ku(altdata.LAT); %+ altdata.COR.inv_bar_cor_01(altdata.LAT)
        altdata.SSHi = altdata.SSHi(:) - corr(:);
        %field = altdata.SSHi(:) - corr(:) - altdata.SSHi_sim(:);
        altdata.SSHi = altdata.SSHi(:);
        
    else
        warning('no altimeter data file with specified cycle found')
        altdata=[];
    end
end

