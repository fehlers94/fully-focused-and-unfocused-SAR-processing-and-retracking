%% check the GEOID model from Cornelis

clear all
LoadCommonSettings
%% will load a structure called "SYN"
%PathDATA = '/home/fehelers/PhD Delft/Projects/FFSAR4ROFI/Data';
load(fullfile(PathDATA,'ROFI','Geoid','NLGEO2018_Netherlands.mat'))

% figure;
% scatter(SYN.LONgrd(:),SYN.LATgrd(:),[],SYN.NLGEO2018(:),'filled');hold on
% colorbar()
% caxis([40 52])

% make gridded Interpolant
Geoid = griddedInterpolant(SYN.LATgrd,SYN.LONgrd,SYN.NLGEO2018)
% lat_g = (47:0.05:58)';
% lon_g = 8*ones(1,numel(lat_g))';
% scatter(lon_g,lat_g,[],Geoid(lat_g,lon_g)','filled')
% colorbar()
% caxis([40 52])
clear SYN
%% #################################################################################
L2_dir = '/home/fehelers/ownCloud/sarsim4rofi/Data/altimetry/ROFI_Level2';


%% plot data and tide gauge positions on map
L2_files   = dir(fullfile(L2_dir,['*_SR_1_SRA_A*']));

figure;
for i=1:5:numel(L2_files)
    load(fullfile(L2_files(i).folder,L2_files(i).name));
    geoplot(CS2.f_080Hz.UF.LAT,CS2.f_080Hz.UF.LON);hold on
end
geobasemap();

%load and plot tide gauge data 
TGdata = load('/home/fehelers/ownCloud/sarsim4rofi/Data/tide_gauges/OBS_NLTGs_wlseries_199701_202101.mat');
TGdata = TGdata.OBS;

geoplot(TGdata.Lat,TGdata.Lon,'b.','LineWidth',3);

% choose tide gauge index and tracks for comparison

% Offshore:
TGidx = find((52.75<TGdata.Lat)&(TGdata.Lat<53)&(4<TGdata.Lon)&(TGdata.Lon<4.5));

% within Rhine Mouth:
%TGidx = find((51.9<TGdata.Lat)&(TGdata.Lat<52)&(4<TGdata.Lon)&(TGdata.Lon<4.25));

%track(1).mission = 'S3B';
%track(1).orbit = '279';
track(1).mission = 'S3A';
track(1).orbit = '370';

geoplot(TGdata.Lat(TGidx),TGdata.Lon(TGidx),'go','LineWidth',3);
%% make correlation analysis

TGdistance_threshold = 10e3; % 8 km distance
sampling_freq = 'f_080Hz';

L2_files   = dir(fullfile(L2_dir,[track.mission,'*_SR_1_SRA_A*',track.orbit,'*']));
%L2_files   = dir(fullfile(L2_dir,['*_SR_1_SRA_A*']));
lat_threshold = 52;

SSH_noise_threshold = 0.2;
movmedian_win_size = 20;

TG_height = [];
ALT_height = [];
D_obs = [];

for i = 1:numel(L2_files)
    % load CS2 data and which array to take
    load(fullfile(L2_files(i).folder,L2_files(i).name));
    ALTdata = CS2.(sampling_freq);
    
    % substract geoid
    ALTdata.UF.Geoid = Geoid(ALTdata.UF.LAT,ALTdata.UF.LON);
    ALTdata.UF.SSHi = ALTdata.UF.SSHi(:) - ALTdata.UF.Geoid(:);
    
    % outlier with respect to moving median of 20 samples
    mask_outlier = (abs(ALTdata.UF.SSHi-movmedian(ALTdata.UF.SSHi,movmedian_win_size,'omitnan')) > SSH_noise_threshold);
    ALTdata.UF.SSHi(mask_outlier) = NaN;
    
    % determine the distance to the tide gauge and make mask
    dist2TG = abs(1e3*deg2km(distance(ALTdata.UF.LAT,ALTdata.UF.LON,TGdata.Lat(TGidx),TGdata.Lon(TGidx))));
    dist2TG = dist2TG(:);
    dist_mask = (dist2TG<TGdistance_threshold);
    
    % determine latitude mask
    mask_lat = (ALTdata.UF.LAT > lat_threshold)';
    
    % mask
    mask = dist_mask(:)&mask_lat(:);
    
    % find closest time index in TG observations
    TGtime_idx = interp1(TGdata.times{TGidx},1:numel(TGdata.times{TGidx}),datenum(CS2.StartTime),'nearest');
	if (~isnan(TGtime_idx))&(sum(mask)>0)
        %plot(TGdata.wl{TGidx}(TGtime_idx)*ones(1,numel(ALTdata.UF.SSHi(mask))),ALTdata.UF.SSHi(mask),'b.');hold on
        ALT_height = [ALT_height(:); ALTdata.UF.SSHi(mask)];
        TG_height = [TG_height(:); TGdata.wl{TGidx}(TGtime_idx)*ones(numel(ALTdata.UF.SSHi(mask)),1)];
        D_obs = [D_obs(:); dist2TG(mask)];
    end
end
nan_mask = (~isnan(TG_height))&(~isnan(ALT_height))&(~isnan(D_obs));
outlier_mask = ones(size(nan_mask));%(abs(ALT_height-TG_height)<0.5)
%outlier_mask = (abs(ALT_height-TG_height)<0.5);
mask = nan_mask&outlier_mask;
TG_height = TG_height(mask);
ALT_height = ALT_height(mask);
D_obs = D_obs(mask);

figure;
subplot(2,1,1)
plot(TG_height,ALT_height,'b.');hold on
plot([-1.5 1.5],[-1.5 1.5],'r--')
P = polyfit(TG_height,ALT_height,1)
plot(TG_height,P(1)*TG_height+P(2),'k-')

xlabel('Tide Gauge')
ylabel('Altimeter')

residuals = ALT_height - (P(1)*TG_height+P(2));
std_residuals = std(residuals);
[RR,PP,RRL,RRU] = corrcoef(TG_height,ALT_height);

title({[track.mission, track.orbit];['mx+b; m=',num2str(P(1)),' m; b=',num2str(P(2)),' m; ','R=',num2str(RR(1,2))];['standard deviation of residuals ',num2str(std_residuals),' m; TG distance threshold ',num2str(TGdistance_threshold),' m']})

grid on

subplot(2,1,2)
plot(D_obs,abs(ALT_height-TG_height),'b.');hold on
grid on

%% check data along one track after geoid substraction, mean substraction
%% track S3B 370
mission = 'S3B'
orbit = '_370_'
L2_files   = dir(fullfile(L2_dir,[mission,'*_SR_1_SRA_A*',orbit,'*']));

lat_threshold = 52.03; % degree
SSH_noise_threshold = 0.2; % m, 20cm noise threshold
movmedian_win_size = 20;

figure;
for i=1:numel(L2_files)
    load(fullfile(L2_files(i).folder,L2_files(i).name))
    %CS2.f_020Hz_matched.EUM.LAT = CS2.f_020Hz_matched.UF.LAT;
    %CS2.f_020Hz_matched.EUM.LON = CS2.f_020Hz_matched.UF.LON;
    data = CS2.f_020Hz_matched.UF;
    
    % ############ for SSHi ###############
    
    % apply different masks
    % latitude
    mask_lat = (data.LAT < lat_threshold)';
    data.SSHi(mask_lat) = NaN;
    data.mask = mask_lat;
    
    % outlier with respect to moving median of 20 samples
    mask_outlier = (abs(data.SSHi-movmedian(data.SSHi,movmedian_win_size,'omitnan')) > SSH_noise_threshold);
    data.SSHi(mask_outlier) = NaN;
    data.mask = mask_outlier|data.mask|isnan(data.SSHi);
    
    % substract geoid
    data.Geoid = Geoid(data.LAT,data.LON);
    data.SSHi = data.SSHi(:) - data.Geoid(:);
    
    % substract linear slope
    P = polyfit(data.LAT(~data.mask),data.SSHi(~data.mask),1);
    %data.SSHi_slope = P(1)*data.LAT.^2 + P(2)*data.LAT+P(3);
    data.SSHi_slope = P(1)*data.LAT + P(2);
    data.SSHi_detrended = data.SSHi - data.SSHi_slope';
    
    % form cumulative sum of estimates
    data.SSHi_cumsum = cumsum(data.SSHi_detrended,'omitnan');

    % ######### instr corrected range #########
    % substract geoid
    data.Geoid = Geoid(data.LAT,data.LON);
    data.r = data.ALT(:) - data.instr_corrected_range(:) - data.Geoid(:);
    
    % substract linear slope
    P = polyfit(data.LAT(~data.mask),data.r(~data.mask),1);
    %data.SSHi_slope = P(1)*data.LAT.^2 + P(2)*data.LAT+P(3);
    data.r_slope = P(1)*data.LAT + P(2);
    data.r_detrended = data.r - data.r_slope';
    
    
    % ######### plotting ###########
%     subplot(3,1,1)
%     title('SSHi minus Geoid')
%     plot(data.LAT,data.SSHi);hold on
%     
%     subplot(3,1,2)
%     title('SSHi minus Geoid; detrended')
%     plot(data.LAT,data.SSHi_detrended);hold on
%     
%     subplot(3,1,3)
%     title('SSHi minus Geoid; detrended; 3 km moving median ')
%     plot(data.LAT,movmedian(data.SSHi_detrended,movmedian_win_size),'.');hold on

    ax(i) = subplot((numel(L2_files) + (4-mod(numel(L2_files),4)))./4,4,i)
    if i==1
        title('SSHi minus Geoid; detrended; 3 km moving median ')
    end
    
    %plot(data.LAT,data.SSHi_detrended-data.r_detrended,'-');hold on % gives you the part of the corrections that is not explained by a linear slope
    plot(data.LAT,data.SSHi_detrended,'.');hold on % detrended SSHi
    plot(data.LAT,movmedian(data.SSHi_detrended,60,'omitnan'),'-');hold on % detrended SSHi
    title(i)
end
linkaxes(ax,'xy')
xlim([52 53])
ylim([-0.15 0.15])

%% along-track noise resolution etc. analysis track S3B 370
clear ax
clear ax2

mission = 'S3B'
orbit = '_370_'
L2_files   = dir(fullfile(L2_dir,[mission,'*_SR_1_SRA_A*',orbit,'*']));

lat_threshold = 52.5; % degree
SSH_noise_threshold = 1; % m, 20cm noise threshold
movmedian_win_size = 240; % choose in accordance with the Hz product 

n=4
L2_files = L2_files(1:n);

figure;
for i=1:numel(L2_files)
    load(fullfile(L2_files(i).folder,L2_files(i).name))
    %CS2.f_020Hz_matched.EUM.LAT = CS2.f_020Hz_matched.UF.LAT;
    %CS2.f_020Hz_matched.EUM.LON = CS2.f_020Hz_matched.UF.LON;
    data = CS2.f_240Hz.UF;
    %data.SSHi = double(data.SWH);
    
    % ############ for SSHi ###############
    
    % apply different masks
    % latitude
    mask_lat = (data.LAT < lat_threshold)';
    data.SSHi(mask_lat) = NaN;
    data.mask = mask_lat;
    
    % outlier with respect to moving median of 20 samples
    mask_outlier = (abs(data.SSHi-movmedian(data.SSHi,movmedian_win_size,'omitnan')) > SSH_noise_threshold);
    data.SSHi(mask_outlier) = NaN;
    data.mask = mask_outlier|data.mask|isnan(data.SSHi);
    
    % substract geoid
    data.Geoid = Geoid(data.LAT,data.LON);
    data.SSHi = data.SSHi(:) - data.Geoid(:);
    
    % substract linear slope
    P = polyfit(data.LAT(~data.mask),data.SSHi(~data.mask),2);
    data.SSHi_slope = P(1)*data.LAT.^2 + P(2)*data.LAT+P(3);
    %data.SSHi_slope = P(1)*data.LAT + P(2);
    data.SSHi_detrended = data.SSHi - data.SSHi_slope';
    
    % form cumulative sum of estimates
    data.SSHi_cumsum = cumsum(data.SSHi_detrended,'omitnan');

    % ######### instr corrected range #########
    % substract geoid
    data.Geoid = Geoid(data.LAT,data.LON);
    data.r = data.ALT(:) - data.instr_corrected_range(:) - data.Geoid(:);
    
    % substract linear slope
    P = polyfit(data.LAT(~data.mask),data.r(~data.mask),1);
    %data.SSHi_slope = P(1)*data.LAT.^2 + P(2)*data.LAT+P(3);
    data.r_slope = P(1)*data.LAT + P(2);
    data.r_detrended = data.r - data.r_slope';
    
    
    % ######### plotting ###########
%     subplot(3,1,1)
%     title('SSHi minus Geoid')
%     plot(data.LAT,data.SSHi);hold on
%     
%     subplot(3,1,2)
%     title('SSHi minus Geoid; detrended')
%     plot(data.LAT,data.SSHi_detrended);hold on
%     
%     subplot(3,1,3)
%     title('SSHi minus Geoid; detrended; 3 km moving median ')
%     plot(data.LAT,movmedian(data.SSHi_detrended,movmedian_win_size),'.');hold on

    ax(i) = subplot(2,numel(L2_files),i)
    if i==1
        title('SSHi minus Geoid; detrended; 3 km moving median ')
    end
    
    along_track_distance = cumsum(deg2km(distance(data.LAT(1:end-1),data.LON(1:end-1),data.LAT(2:end),data.LON(2:end))));
    %plot(data.LAT,data.SSHi_detrended-data.r_detrended,'-');hold on % gives you the part of the corrections that is not explained by a linear slope
    plot(along_track_distance,data.SSHi_detrended(1:end-1),'k-');hold on % detrended SSHi
    plot(along_track_distance,data.SSHi_detrended(1:end-1),'k.');hold on % detrended SSHi
    %plot(data.LAT,movmedian(data.SSHi_detrended,60,'omitnan'),'-');hold on % detrended SSHi
    title(i)
    xlabel('km')
    
    ax2(i) = subplot(2,numel(L2_files),i+numel(L2_files))
    [acf,lags] = autocorr(data.SSHi_detrended,'NumLags',100);
    azimuth_res = 1e3*median(deg2km(distance(data.LAT(1:end-1),data.LON(1:end-1),data.LAT(2:end),data.LON(2:end))));
    plot(-600:600,sinc((-600:600)./300).^2,'r--','LineWidth',2);hold on
    plot( azimuth_res*[-flip(lags(2:end));lags(1:end)] , [ flip(acf(2:end));acf(1:end) ],'k-');hold on
    grid on
    if i==1
        xlabel('lag in m')
    end
    
end
linkaxes(ax,'xy')
ax(1).YLim = [-0.15,0.15];

linkaxes(ax2,'xy')
ax2(1).XLim = [-600,600];

%% check along track ACF for threshold retracker result:
i = 2;
clear CS1b
%load(fullfile(L2_files(i).folder,L2_files(i).name))
CS1b = load(fullfile(PathDATA,'ROFI','L1b',L2_files(i).name));
data = CS1b.CS1b.SAR.data_pseudoDD(:,50000:100000);

%%
azimuth_res = 1e3*median(deg2km(distance(CS1b.CS1b.GEO.LAT(1:end-1),CS1b.CS1b.GEO.LON(1:end-1),CS1b.CS1b.GEO.LAT(2:end),CS1b.CS1b.GEO.LON(2:end))));
r = threshold_retracker(data,87,0.23);

figure;
[acf,lags] = autocorr(detrend(r),'NumLags',1200);
plot(-1200:1200,sinc((-1200:1200)./300).^2,'r--','LineWidth',2);hold on
plot( azimuth_res*[-flip(lags(2:end));lags(1:end)] , [ flip(acf(2:end));acf(1:end) ],'k-');hold on

%% plotting the waveforms of a dedicated file index
i = 2;
clear CS1b
%load(fullfile(L2_files(i).folder,L2_files(i).name))
CS1b = load(fullfile(PathDATA,'ROFI','L1b',L2_files(i).name));
%%

figure;
imagesc(CS1b.CS1b.GEO.LAT,1:256,log10(CS1b.CS1b.SAR.data_pseudoDD))
colorbar()
colormap('pink')
caxis([10,15])

%%
figure;
imagesc(CS1b.CS1b.GEO.LAT,1:256,log10(movmean(CS1b.CS1b.SAR.data,50,2)))
colorbar()
colormap('pink')
caxis([11,15])

%% as an overview to detrended SSHi, plot also UF SAR waveforms
figure;
for i=1:numel(L2_files)
    L2_files(i).name
    %load(fullfile(L2_files(i).folder,L2_files(i).name))
end
%% testing corrections values:
x = data.LAT
dummy = CS2.COR.mod_dry_tropo_cor_meas_altitude_01(x) + ...
CS2.COR.mod_wet_tropo_cor_meas_altitude_01(x) + ...
CS2.COR.load_tide_sol2_01(x) + ...
CS2.COR.iono_cor_gim_01_ku(x) + ...
CS2.COR.solid_earth_tide_01(x) + ...
0.468*CS2.COR.pole_tide_01(x)

%% read and write the sensing times for Lennart
L2_files   = dir(fullfile(L2_dir,'*_SR_1_SRA_A*'));

for i =1:numel(L2_files)
    load(fullfile(L2_files(i).folder,L2_files(i).name))
    fprintf([CS2.StartTime,'\n'])
end

%% check per track the bias against eumetsat data north of 52 deg latitude

mission = 'S3A'%S3B

orbit = '_370_'
%orbit = '_279_'
lat_threshold = 52.5;

L2_files   = dir(fullfile(L2_dir,[mission,'*_SR_1_SRA_A*',orbit,'*']));

bias = zeros(numel(L2_files),1);
std = zeros(numel(L2_files),1);

figure;
for i=1:numel(L2_files)
    load(fullfile(L2_files(i).folder,L2_files(i).name))
    
    mask_lat = (CS2.f_020Hz_matched.UF.LAT > lat_threshold)';
    
    mask_outlier = (abs(diff(CS2.f_020Hz_matched.UF.SSHi)) < 0.2)&(abs(diff(CS2.f_020Hz_matched.EUM.SSHi)) < 0.2); % deletes one neighbouring value as well, doesn matter for now. 20cm noise threshold
    mask_outlier(end+1) = 0;
    
    mask = mask_lat & mask_outlier;
    
    plot(CS2.f_020Hz_matched.UF.SSHi(mask));hold on;plot(CS2.f_020Hz_matched.EUM.SSHi(mask),'o')
    
    bias(i) = nanmean(CS2.f_020Hz_matched.UF.SSHi(mask) - CS2.f_020Hz_matched.EUM.SSHi(mask));
    std(i) = nanstd(CS2.f_020Hz_matched.UF.SSHi(mask) - CS2.f_020Hz_matched.EUM.SSHi(mask));
    
    clear CS2
end

figure;
errorbar(bias,std)
% looking reasonable, track based biases of less than 1 cm

%% comparison to tide gauge data from Yosra



