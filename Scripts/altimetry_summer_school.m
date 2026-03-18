clear all;

%% ##########################################
%               coastal case
% ##########################################

filename = 'S3A_SR_1_SRA_A__20211210T101721_20211210T110750_20220105T012443_3029_079_279______MAR_O_NT_004.SEN3.mat';

load(fullfile('/home/fehelers/PhD Delft/Documents/Courses/Teaching Assistant/Summer school SAR altimetry/coastal',filename));

%%
% make default plot backgroudn white
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);

%% plot track on map
figure;
geobasemap();
geoplot(CS1b.GEO.LAT,CS1b.GEO.LON,'LineWidth',3)

%% example for a dB waveform plot (UF-SAR)
% here the tail of the ocean waveform drops by about 50 % because only one
% side is land, see subplot 2

id1 = 70000;
id2 = 96000;
id3 = 200000;

figure;
subplot(2,1,1)
imagesc(CS1b.GEO.LAT,1:256,10*log10(CS1b.SAR.data_pseudoDD(:,1:end)));hold on
plot(repmat(CS1b.GEO.LAT(id1),256),1:256,'k-')
plot(repmat(CS1b.GEO.LAT(id2),256),1:256,'k--')
plot(repmat(CS1b.GEO.LAT(id3),256),1:256,'k-.')
s = colorbar()
ylabel(s,'power in dB')
caxis([90 130])
ylabel('range gate')
xlabel('latitude (deg)')
set(gca, 'YDir','reverse')
colormap('pink')
subplot(2,2,4)
plot(CS1b.SAR.data_pseudoDD(:,id1)./1e6,1:256);hold on
plot(CS1b.SAR.data_pseudoDD(:,id2)./1e6,1:256)
ylim([1,256])
set(gca, 'YDir','reverse')
xlabel('power')
subplot(2,2,3)

plot(CS1b.SAR.data_pseudoDD(:,id3)./1e6,1:256);hold on
set(gca, 'YDir','reverse')
ylim([1,256])
ylabel('range gate')
xlabel('power')

%% running standard deviation over ~80 Hz data as an example for how to detect signal contamination
data = log10(CS1b.SAR.data_pseudoDD(:,1:150:end));
std_data = movstd(data,[5 5],0,2);
%figure;
%imagesc(std_data)
%colorbar()

figure;
ax1 = subplot(3,1,1)
imagesc(data)
colormap('pink')
title("log10 power radargram")

ax2 = subplot(3,1,2)
imagesc(std_data)
title("moving std in log10 power")
colormap('pink')

ax3 = subplot(3,1,3)
masked_data = data;
masked_data(std_data>0.1)=NaN;
imagesc(masked_data)
colormap('pink')
title("masked log10 power radargram")
% think about measures how to extract a quality flag from such analysis
% e.g.

linkaxes([ax1 ax2 ax3],'xy')


%% 

clear all


%% ##########################################
% part on lead detection and radargram
%% ##########################################

% make default plot backgroudn white
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);

%% plot the olci geotiff from sentinel hub

[OLCI,R] = readgeoraster('/home/fehelers/PhD Delft/Documents/Courses/Teaching Assistant/Summer school SAR altimetry/leads/2017-04-29-00 00_2017-04-29-23 59_Sentinel-3_OLCI_RGB_(17,_6,_3)(1).tiff');
%[EUML2,COR] = S3_read_L2('/home/fehelers/PhD Delft/Documents/Courses/Teaching Assistant/Summer school SAR altimetry/leads/S3A_SR_2_WAT____20170429T203437_20170429T211949_20191212T230550_2711_017_114______MR1_R_NT_004.SEN3/enhanced_measurement.nc',[71,73.3;-141,-123]);
load('/home/fehelers/PhD Delft/Documents/Courses/Teaching Assistant/Summer school SAR altimetry/leads/S3A_SR_1_SRA_A__20170429T203415_20170429T212444_20190901T010333_3029_017_114______MR1_R_NT_004.SEN3.mat')
%%
CS1b.AVG = [];
CS1b.COR = [];
CS1b.SAR.data_pseudoDD_zeroDop = [];
%%
figure;
worldmap([71,73.3],[-134,-128])
geoshow(OLCI,R);
geoshow(CS1b.GEO.LAT(1:1000:end),CS1b.GEO.LON(1:1000:end),'LineWidth',3)


%% plot it UF-SAR at roughly 80 Hz
figure;
%imagesc(1:256,CS1b.GEO.LAT(:,1:150:end),10*log10(CS1b.SAR.data_pseudoDD(:,1:150:end))')
imagesc(10*log10(CS1b.SAR.data_pseudoDD(:,1:150:end)'))
set(gca, 'YDir','reverse')
colormap('pink')
s = colorbar()
ylabel(s,'power in dB')
caxis([90,150])
set(gca, 'YDir','normal')

%% plot an FF-SAR excerpt
figure;
%imagesc(1:256,CS1b.GEO.LAT(:,1:150:end),10*log10(CS1b.SAR.data_pseudoDD(:,1:150:end))')
imagesc(10*log10(CS1b.SAR.data(:,1:250*150)'))
set(gca, 'YDir','reverse')
colormap('pink')
s = colorbar()
ylabel(s,'power in dB')
caxis([90,150])
set(gca, 'YDir','normal')


%% plot it FF-SAR (zoom in to see leads with ghosts and estimate/visualise width of the lead)
%data = movmean(CS1b.SAR.data,50,2);
%data = data(:,1:end);
%data = CS1b.SAR.data(:,1:2000);
figure;
ax1 = subplot(1,2,1)
imagesc(1:256,CS1b.GEO.LAT,log10(CS1b.SAR.data)')
caxis([11,15])
set(gca, 'YDir','normal')
ax2 = subplot(1,2,2)
plot(sum(CS1b.SAR.data),CS1b.GEO.LAT)
colormap('pink')


linkaxes([ax1 ax2],'y')

%% ##########################################
% some plotting of files
%% ##########################################
% load transponder image
load('/home/fehelers/PhD_Delft/pseudoDDprocessing/paperIRF_S3B.mat')
%%
CS1b.SAR.data_pseudoDD = CS1b.SAR.data_pseudoDD./max(CS1b.SAR.data_pseudoDD(:));
CS1b.SAR.data = CS1b.SAR.data./max(CS1b.SAR.data(:));
%%
figure;
ax1 = subplot(2,1,1)
imagesc(10*log10(CS1b.SAR.data))
colormap('pink')
colorbar()
caxis([-30,10])

ax2 = subplot(2,1,2)
imagesc((CS1b.SAR.data))
colormap('pink')
colorbar()
%caxis([120,160])

linkaxes([ax1 ax2],'xy')
%% downsample to 20 Hz
dt = mean(diff(CS1b.GEO.Elapsed_Time))
n = 1;%round(1/20/dt);
Ntot = size(CS1b.SAR.data,2);

figure;
ax1 = subplot(2,1,1)
data = CS1b.SAR.data(:,1:n:end);
imagesc(1:numel(data)*n,1:256,10*log10(data))
colormap('pink')
s = colorbar()
ylabel(s,'power (dB)')
caxis([-30,5])
ylabel('range gate')


ax1.XLim = [1.5e6,4.5e6]
ax1.YLim = [20,90]