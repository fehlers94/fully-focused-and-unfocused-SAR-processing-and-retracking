clear all
close all
clc
%% S3B case

mission = 'S3B'
data_dir = fullfile([getenv('HOME') '/TUDTUM/ffsar_data/']);
%fname ='s3a/l1a/crete/S3A_SR_1_SRA_A__20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement_l1a.nc'; 
%fname_l1b = '/home/fehelers/TUDTUM/ffsar_data/s3a/l1a-l1b-l2_sets/l1b/S3A_SR_1_SRA____20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement.nc';
fname = 's3b/l1a-l1b-l2_sets/l1a/S3B_SR_1_SRA_A__20190620T083403_20190620T092432_20191107T073220_3029_026_335______MR1_R_NT_004.SEN3/measurement_l1a.nc'
fname_l1b = '/home/fehelers/TUDTUM/ffsar_data/s3b/l1a-l1b-l2_sets/l1b/S3B_SR_1_SRA____20190620T083403_20190620T092432_20191107T073220_3029_026_335______MR1_R_NT_004.SEN3/measurement.nc'
%DOM = [(35.34-0.005-0.07) (35.34+0.0051+0.07); 0 360]
DOM = [35.33 35.37; 0 360]

%%
mission = 'S6A'
data_dir = fullfile([getenv('HOME') '/TUDTUM/ffsar_data/']);
fname = ['s6a/l1a-l1b-l2/crete/S6A_P4_1A_HR______20210901T212456_20210901T222027_20210902T134010_3331_030_018_009_EUM__OPE_ST_F03.nc'];
fname_l1b = ['s6a/l1a-l1b-l2/crete/S6A_P4_1B_HR______20210901T212456_20210901T222025_20210902T134009_3329_030_018_009_EUM__OPE_ST_F03.nc'];
%L2_file = [getenv('HOME') '/TUDTUM/ffsar_data/s6a/l1a-l1b-l2/crete/S6A_P4_2__HR_STD__ST_030_018_20210901T212456_20210901T222025_F03.nc'];
DOM = [35.33 35.35; 0 360]

%%
DDAcf = DDA_ConfigFile(mission,'SAR');
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
FFSAR_processing_settings.combine_n_look_locations = 50;
FFSAR_processing_settings.output_pseudo_delay_doppler_processing = true;
FFSAR_processing_settings.num_coherent_bursts = 1;
FFSAR_processing_settings.simulate_range_walk = false;
FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant';
%%

ffproc = FFSAR_Processor([data_dir fname], DOM, FFSAR_processing_settings);
ffproc.setup_proc();
ffproc.proc();

CS1b = ffproc.CS1b;

%%
% compare if total power is approximately equal:
sum(sum(CS1b.SAR.data_pseudoDD))
sum(sum(CS1b.SAR.data))

%% visualize the different point target responses and their relation
%% is the diffusion in range reasonable!? What is delay doppler response over transponder!? look at S3 overpass L1b file!!!
figure; plot(sum(CS1b.SAR.data_pseudoDD,2)); hold on;
plot(sum(CS1b.SAR.data,2)); hold on;

xlabel('range bin')

sum(sum(CS1b.SAR.data_pseudoDD,2))
sum(sum(CS1b.SAR.data,2))
DDkernel = sum(CS1b.SAR.data_pseudoDD,2);
FFkernel = sum(CS1b.SAR.data,2);

%% with IRF correction
figure; plot(sum(CS1b.SAR.data_pseudoDD,2)); hold on;
FFSAR_processing_settings.integration_time = 3.6;
kernel = IRF_corr(CS1b, FFSAR_processing_settings);
FF_corrected = conv2(CS1b.SAR.data,kernel,'same');

plot(sum(FF_corrected,2)); hold on;

xlabel('range bin')

sum(sum(CS1b.SAR.data_pseudoDD,2))
sum(sum(FF_corrected,2))

%% 
[val, ind] = max(FFkernel);
FFkernel = FFkernel(ind-10:ind+10)/sum(FFkernel(ind-10:ind+10));
DDkernel = DDkernel(ind-10:ind+10)/sum(DDkernel(ind-10:ind+10));
plot(FFkernel); hold on;plot(DDkernel);hold on

FTInput = fft(FFkernel);
FtOutput = fft(DDkernel);
kernel = fftshift(ifft(FtOutput./FTInput));

plot(kernel); hold on
plot(conv2(FFkernel,kernel,'same'),'ro'); hold on
%%
figure;
imagesc(CS1b.GEO.LAT,1:256,log10(CS1b.SAR.data_pseudoDD))
cax = colorbar()
ylabel(cax,'log10(power)')
caxis([11,17])
colormap('pink')
%%
DD_scale = (DDAcf.Nb*DDAcf.PRI)*2*DDAcf.fc*median(CS1b.GEO.V.V)/(DDAcf.c*median(CS1b.MEA.tracker_range)); %now, multiplied with d (y_a in Egido), we get along track the since argument
m_per_look = 1e3*deg2km(distance(CS1b.GEO.LAT(1),CS1b.GEO.LON(1),CS1b.GEO.LAT(end),CS1b.GEO.LON(end)))/numel(CS1b.GEO.LON);
N_looks = 1:size(CS1b.SAR.data,2);

figure;
ax1 = subplot(3,1,1);
%imagesc(CS1b.GEO.LAT,1:256,log10(CS1b.SAR.data_pseudoDD))
imagesc(log10(CS1b.SAR.data_pseudoDD))
cax = colorbar()
ylabel(cax,'log10(power)')
caxis([11,17])
colormap('pink')
ax2 = subplot(3,1,2);
%imagesc(CS1b.GEO.LAT,1:256,log10(CS1b.SAR.data))
imagesc(log10(CS1b.SAR.data))
cax = colorbar()
ylabel(cax,'log10(power)')
caxis([11,17])
ax3 = subplot(3,1,3);

% sum along range bin dim to overcome stretching
plot((sum(CS1b.SAR.data_pseudoDD) / max(sum(CS1b.SAR.data_pseudoDD)))); hold on
% plot theoretical sinc over it
[val,ind] = max(sum(CS1b.SAR.data_pseudoDD));
plot(sinc(DD_scale*m_per_look*(N_looks-ind)).^2,'--')

plot(   (sum(CS1b.SAR.data) / max(sum(CS1b.SAR.data)))); hold on

% try to account for width of the peaks as well to correct the misfit with the envelope
plot(   (movsum(sum(CS1b.SAR.data),50) / max(movsum(sum(CS1b.SAR.data),50))))
colorbar()
legend('sum(pseudo DD)','sinc model', 'sum(FF SAR)', 'movsum(sum(FF SAR))')

linkaxes([ax1,ax2,ax3],'x')
linkaxes([ax1,ax2],'y')
%TO-DO: put correct axes

%% look into corresponding L1b file for transponder signal
% correct for the weird scaling / unit mismatch!

[S3,CS1b_eum] = S3_L1b_read(fname_l1b);
idx = (CS1b_eum.GEO.LAT>35.33)&(CS1b_eum.GEO.LAT<35.35);

% find colocations with FF SAR
[minValue,closestIndex] = min(abs(CS1b_eum.GEO.LAT(idx) - CS1b.GEO.LAT),[],2);

figure;
ax1 = subplot(2,1,1);
dummy1 = CS1b.SAR.data_pseudoDD(:,closestIndex);
dummy1_plot = 10*log10(dummy1./max(dummy1(:)));
imagesc(dummy1_plot);
title('pseudoDD')
cax = colorbar()
ylabel(cax,'power (dB)')
caxis([-60,0])
colormap('pink')

ax2 = subplot(2,1,2);
%imagesc(log10(CS1b_eum.SAR.data(:,idx)));
a = CS1b_eum.SAR.agc_ku(idx); % dB
b = CS1b_eum.SAR.scale_factor_ku(idx); % dB

% THERE IS REALLY SOMETHING ODD WITH THE AGC CORRECTION, maybe the relict
% at the border is from using less data?
%dummy = CS1b_eum.SAR.data(:,idx).*10.^(a'/10).*10.^(b'/10)*180;
dummy = CS1b_eum.SAR.data(:,idx).*220.*10.^(-b'/20); % under approximation of ~ 180 doppler beams when taking average instead of sum
dummy_plot = 10*log10(dummy./max(dummy(:)));
imagesc(dummy_plot);

title('Eumetsat L1b')
cax = colorbar()
ylabel(cax,'power (dB)')
caxis([-60,0])
colormap('pink')

linkaxes([ax1,ax2],'x')
linkaxes([ax1,ax2],'y')

% plot waveforms over the domain
figure;
n = size(dummy1_plot,2);
for j = 1:n
    subplot(1,n,j)
    plot(dummy1_plot(:,j),1:128);hold on
    plot(dummy_plot(:,j),1:128);hold on
    ylim([0.5,128.5])
    set(gca, 'YDir','reverse')
end

figure;
%subplot(2,1,1);
plot(log10(sum((CS1b.SAR.data_pseudoDD(:,closestIndex))))); hold on
plot(log10(sum(dummy)));
ylabel('log10(power)')
xlabel('# waveform')
title(['pseudoDD, mismatch factor = ' num2str(median(sum((CS1b.SAR.data_pseudoDD(:,closestIndex)))./sum(dummy)))])



%% example: simulate a multilook of some target and its response
% multilooking is basically adding up point target responses along track
% can be approximated by convolution along track

n_looks = 1200
multilook_PTR = conv2(CS1b.SAR.data,ones(1,n_looks),'same');

figure;
ax1 = subplot(3,1,1);
imagesc(   (CS1b.SAR.data_pseudoDD))
colormap('pink')
ax2 = subplot(3,1,2);
imagesc(   (multilook_PTR))
ax3 = subplot(3,1,3);
plot(   sum(multilook_PTR,1) / max(sum(multilook_PTR,1))); hold on
plot(   sum(CS1b.SAR.data_pseudoDD,1) / max(sum(CS1b.SAR.data_pseudoDD,1))); hold on

linkaxes([ax1,ax2,ax3],'x')
linkaxes([ax1,ax2],'y')

%% actually process a water target south of crete

FFSAR_processing_settings.combine_n_look_locations = 1;
%DOM = [(34.0-0.07) (34.05+0.07); 0 360] for a longer segment and D/D autocorrelation
DOM = [(34.0-0.07) (34.0+0.07); 0 360]

ffproc = FFSAR_Processor([data_dir fname], DOM, FFSAR_processing_settings);
ffproc.setup_proc();
ffproc.proc();

CS1b_wat = ffproc.CS1b;

%% plot the averaged D/D SAR and FF SAR waveforms

figure;
ax1=subplot(1,3,1);
imagesc(CS1b_wat.SAR.data_pseudoDD)
colorbar()
colormap(pink)
title('pseudo DD SAR')
ax2=subplot(1,3,2);
imagesc(movmean(CS1b_wat.SAR.data,2500/4,2))
colorbar()
colormap(pink)
title('FF SAR multilooked')
ax3 = subplot(1,3,3);

%plot((mean(CS1b_wat.SAR.data_pseudoDD(:,1:500:end),2)),1:128','LineWidth',2);hold on
plot((mean(CS1b_wat.SAR.data,2)),1:128','LineWidth',2);hold on
plot((mean(CS1b_wat.SAR.data_pseudoDD,2)),1:128','LineWidth',2); hold on
set(gca, 'YDir','reverse')
legend('FFSAR','pseudo DD SAR')


linkaxes([ax1,ax2],'x')
linkaxes([ax1,ax2],'y')

%% zoomed version
figure;
ax1=subplot(1,2,1);
dummy_plot = movmean(CS1b_wat.SAR.data_pseudoDD,2500/4,2);
imagesc(dummy_plot(:,1:6000))
colorbar()
colormap(pink)
title('pseudo DD SAR')
ax2=subplot(1,2,2);
dummy_plot = movmean(CS1b_wat.SAR.data,2500/4,2);
imagesc(dummy_plot(:,1:6000))
colorbar()
colormap(pink)
title('FF SAR multilooked')



%% calculate estimated number of looks by mean power squared divided by variance according to Quartly et al. 2001 and Egido et al. 2020
%% this is only approximate because the waveforms are not perfectly aligned, but it should be good enough for an impression, if the data along a range bin is divided by the moving mean
%% also compare the so-obtained correlation length with Egido et al. ~160 m ~ 320 looks
%plot(CS1b_wat.SAR.data(40,:)); hold on
%plot(movmean(CS1b_wat.SAR.data(40,:),1000))

% visualise the DD data
rbin = 60;
multilook = 2500
%figure;
%plot(CS1b_wat.SAR.data_pseudoDD(rbin,:)); hold on
%plot(movmean(CS1b_wat.SAR.data(rbin,:),multilook)); hold on
%plot(movmean(CS1b_wat.SAR.data_pseudoDD(rbin,:),multilook))


DDnorm = 1;%movmean(CS1b_wat.SAR.data_pseudoDD(rbin,:),2000);
DDdata = CS1b_wat.SAR.data_pseudoDD(rbin,:)./DDnorm;
FFdata = CS1b_wat.SAR.data(rbin,1:1000);

% after normalization
%figure;
%plot(DDdata)

% calculate autocorrelation D/D SAR:
figure;
title('DD SAR autocorrelation')
plot(0.56/10*(1:3601),autocorr(DDdata,'NumLags',3600))
xlabel('distance (m)')
ylabel('Autocorrelation function')
grid

% calculate autocorrelation FFSAR:
figure;
title('DD SAR autocorrelation')
plot(0:10,autocorr(FFdata(1:1:end),'NumLags',10))
xlabel('look')
ylabel('Autocorrelation function')
grid

%%

% now for DD SAR and FF SAR take corresponding an average and std corresponding to
% almost 20 Hz data after normalization

figure;
FF = movmean(CS1b_wat.SAR.data(rbin,:),multilook);
DD = movmean(CS1b_wat.SAR.data_pseudoDD(rbin,:),multilook);

title('moving mean with same multilooking applied')
plot(FF); hold on
plot(DD)
legend('FF SAR','DD SAR')

ENL_FF = mean(FF).^2/var(FF)
ENL_DD = mean(DD).^2/var(DD)


%% compare the psuedo DD waveforms to L1b SAR waveforms and single FFSAR multilooked waveforms

[S3,CS1b_eum] = S3_L1b_read(fname_l1b);
idx = (CS1b_eum.GEO.LAT>33.99)&(CS1b_eum.GEO.LAT<34.01);

% find colocations with FF SAR
[minValue,closestIndex] = min(abs(CS1b_eum.GEO.LAT(idx) - CS1b_wat.GEO.LAT),[],2);


figure;
ax1 = subplot(2,1,1);
dummy1 = CS1b_wat.SAR.data_pseudoDD(:,closestIndex);
dummy1_plot = 10*log10(dummy1./max(dummy1(:)));
imagesc(dummy1_plot);
title('pseudoDD')
cax = colorbar()
ylabel(cax,'power (dB)')
caxis([-30,0])
colormap('pink')

ax2 = subplot(2,1,2);
%imagesc(log10(CS1b_eum.SAR.data(:,idx)));
a = CS1b_eum.SAR.agc_ku(idx); % dB
b = CS1b_eum.SAR.scale_factor_ku(idx); % dB

%dummy = CS1b_eum.SAR.data(:,idx).*10.^(a'/10).*10.^(b'/10)*180;
dummy = CS1b_eum.SAR.data(:,idx).*10.^(a'/10).*180; % under approximation of ~ 180 doppler beams when taking average instead of sum
dummy_plot = 10*log10(dummy./max(dummy(:)));
imagesc(dummy_plot);

title('Eumetsat L1b')
cax = colorbar()
ylabel(cax,'power (dB)')
caxis([-30,0])
colormap('pink')

linkaxes([ax1,ax2],'x')
linkaxes([ax1,ax2],'y')

% plot waveforms over the domain
figure;
n = size(dummy1_plot,2);
for j = 1:n
    subplot(1,n,j)
    plot(10.^(dummy1_plot(:,j)./10),1:128,'LineWidth',2);hold on
    plot(10.^(dummy_plot(:,j)./10),1:128,'LineWidth',2);hold on
    ylim([0.5,128.5])
    set(gca, 'YDir','reverse')
end

%%
figure;
n = size(dummy1_plot,2);
for j = 1:1
    subplot(1,3,j)
    plot(10.^(dummy1_plot(:,j)./10),1:128,'LineWidth',2);hold on
    plot(10.^(dummy_plot(:,j)./10),1:128,'LineWidth',2);hold on
    ylim([0.5,128.5])
    set(gca, 'YDir','reverse');    
    legend('pseudo D/D','L1b D/D');
end
subplot(1,3,3)
plot((10.^(dummy1_plot(:,1)./10)-10.^(dummy_plot(:,1)./10))./(10.^(dummy_plot(:,1)./10)),1:128)
set(gca, 'YDir','reverse');
title('relative difference')

