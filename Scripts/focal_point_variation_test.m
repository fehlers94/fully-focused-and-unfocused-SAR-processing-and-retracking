clc
clear all
close all

%%
LoadCommonSettings
DDAcf   = DDA_ConfigFile('S3A','SAR');
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
%DOM = [52.197 52.51; -30 30] % actual DOM
DOM = [52.197 52.31; -30 30] % actual DOM

L1a_file = '/home/fehelers/PhD Delft/Projects/FFSAR4ROFI/Data/ROFI/L1a/S3A_SR_1_SRA_A__20181129T101709_20181129T110739_20181225T013537_3029_038_279______MAR_O_NT_003.SEN3/measurement_l1a.nc'

%%

tic
FFSAR_processing_settings.combine_n_look_locations = 50;
FFSAR_processing_settings.output_pseudo_delay_doppler_processing = true;
FFSAR_processing_settings.num_coherent_bursts = 1;
FFSAR_processing_settings.integration_time = DDAcf.BRI*180.01; % corresponding to 180 bursts
FFSAR_processing_settings.along_track_oversampling = 1; % corresponding to ~80 Hz
FFSAR_processing_settings.simulate_range_walk = false; % corresponding to ~80 Hz
FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant';

ffproc = FFSAR_Processor(L1a_file, DOM, FFSAR_processing_settings);
ffproc.setup_proc();
ffproc.proc();
toc
%CS1b = ffproc.CS1b;
%%
%CS1b_shifted_rev = ffproc.CS1b;
%%
%CS1b_shifted = ffproc.CS1b;
%%
%CS1b_normal = ffproc.CS1b;
%% plotting UF-SAR
cmin = 8;
cmax = 13;

figure;
ax1 = subplot(1,3,1)
imagesc(log10(CS1b_shifted.SAR.data_pseudoDD));
colorbar()
caxis([cmin cmax])

ax2 = subplot(1,3,2)
imagesc(log10(CS1b_shifted_rev.SAR.data_pseudoDD));
colorbar()
caxis([cmin cmax])

ax3 = subplot(1,3,3)
imagesc(log10(abs(CS1b_shifted_rev.SAR.data_pseudoDD-CS1b_shifted.SAR.data_pseudoDD)));
colorbar()
caxis([cmin-2 cmax-2])
linkaxes([ax1 ax2 ax3],'xy')

%% plotting FF-SAR
cmin = -9;
cmax = -5;
n_mean = 50;

figure;

data1 = movmean(CS1b_shifted.SAR.data,n_mean,2);
%data1 = CS1b_shifted.SAR.data;
ax1 = subplot(1,3,1)
imagesc(log10(data1./sum(data1(:))));
title('left side')
colorbar()
colormap('pink')
caxis([cmin cmax])


data2 = movmean(CS1b_shifted_rev.SAR.data,n_mean,2);
%data2 = CS1b_shifted_rev.SAR.data;
ax2 = subplot(1,3,2)
imagesc(log10(data2./sum(data2(:))));
title('right side')
colorbar()
colormap('pink')
caxis([cmin cmax])

data3 = abs(CS1b_shifted_rev.SAR.data + CS1b_shifted.SAR.data);
%data3 = abs(data1-data2);

%data3 = 10.^(log10(data3)./log10(data1+data2));
% data3= abs(data3./(data2+data1));
data3 = movmean(data3,n_mean,2);
data3 = log10(data3./sum(data3(:)));
%data = angle(CS1b_normal.SAR.dataIQ)-angle(CS1b_shifted.SAR.dataIQ);



ax3 = subplot(1,3,3)
imagesc(data3);
title('difference, normalized')
colorbar()
colormap('pink')
caxis([cmin cmax])
linkaxes([ax1 ax2 ax3],'xy')

%% plot correlation of radargram entries
hist3([log10(data1(:)), log10(data2(:))],'CDataMode','auto','Nbins',[1 1]*100)
view(2)
%hist3([log10(data1(:)), log10(data2(:))],'Nbins',[1 1]*200)


%plot(log10(data1(:)),log10(data2(:)),'b.')
