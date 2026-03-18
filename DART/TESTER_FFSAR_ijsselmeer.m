clc;clear all;
close all;

% test FFSAR_Processor
LoadCommonSettings
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
data_dir = fullfile([getenv('HOME') '/TUDTUM/ffsar_data/']);

% preconditions
FName = [data_dir 's3b/l1a/ijsselmeer/S3B_SR_1_SRA_A__20181223T201137_20181223T210201_20190801T130555_3024_020_099______MR1_R_NT_004.SEN3/measurement_l1a.nc'];
DOM = [52.6 52.9; 0 180]; %crete transponder, lat 35.3379, lon 23.7795

%%
ffproc = FFSAR_Processor(FName,DOM,FFSAR_processing_settings);
ffproc.setup_proc()
ffproc.proc()

%%

%imagesc(1:256,ffproc.CS1b.GEO.LAT,ffproc.CS1b.SAR.data(:,1:end)'.^2)

CS1ba = FF_SAR_Average_Waveforms(ffproc.CS1b,50,'S3A');

%%

image = CS1ba.SAR.data(60:end,:)'.^2;

figure;
imagesc(1:256,CS1ba.GEO.LAT,log(image))

alt = mean(CS1ba.GEO.H);
[X,Y] = meshgrid(1:size(image,1),1:size(image,2));
Y = sqrt( (alt+Y.*0.47/2).^2 - alt.^2 );

figure;
s = pcolor(X,Y,log(image'))
grid off
s.EdgeColor = 'None';
colormap('pink');
%s.LineWidth = 6;

%% read S3B SAR L1B file
L1BFName = '/home/fehelers/PhD Delft/Documents/Meetings/PSD Section meeting/S3B_SR_1_SRA____20181223T201137_20181223T210201_20190801T130555_3024_020_099______MR1_R_NT_004.SEN3/measurement.nc';
[S3,CS1bSAR] = S3_L1b_read(L1BFName);

IDXsar = ((DOM(1)+0.06) < S3.lat_l1b_echo_sar_ku)&(S3.lat_l1b_echo_sar_ku < (DOM(3)-0.06));
imageSAR = S3.i2q2_meas_ku_l1b_echo_sar_ku(40:end-5,IDXsar);

alt = mean(CS1ba.GEO.H);
[X,Y] = meshgrid(1:size(imageSAR',1),1:size(imageSAR',2));
Y = sqrt( (alt+Y.*0.47).^2 - alt.^2 );

figure;
s = pcolor(X,Y,log(imageSAR))
grid off
s.EdgeColor = 'None';
colormap('pink');


%% block averaging
wav = ffproc.CS1b.SAR.data'.^2;
% Define the mean function to be used with blockproc.
myMeanFunction = @(block_struct) mean(block_struct.data);
% Get the means by each 2-by-2 block.
wav_avg = blockproc(wav, [50, 1], myMeanFunction);
%%
imagesc(1:256,ffproc.CS1b.GEO.LAT,wav_avg)
