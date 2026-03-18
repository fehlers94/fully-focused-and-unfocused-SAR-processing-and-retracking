clear all
close all
clc

mission = 'S3B'
%fname ='/home/fehelers/PhD_Delft/Hoek_van_Holland/S3B_SR_1_SRA_A__20191228T201909_20191228T210938_20200210T135810_3029_033_370______MR1_R_NT_004.SEN3/measurement_l1a.nc'; 
%fname ='/home/fehelers/PhD_Delft/Hoek_van_Holland/S3B_SR_1_SRA_A__20191008T201906_20191008T210934_20200209T171316_3028_030_370______MR1_R_NT_004.SEN3/measurement_l1a.nc';
fname = '/home/fehelers/PhD_Delft/Hoek_van_Holland/S3B_SR_1_SRA_A__20190815T201909_20190815T210937_20191220T220744_3028_028_370______MR1_R_NT_004.SEN3/measurement_l1a.nc';
DOM = [(51.99-0.07) (52.14+0.07); 0 360]
%DOM = [(75-0.07) (75+0.07); 0 360]

%data_dir = [getenv('HOME')];
%fname ='/home/fehelers/TUDTUM/ffsar_data/s3a/l1a/crete/S3A_SR_1_SRA_A__20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement_l1a.nc'; 


%%
DDAcf = DDA_ConfigFile(mission,'SAR');
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
FFSAR_processing_settings.combine_n_look_locations = 50;
%%

ffproc = FFSAR_Processor(fname, DOM, FFSAR_processing_settings);
ffproc.setup_proc();
ffproc.proc();

CS1b = ffproc.CS1b;

%%
save('/home/fehelers/PhD_Delft/for Marcel/ice_waves.mat','CS1b','-v7.3')

%%
figure;
dummy_plot = log10(movmean(CS1b.SAR.data,50,2));
imagesc(dummy_plot)
colorbar()
caxis([6,9.5])
colormap(pink)

figure; plot(dummy_plot(47,:))
figure; plot(10.^dummy_plot(:,9670))