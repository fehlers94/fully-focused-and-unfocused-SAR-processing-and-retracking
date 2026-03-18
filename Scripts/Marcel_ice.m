clear all
close all
clc

mission = 'S3B'
%fname ='/home/fehelers/PhD_Delft/for Marcel/S3B_SR_1_SRA_A__20210312T105625_20210312T114655_20210407T021254_3029_050_094______MAR_O_NT_004.SEN3/measurement_l1a.nc'; 
%fname ='/home/fehelers/PhD_Delft/for Marcel/S3B_SR_1_SRA_A__20210312T114655_20210312T123724_20210407T025921_3029_050_094______MAR_O_NT_004.SEN3/measurement_l1a.nc';
%fname ='/home/fehelers/PhD_Delft/for Marcel/S3B_SR_1_SRA_A__20220122T225052_20220122T234122_20220217T140219_3029_061_372______MAR_O_NT_004.SEN3/measurement_l1a.nc';
fname = '/home/fehelers/PhD_Delft/for Marcel/S3B_SR_1_SRA_A__20220122T220023_20220122T225052_20220217T132952_3029_061_371______MAR_O_NT_004.SEN3/measurement_l1a.nc'

DOM = [65.6 67.6; -180 180]
%DOM = [65.6 65.61; -180 180]
parts = 4;

%data_dir = [getenv('HOME')];
%fname ='/home/fehelers/TUDTUM/ffsar_data/s3a/l1a/crete/S3A_SR_1_SRA_A__20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement_l1a.nc'; 


%%
DDAcf = DDA_ConfigFile(mission,'SAR');
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
FFSAR_processing_settings.combine_n_look_locations = 50;
FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant';
FFSAR_processing_settings.n_cores_parfor = 2;
%FFSAR_processing_settings.along_track_oversampling = 0.01;
%FFSAR_processing_settings.combine_n_look_locations = 1;
%FFSAR_processing_settings.combine_n_look_locations = 50;
%FFSAR_processing_settings.output_pseudo_delay_doppler_processing = true;
%FFSAR_processing_settings.num_coherent_bursts = 1;
%FFSAR_processing_settings.integration_time = DDAcf.BRI*180.01; % corresponding to 180 bursts
%FFSAR_processing_settings.along_track_oversampling = 1/(0.25*585); % corresponding to ~80 Hz
%FFSAR_processing_settings.along_track_oversampling = 1; % corresponding to ~80 Hz
%FFSAR_processing_settings.simulate_range_walk = false; % corresponding to ~80 Hz
%FFSAR_processing_settings.integration_time = 2;
%FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant'; % corresponding to ~80 Hz
%DOM(3) = DOM(1)+0.025;
%%
lat_max = DOM(3);
lat_min = DOM(1);
d_lat = (DOM(3)-DOM(1))/4;

%%
for part = 3:parts
    DOM_proc = DOM;
    DOM_proc(1) = lat_min + (part-1)*d_lat
    DOM_proc(3) = lat_min + (part)*d_lat

    ffproc = FFSAR_Processor(fname, DOM_proc, FFSAR_processing_settings);
    ffproc.setup_proc();
    ffproc.proc();

    CS1b = ffproc.CS1b;


    save(['/home/fehelers/ownCloud/ice_waves_4-' num2str(part) '.mat'],'CS1b','-v7.3')
end
%%
figure
imagesc(CS1b.GEO.LAT,1:256,log10(movmean(CS1b.SAR.data,50,2)))
%figure
%imagesc(log10(movmean(CS1b.SAR.data,1,2)))
colorbar()
colormap(pink)

% plot(sum(CS1b.SAR.data,2))
