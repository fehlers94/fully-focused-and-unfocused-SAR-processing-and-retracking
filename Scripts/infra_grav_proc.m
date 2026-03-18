clear all
LoadCommonSettings

PathDATA = '/home/fehelers/PhD Delft/Projects/infra_grav/';
DOM = [53.5 58; -30 30] % for 14-01 pass
%DOM = [51.5 58; -30 30] % for 15-01 pass

% read the L1a data directory
L1A_data   = dir([PathDATA,'L1a/*_SR_1_SRA_A*']);

i = 1;
    
L1b_filepath = fullfile(PathDATA,'L1b',[L1A_data(i).name '.mat'])
%%
% Load parameters and settings
DDAcf   = DDA_ConfigFile('S3A','SAR');

% determine FFSAR processing settings
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
FFSAR_processing_settings.combine_n_look_locations = 50;
FFSAR_processing_settings.output_pseudo_delay_doppler_processing = true;
FFSAR_processing_settings.num_coherent_bursts = 1;
FFSAR_default_processing_settings.along_track_oversampling = 1;
FFSAR_processing_settings.integration_time = DDAcf.BRI*180.01%DDAcf.BRI*180.01; % corresponding to 180 bursts as in L1b products
%FFSAR_processing_settings.integration_time = 2;
%FFSAR_processing_settings.along_track_oversampling = 1/(0.25*585); % corresponding to ~80 Hz
FFSAR_processing_settings.along_track_oversampling = 1; % corresponding to FFSAR IRF speckle noise resolution
FFSAR_processing_settings.simulate_range_walk = false;
FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant';
%FFSAR_processing_settings.n_cores_parfor = num_workers;

L1a_filepath = fullfile(PathDATA,'L1a',L1A_data(i).name,'measurement_l1a.nc');
ffproc = FFSAR_Processor(L1a_filepath, DOM, FFSAR_processing_settings);
%ffproc = FFSAR_Processor(L1a_file, DOM, FFSAR_processing_settings);
ffproc.setup_proc();
%%

tic
ffproc.proc();
toc

%%
CS1b = ffproc.CS1b;

save(L1b_filepath,'CS1b','-v7.3')

%%
clear ffproc

%%