function [CONST, TR, FFSAR_default_processing_settings] = FFSAR_LoadCommonSettings()
% CONST
% Geophysical parameters
CONST.Re=6371008.8; % radius Earth [m]
CONST.c=299792458; % speed-of-light [m]
CONST.a=6378137.0; % semi-major axis [m]
CONST.b=6356752.31425; % semi-minor axis [m] 
%a=6378136.3; % TOPEX
%f=1/298.257; % TOPEX
%b=(1-f)*a; % TOPEX
CONST.e=sqrt(1-CONST.b^2/CONST.a^2); % eccentricity 

% TR (transponder locations)
% crete
TR.crete.lon=23+46/60+46.266892/3600; %23.7795
TR.crete.lat=35+20/60+16.547930/3600; %35.3379
TR.crete.h=1048.8258; %1046

% svalbard
TR.svalbard.lon=15.39376997;
TR.svalbard.lat=78.23052306;
TR.svalbard.h=492.772-9.88/2;

% FFSAR_default_processing_settings
% default settings for FFSAR_Processor
FFSAR_default_processing_settings = {};
FFSAR_default_processing_settings.combine_n_look_locations = 50;
FFSAR_default_processing_settings.transponder_test_mode = false;
FFSAR_default_processing_settings.transponder_latlonalt = [TR.crete.lat, TR.crete.lon, TR.crete.h];
FFSAR_default_processing_settings.output_pseudo_delay_doppler_processing = false;
FFSAR_default_processing_settings.num_coherent_bursts = 1; % experimental setting for pseudo delay doppler, determines the number of bursts that are coherently summed
FFSAR_default_processing_settings.integration_time = []; % if set: number of seconds
FFSAR_default_processing_settings.along_track_oversampling = 1; %1 corresponds to the locally calculated azimuth resolution, but can be chosen as any floating point value, due to time interpolation (undersampling possible)
FFSAR_default_processing_settings.simulate_range_walk = false; %should remain false, but for comparison with outdated ESA L1b data
FFSAR_default_processing_settings.focal_point_height_interpolation = 'pchip_smooth'; % either 'piecewise_constant' or 'pchip_smooth', the latter might lead to extreme waveheight biases, whenever the tracker is not well behaved. The former is more robust, as it references waveforms (piecewise) along isocontours determined by the erratic tracker settings. This setting is not used / inactive, whenever an interpolator is handed over to the processor directly.
FFSAR_default_processing_settings.do_correct_IRF = false; % correct based on our IRF correction model
FFSAR_default_processing_settings.n_cores_parfor = 0;%feature('numcores'); % 0 means no parallel execution, feature('numcores') means to use the number of physical cores (not logical, which is the default setting)  
FFSAR_default_processing_settings.split_aperture = false;
end
