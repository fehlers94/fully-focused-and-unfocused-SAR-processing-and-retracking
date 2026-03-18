clc;clear all;
close all;

% test FFSAR_Processor
LoadCommonSettings
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
data_dir = fullfile([getenv('HOME') '/TUDTUM/ffsar_data/']);

% preconditions
FName = [data_dir 's6a/S6A_GPP__P4__FR__1A_20210116T201350_20210116T201413_0001.NC'];
DOM = [35.15 35.4; 23 24]; %crete transponder, lat 35.3379, lon 23.7795
% DOM_crete = [35.15 35.4; 23 24]; %crete transponder, lat 35.3379, lon 23.7795

%% test basic class setup
ffproc = FFSAR_Processor(FName,DOM,FFSAR_processing_settings);
assert(strcmp(ffproc.mission,'S6A'));

%% setup processing
ffproc = FFSAR_Processor(FName,DOM,FFSAR_processing_settings);
ffproc.setup_proc()
assert(~isempty(ffproc.wav));
assert(~isempty(ffproc.CS1a));

%% custom FFSAR_Processor setting: number of look loc inds
FFSAR_processing_settings.combine_n_look_locations = 50;

% each combine_n_look_locations look loc is processed
ffproc = FFSAR_Processor(FName,DOM,FFSAR_processing_settings);
ffproc.setup_proc()
assert(length(ffproc.look_loc_inds),length(ffproc.h_l)/FFSAR_processing_settings.combine_n_look_locations)

% only transponder loc is processed
FFSAR_processing_settings.transponder_test_mode = true;
ffproc = FFSAR_Processor(FName,DOM,FFSAR_processing_settings);
ffproc.setup_proc()
assert(length(ffproc.look_loc_inds),1)

%% flat-phase-transponder scenario test
runtests('TESTER_FFSAR_Processor_TransponderTest')
