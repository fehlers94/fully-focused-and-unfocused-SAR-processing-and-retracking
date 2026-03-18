clear all
close all
clc

data_dir = [getenv('HOME') '/TUDTUM/ffsar_data/'];
fname ='s3a/l1a/crete/S3A_SR_1_SRA_A__20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement_l1a.nc'; 
DOM = [35.15 35.4; 23 24]

DDAcf = DDA_ConfigFile(mission,'SAR');
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
FFSAR_default_processing_settings.transponder_test_mode = true;



%%

kernel = [0.3 0 0 0 0.5 0 0 0 1 0 0 0 0.5 0 0 0 0.3];
multilook = [1 1 1 1 1 1];

res = conv(multilook, kernel);

plot(res,'b')