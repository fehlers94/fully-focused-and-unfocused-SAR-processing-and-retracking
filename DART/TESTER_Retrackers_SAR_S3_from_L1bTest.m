clc;clear all;
close all;

LoadCommonSettings

data_dir = fullfile([getenv('HOME') '/TUDTUM/ffsar_data/']);

% S3
FName_l1b = [data_dir '/s3a/l1a-l1b-l2_sets/l1b/S3A_SR_1_SRA____20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement.nc'];
FName_l2 = [data_dir '/s3a/l1a-l1b-l2_sets/l2/S3A_SR_2_WAT____20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004_concat.nc'];
DOM = [33.45 33.60; 0, 360];

% % CS2
% FName_l1b = [data_dir 'cs/l1a-l1b-l2_sets/l1b/CS_LTA__SIR_SAR_1B_20150503T160800_20150503T161329_D001.nc'];
% FName_l2 = [data_dir 'cs/l1a-l1b-l2_sets/l2/CS_LTA__SIR_SAR_2__20150503T160800_20150503T161329_D001_realigned.nc'];
% DOM = [-15.571 -15.491; 0, 360];

mission = mission_from_fname(FName_l1b);
DDAcf = DDA_ConfigFile(mission,'SAR');

do_plot = true;

%% Test: Validate retracking result against EUM L2

SAMOSAimpl = {...
    "'DPM_v2_5_2'"...
%     "'S+'"...
    };

for samimpl=SAMOSAimpl
samimpl = samimpl{1};
    
% retrack data
[L2_retracked,CS1b_eum] = SAR_L1b_to_L2(mission,FName_l1b,DOM,'SAMOSA2',...
    {'SAMOSAimpl', samimpl;'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'});

% read original L2 for comparison
startLoc = L2_retracked.IDXpoi(1);
count = length(L2_retracked.IDXpoi);
L2.SWH = ncread(FName_l2,DDAcf.var_mapping_l2.swh,startLoc,count);
L2.LON = ncread(FName_l2,DDAcf.var_mapping_l2.lon,startLoc,count);
L2.LAT = ncread(FName_l2,DDAcf.var_mapping_l2.lat,startLoc,count);

%reverse, that cog_cor has been added to range already to reproduce approximately 'raw' retracking result
reverse_cog = 0.0;
if contains(mission, 'S3')
    reverse_cog = 0.55590;
end
L2.HEI = ncread(FName_l2,DDAcf.var_mapping_l2.alt,startLoc,count) - (ncread(FName_l2,DDAcf.var_mapping_l2.range,startLoc,count)) + reverse_cog;

% evaluate stuff
rmse_swh = sqrt(mean((L2_retracked.SWH - L2.SWH).^2));
rmse_ssh = sqrt(mean((L2_retracked.HEI - L2.HEI).^2));

% assert(rmse_swh<0.1)
% assert(rmse_ssh<0.05)

% plot stuff
if do_plot
figure;
subplot(2,1,1)
plot(L2_retracked.SWH);hold on
plot(L2.SWH);hold on
ylabel('SWH [m]')
xlabel('i')
legend(samimpl + ':  L1b->L2','EUM L2','Interpreter', 'none')
title(samimpl+ ", RMSE=" + num2str(rmse_swh),'Interpreter', 'none')
grid on

subplot(2,1,2)
plot(L2_retracked.HEI);hold on
plot(L2.HEI);hold on
ylabel('uncorrected SSH [m]')
xlabel('i')
legend(samimpl + ':  L1b->L2','EUM L2','Interpreter', 'none')
title(samimpl + ", RMSE=" + num2str(rmse_ssh),'Interpreter', 'none')
grid on
end

end