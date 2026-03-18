%%
clc
clear all
close all
%%
%set = 'swh4.5m'
%set = 'swh1.5m'
%set = 'test'
%set = 'S6test_FF'
%dataset = 'S6test_UF'
dataset = 'S6test_FF'
%mission = 'S3B'
mission = 'S6A'
LoadCommonSettings
DDAcf   = DDA_ConfigFile(mission,'SAR');
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;

%% define (L1A and) L2 directories
L1a_file = [getenv('HOME') '/TUDTUM/ffsar_data/s6a/l1a-l1b-l2/crete/S6A_P4_1A_HR______20210901T212456_20210901T222027_20210902T134010_3331_030_018_009_EUM__OPE_ST_F03.nc'];
L1b_file = [getenv('HOME') '/TUDTUM/ffsar_data/s6a/l1a-l1b-l2/crete/S6A_P4_1B_HR______20210901T212456_20210901T222025_20210902T134009_3329_030_018_009_EUM__OPE_ST_F03.nc'];
L2_file = [getenv('HOME') '/TUDTUM/ffsar_data/s6a/l1a-l1b-l2/crete/S6A_P4_2__HR_STD__ST_030_018_20210901T212456_20210901T222025_F03.nc'];

[CS1b_eum] = S6_L1bs_read(L1b_file,[],true);
%S6.dop_beam_stack = complex(S3.i_echoes_ku_l1bs_echo_sar_ku, S3.q_echoes_ku_l1bs_echo_sar_ku);

% %% illustrate the data
% % figure; imagesc(log10(CS1b_eum.SAR.data));
% 
% S6ind = 7002;%14000
% S6num_bursts = CS1b_eum.SAR.N_averaged_echoes(S6ind)
% %figure; imagesc(CS1b_eum.GEO.LAT(S6ind-200:S6ind+200),1:256,log10(CS1b_eum.SAR.data(:,S6ind-200:S6ind+200)));
% 
% S6Lat = CS1b_eum.GEO.LAT(S6ind);
% stack_mask = CS1b_eum.SAR.stack_mask_start_stop(:,S6ind);
% binary_stack_mask_eum = ones(size(stack_mask,1),DDAcf.Np*DDAcf.os_ZP);
% 
% for i = 1:numel(stack_mask)
%     binary_stack_mask_eum(i,1:stack_mask(i)) = 0;
% end
% 
% %imagesc(~binary_stack_mask_eum)

%% read L2 data and illustrate data
gr = 'data_20/ku/' 
% 'sig0_ocean'
% 'swh_ocean'
% 'longitude'
% 'latitude'

L2.LAT_20                    = ncread(L2_file,[gr 'latitude']);
L2.LON_20                    = ncread(L2_file,[gr 'longitude']);
L2.SWH_20                    = ncread(L2_file,[gr 'swh_ocean']);
L2.SIGMA_20                    = ncread(L2_file,[gr 'sig0_ocean']);

S6_snippet = 7000:7500;
%S6_snippet = 7000:7250;

figure; 
ax1 = subplot(3,1,1)
imagesc(log10(CS1b_eum.SAR.data(:,S6_snippet)));
ax2 = subplot(3,1,2)
plot(L2.SWH_20(S6_snippet));hold on
plot(movmean(L2.SWH_20(S6_snippet),20));hold on
ax3 = subplot(3,1,3)
plot(L2.SIGMA_20(S6_snippet))
plot(movmean(L2.SIGMA_20(S6_snippet),20));hold on

linkaxes([ax1,ax2,ax3],'x')
% L2.ALT_20                    = ncread(L2_file,'alt_20_ku',startLoc,count);
% L2.range_20                  = ncread(L2_file,'range_ocean_20_ku',startLoc,count);
% L2.mean_ssh_20               = ncread(L2_file,'mean_sea_surf_sol1_20_ku',startLoc,count);
% L2.tracker_range_20          = ncread(L2_file,'tracker_range_20_ku',startLoc,count);

%% Define DOM for processing 
DOM = [min(L2.LAT_20(S6_snippet)) max(L2.LAT_20(S6_snippet)); 0, 360] % open ocean

% %%
% plot(CS1b_eum.SAR.data(:,S6ind))
% S6num_bursts = size(CS1b_eum.SAR.BeamAngle(:,S6ind),1) - sum(isnan(CS1b_eum.SAR.BeamAngle(:,S6ind)))
% % now, what integration time does this correpsond to?
% 
% % obviously, range zero padding is = 2

%% make h_l interpolator based on the local tracker range:
[~,ind] = unique(CS1b_eum.GEO.LAT);
h_interpolator = griddedInterpolant(CS1b_eum.GEO.LAT(ind),CS1b_eum.GEO.H(ind) - CS1b_eum.MEA.ref_range(ind),'nearest','nearest');

%% make some pseudo DD processing
if strcmp(dataset,'S6test_UF')
    %FFSAR_processing_settings.combine_n_look_locations = 1;
    FFSAR_processing_settings.combine_n_look_locations = 1;
    FFSAR_processing_settings.output_pseudo_delay_doppler_processing = true;
    FFSAR_processing_settings.num_coherent_bursts = 1;
    %FFSAR_processing_settings.integration_time = DDAcf.BRI*180.01; % corresponding to 180 bursts
    FFSAR_processing_settings.along_track_oversampling = 1/(0.25*585); % corresponding to ~80 Hz
    %FFSAR_processing_settings.along_track_oversampling = 1; % corresponding to ~80 Hz
    FFSAR_processing_settings.simulate_range_walk = true; % corresponding to ~80 Hz
    %FFSAR_processing_settings.integration_time = 2;
    %FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant'; % corresponding to ~80 Hz
    %DOM(3) = DOM(1)+0.025;
    FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant';
elseif strcmp(dataset,'S6test_FF')
    tic
    %FFSAR_processing_settings.combine_n_look_locations = 1;
    FFSAR_processing_settings.combine_n_look_locations = 50;
    FFSAR_processing_settings.output_pseudo_delay_doppler_processing = true;
    FFSAR_processing_settings.num_coherent_bursts = 1;
    %FFSAR_processing_settings.integration_time = DDAcf.BRI*180.01; % corresponding to 180 bursts
    %FFSAR_processing_settings.along_track_oversampling = 1/(0.25*585); % corresponding to ~80 Hz
    FFSAR_processing_settings.along_track_oversampling = 1; % corresponding to ~80 Hz
    FFSAR_processing_settings.simulate_range_walk = false; % corresponding to ~80 Hz
    FFSAR_processing_settings.integration_time = 2;
    %FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant'; % corresponding to ~80 Hz
    %DOM(3) = DOM(1)+0.025;
    FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant';
end
%%
tic
ffproc = FFSAR_Processor(L1a_file, DOM, FFSAR_processing_settings, h_interpolator);
%ffproc = FFSAR_Processor(L1a_file, DOM, FFSAR_processing_settings);
ffproc.setup_proc();
ffproc.proc();
CS1b = ffproc.CS1b;
toc

save(['/home/fehelers/PhD_Delft/pseudoDDprocessing/' dataset '.mat'],'ffproc','-v7.3' )
%%
load(['/home/fehelers/PhD_Delft/pseudoDDprocessing/' dataset '.mat'],'ffproc')
CS1b = ffproc.CS1b;
clear ffproc
%%
load(['/home/fehelers/PhD_Delft/pseudoDDprocessing/' 'S6test_UF' '.mat'],'ffproc')
CS1b_UF = ffproc.CS1b;
clear ffproc
%%
% make default plot backgroudn white
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);
%%
% choose color palette
cFF = [207, 0, 110]/255;
cDD = [0, 185, 227]/255;
c1b = [252, 157, 3]/255;
cT = [214, 214, 214]/255;
%% do some retracking here after averaging on about 20 Hz
% average on 80 Hz
N_avg = round((1/20)/median(diff(CS1b.GEO.Elapsed_Time)));
CS1b_avg = FF_SAR_Average_Waveforms(CS1b,N_avg,mission);
%%

CS1b_avg.SAR.data_FF = CS1b_avg.SAR.data;
CS1b_avg.MEA.ref_range = CS1b_avg.MEA.tracker_range;
CS1b_avg.SAR.scale_factor_ku = CS1b_avg.SAR.echo_scale_power;

%% 
% FFSAR with SAMOSA zero Doppler beam
CS1b_avg.SAR.data = CS1b_avg.SAR.data_FF
%[L2_FF,~] = SAR_L1b_to_L2(mission,CS1b_avg,DOM,'SAMOSA2',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'})
[L2_FF,~] = SAR_L1b_to_L2(mission,CS1b_avg,DOM,'SAMOSA2FF',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'});
%%
% UFSAR
CS1b_avg.SAR.data = CS1b_avg.SAR.data_pseudoDD
[L2_UF,~] = SAR_L1b_to_L2(mission,CS1b_avg,DOM,'SAMOSA2',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'})


%% #####################   plotting retracking results   #######################
% without moving mean
%L2.SSH_20 = L2.ALT_20-L2.range_20;
%maskL2 = (L2.LAT_20>DOM(1,1))&(L2.LAT_20<DOM(1,2));

figure('units','inch','position',[0,0,6.5,5/3*2]);
lw=1;

mask = ones(size(L2_FF.HEI));
mask(793) = 0; 
mask = logical(mask);


ax1= subplot(2,1,1)
dataFF = L2_FF.HEI(mask);
dataUF = L2_UF.HEI(mask);

% FF data
plot(L2_FF.LAT(mask),dataFF,'Color',cFF,'LineStyle','-','LineWidth',lw);hold on
% UF data
plot(L2_UF.LAT(mask),dataUF,'Color',cDD,'LineStyle','-','LineWidth',lw);hold on
grid on

R = corrcoef(dataFF,dataUF)
std_error = std(dataFF-dataUF)/sqrt(numel(dataFF)/1)
text(56.301,41.25,{['bias =(',num2str(1000*mean(dataFF-dataUF),2),'\pm',num2str(1.96*1000*std_error,2),') mm',', R = ',num2str(R(1,2),4)]})
%title({['bias =',num2str(1000*mean(dataFF-dataUF),2),' mm',', R = ',num2str(R(1,2),2)]})
ylabel('uncorrected SSH (m)')
h = legend('FF-SAR; SAMOSA zero-Doppler beam fit (2 s)','UF-SAR; Full SAMOSA fit (2 s)')
set(h,'Location','northwest')


%plot(L2.LAT_20(maskL2),L2.SSH_20(maskL2))

ax2= subplot(2,1,2)
dataFF = L2_FF.SWH(mask);
dataUF = L2_UF.SWH(mask);
% FF data
plot(L2_FF.LAT(mask),dataFF,'Color',cFF,'LineStyle','-','LineWidth',lw);hold on
% UF data
plot(L2_UF.LAT(mask),dataUF,'Color',cDD,'LineStyle','-','LineWidth',lw);hold on
%plot(L2.LAT_20,L2.SWH_20);hold on
grid on
R = corrcoef(dataFF,dataUF)
std_error = std(dataFF-dataUF)/sqrt(numel(dataFF)/1)
text(56.301,2.85,{['bias =(',num2str(100*mean(dataFF-dataUF),2),'\pm',num2str(1.96*100*std_error,2),') cm',', R = ',num2str(R(1,2),2)]})
%title({['bias =',num2str(100*mean(dataFF-dataUF),2),' cm',', R = ',num2str(R(1,2),2)]})
ylabel('SWH (m)')

%plot(L2.LAT_20(maskL2),L2.SWH_20(maskL2))

% ax3= subplot(3,1,3)
% dataFF = L2_FF.sigma0(mask);
% dataUF = L2_UF.sigma0(mask);
% sig0_bias=92.5;
% % FF data
% plot(L2_FF.LAT(mask),dataFF-sig0_bias,'Color',cFF,'LineStyle','-','LineWidth',lw);hold on
% % UF data
% plot(L2_UF.LAT(mask),dataUF-sig0_bias,'Color',cDD,'LineStyle','-','LineWidth',lw);hold on
% plot(L2.LAT_20,L2.SIGMA_20);hold on
% grid on
% R = corrcoef(dataFF,dataUF)
% %text(54.371,10.2,{['bias =',num2str(mean(dataFF-dataUF),2),' dB',', R = ',num2str(R(1,2),2)]})
% title({['bias =',num2str(mean(dataFF-dataUF),2),' dB',', R = ',num2str(R(1,2),2)]})
% ylabel('sigma0 (dB)')
xlabel('latitude (°N)')

linkaxes([ax1 ax2 ax3],'x')
ax1.XLim = [56.2,56.5]
ax2.YLim = [1.5,3.0]
exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/paper_S6_retracking_results.pdf'],'Resolution',300)
