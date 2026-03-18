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

%% check variability of the kernel input
v_max = max(CS1b.GEO.V.V)
v_min = min(CS1b.GEO.V.V)
h_max = max(CS1b.GEO.H)
h_min = min(CS1b.GEO.H)
hdot_max = max(CS1b.GEO.H_rate)
hdot_min = min(CS1b.GEO.H_rate)
%% optional: convolute waveform with kernel
% make sure its centered symmetrically on odd support, so not to introduce a shift
kernel;
figure;plot(log10(abs(kernel)));

%%
%dummy = conv2(CS1b.SAR.data,kernel,'same');
% figure;
% ax1 = subplot(1,2,2)
% imagesc(log(dummy))
% ax1 = subplot(1,2,1)
% imagesc(log(CS1b.SAR.data))

CS1b.SAR.data_conv = conv2(CS1b.SAR.data,kernel,'same');

%% ##################### average waveforms (not pseudo_DD) ############
N_avg = 100;

CS1b_avg = FF_SAR_Average_Waveforms(CS1b,N_avg,mission);
CS1b_avg.SAR.data_FF = CS1b_avg.SAR.data;
CS1b_avg.MEA.ref_range = CS1b_avg.MEA.tracker_range;
CS1b_avg.SAR.scale_factor_ku = CS1b_avg.SAR.echo_scale_power;

CS1b_UF.MEA.ref_range = CS1b_UF.MEA.tracker_range;
CS1b_UF.SAR.scale_factor_ku = CS1b_UF.SAR.echo_scale_power;
CS1b_UF.SAR.data = CS1b_UF.SAR.data_pseudoDD;

%%
figure;
subplot(1,3,1)
imagesc(CS1b_avg.SAR.data)
subplot(1,3,2)
imagesc(CS1b_avg.SAR.data_conv)
subplot(1,3,3)
imagesc(CS1b_avg.SAR.data_pseudoDD)

%% #####################   do some retracking   #######################

% FFSAR without convolution
CS1b_avg.SAR.data = CS1b_avg.SAR.data_FF
[L2_FF,~] = SAR_L1b_to_L2(mission,CS1b_avg,DOM,'SAMOSA2FF',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'});

%% 
% FFSAR with full SAMOSA
CS1b_avg.SAR.data = CS1b_avg.SAR.data_FF
[L2_FF_wrong,~] = SAR_L1b_to_L2(mission,CS1b_avg,DOM,'SAMOSA2',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'});

%%
% FFSAR with convolution
CS1b_avg.SAR.data = CS1b_avg.SAR.data_conv
[L2_FF_conv,~] = SAR_L1b_to_L2(mission,CS1b_avg,DOM,'SAMOSA2',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'});

%%
% FFSAR with convolution
CS1b_avg.SAR.data = CS1b_avg.SAR.data_conv
[L2_FF_conv_wrong,~] = SAR_L1b_to_L2(mission,CS1b_avg,DOM,'SAMOSA2FF',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'});

%%
% UFSAR
CS1b_avg.SAR.data = CS1b_avg.SAR.data_pseudoDD
[L2_UF,~] = SAR_L1b_to_L2(mission,CS1b_avg,DOM,'SAMOSA2',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'});

%%
% UF-SAR with range walk and 3.4 s integration time
% UFSAR
[L2_UF_rw,~] = SAR_L1b_to_L2(mission,CS1b_UF,DOM,'SAMOSA2',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'});
%% #####################   plotting retracking results   #######################
% without moving mean
figure;
% L2 data
plot(L2.LAT_20(S6_snippet),L2.SWH_20(S6_snippet),'k');hold on
% FF data
plot(L2_FF.LAT,L2_FF.SWH,'Color',cFF);hold on
% FF data
plot(L2_FF_conv.LAT,L2_FF_conv.SWH,'Color',cFF);hold on
% UF data
plot(L2_UF.LAT,L2_UF.SWH,'Color',cDD);hold on
% UF data with rw and 3.41 s integration time
plot(L2_UF_rw.LAT,L2_UF_rw.SWH,'Color','cyan');hold on



%% with moving mean
numel(L2.SWH_20(S6_snippet))
numel(L2_FF.SWH)
numel(L2_UF.SWH)
numel(L2_UF_rw.SWH)
%%
n_mov_eum = 20
n_mov_FF = round(n_mov_eum*numel(L2_FF.SWH)/numel(L2.SWH_20(S6_snippet)))
n_mov_UF_rw = round(n_mov_eum*numel(L2_UF_rw.SWH)/numel(L2.SWH_20(S6_snippet)))

%%
figure;
lw=2;
% L2 data
plot(L2.LAT_20(S6_snippet),movmean(L2.SWH_20(S6_snippet),n_mov_eum),'k','LineWidth',lw);hold on
% FF data
plot(L2_FF.LAT,movmean(L2_FF.SWH,n_mov_FF),'Color',cFF,'LineStyle','--','LineWidth',lw);hold on
% FF data
plot(L2_FF_conv.LAT,movmean(L2_FF_conv.SWH,n_mov_FF),'Color',cFF,'LineStyle','-','LineWidth',lw);hold on
% UF data
plot(L2_UF.LAT,movmean(L2_UF.SWH,n_mov_FF),'Color',cDD,'LineStyle','-','LineWidth',lw);hold on
% UF data with rw and 3.41 s integration time
plot(L2_UF_rw.LAT,movmean(L2_UF_rw.SWH,n_mov_UF_rw),'Color',cDD,'LineStyle','--','LineWidth',lw);hold on
grid on
ylabel('significant waveheight SWH (m)')
xlabel('latitude (°N)')

legend('UF-SAR L2 estimate (3.4 s)','FF-SAR with SAMOSA zero-Doppler beam (2 s)','convoluted FF-SAR with Full SAMOSA (2 s)','UF-SAR with Full SAMOSA (2 s)','UF-SAR inclusive simulated range-walk with Full SAMOSA (3.4 s)')

%% plotted at ~1 Hz only 
figure('units','inch','position',[0,0,6.5,5]);
lw=2;

% L2 data
LAT = L2.LAT_20(S6_snippet)
LAT = LAT(1:n_mov_eum:end)
SWH = movmean(L2.SWH_20(S6_snippet),n_mov_eum)
SWH = SWH(1:n_mov_eum:end)
plot(LAT,SWH,'Color','k','LineStyle','-','LineWidth',lw);hold on

% FF data
LAT = L2_FF.LAT
LAT = LAT(1:n_mov_FF:end)
SWH = movmean(L2_FF.SWH,n_mov_FF)
SWH = SWH(1:n_mov_FF:end)
plot(LAT,SWH,'Color',cFF,'LineStyle','--','LineWidth',lw);hold on

% FF data wrong (full SAMOSA)
LAT = L2_FF.LAT
LAT = LAT(1:n_mov_FF:end)
SWH = movmean(L2_FF_wrong.SWH,n_mov_FF)
SWH = SWH(1:n_mov_FF:end)
plot(LAT,SWH,'Color',cFF,'LineStyle','-.','LineWidth',lw);hold on

% FF data conv (full SAMOSA)
LAT = L2_FF_conv.LAT
LAT = LAT(1:n_mov_FF:end)
SWH = movmean(L2_FF_conv.SWH,n_mov_FF)
SWH = SWH(1:n_mov_FF:end)
plot(LAT,SWH,'Color',cFF,'LineStyle','-','LineWidth',lw);hold on

% % FF data conv wrong (zero-Doppler beam SAMOSA)
% LAT = L2_FF_conv.LAT
% LAT = LAT(1:n_mov_FF:end)
% SWH = movmean(L2_FF_conv_wrong.SWH,n_mov_FF)
% SWH = SWH(1:n_mov_FF:end)
% plot(LAT,SWH,'Color',cFF,'LineStyle',':','LineWidth',lw);hold on


% UF data
LAT = L2_UF.LAT
LAT = LAT(1:n_mov_FF:end)
SWH = movmean(L2_UF.SWH,n_mov_FF)
SWH = SWH(1:n_mov_FF:end)
plot(LAT,SWH,'Color',cDD,'LineStyle','-','LineWidth',lw);hold on

% UF data with rw and 3.41 s integration time
LAT = L2_UF_rw.LAT
LAT = LAT(1:n_mov_UF_rw:end)
SWH = movmean(L2_UF_rw.SWH,n_mov_UF_rw)
SWH = SWH(1:n_mov_UF_rw:end)

plot(LAT,SWH,'Color',cDD,'LineStyle','--','LineWidth',lw);hold on
grid on
ylabel('significant waveheight SWH (m, 1 Hz average)')
xlabel('latitude (°N)')

ylim([1.,2.6])
xlim([55.65,56.68])

h = legend('official UF-SAR Eumetsat estimate (3.4 s)','FF-SAR; SAMOSA zero-Doppler beam fit (2 s)','FF-SAR; Full SAMOSA fit (2 s)','convoluted FF-SAR; Full SAMOSA fit (2 s)','UF-SAR; Full SAMOSA fit (2 s)','UF-SAR including range walk error; Full SAMOSA fit (3.4 s)')
set(h,'Location','best')
exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/paper_SWH_comparison_new_3_after_discussion_flo.pdf'],'Resolution',300)

%% ############## compare pseudo D/D to FF-SAR averaged waveform #########
% for ENL calculation we need to first average the waveforms onto 20 Hz

time_tot = CS1b.GEO.Elapsed_Time(end)-CS1b.GEO.Elapsed_Time(1);
num_samples = time_tot*20;
averaging_size = round(size(CS1b.SAR.data,2)/num_samples)

CS1b.SAR.data_pseudoDD = movmean(CS1b.SAR.data_pseudoDD,averaging_size,2)
CS1b.SAR.data = movmean(CS1b.SAR.data,averaging_size,2)

%%
figure;
ax1 = subplot(2,1,1)
imagesc(CS1b.SAR.data)
ax2 = subplot(2,1,2)
imagesc(CS1b.SAR.data_pseudoDD)

linkaxes([ax1,ax2],'xy')
%linkaxes([ax1,ax2],'y')

%%
figure;
plot(CS1b.SAR.data(120,:));hold on
plot(CS1b.SAR.data_pseudoDD(120,:))

%%

comp.avg_pseudo_DD = nanmean(CS1b.SAR.data_pseudoDD,2);
comp.avg_FF = nanmean(CS1b.SAR.data,2);

comp.ENL_pseudo_DD = (nanvar(CS1b.SAR.data_pseudoDD,[],2)./nanmean(CS1b.SAR.data_pseudoDD,2).^2).^(-1);
comp.ENL_FF = (nanvar(CS1b.SAR.data,[],2)./nanmean(CS1b.SAR.data,2).^2).^(-1);


%% ############## compare pseudo D/D to FF-SAR averaged waveform #########

figure('units','inch','position',[0,0,4,4]);
ax1 = subplot(2,1,1)
lw = 1;

yyaxis left
plot(comp.avg_FF./max(comp.avg_FF),'Color',cFF,'LineWidth',lw,'LineStyle','-'); hold on
plot(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD),'Color',cDD,'LineWidth',lw,'LineStyle','--'); hold on
legend('FF-SAR','UF-SAR')
%xlabel('range bin')
ylabel('power')
ylim([0,1.05])

yyaxis right
p1 = plot(comp.ENL_FF,'Color',cFF,'LineWidth',lw,'LineStyle','-'); hold on
p2 = plot(comp.ENL_pseudo_DD,'Color',cDD,'LineWidth',lw,'LineStyle','-'); hold on
p1.Color(4) = 0.2;
p2.Color(4) = 0.2;
ylabel('20 Hz ENL')
legend('FF-SAR','UF-SAR')
ylim([0,1000])

ax2 = subplot(4,1,3)
plot(100*(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD)-comp.avg_FF./max(comp.avg_FF))./(comp.avg_FF./max(comp.avg_FF)),'Color',[0.5,0.5,0.5]); hold on
%plot(100*(comp.avg_pseudo_DD_lin_scaled-comp.avg_L1b)./comp.avg_L1b); hold on
%title('(pseudo D/D - L1b) / L1b')
ylabel({'relative'; 'difference (%)'})
%xlabel('range bin')
%ylabel('%')
grid on
ylim([-7,10])

ax3 = subplot(4,1,4)
plot(100*(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD) - comp.avg_FF./max(comp.avg_FF)),'Color',[0.5,0.5,0.5]); hold on
%plot(100*(comp.avg_pseudo_DD_lin_scaled - comp.avg_L1b)/max(comp.avg_L1b)); hold on
%title('difference relative to maximum: pseudo D/D - L1b')
ylabel({'absolute'; 'difference (%)'})
xlabel('range bin')

grid on
ylim([-1,1])


linkaxes([ax1,ax2,ax3],'x')
xlim([0,numel(comp.avg_pseudo_DD)])

set(ax1,'xticklabel',[])
set(ax2,'xticklabel',[])

ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = 'k';

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
%print -dpdf -painters epsFig

%exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/' dataset '.pdf'],'Resolution',300)
exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/convoluted_' dataset '.pdf'],'Resolution',300)
%saveas(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/' dataset '.pdf'])

%% now compare to L1b file
%% filter for the overlapping waveforms
[minval,idx] = min(abs(CS1b.GEO.LAT - CS1b_eum.GEO.LAT),[],2);
%size(idx)

% get indices where lat projected distance is smaller than 50 m:
L1bidx = find(111000*minval < 50);
FFidx = idx(111000*minval < 50);

%% plot radargram
figure;
subplot(2,1,1)
imagesc(CS1b_eum.GEO.LAT, 1:512, CS1b_eum.SAR.data(:,L1bidx))
title('L1b SAR')
subplot(2,1,2)
imagesc(CS1b.GEO.LAT(FFidx), 1:512, CS1b.SAR.data_pseudoDD(:,FFidx))
title('pseudo D/D SAR')
colormap('pink')

%% compare integration time
figure;
plot(CS1b_eum.SAR.N_averaged_echoes(L1bidx));hold on
plot(CS1b.SAR.num_bursts(FFidx))

%% plot single look waveforms
N = numel(L1bidx)
figure;
for i = 1:min(N,5)
    subplot(min(N,5),1,i)
    
    % scaled with maximum
    plot(CS1b_eum.SAR.data(:,L1bidx(i)) / max(CS1b_eum.SAR.data(:,L1bidx(i)))  );hold on
    plot(CS1b.SAR.data_pseudoDD(:,FFidx(i))    /   max(CS1b.SAR.data_pseudoDD(:,FFidx(i)))   );hold on
    
    legend('Level-1b waveform', 'SAR pseudo D/D waveform')
end


%% based on the matching, calculate average waveforms for L1b and pseudo D/D
comp.avg_pseudo_DD = nanmean(CS1b.SAR.data_pseudoDD(:,FFidx),2);
offset_corr = 0.0015%0.0012;
comp.avg_pseudo_DD = max(comp.avg_pseudo_DD)*(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD)-offset_corr)./(1-offset_corr);
comp.avg_L1b = nanmean(CS1b_eum.SAR.data(:,L1bidx),2);

comp.ENL_pseudo_DD = (nanvar(CS1b.SAR.data_pseudoDD(:,FFidx),[],2)./nanmean(CS1b.SAR.data_pseudoDD(:,FFidx),2).^2).^(-1)
comp.ENL_L1b = (nanvar(CS1b_eum.SAR.data(:,L1bidx),[],2)./nanmean(CS1b_eum.SAR.data(:,L1bidx),2).^2).^(-1)

% for matlab plotting of several y axes see https://uk.mathworks.com/help/matlab/creating_plots/plotting-with-two-y-axes.html
%% plot comparison with max scale

figure('units','inch','position',[0,0,4,4]);
ax1 = subplot(2,1,1)
lw = 1;

yyaxis left
plot(comp.avg_L1b./max(comp.avg_L1b),'Color',c1b,'LineWidth',lw,'LineStyle','-'); hold on
plot(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD),'Color',cDD,'LineWidth',lw,'LineStyle','--'); hold on
legend('L1b D/D','UF-SAR')
%xlabel('range bin')
ylabel('power')
ylim([0,1.05])

yyaxis right
p1 = plot(comp.ENL_L1b,'Color',c1b,'LineWidth',lw,'LineStyle','-'); hold on
p2 = plot(comp.ENL_pseudo_DD,'Color',cDD,'LineWidth',lw,'LineStyle','-'); hold on
p1.Color(4) = 0.2;
p2.Color(4) = 0.2;
ylabel('20 Hz ENL')
legend('L1b D/D','UF-SAR')
ylim([0,1000])

ax2 = subplot(4,1,3)
plot(100*(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD)-comp.avg_L1b./max(comp.avg_L1b))./(comp.avg_L1b./max(comp.avg_L1b)),'Color',[0.5,0.5,0.5]); hold on
%plot(100*(comp.avg_pseudo_DD_lin_scaled-comp.avg_L1b)./comp.avg_L1b); hold on
%title('(pseudo D/D - L1b) / L1b')
ylabel({'relative'; 'difference (%)'})
%xlabel('range bin')
%ylabel('%')
grid on
ylim([-5,5])

ax3 = subplot(4,1,4)
plot(100*(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD) - comp.avg_L1b./max(comp.avg_L1b)),'Color',[0.5,0.5,0.5]); hold on
%plot(100*(comp.avg_pseudo_DD_lin_scaled - comp.avg_L1b)/max(comp.avg_L1b)); hold on
%title('difference relative to maximum: pseudo D/D - L1b')
ylabel({'absolute'; 'difference (%)'})
xlabel('range bin')

grid on
ylim([-1,1])


linkaxes([ax1,ax2,ax3],'x')
xlim([0,numel(comp.avg_pseudo_DD)])

set(ax1,'xticklabel',[])
set(ax2,'xticklabel',[])

ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = 'k';

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
%print -dpdf -painters epsFig

exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/' dataset '.pdf'],'Resolution',300)
%saveas(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/' dataset '.pdf'])

%% make some hacky retracking of L1b, pseudoDD and FFSAR
%% bare in mind that FFSAR waveforms are not perfectly range aligned!
plot(CS1b.GEO.H - CS1b.MEA.tracker_range)



% 
% %% plot CS1b waveforms for comparison:
% normL1b = CS1b_eum.SAR.data(:,S6ind)./max(CS1b_eum.SAR.data(:,S6ind));
% normDD = CS1b.SAR.data_pseudoDD(:,ffind)./max(CS1b.SAR.data_pseudoDD(:,ffind));
% figure;
% plot(normL1b);hold on
% plot(normDD);hold on
% legend('L1b waveform','pseudo D/D waveform')
% 
% % plot differences
% figure;
% plot(normL1b-normDD);hold on
% %plot(movmedian(normL1b-normDD,10));hold on
% legend('pseudo D/D minus L1b')
% 
% figure;
% plot(-CS1b.SAR.BeamAngle(:,ffind),'r.');hold on
% plot(CS1b_eum.SAR.BeamAngle(:,S6ind),'b-');hold on
% 
% %% plot the stack masks for comparison:
% figure;
% subplot(1,3,1)
% stack_mask = CS1b_eum.SAR.stack_mask_start_stop(:,S6ind);
% binary_stack_mask_eum = ones(size(stack_mask,1),DDAcf.Np*DDAcf.os_ZP);
% 
% for i = 1:numel(stack_mask)
%     binary_stack_mask_eum(i,1:stack_mask(i)) = 0;
% end
% 
% imagesc(~binary_stack_mask_eum)
% title('L1b range mask')
% 
% subplot(1,3,2)
% 
% stack_mask = CS1b.SAR.stack_mask_start_stop(:,ffind);
% binary_stack_mask = ones(size(stack_mask,1),DDAcf.Np*DDAcf.os_ZP);
% 
% for i = 1:numel(stack_mask)
%     binary_stack_mask(i,1:stack_mask(i)) = 0;
% end
% 
% imagesc(~binary_stack_mask)
% title('Our range mask')
% 
% subplot(1,3,3)
% imagesc(~binary_stack_mask_eum - ~binary_stack_mask)
% title('residual')
% colorbar()

