% This script tries to rigorously compare the D/D L1b waveforms with FF-SAR
% and unfcused SAR output from put processor. We start from the L2 product
% SSH, which will be taken as reference range within our processing, and
% to which the L1b waveform Doppler beam (IQ) stacks will be shifted. This
% way, the effects of the reference range on the multilooking will be
% minimized and the waveforms become comparable due to identical
% references. Furthermore, we try to identify the bursts that are being
% used by our processor and remove any additional bursts from the L1b
% doppler beam stack as well, to produce the most possibly aligned product.
%%
clc
clear all
close all
%%
%dataset = 'swh4.5m'
dataset = 'swh1.5m'
mode = 'FF'
%set = 'test'
%set = 'S6test'
mission = 'S3B'
%mission = 'S6A'
LoadCommonSettings
DDAcf   = DDA_ConfigFile(mission,'SAR');
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
%lat_buffer = 0.07273; % necessary to account for the whole DOM on ground
%DOM = [52 60; -30 30] % actual DOM

% define excerpt for examples and debugging
%DOM_proc = [55.1 55.1; -180 180]
%DOM_proc = [58 60; -180 180]


%% define (L1A and) L2 directories
switch dataset
    case 'test'
        L1a_file = [getenv('HOME') '/TUDTUM/pseudo_DD_analysis/S3B_SR_1_SRA_A__20211102T201926_20211102T210955_20211128T112328_3029_058_370______MAR_O_NT_004.SEN3/measurement_l1a.nc'];
        L1b_file = [getenv('HOME') '/TUDTUM/pseudo_DD_analysis/S3B_SR_1_SRA_BS_20211102T201926_20211102T210955_20211128T112328_3029_058_370______MAR_O_NT_004.SEN3/measurement_l1bs.nc'];
        L2_file = [getenv('HOME') '/TUDTUM/pseudo_DD_analysis/S3B_SR_2_WAT____20211102T202512_20211102T210955_20211128T121218_2683_058_370______MAR_O_NT_004.SEN3/enhanced_measurement.nc'];
        DOM = [52 60; -30 30] % actual DOM
        %DOM_proc = DOM + [-lat_buffer lat_buffer; 0 0]
    case 'swh1.5m'
        L1a_file = [getenv('HOME') '/TUDTUM/pseudo_DD_analysis/S3A_SR_1_SRA_A__20211017T101723_20211017T110752_20211112T011946_3029_077_279______MAR_O_NT_004.SEN3/measurement_l1a.nc'];
        L1b_file = [getenv('HOME') '/TUDTUM/pseudo_DD_analysis/S3A_SR_1_SRA_BS_20211017T101723_20211017T110752_20211112T011946_3029_077_279______MAR_O_NT_004.SEN3/measurement_l1bs.nc'];
        L2_file = [getenv('HOME') '/TUDTUM/pseudo_DD_analysis/S3A_SR_2_WAT____20211017T101723_20211017T110458_20211112T021738_2855_077_279______MAR_O_NT_004.SEN3/enhanced_measurement.nc'];
        DOM = [53.5 57; -30 30] % actual DOM
        %DOM = [54 54.42; -30 30] % actual DOM
        %DOM_proc = DOM + [-lat_buffer lat_buffer; 0 0]
    case 'swh4.5m'
        L1a_file = [getenv('HOME') '/TUDTUM/pseudo_DD_analysis/S3A_SR_1_SRA_A__20211101T102835_20211101T111904_20211127T013420_3029_078_108______MAR_O_NT_004.SEN3/measurement_l1a.nc'];
        L1b_file = [getenv('HOME') '/TUDTUM/pseudo_DD_analysis/S3A_SR_1_SRA_BS_20211101T102835_20211101T111904_20211127T013420_3029_078_108______MAR_O_NT_004.SEN3/measurement_l1bs.nc'];
        L2_file = [getenv('HOME') '/TUDTUM/pseudo_DD_analysis/S3A_SR_2_WAT____20211101T102835_20211101T111602_20211127T024632_2847_078_108______MAR_O_NT_004.SEN3/enhanced_measurement.nc'];
        DOM = [54.5 58; -30 30] % actual DOM
        %DOM_proc = DOM + [-lat_buffer lat_buffer; 0 0]
    case 'S6test'
        L1a_file = [getenv('HOME') '/TUDTUM/ffsar_data/s6a/l1a-l1b-l2/crete/S6A_P4_1A_HR______20210901T212456_20210901T222027_20210902T134010_3331_030_018_009_EUM__OPE_ST_F03.nc'];
        L1b_file = [getenv('HOME') '/TUDTUM/ffsar_data/s6a/l1a-l1b-l2/crete/S6A_P4_1B_HR______20210901T212456_20210901T222025_20210902T134009_3329_030_018_009_EUM__OPE_ST_F03.nc'];
        L2_file = [getenv('HOME') '/TUDTUM/ffsar_data/s6a/l1a-l1b-l2/crete/S6A_P4_2__HR_STD__ST_030_018_20210901T212456_20210901T222025_F03.nc'];
        DOM = [33.55 33.60; 0, 360]; % open ocean
end

%% ###### read the L2 file data, make interpolator of sea surface #########

% 20 Hz
LAT                    = ncread(L2_file,'lat_20_ku');
LON                    = ncread(L2_file,'lon_20_ku');
IDX_mask = ingeoquad(LAT,LON,DOM(1,:),DOM(2,:));
IDX = find(IDX_mask);
startLoc=min(IDX);
count=max(IDX)-startLoc+1;

L2.LAT_20                    = ncread(L2_file,'lat_20_ku',startLoc,count);
L2.LON_20                    = ncread(L2_file,'lon_20_ku',startLoc,count);
L2.SWH_20                    = ncread(L2_file,'swh_ocean_20_ku',startLoc,count);
L2.SIG0_20                    = ncread(L2_file,'sig0_ocean_20_ku',startLoc,count);
L2.ALT_20                    = ncread(L2_file,'alt_20_ku',startLoc,count);
L2.range_20                  = ncread(L2_file,'range_ocean_20_ku',startLoc,count); %SAR mode : ocean/coastal retracking. Instrumental corrections included : USO drift correction (uso_cor_20_ku), internal path correction (int_path_cor_20_ku), distance antenna-COG (cog_cor_01), Doppler correction (dop_cor_20_ku), modeled instrumental errors correction (mod_instr_cor_range_01_ku) and system bias
L2.mean_ssh_20               = ncread(L2_file,'mean_sea_surf_sol1_20_ku',startLoc,count);
L2.tracker_range_20          = ncread(L2_file,'tracker_range_20_ku',startLoc,count); %SAR mode : reference range corrected for USO frequency drift (uso_cor_20_ku) and internal path corrections (int_path_cor_20_ku)

% range corrections to transform into l1a data:
L2.dop_cor_20                = ncread(L2_file,'dop_cor_20_ku',startLoc,count);

% 01 Hz
LAT                    = ncread(L2_file,'lat_01');
LON                    = ncread(L2_file,'lon_01');
IDX_mask = ingeoquad(LAT,LON,DOM(1,:),DOM(2,:));
IDX = find(IDX_mask);
startLoc=min(IDX);
count=max(IDX)-startLoc+1;

L2.geoid_01                  = ncread(L2_file,'geoid_01',startLoc,count);
L2.LAT_01                    = ncread(L2_file,'lat_01',startLoc,count);
L2.LON_01                    = ncread(L2_file,'lon_01',startLoc,count);
L2.ALT_01                    = ncread(L2_file,'alt_01',startLoc,count);

% range corrections to transform into l1a data:
L2.cog_cor_01                   = ncread(L2_file,'cog_cor_01',startLoc,count);
L2.mod_instr_cor_range_01    = ncread(L2_file,'mod_instr_cor_range_01_ku',startLoc,count);

%% remove cog_cor_01, dop_cor_20_ku, mod_instr_cor_range_01_ku from retracked range to align with tracker range
L2.range_20         =   L2.range_20 - L2.dop_cor_20 ...
                                    - interp1(L2.LAT_01, L2.cog_cor_01, L2.LAT_20,'linear') ...
                                    - interp1(L2.LAT_01, L2.mod_instr_cor_range_01, L2.LAT_20, 'linear');
%% some plotting (sanity check of retracked range versus geoid)
mov_mean_window = 20;

figure; 
plot(L2.LAT_20,L2.ALT_20 - L2.tracker_range_20); hold on
scatter(L2.LAT_20,L2.ALT_20 - L2.range_20); hold on
plot(L2.LAT_20,movmean(L2.ALT_20 - L2.range_20,mov_mean_window)); hold on
scatter(L2.LAT_20,L2.mean_ssh_20); hold on
scatter(L2.LAT_01,L2.geoid_01); hold on
legend('tracker range','retracked height','moving mean of retracked height','mean ssh', 'geoid')
ylabel('height over reference ellipsoid WGS84 (m)')
xlabel('altitude (deg)')

%% some plotting of SWH to define areas of interest
mov_mean_window = 20;

figure; 
subplot(1,2,1)
scatter(L2.LAT_20,L2.SWH_20); hold on
plot(L2.LAT_20,movmedian(L2.SWH_20,mov_mean_window),'LineWidth',3); hold on
legend('SWH','moving mean of SWH')
ylabel('SWH (m)')
xlabel('altitude (deg)')

subplot(1,2,2)
scatter(L2.LAT_20,L2.SIG0_20); hold on
plot(L2.LAT_20,movmedian(L2.SIG0_20,mov_mean_window),'LineWidth',3); hold on
legend('SIG0','moving mean of SIG0')
ylabel('SIG0 (/)')
xlabel('altitude (deg)')
ylim([0 20])

% area of interest for the case
if strcmp(dataset,'swh1.5m')
    DOM = [53.5 54.5; -30 30]
end

if strcmp(dataset,'swh4.5m')
    DOM = [55 56; -30 30]
end



%% build interpolator of geoid and retracked range
% for interpolation only ascending x-vectors are accepted
if sum(diff(L2.LAT_20))>0
    F_retracked = griddedInterpolant(L2.LAT_20,movmean(L2.ALT_20 - L2.range_20,mov_mean_window),'pchip','none');
    F_geoid = griddedInterpolant(L2.LAT_01,L2.geoid_01,'pchip','none');
else
    F_retracked = griddedInterpolant(flip(L2.LAT_20),flip(movmean(L2.ALT_20 - L2.range_20,mov_mean_window)),'pchip','none');
    F_geoid = griddedInterpolant(flip(L2.LAT_01),flip(L2.geoid_01),'pchip','none');
end
%% plot interpolators
figure;
plot(L2.LAT_20,F_retracked(L2.LAT_20)); hold on
plot(L2.LAT_20,F_geoid(L2.LAT_20)); hold on

%% plot interpolator differences
figure;
km_per_deg_lat = 111; % 1 deg altitude ~ 111 km, minor differences on 300 m
plot(111.*L2.LAT_20,F_retracked(L2.LAT_20)-F_geoid(L2.LAT_20)); hold on
ylabel('difference height (m)')
xlabel('latitude projected distance (km)')

%% ############## do some FF SAR processing ###############################
%% make some pseudo DD processing
if strcmp(mode,'UF')
    tic
    %FFSAR_processing_settings.combine_n_look_locations = 1;
    FFSAR_processing_settings.combine_n_look_locations = 1;
    FFSAR_processing_settings.output_pseudo_delay_doppler_processing = true;
    FFSAR_processing_settings.num_coherent_bursts = 1;
    FFSAR_processing_settings.integration_time = DDAcf.BRI*180.01; % corresponding to 180 bursts
    FFSAR_processing_settings.along_track_oversampling = 1/(0.25*585); % corresponding to ~80 Hz
    %FFSAR_processing_settings.along_track_oversampling = 1; % corresponding to ~80 Hz
    FFSAR_processing_settings.simulate_range_walk = true; % corresponding to ~80 Hz
    %FFSAR_processing_settings.integration_time = 2;
    %FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant'; % corresponding to ~80 Hz
    %DOM(3) = DOM(1)+0.025;
    FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant';

    ffproc = FFSAR_Processor(L1a_file, DOM, FFSAR_processing_settings, F_retracked);
    %ffproc = FFSAR_Processor(L1a_file, DOM, FFSAR_processing_settings);
    ffproc.setup_proc();
    ffproc.proc();
    CS1b = ffproc.CS1b;
    toc
elseif strcmp(mode,'FF')
    tic
    %FFSAR_processing_settings.combine_n_look_locations = 1;
    FFSAR_processing_settings.combine_n_look_locations = 50;
    FFSAR_processing_settings.output_pseudo_delay_doppler_processing = true;
    FFSAR_processing_settings.num_coherent_bursts = 1;
    FFSAR_processing_settings.integration_time = DDAcf.BRI*180.01; % corresponding to 180 bursts
    %FFSAR_processing_settings.along_track_oversampling = 1/(0.25*585); % corresponding to ~80 Hz
    FFSAR_processing_settings.along_track_oversampling = 1; % corresponding to ~80 Hz
    FFSAR_processing_settings.simulate_range_walk = false; % corresponding to ~80 Hz
    %FFSAR_processing_settings.integration_time = 2;
    %FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant'; % corresponding to ~80 Hz
    %DOM(3) = DOM(1)+0.025;
    FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant';

    ffproc = FFSAR_Processor(L1a_file, DOM, FFSAR_processing_settings, F_retracked);
    %ffproc = FFSAR_Processor(L1a_file, DOM, FFSAR_processing_settings);
    ffproc.setup_proc();
    ffproc.proc();
    CS1b = ffproc.CS1b;
    toc
end

save(['/home/fehelers/PhD_Delft/pseudoDDprocessing/' dataset '_' mode '.mat'],'ffproc','-v7.3' )
%%
load(['/home/fehelers/PhD_Delft/pseudoDDprocessing/' dataset '_' mode '.mat'],'ffproc')
CS1b = ffproc.CS1b;
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
%%
% tic
% %FFSAR_processing_settings.combine_n_look_locations = 1;
% FFSAR_processing_settings.combine_n_look_locations = 50;
% FFSAR_processing_settings.output_pseudo_delay_doppler_processing = true;
% FFSAR_processing_settings.num_coherent_bursts = 1;
% %FFSAR_processing_settings.integration_time = DDAcf.BRI*180.01; % corresponding to 180 bursts
% %FFSAR_processing_settings.along_track_oversampling = 1/(0.05*585); % corresponding to ~80 Hz
% FFSAR_processing_settings.along_track_oversampling = 1; % corresponding to ~80 Hz
% FFSAR_processing_settings.simulate_range_walk = false; % corresponding to ~80 Hz
% %FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant'; % corresponding to ~80 Hz
% %DOM = [55. 56.; -30 30]
% FFSAR_default_processing_settings.focal_point_height_interpolation = 'piecewise_constant';
% 
% %ffproc = FFSAR_Processor(L1a_file, DOM, FFSAR_processing_settings, F_retracked);
% ffproc = FFSAR_Processor(L1a_file, DOM, FFSAR_processing_settings);
% ffproc.setup_proc();
% ffproc.proc();
% CS1b = ffproc.CS1b;
% toc
%%
%save(['/home/fehelers/PhD_Delft/pseudoDDprocessing/' dataset '.mat'],'CS1b','-v7.3' )
%%
num_bursts_used = round(mean(CS1b.SAR.num_bursts));
%% plot FF SAR ref_range against interpolator
figure;
plot(L2.LAT_20,F_retracked(L2.LAT_20)); hold on
plot(CS1b.GEO.LAT,CS1b.GEO.H - CS1b.MEA.tracker_range,'ro'); hold on
plot(CS1b.GEO.LAT,CS1b.GEO.H - CS1b.MEA.win_delay*(DDAcf.c/2),'bo'); hold on

%% optional: convolute waveform with kernel
% make sure its centered symmetrically on odd support, so not to introduce a shift
kernel;
figure;plot(kernel);

%%
dummy = conv2(CS1b.SAR.data,kernel,'same');
figure;
ax1 = subplot(1,2,2)
imagesc(log(dummy))
ax1 = subplot(1,2,1)
imagesc(log(CS1b.SAR.data))

CS1b.SAR.data = dummy;
% %% plot waveforms (leading edge should coincide with bin 44)
% averaging = 585;
% figure;
% ax1 = subplot(1,2,1)
% imagesc(CS1b.GEO.LAT,1:128,movmean(CS1b.SAR.data,averaging,2)); hold on
% xlabel('latitude (deg)')
% ylabel('range bin')
% title('FF-SAR multilooked 20 Hz')
% ax2 = subplot(1,2,2)
% imagesc(CS1b.GEO.LAT,1:128,movmean(CS1b.SAR.data_pseudoDD,averaging,2)); hold on
% %imagesc(CS1b.GEO.LAT,1:128,movmean(CS1b.SAR.data_pseudoDD,1,2)); hold on
% title('emulated SAR multilooked 20 Hz')
% colormap('pink')
% xlabel('latitude (deg)')
% 
% linkaxes([ax1 ax2],'xy')
% %%
% figure;
% plot(CS1b.GEO.LAT(1:round(end/2)),movmean(CS1b.SAR.data(100,1:round(end/2)),averaging)); hold on
% xlabel('latitude (deg)')
% ylabel('power (arbitrary)')
% title('power along range bin 100')
% %plot(CS1b.GEO.LAT(1:round(end/2)),movmean(CS1b.SAR.data_pseudoDD(50,1:round(end/2)),averaging),'LineWidth',2); hold on
% plot(CS1b.GEO.LAT(1:round(end/2)),movmean(CS1b.SAR.data_pseudoDD(100,1:round(end/2)),averaging),'LineWidth',2); hold on
% legend('FF-SAR multilooked 20 Hz', 'emulated SAR multilooked 20 Hz')
% grid()


%% ############## compare pseudo D/D to FF-SAR averaged waveform #########
% for ENL calculation we need to first average the waveforms onto 20 Hz

time_tot = CS1b.GEO.Elapsed_Time(end)-CS1b.GEO.Elapsed_Time(1)
num_samples = time_tot*20
averaging_size = round(size(CS1b.SAR.data,2)/num_samples)
%%
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

% yyaxis left
plot(comp.avg_FF./max(comp.avg_FF),'Color',cFF,'LineWidth',lw,'LineStyle','-'); hold on
plot(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD),'Color',cDD,'LineWidth',lw,'LineStyle','--'); hold on
legend('FF-SAR','UF-SAR')
%xlabel('range bin')
ylabel('power')
ylim([0,1.05])

% yyaxis right
% p1 = plot(comp.ENL_FF,'Color',cFF,'LineWidth',lw,'LineStyle','-'); hold on
% p2 = plot(comp.ENL_pseudo_DD,'Color',cDD,'LineWidth',lw,'LineStyle','-'); hold on
% p1.Color(4) = 0.2;
% p2.Color(4) = 0.2;
% ylabel('20 Hz ENL')
% legend('FF-SAR','UF-SAR')
% ylim([0,500])

ax2 = subplot(4,1,3)
plot(100*(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD)-comp.avg_FF./max(comp.avg_FF))./(comp.avg_FF./max(comp.avg_FF)),'Color',[0.5,0.5,0.5]); hold on
%plot(100*(comp.avg_pseudo_DD_lin_scaled-comp.avg_L1b)./comp.avg_L1b); hold on
%title('(pseudo D/D - L1b) / L1b')
ylabel({'relative'; 'difference (%)'})
%xlabel('range bin')
%ylabel('%')
grid on
ylim([-2,2])

ax3 = subplot(4,1,4)
plot(100*(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD) - comp.avg_FF./max(comp.avg_FF)),'Color',[0.5,0.5,0.5]); hold on
%plot(100*(comp.avg_pseudo_DD_lin_scaled - comp.avg_L1b)/max(comp.avg_L1b)); hold on
%title('difference relative to maximum: pseudo D/D - L1b')
ylabel({'absolute'; 'difference (%)'})
xlabel('range bin')

grid on
ylim([-0.15,0.15])


linkaxes([ax1,ax2,ax3],'x')
xlim([0,numel(comp.avg_pseudo_DD)])

set(ax1,'xticklabel',[])
set(ax2,'xticklabel',[])

ax1.YAxis(1).Color = 'k';
% ax1.YAxis(2).Color = 'k';

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
%print -dpdf -painters epsFig

%exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/' dataset '.pdf'],'Resolution',300)
exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/' dataset '_' mode 'new.pdf'],'Resolution',300)
%saveas(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/' dataset '.pdf'])

%% ############## load L1b file and compare stack to waveform #############

[S3,CS1b_eum] = S3_L1bs_read(L1b_file,DOM,true);
S3.dop_beam_stack = complex(S3.i_echoes_ku_l1bs_echo_sar_ku, S3.q_echoes_ku_l1bs_echo_sar_ku);

%% apply the scaling factor to all doppler beam stacks
dummy = [];
dummy(1,:,:) = S3.iq_scale_factor_l1bs_echo_sar_ku;
S3.dop_beam_stack = S3.dop_beam_stack.*dummy;
clear dummy

%% apply agc scaling factor to L1b multilooked waveforms
S3.i2q2_meas_ku_l1bs_echo_sar_ku = S3.i2q2_meas_ku_l1bs_echo_sar_ku.*10.^(S3.agc_ku_l1bs_echo_sar_ku'./10);

%% compare i2q2_meas_ku_l1bs_echo_sar_ku with summed dop_beam_stack
num_wav = 1;

figure;
plot(S3.i2q2_meas_ku_l1bs_echo_sar_ku(:,num_wav));hold on
plot(nansum(abs(S3.dop_beam_stack(:,:,num_wav)).^2,2));
xlabel('range bin')
ylabel('waveform power after agc application')

% the result agrees well now
figure; plot(nansum(abs(S3.dop_beam_stack(:,:,num_wav)).^2,2)./(S3.i2q2_meas_ku_l1bs_echo_sar_ku(:,num_wav)))

% good agreement down to 0.2% within the waveform and 2% within the noise floor. Where does the remaining
% discrepancy come from?

%% ############## re-reference L1bs stack ##################################
% for this we use the difference between the ref_range interpolator F_retracked
% and the reference range in the L1bs product
S3.range_ku_l1bs_echo_sar_ku; % "Distance between the altimeter reference point and the surface height associated to a range gate used as reference inside the tracking window (reference tracking point), corrected for USO frequency drift and internal path correction"
S3.alt_l1bs_echo_sar_ku;
S3.new_reference_range = S3.alt_l1bs_echo_sar_ku - F_retracked(S3.lat_l1bs_echo_sar_ku);

S3.range_shift = S3.range_ku_l1bs_echo_sar_ku - S3.new_reference_range;

% plot range shift
figure;
subplot(1,3,1)
plot(S3.range_shift)
subplot(1,3,2)
% verify by plotting L1b waveforms
imagesc(S3.i2q2_meas_ku_l1bs_echo_sar_ku)
subplot(1,3,3)
imagesc(squeeze(nansum(abs(S3.dop_beam_stack).^2,2)))

%% now apply a fourier transform shift by the required range bin number for each waveform
dt       = 1/(DDAcf.B);          %Waveform sampling interval [s] % valid for no zero padding
dr       = dt*CONST.c/2;
N_samples = 128

wav = S3.dop_beam_stack;
wdims = size(S3.dop_beam_stack);
wav = reshape(wav, [wdims(1) wdims(2)*wdims(3)] );

% %% plotting to see the right dimensions
% figure;imagesc(abs(wav))
range_shift_upscaled = repmat(S3.range_shift',256,1);
range_shift_upscaled = range_shift_upscaled(:)';
% figure;plot(range_shift_upscaled)

%% perform range bin shift on all stacks simultaneously
shift_phasor = exp(-2*pi*1i*((1:N_samples)-1)'/N_samples*range_shift_upscaled/dr);

wav = fft(wav,[],1);
wav = wav.*shift_phasor;
wav = ifft(wav,[],1);

% reshape into stack and multilook, but beware that 180 bursts are 
% typically being used for multilooking, whereas we have num_bursts_used in
% the FF-SAR processing. Hence, we want to cut these away along the doppler
% beam dimension (more or less accurately, one more or less bursts is not expected to make a huge difference)

%b_edge = round((180-num_bursts_used)/2);
wav = reshape(wav,size(S3.dop_beam_stack));

%wav = nansum(abs(wav(:,b_edge:180-b_edge,:)).^2,2);
wav = nansum(abs(wav(:,:,:)).^2,2);
wav = squeeze(wav);
S3.wav_aligned = wav;
clear wav
%% illustrate the alignment against the original data
figure; 
subplot(2,1,1);
imagesc(S3.lat_l1bs_echo_sar_ku, 1:128, S3.wav_aligned)
title('L1b after alignment')
subplot(2,1,2);
imagesc(S3.lat_l1bs_echo_sar_ku, 1:128, S3.i2q2_meas_ku_l1bs_echo_sar_ku)
title('L1b before alignment')
colormap('pink')
% the borders can be nan-padded due to the sea surface interpolation, which does not
% extrapolate, just take care in computing averages


%% ########### plot the pseudo D/D results vs L1b waveforms ##############
[minval,idx] = min(abs(CS1b.GEO.LAT - S3.lat_l1bs_echo_sar_ku),[],2);
%size(idx)

% get indices where lat projected distance is smaller than 100 m:
L1bidx = find(111000*minval < 50);
FFidx = idx(111000*minval < 50);

% qualitative impression
figure;
subplot(2,1,1)
imagesc(S3.lat_l1bs_echo_sar_ku(L1bidx), 1:128, S3.wav_aligned(:,L1bidx))
title('L1b SAR')
subplot(2,1,2)
imagesc(CS1b.GEO.LAT(FFidx), 1:128, CS1b.SAR.data_pseudoDD(:,FFidx))
title('pseudo D/D SAR')
colormap('pink')

% sum over waveforms
figure;
subplot(4,1,1)
plot(S3.lat_l1bs_echo_sar_ku(L1bidx), sum(S3.wav_aligned(:,L1bidx),1))
title('cum. power L1b SAR')
subplot(4,1,2)
plot(CS1b.GEO.LAT(FFidx), sum(CS1b.SAR.data_pseudoDD(:,FFidx),1))
title('cum. power pseudo D/D SAR')
subplot(4,1,3)
plot(CS1b.GEO.LAT(FFidx), sum(CS1b.SAR.data_pseudoDD(:,FFidx),1)./sum(S3.wav_aligned(:,L1bidx),1))
title('cum. power ratio (pseudo D/D SAR  to  L1b SAR)')
grid()
subplot(4,1,4)
plot(CS1b.GEO.LAT(FFidx), max(CS1b.SAR.data_pseudoDD(:,FFidx),[],1)./max(S3.wav_aligned(:,L1bidx),[],1))
title('max power ratio (pseudo D/D SAR  to  L1b SAR)')
grid()
xlabel('latitude')

%% based on the matching, calculate average waveforms for L1b and pseudo D/D
% comp.avg_pseudo_DD = nanmean(CS1b.SAR.data_pseudoDD(:,FFidx),2)
% comp.avg_L1b = nanmean(S3.wav_aligned(:,L1bidx),2)

dummy = S3.wav_aligned(:,L1bidx);
dummy(isnan(CS1b.SAR.data_pseudoDD(:,FFidx))) = NaN;
S3.wav_aligned(:,L1bidx) = dummy;

%% based on the matching, calculate average waveforms for L1b and pseudo D/D
comp.avg_pseudo_DD = nanmean(CS1b.SAR.data_pseudoDD(:,FFidx),2);
offset_corr = -0.0008%12;
comp.avg_pseudo_DD = max(comp.avg_pseudo_DD)*(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD)-offset_corr)./(1-offset_corr);
comp.avg_L1b = nanmean(S3.wav_aligned(:,L1bidx),2);

comp.ENL_pseudo_DD = (nanvar(CS1b.SAR.data_pseudoDD(:,FFidx),[],2)./nanmean(CS1b.SAR.data_pseudoDD(:,FFidx),2).^2).^(-1)
comp.ENL_L1b = (nanvar(S3.wav_aligned(:,L1bidx),[],2)./nanmean(S3.wav_aligned(:,L1bidx),2).^2).^(-1)

%% figures show same behaviour, but in L1b case the array is filled with zeros
figure;subplot(1,2,1);imagesc(S3.wav_aligned(:,L1bidx));subplot(1,2,2);imagesc(CS1b.SAR.data_pseudoDD(:,FFidx))


% for matlab plotting of several y axes see https://uk.mathworks.com/help/matlab/creating_plots/plotting-with-two-y-axes.html
% for exportgraphics : https://uk.mathworks.com/help/matlab/creating_plots/save-figure-at-specific-size-and-resolution.html
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
ylim([0,500])

ax2 = subplot(4,1,3)
plot(100*(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD)-comp.avg_L1b./max(comp.avg_L1b))./(comp.avg_L1b./max(comp.avg_L1b)),'Color',[0.5,0.5,0.5],'LineWidth',lw); hold on
%plot(100*(comp.avg_pseudo_DD_lin_scaled-comp.avg_L1b)./comp.avg_L1b); hold on
%title('(pseudo D/D - L1b) / L1b')
ylabel({'relative'; 'difference (%)'})
%xlabel('range bin')
%ylabel('%')
grid on
ylim([-8,8])

ax3 = subplot(4,1,4)
plot(100*(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD) - comp.avg_L1b./max(comp.avg_L1b)),'Color',[0.5,0.5,0.5],'LineWidth',lw); hold on
%plot(100*(comp.avg_pseudo_DD_lin_scaled - comp.avg_L1b)/max(comp.avg_L1b)); hold on
%title('difference relative to maximum: pseudo D/D - L1b')
ylabel({'absolute'; 'difference (%)'})
xlabel('range bin')

grid on
ylim([-1.5,1.5])


linkaxes([ax1,ax2,ax3],'x')
xlim([0,numel(comp.avg_pseudo_DD)])

set(ax1,'xticklabel',[])
set(ax2,'xticklabel',[])

ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = 'k';

exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/' dataset '.pdf'],'Resolution',300)


%% ########### compare psuedo D/D results with L1b waveforms #############
%CS1b.SAR.data_pseudoDD
%size(CS1b.GEO.LAT)
%size(S3.lat_l1bs_echo_sar_ku)
%size(abs(CS1b.GEO.LAT - S3.lat_l1bs_echo_sar_ku))

[minval,idx] = min(abs(CS1b.GEO.LAT - S3.lat_l1bs_echo_sar_ku),[],2);
%size(idx)

% get indices where lat projected distance is smaller than 5 m:
L1bidx = find(111000*minval < 2);
FFidx = idx(111000*minval < 2);

N = numel(L1bidx)
figure;
for i = 1:min(N,5)
    subplot(min(N,5),1,i)
    
    % scaled with maximum
    plot(S3.wav_aligned(:,L1bidx(i)) / max(S3.wav_aligned(:,L1bidx(i)))  );hold on
    plot(CS1b.SAR.data_pseudoDD(:,FFidx(i))    /   max(CS1b.SAR.data_pseudoDD(:,FFidx(i)))   );hold on
    
    % divided by 10^6
    %plot(S3.wav_aligned(:,L1bidx(i)) );hold on
    %plot(CS1b.SAR.data_pseudoDD(:,FFidx(i))    /   1e6 )   ;hold on
    
    % without scale
    %plot(S3.wav_aligned(:,L1bidx(i)));hold on
    %plot(CS1b.SAR.data_pseudoDD(:,FFidx(i)));hold on
    
    % with hamming windowing correction
    %plot(S3.wav_aligned(:,L1bidx(i)));hold on
    %plot(1/1.5860^2*CS1b.SAR.data_pseudoDD(:,FFidx(i)));hold on
    
    % with hanning windowing correction
    %plot(    (  S3.wav_aligned(:,L1bidx(i))));hold on
    %plot(    (  1/1.6330^2*CS1b.SAR.data_pseudoDD(:,FFidx(i))));hold on
    
    % scaled with range bin sum
    %plot(S3.wav_aligned(:,L1bidx(i)) / sum(S3.wav_aligned(:,L1bidx(i)))  );hold on
    %plot(CS1b.SAR.data_pseudoDD(:,FFidx(i))    /   sum(CS1b.SAR.data_pseudoDD(:,FFidx(i)))   );hold on
    
    legend('Level-1b waveform', 'SAR pseudo D/D waveform')
    
end

figure;
for i = 1:min(N,5)
    subplot(min(N,5),1,i)
    
    % scaled with maximum
    %plot(S3.wav_aligned(:,L1bidx(i)) / max(S3.wav_aligned(:,L1bidx(i)))  );hold on
    %plot(CS1b.SAR.data_pseudoDD(:,FFidx(i))    /   max(CS1b.SAR.data_pseudoDD(:,FFidx(i)))   );hold on
    
    % without scale
    plot(CS1b.SAR.data_pseudoDD(:,FFidx(i))./1e6 ./ S3.wav_aligned(:,L1bidx(i)));hold on
end


figure;
for i = 1:min(N,5)
    subplot(min(N,5),1,i)
    
     % plot relative differences
    plot(  (  CS1b.SAR.data_pseudoDD(:,FFidx(i))/1e6 - S3.wav_aligned(:,L1bidx(i))  ) ./ S3.wav_aligned(:,L1bidx(i)));hold on
end

% online research: Hanning window does reduce the received power by a
% factor of 1/1.6330^2 = 1/2.6667. That is almost exactly the mismatch that
% I see between the pseudo D/D waveforms and the L1b waveforms.
% https://www.physicsforums.com/threads/energy-loss-of-spectra-associated-with-window-functions.202298/
% Moreover, this scaling seems to be applied to the L1b waveforms, although
% the result fits best with ours, when no windowing in azimuth is applied!
% This might be a bug on their side, or the factor is applied in order to
% keep untouched the waveform scale variable for sigma_0 or something...

% The above is outdated, other misalignments can be present on top of this.


%% ############## make some hacky retracking of both ######################
CS1b_eum;
CS1b_eum_shift = CS1b_eum;
CS1b_own_shift = CS1b_eum;

% now change only the data in CS1b.SAR.data and CS1b.MEA.ref_range
CS1b_eum_shift.SAR.data = S3.wav_aligned;
CS1b_own_shift.SAR.data = CS1b.SAR.data_pseudoDD(:,FFidx)/1e6;

CS1b_eum_shift.MEA.ref_range = CS1b.MEA.tracker_range(FFidx);
CS1b_own_shift.MEA.ref_range = CS1b.MEA.tracker_range(FFidx);
%%

[DATA_eum,CS]                           = SAR_L1b_to_L2(mission,CS1b_eum_shift,DOM,'SAMOSA2',{'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'});
[DATA_own,CS]                           = SAR_L1b_to_L2(mission,CS1b_own_shift,DOM,'SAMOSA2',{'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'});

%% plot results
figure;
subplot(2,1,1)
plot(DATA_eum.HEI);hold on
plot(DATA_own.HEI);hold on
title('SSH')
legend('L1b D/D','pseudo D/D')

subplot(2,1,2)
plot(DATA_eum.SWH);hold on
plot(DATA_own.SWH);hold on
title('SWH')
legend('L1b D/D','pseudo D/D')