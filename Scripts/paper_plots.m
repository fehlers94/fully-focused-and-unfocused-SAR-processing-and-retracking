clear all
clc
%%
% make default plot backgroudn white
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);

% choose color palette
cFF = [207, 0, 110]/255;
cDD = [0, 185, 227]/255;
cT = [214, 214, 214]/255;
%% transponder overpass impulse response functions for Sentinel-3B and S6
%% S3B case
mission = 'S3B'
data_dir = [getenv('HOME') '/TUDTUM/ffsar_data/'];
%fname ='s3a/l1a/crete/S3A_SR_1_SRA_A__20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement_l1a.nc'; 
%fname_l1b = '/home/fehelers/TUDTUM/ffsar_data/s3a/l1a-l1b-l2_sets/l1b/S3A_SR_1_SRA____20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement.nc';
fname = 's3b/l1a-l1b-l2_sets/l1a/S3B_SR_1_SRA_A__20190620T083403_20190620T092432_20191107T073220_3029_026_335______MR1_R_NT_004.SEN3/measurement_l1a.nc'
fname_l1b = '/home/fehelers/TUDTUM/ffsar_data/s3b/l1a-l1b-l2_sets/l1b/S3B_SR_1_SRA____20190620T083403_20190620T092432_20191107T073220_3029_026_335______MR1_R_NT_004.SEN3/measurement.nc'
%DOM = [35.33 35.37; 0 360]

%% ########### S3B case SVALBARD ###########
mission = 'S3B'
data_dir = [getenv('HOME') '/TUDTUM/ffsar_data/'];
%fname ='s3a/l1a/crete/S3A_SR_1_SRA_A__20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement_l1a.nc'; 
%fname_l1b = '/home/fehelers/TUDTUM/ffsar_data/s3a/l1a-l1b-l2_sets/l1b/S3A_SR_1_SRA____20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement.nc';
fname = 's3b/l1a/svalbard/S3B_SR_1_SRA_A__20191220T170439_20191220T175507_20200210T115323_3028_033_254______MR1_R_NT_004.SEN3/measurement_l1a.nc'
%fname_l1b = '/home/fehelers/TUDTUM/ffsar_data/s3b/l1a-l1b-l2_sets/l1b/S3B_SR_1_SRA____20190620T083403_20190620T092432_20191107T073220_3029_026_335______MR1_R_NT_004.SEN3/measurement.nc'
%DOM = [78.1 78.3; 15.3 15.5]

%% S6A case
mission = 'S6A'
data_dir = [getenv('HOME') '/TUDTUM/ffsar_data/'];
fname = ['s6a/l1a-l1b-l2/crete/S6A_P4_1A_HR______20210901T212456_20210901T222027_20210902T134010_3331_030_018_009_EUM__OPE_ST_F03.nc'];
fname_l1b = ['s6a/l1a-l1b-l2/crete/S6A_P4_1B_HR______20210901T212456_20210901T222025_20210902T134009_3329_030_018_009_EUM__OPE_ST_F03.nc'];
%L2_file = [getenv('HOME') '/TUDTUM/ffsar_data/s6a/l1a-l1b-l2/crete/S6A_P4_2__HR_STD__ST_030_018_20210901T212456_20210901T222025_F03.nc'];
%DOM = [35.33 35.35; 0 360]

%% set FF-SAR processing settings

DDAcf = DDA_ConfigFile(mission,'SAR');
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
FFSAR_processing_settings.combine_n_look_locations = 50;
FFSAR_processing_settings.output_pseudo_delay_doppler_processing = true;
FFSAR_processing_settings.num_coherent_bursts = 1;
FFSAR_processing_settings.simulate_range_walk = false;
FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant';
FFSAR_processing_settings.along_track_oversampling = 1;
FFSAR_processing_settings.split_aperture = true;
FFSAR_processing_settings.integration_time = 2;

% choose DOM symmetrically around the transponder
DOM = [TR.crete.lat-0.015 TR.crete.lat+0.015; -180 180];
%DOM = [TR.svalbard.lat-0.01 TR.svalbard.lat+0.01; -180 180];

%%

ffproc = FFSAR_Processor_sa([data_dir fname], DOM, FFSAR_processing_settings);
ffproc.setup_proc();
ffproc.proc();

CS1b = ffproc.CS1b;
%%
save(['/home/fehelers/PhD_Delft/pseudoDDprocessing/paperIRF_' mission '_T' num2str(FFSAR_processing_settings.integration_time) '.mat'],'CS1b','-v7.3' )
%% make IRF plot
%% double check factorization/normalization and dB recalculation

if strcmp(mission,'S6A')
    DDAcf.PRI = median(CS1b.MEA.PRI);
end

irfDD = CS1b.SAR.data_pseudoDD;
irfFF = CS1b.SAR.data;

PtotDD = sum(sum(irfDD))
PtotFF = sum(sum(irfFF))

irfDD = CS1b.SAR.data_pseudoDD/PtotDD;
irfFF = CS1b.SAR.data/PtotFF;

% calculate scale of DD along track PTR sinc function
DD_scale = (DDAcf.Nb*DDAcf.PRI)*2*DDAcf.fc*median(CS1b.GEO.V.V)/(DDAcf.c*median(CS1b.MEA.tracker_range)); %now, multiplied with d (y_a in Egido), we get along track the sinc argument
m_per_look = 1e3*deg2km(distance(CS1b.GEO.LAT(1),CS1b.GEO.LON(1),CS1b.GEO.LAT(end),CS1b.GEO.LON(end)))/numel(CS1b.GEO.LON);
N_looks = 1:size(CS1b.SAR.data,2);

% range sinc parameter
Range_scale = 1/DDAcf.os_ZP; % correct for S3 but also for S6!?
if strcmp(mission,'S6A')
    Range_scale = DDAcf.Bt/DDAcf.B/DDAcf.os_ZP
end


% calculate distance to transponder projected in along track direction
dist2Tr = 1e3*deg2km(distance(CS1b.GEO.LAT(:),CS1b.GEO.LON(:),TR.crete.lat,TR.crete.lon));
% this distance is an hyperpbolic, because of the remaining cross-track
% distance of the transponder, correct for this:
dist2Tr = sqrt(dist2Tr.^2 - min(dist2Tr.^2));
[val,ind_x] = min(dist2Tr);
dist2Tr(1:ind_x) = -dist2Tr(1:ind_x);

% position of transponder in range:
[val,ind_r] = max(irfDD(:,ind_x))

DD_theo = sinc(DD_scale*dist2Tr).^2/sum(sinc(DD_scale*dist2Tr).^2);
dr = 0.1*(1:DDAcf.Np*DDAcf.os_ZP*10);
if strcmp(mission,'S6A')
    offset = 0.5;
elseif strcmp(mission,'S3A')|strcmp(mission,'S3B')
    offset = 0;
end
Range_theo = 10*sinc(0.1*Range_scale*(((1:DDAcf.Np*DDAcf.os_ZP*10) - 10*(ind_r+offset)))).^2/ sum(sinc(0.1*Range_scale*(((1:DDAcf.Np*DDAcf.os_ZP*10) - 10*ind_r))).^2);


fig = figure('units','inch','position',[0,0,10,5]);
% first plot the integrated and cumulated power for DD:
ax1 = subplot(6,2,1);
plot(dist2Tr,cumsum(DD_theo),'LineWidth',4,'Color',cT); hold on
plot(dist2Tr,cumsum(sum(irfDD)),'Color',cDD); hold on
plot(dist2Tr,cumsum(sum(irfFF)),'Color',cFF); hold on
%legend('UF-SAR theoretical','UF-SAR','FF-SAR')
ylabel('\int P(x)/P_{tot}')%,'Interpreter','latex')
grid()
colorbar() % placeholder

ax2 = subplot(6,2,3);
plot(dist2Tr,10*log10(DD_theo),'LineWidth',4,'Color',cT); hold on
plot(dist2Tr,10*log10(sum(irfDD)),'Color',cDD); hold on
plot(dist2Tr,10*log10(sum(irfFF)),'Color',cFF); hold on
int_irfFF = movsum(sum(irfFF),300);
[peak_val,peak_ind] = findpeaks(int_irfFF,'MinPeakProminence',0.00001);
%plot(dist2Tr,10*log10(int_irfFF),'Color',cFF); hold on
s = scatter(dist2Tr(peak_ind),10*log10(int_irfFF(peak_ind)/max(int_irfFF(peak_ind))*max(sum(irfDD))),50,'filled','MarkerEdgeColor',cFF,'MarkerFaceColor',cFF); hold on
s.MarkerEdgeAlpha = 0.5;
s.MarkerFaceAlpha = .5;

legend('UF-SAR theoretical','UF-SAR','FF-SAR','FF-SAR integrated peak energy')
ylabel('P(x)/P_{tot} (dB)')
ylim([-80,0])
grid()
colorbar() % placeholder

% now along range axis
axr1 = subplot(3,6,10);
plot(10*log10(Range_theo),dr,'k','LineWidth',2,'Color',cT); hold on
plot(10*log10(sum(irfDD,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cDD); hold on
xlim([-45,0])
set(gca, 'YDir','reverse');
%set(gca, 'XDir','reverse');
%legend('UF-SAR theoretical','UF-SAR')
xlabel('P(r)/P_{tot} (dB)')
grid()

axr2 = subplot(3,6,16);
plot(10*log10(Range_theo),dr,'k','LineWidth',2,'Color',cT); hold on
plot(10*log10(sum(irfFF,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cFF); hold on
%plot(10*log10(irfFF(:,11650)/max(irfFF(:,11650))),1:DDAcf.Np*DDAcf.os_ZP,'Color',cFF); hold on
set(gca, 'YDir','reverse');
xlim([-45,0])
xlabel('P(r)/P_{tot} (dB)')
%set(gca, 'XDir','reverse');
%legend('FF-SAR theoretical','FF-SAR')
grid()

% now along range axis
axr3 = subplot(3,6,11);
plot(cumsum(Range_theo/10),dr,'k','LineWidth',2,'Color',cT); hold on
plot(cumsum(sum(irfDD,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cDD); hold on
set(gca, 'YDir','reverse');
%set(gca, 'XDir','reverse');
%legend('UF-SAR')
xlabel('\int P(r)/P_{tot}')
grid()

axr4 = subplot(3,6,17);
plot(cumsum(Range_theo/10),dr,'k','LineWidth',2,'Color',cT); hold on
plot(cumsum(sum(irfFF,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cFF); hold on
set(gca, 'YDir','reverse');
%set(gca, 'XDir','reverse');
%legend('FF-SAR')
xlabel('\int P(r)/P_{tot}')
grid()

% now plot irf's

ax3 = subplot(3,2,3);
imagesc(dist2Tr,1:DDAcf.Np*DDAcf.os_ZP,10*log10(irfDD))
cax = colorbar()
ylabel(cax,'P(r,x)/P_{tot} (dB)')
caxis([-80,0])
colormap('pink')
ylabel('range bin index')

ax4 = subplot(3,2,5);
imagesc(dist2Tr,1:DDAcf.Np*DDAcf.os_ZP,10*log10(irfFF))
cax = colorbar()
ylabel(cax,'P(r,x)/P_{tot} (dB)')
caxis([-80,0])
xlabel('transponder distance (m)')
ylabel('range bin index')

linkaxes([ax1,ax2,ax3,ax4],'x')
linkaxes([ax3,ax4,axr1,axr2,axr3,axr4],'y')
linkaxes([axr1,axr2],'x')

if strcmp(mission,'S6A')
    ax3.XLim = [-600,600]
    ax3.YLim = [30,90]%60*DDAcf.os_ZP]
end


if strcmp(mission,'S3B')
    ax3.XLim = [-650,650]
    ax3.YLim = [52-30,52+30]%60*DDAcf.os_ZP]
end

ax2.YTick = [-80 -60 -40 -20 0]

fig.Renderer='Painters';
%exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/paperIRF_' mission '.pdf'],'Resolution',300)
%TO-DO: put correct axes

%% plot only range responses P(r) and compare to model output
%run IRF_simulator before to get FFkernel, DDkernel
os_ZP_model = 8;
DDAcf.os_ZP = 2;

% now along range axis
figure('units','inch','position',[0,0,3,5]);
%axr1 = subplot(1,2,1);
plot( (Range_theo),dr,'k','LineWidth',2,'Color',cT); hold on
plot( (sum(irfDD,2)),1:DDAcf.Np*DDAcf.os_ZP,'o','Color',cDD); hold on
plot( (sum(irfFF,2)),1:DDAcf.Np*DDAcf.os_ZP,'+','Color',cFF); hold on


if strcmp(mission,'S3B')
    SHIFT = 0;
    N=length(DDkernel);
    mult = exp(-1i*2*pi/N*(0:N-1)*SHIFT);
    DDkernel2=abs(ifft((fft(DDkernel).*mult' ),'symmetric'));
    FFkernel2=abs(ifft((fft(FFkernel).*mult' ),'symmetric'));
    
    center = 52*os_ZP_model/DDAcf.os_ZP;
    offset = -128*os_ZP_model/DDAcf.os_ZP+center;
    
else
    SHIFT = -0.475*os_ZP_model/DDAcf.os_ZP;
    N=length(DDkernel);
    mult = exp(-1i*2*pi/N*(0:N-1)*SHIFT);
    DDkernel2=abs(ifft((fft(DDkernel).*mult' ),'symmetric'));
    FFkernel2=abs(ifft((fft(FFkernel).*mult' ),'symmetric'));
    
    center = 59*os_ZP_model/DDAcf.os_ZP;
    offset = -256*os_ZP_model/DDAcf.os_ZP+center;
end


scaleDD = 1.1*max(sum(irfDD,2))/max(DDkernel2)
scaleFF = 1.19*max(sum(irfFF,2))/max(FFkernel2)

plot( ((circshift(scaleDD*DDkernel2,offset))),(1:DDAcf.Np*os_ZP_model)*DDAcf.os_ZP/os_ZP_model,'Color',cDD); hold on
plot( ((circshift(scaleFF*FFkernel2,offset))),(1:DDAcf.Np*os_ZP_model)*DDAcf.os_ZP/os_ZP_model,'Color',cFF); hold on
%xlim([0,1])
set(gca, 'YDir','reverse');
%set(gca, 'XDir','reverse');
%legend('UF-SAR theoretical','UF-SAR')
ylabel('range bin index')
xlabel('P(r)/P_{tot} (no unit)')
%grid()
%plot( (Range_theo),dr,'k','LineWidth',2,'Color',cT); hold on

%plot(10*log10(irfFF(:,11650)/max(irfFF(:,11650))),1:DDAcf.Np*DDAcf.os_ZP,'Color',cFF); hold on
set(gca, 'YDir','reverse');
%xlim([-45,0])
%xlabel('P(r)/P_{tot} (dB)')
%set(gca, 'XDir','reverse');
%legend('FF-SAR theoretical','FF-SAR')
grid on

radius = 10
ylim([center*DDAcf.os_ZP/os_ZP_model-radius,center*DDAcf.os_ZP/os_ZP_model+radius])
xlim([-0.03 0.53])
legend('UF-SAR theoretical','UF-SAR','FF-SAR','UF-SAR model','FF-SAR model')
%xlim([-25,0])

% axr2 = subplot(1,2,2);
% 
% plot(cumsum(Range_theo/10),dr,'k','LineWidth',2,'Color',cT); hold on
% plot(cumsum(sum(irfDD,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cDD); hold on
% plot( cumsum(circshift(scaleDD*DDkernel,-128+52)),1:DDAcf.Np*DDAcf.os_ZP,'o','Color',cDD); hold on
% plot( cumsum(circshift(scaleFF*FFkernel,-128+52)),1:DDAcf.Np*DDAcf.os_ZP,'o','Color',cFF); hold on
% 
% set(gca, 'YDir','reverse');
% %set(gca, 'XDir','reverse');
% %legend('UF-SAR')
% xlabel('\int P(r)/P_{tot}')
% grid()
% 
% 
% plot(cumsum(Range_theo/10),dr,'k','LineWidth',2,'Color',cT); hold on
% plot(cumsum(sum(irfFF,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cFF); hold on
% set(gca, 'YDir','reverse');
% %set(gca, 'XDir','reverse');
% %legend('FF-SAR')
% xlabel('\int P(r)/P_{tot}')
% grid()
% 
% linkaxes([axr1 axr2],'y')

exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/paperIRF_P(r)_' mission '.pdf'],'Resolution',300)



%% plot IRF_DD slices to see whether there is a range shift or anything alike!
figure;
[val,ind] = max(irfDD,[],2);
ind = median(ind);
plot(log10(irfDD(:,ind)));hold on
plot(log10(irfDD(:,ind+1000)));hold on
plot(log10(irfDD(:,ind-1000)));hold on
plot(log10(irfDD(:,ind+3000)));hold on
plot(log10(irfDD(:,ind-3000)));hold on
xlim([0 120])
grid on

%% integrated IRF
figure;
plot((sum(irfDD,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cDD); hold on
plot((sum(irfFF,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cFF); hold on

%% plot ghosts of IRF in a nice way (zoomed)
figure('units','inch','position',[0,0,10,2]);
%ind = peak_ind(15:22)
scale_d =1.01;
off_d = 7.5;
spacing = round((1.23e4-1.1e4)/2);
center = 1.165e4;

for i=-3:3
    subplot(1,7,i+4)
    imagesc(dist2Tr+1.5,1:256,10*log10(irfFF))
    xlim([scale_d*(dist2Tr(center+i*spacing)-off_d+1.5) scale_d*(dist2Tr(center+i*spacing)+off_d+1.5)])
    ylim([52-30,52+30])
    colormap('pink')
    caxis([-70 0])
    if i==0
        xlabel('transponder distance (m)')
    end
    if i==-3
        ylabel('range bin index')
    end
end

exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/paperIRF_ghosts_' mission '.pdf'],'Resolution',300)

%% plot example of a contamination 
load(['/home/fehelers/PhD_Delft/pseudoDDprocessing/' 'S6test_FF.mat'],'ffproc')
T = ffproc.proc_sets.integration_time;
CS1b = ffproc.CS1b;
clear ffproc
%% convolve waveforms

dummy = conv2(CS1b.SAR.data,kernel,'same');
CS1b.SAR.data_conv = dummy;

%% differences of P(r)
figure;
plot(sum(irfDD,2)); hold on
plot(sum(irfFF,2)); hold on
figure;
plot(  (sum(irfFF,2) - sum(irfDD,2))./sum(irfDD,2)  );




%%

crange = [0,15];

figure('units','inch','position',[0,0,5,3]);
ax1 = subplot(1,3,1)
FFdat = movsum(CS1b.SAR.data,20,2);
imagesc(CS1b.GEO.LAT,1:512, (FFdat./mean(FFdat(:))))
caxis(crange)
title('FF-SAR')
colormap('pink')
ylabel('range bin index')

ax2 = subplot(1,3,2)
FFdat = movsum(CS1b.SAR.data_conv,20,2);
imagesc(CS1b.GEO.LAT,1:512,   (FFdat./mean(FFdat(:))))
caxis(crange)
title('convoluted FF-SAR')
colormap('pink')
xlabel('latitude (deg)')

ax3 = subplot(1,3,3)
DDdat = CS1b.SAR.data_pseudoDD;
imagesc(CS1b.GEO.LAT,1:512,   (DDdat./mean(DDdat(:))))
title('UF-SAR')
caxis(crange)
colormap('pink')

linkaxes([ax1,ax2,ax3],'xy')
xlim([56.278, 56.287])
ylim([30,160])

exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/paper_contamination_S6_conv.pdf'],'Resolution',300)