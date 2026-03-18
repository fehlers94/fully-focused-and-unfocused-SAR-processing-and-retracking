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
%% file for comparison
load('/home/fehelers/PhD_Delft/pseudoDDprocessing/S6test.mat','CS1b')

%%
figure;
ax1 = subplot(2,2,1);
ax2 = subplot(2,2,2);
ax3 = subplot(4,2,5);
ax4 = subplot(4,2,6);
%% now compare to L1b file
%% filter for the overlapping waveforms
[minval,idx] = min(abs(CS1b.GEO.LAT - CS1b_eum.GEO.LAT),[],2);
%size(idx)

% get indices where lat projected distance is smaller than 100 m:
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
comp.avg_pseudo_DD = nanmean(CS1b.SAR.data_pseudoDD(:,FFidx),2)
offset_corr = 0.0012;
comp.avg_pseudo_DD = max(comp.avg_pseudo_DD)*(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD)-offset_corr)./(1-offset_corr);
comp.avg_L1b = nanmean(CS1b_eum.SAR.data(:,L1bidx),2)

comp.ENL_pseudo_DD = (nanvar(CS1b.SAR.data_pseudoDD(:,FFidx),[],2)./nanmean(CS1b.SAR.data_pseudoDD(:,FFidx),2).^2).^(-1)
comp.ENL_L1b = (nanvar(CS1b_eum.SAR.data(:,L1bidx),[],2)./nanmean(CS1b_eum.SAR.data(:,L1bidx),2).^2).^(-1)

% for matlab plotting of several y axes see https://uk.mathworks.com/help/matlab/creating_plots/plotting-with-two-y-axes.html
%% plot comparison with max scale
figure;
ax1 = subplot(3,1,1)
yyaxis left
plot(comp.avg_L1b./max(comp.avg_L1b),'LineWidth',2); hold on
plot(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD),'r','LineWidth',2); hold on
legend('L1b','pseudo D/D')
xlabel('range bin')
ylabel('power')
ylim([0,1.05])

yyaxis right
plot(comp.ENL_L1b); hold on
plot(comp.ENL_pseudo_DD); hold on
ylabel('20 Hz ENL')
legend('L1b','pseudo D/D')

ax2 = subplot(3,1,2)
plot(100*(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD)-comp.avg_L1b./max(comp.avg_L1b))./(comp.avg_L1b./max(comp.avg_L1b))); hold on
%plot(100*(comp.avg_pseudo_DD_lin_scaled-comp.avg_L1b)./comp.avg_L1b); hold on
title('(pseudo D/D - L1b) / L1b')
legend('pseudo DD')
xlabel('range bin')
ylabel('%')
grid()
ylim([-10,10])

ax3 = subplot(3,1,3)
plot(100*(comp.avg_pseudo_DD./max(comp.avg_pseudo_DD) - comp.avg_L1b./max(comp.avg_L1b))); hold on
%plot(100*(comp.avg_pseudo_DD_lin_scaled - comp.avg_L1b)/max(comp.avg_L1b)); hold on
title('difference relative to maximum: pseudo D/D - L1b')
legend('pseudo DD')
xlabel('range bin')
ylabel('%')
grid()
ylim([-2,2])

linkaxes([ax1,ax2,ax3],'x')


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

figure;
% first plot the integrated and cumulated power for DD:
ax1 = subplot(6,2,1);
plot(dist2Tr,cumsum(DD_theo),'LineWidth',5,'Color',cT); hold on
plot(dist2Tr,cumsum(sum(irfDD)),'Color',cDD); hold on
plot(dist2Tr,cumsum(sum(irfFF)),'Color',cFF); hold on
legend('UF-SAR theoretical','UF-SAR','FF-SAR')
ylabel('P/P_{tot} (no unit)')
grid()
colorbar() % placeholder

ax2 = subplot(6,2,3);
plot(dist2Tr,10*log10(DD_theo),'LineWidth',5,'Color',cT); hold on
plot(dist2Tr,10*log10(sum(irfDD)),'Color',cDD); hold on
plot(dist2Tr,10*log10(sum(irfFF)),'Color',cFF); hold on
int_irfFF = movsum(sum(irfFF),300);
[peak_val,peak_ind] = findpeaks(int_irfFF,'MinPeakProminence',0.00001);
%plot(dist2Tr,10*log10(int_irfFF),'Color',cFF); hold on
s = scatter(dist2Tr(peak_ind),10*log10(int_irfFF(peak_ind)/max(int_irfFF(peak_ind))*max(sum(irfDD))),50,'filled','MarkerEdgeColor',cFF,'MarkerFaceColor',cFF); hold on
s.MarkerEdgeAlpha = 0.5;
s.MarkerFaceAlpha = .5;

legend('UF-SAR theoretical','UF-SAR','FF-SAR','FF-SAR integrated peak energy')
ylabel('P/P_{tot} (dB)')
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
legend('UF-SAR theoretical','UF-SAR')
xlabel('P/P_{tot} (dB)')
grid()

axr2 = subplot(3,6,16);
plot(10*log10(Range_theo),dr,'k','LineWidth',2,'Color',cT); hold on
plot(10*log10(sum(irfFF,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cFF); hold on
%plot(10*log10(irfFF(:,11650)/max(irfFF(:,11650))),1:DDAcf.Np*DDAcf.os_ZP,'Color',cFF); hold on
set(gca, 'YDir','reverse');
xlim([-45,0])
xlabel('P/P_{tot} (dB)')
%set(gca, 'XDir','reverse');
legend('FF-SAR theoretical','FF-SAR')
grid()

% now along range axis
axr3 = subplot(3,6,11);
plot(cumsum(Range_theo/10),dr,'k','LineWidth',2,'Color',cT); hold on
plot(cumsum(sum(irfDD,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cDD); hold on
set(gca, 'YDir','reverse');
%set(gca, 'XDir','reverse');
legend('UF-SAR')
xlabel('P/P_{tot} (no unit)')
grid()

axr4 = subplot(3,6,17);
plot(cumsum(Range_theo/10),dr,'k','LineWidth',2,'Color',cT); hold on
plot(cumsum(sum(irfFF,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cFF); hold on
set(gca, 'YDir','reverse');
%set(gca, 'XDir','reverse');
legend('FF-SAR')
xlabel('P/P_{tot} (no unit)')
grid()

% now plot irf's

ax3 = subplot(3,2,3);
imagesc(dist2Tr,1:DDAcf.Np*DDAcf.os_ZP,10*log10(irfDD))
cax = colorbar()
ylabel(cax,'P/P_{tot} (dB)')
caxis([-80,0])
colormap('pink')
ylabel('range bin index')

ax4 = subplot(3,2,5);
imagesc(dist2Tr,1:DDAcf.Np*DDAcf.os_ZP,10*log10(irfFF))
cax = colorbar()
ylabel(cax,'P/P_{tot} (dB)')
caxis([-80,0])
xlabel('transponder distance (m)')
ylabel('range bin index')

linkaxes([ax1,ax2,ax3,ax4],'x')
linkaxes([ax3,ax4,axr1,axr2,axr3,axr4],'y')
linkaxes([axr1,axr2],'x')

ax3.XLim = [-1320,1320]
ax3.YLim = [0,60*DDAcf.os_ZP]
%TO-DO: put correct axes

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


