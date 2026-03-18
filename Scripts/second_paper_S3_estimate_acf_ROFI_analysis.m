clear all
LoadCommonSettings
%%
L2_dir = '/home/fehelers/ownCloud/sarsim4rofi/Data/altimetry/ROFI_Level2';
L1_dir = fullfile(PathDATA,'ROFI','L1b');
mission = 'S3B'
orbit = '_370_'

%%
L2_files   = dir(fullfile(L2_dir,[mission,'*_SR_1_SRA_A*',orbit,'*']));
L1_files   = dir(fullfile(L1_dir,[mission,'*_SR_1_SRA_A*',orbit,'*']));

lat_threshold = 52.4; % degree
noise_threshold = 0.2; % m, 20cm noise threshold
movmedian_win_size = 240; % choose in accordance with the Hz product 
num_lags = 100;

%%
var = 'SWH'
[ACF_swh,RES_swh,SWH] = compute_acfs(L2_files(1:end),var,num_lags,lat_threshold,3,movmedian_win_size);
%%
var = 'HEI'
[ACF,RES,SWH] = compute_acfs(L2_files(1:end),var,num_lags,lat_threshold,noise_threshold,movmedian_win_size);
%plot(SWH,ACF(5,:),'.') % clearly showing a dependence on significant waveheight!
[ACF_thr,ACF_sig,ACF_r,RES_thr] = compute_acfs_thr(L1_files(1:end),num_lags,lat_threshold,2*noise_threshold,movmedian_win_size)


%% save results things in a struct:
acfdata.ACF = ACF;
acfdata.RES = RES;
acfdata.SWH = SWH;
acfdata.ACF_thr = ACF_thr;
acfdata.ACF_sig = ACF_sig;
acfdata.ACF_r = ACF_r;
acfdata.RES_thr = RES_thr;

save(['/home/fehelers/PhD_Delft/2nd_paper/acf_S3_ROFI.mat'],'acfdata','-v7.3' )

%%
load(['/home/fehelers/PhD_Delft/2nd_paper/acf_S3_ROFI.mat'],'acfdata')

ACF = acfdata.ACF;
RES = acfdata.RES;
SWH = acfdata.SWH;
ACF_thr = acfdata.ACF_thr;
ACF_sig = acfdata.ACF_sig;
ACF_r = acfdata.ACF_r;
RES_thr = acfdata.RES_thr;

%% visual inspection: throw out the outliers:

if numel(RES)>=42
    n = [36,32,30,24,19,13] % highest to lowest

    ACF(:,n) = [];
    ACF_sig(:,n) = [];
    ACF_swh(:,n) = [];
    ACF_thr(:,n) = [];
    ACF_r(:,n) = [];
    SWH(n) = [];
    RES(n) = [];
    RES_thr(n) = [];
end

%% make custom colormap
c2_rgb = str2rgb(c2);
c1_rgb = str2rgb(c1);

n_levels = 100;
level = (0:n_levels)'./n_levels;

custommap = level*c1_rgb + (1-level)*c2_rgb;

%%
% var = 'MF'
% [ACF_m,RES_m] = compute_acfs(L2_files(1:4),var,num_lags,lat_threshold,noise_threshold,movmedian_win_size);
%% analyse resulting autocorrelation function
figure;
imagesc(ACF)

figure;
imagesc(ACF_thr)

figure;
imagesc(ACF_sig)

figure;
imagesc(ACF_r)



%% load CS1b file to get calculate the impulse responses:
i = 1;
clear CS1b
CS1b = load(fullfile(L1_dir,L1_files(i).name));
CS1b = CS1b.CS1b;


%% then get the IRF_DD from script

%IRF_simulator

%% IRF_DD
DD_x = x;
DD_envelope = sum(IRF_DD);
DD_envelope = DD_envelope./max(DD_envelope);
DD_cut = IRF_DD(128,:);
DD_cut = DD_cut./max(DD_cut);

%%
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);

c1 = '#D81B60';
c2 = '#1E88E5';
c3 = '#FFC107';

%% plot ACF mean and std
figure('units','inch','position',[0,0,4.5,3.5]);

t = DD_x(8001:end);
s1 = DD_envelope(8001:end);
s2 = DD_cut(8001:end); 
p = patch([t fliplr(t)], [s1 fliplr(s2)], [0.8 0.8 0.8]);hold on
p.FaceAlpha = 0.5;
p.EdgeAlpha = 0;

grid

%errorbar(mean(RES)*(0:num_lags),(median(ACF,2,'omitnan')),1.96*1/sqrt(numel(RES))*std(ACF,[],2,'omitnan'));hold on
errorbar(mean(RES)*(0:num_lags),(median(ACF,2,'omitnan')),mad(ACF,1,2),'Color',c1);hold on
%plot(mean(RES)*(0:num_lags),ACF,'r');hold on


% downsample ACF
downsampling_f = 50;
data = ACF_thr(1:downsampling_f:end,:);
%errorbar(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data,2,'omitnan'),1.96*1/sqrt(numel(RES))*std(data,[],2,'omitnan'));hold on
errorbar(downsampling_f*median(RES_thr)*(0:size(data,1)-1),median(data,2,'omitnan'),mad(data,1,2),'Color',c2);hold on
%plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),data,'g');hold on

data = ACF_sig(1:downsampling_f:end,:);
%errorbar(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data,2,'omitnan'),1.96*1/sqrt(numel(RES))*std(data,[],2,'omitnan'));hold on
errorbar(downsampling_f*median(RES_thr)*(0:size(data,1)-1),median(data,2,'omitnan'),mad(data,1,2),'Color',c3);hold on
%plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),data,'b');hold on

plot(DD_x(8001:end),DD_envelope(8001:end),'k--');hold on
plot(DD_x(8001:end),DD_cut(8001:end),'k-.');hold on
%plot(DD_x(8001:end),0.5*DD_cut(8001:end)+0.5*DD_envelope(8001:end),'k.');hold on

xlim([0,700]);
xlabel('along-track lag (m)')
ylabel('Autocorrelation function')

legend('','detrended SSH (SAMOSA2)','detrended SSH (0.75 threshold)','cumulated waveform power','sinc^2(x/L_x)','along-track PTR cut')

%exportgraphics(gcf,['/home/fehelers/PhD_Delft/2nd_paper/ACF_ssh_estimates_SAMOSA2_threshold_power.pdf'],'Resolution',300)

%% plot ACF single realisations in neighbouring panels

% first get ACF_mod from SAM_model_gradients.m
% ACF_swh 
% index 1: Pu, index 2: epoch, index 3: SWH
% build interpolator based on HSs and ACF_swh
% first dimension is along-track, second SWH

[X1,X2] = ndgrid(DD_x,HSs);

C_Pu = griddedInterpolant(X1,X2,ACF_mod(:,:,1));
C_epoch = griddedInterpolant(X1,X2,ACF_mod(:,:,2));
C_SWH = griddedInterpolant(X1,X2,ACF_mod(:,:,3));

%% then plot the data with average and range

str2rgb=@(x)get(line('color',x),'color');
t = DD_x(8001:end);
s1 = DD_envelope(8001:end);
s2 = DD_cut(8001:end); 

figure('units','inch','position',[0,0,12,3.2]);

ax4 = subplot(1,4,4)
plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',5);hold on

ACF_SAMOSA_mean = C_epoch(t,mean(SWH)*ones(size(t)));
ACF_SAMOSA_min = C_epoch(t,min(SWH)*ones(size(t)));
ACF_SAMOSA_max = C_epoch(t,max(SWH)*ones(size(t)));

% t = DD_x(8001:end);
% s1 = DD_envelope(8001:end);
% s2 = 0*ones(size(DD_cut(8001:end))); 
% p = patch([t fliplr(t)], [s1 fliplr(s2)], [0.8 0.8 0.8]);hold on
% p.FaceAlpha = 0.5;
% p.EdgeAlpha = 0;

grid

%errorbar(mean(RES)*(0:num_lags),(median(ACF,2,'omitnan')),1.96*1/sqrt(numel(RES))*std(ACF,[],2,'omitnan'));hold on
%errorbar(mean(RES)*(0:num_lags),(median(ACF,2,'omitnan')),mad(ACF,1,2));hold on
plot(mean(RES)*(0:num_lags),mean(ACF,2),'Color',c1,'Marker','.','MarkerSize',10);hold on

tau = mean(RES)*(0:num_lags);
dum=ACF;
g1 = mean(dum,2)+std(dum,[],2);
g2 = mean(dum,2)-std(dum,[],2);
p = patch([tau fliplr(tau)], [g1' fliplr(g2')], str2rgb(c1));hold on
p.FaceAlpha = 0.2;
p.EdgeAlpha = 0;


% p = patch([t fliplr(t)], [g1 fliplr(g2)], [0. 0. 0.]);hold on
% p.FaceAlpha = 0.2;
% p.EdgeAlpha = 0;

g1 = ACF_SAMOSA_min;
g2 = ACF_SAMOSA_max; 
plot(t,[g1' g2']','k:','LineWidth',1.4);hold on

%plot(mean(RES)*(0:num_lags),prctile(ACF,75,2),'Color',c2);hold on
%plot(mean(RES)*(0:num_lags),prctile(ACF,25,2),'Color',c2);hold on


%errorbar(mean(RES)*(0:num_lags),mean(ACF,2),std(ACF,[],2),'Color',c2,'Marker','.','MarkerSize',8);hold on
%p = plot(t,[s1],'k-','LineWidth',1.5);
%p = plot(t,[s2],'k--','LineWidth',1.5);
p = plot(t,ACF_SAMOSA_mean,'k--','LineWidth',1);


grid on
%' ','','data: SSH(x)',
legend('sinc^2(x/L_x)','data','','','','','model','Location','northeast')

title({'SSH' '(SAMOSA2)'})

ax3 = subplot(1,4,3)
plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',5);hold on

% 
% t = DD_x(8001:end);
% s1 = DD_envelope(8001:end);
% s2 = DD_cut(8001:end); 
% p = patch([t fliplr(t)], [s1 fliplr(s2)], [0.8 0.8 0.8]);hold on
% p.FaceAlpha = 0.5;
% p.EdgeAlpha = 0;

% downsample ACF
downsampling_f = 1;
data = ACF_thr(1:downsampling_f:end,:);
%errorbar(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data,2,'omitnan'),1.96*1/sqrt(numel(RES))*std(data,[],2,'omitnan'));hold on
%errorbar(downsampling_f*median(RES_thr)*(0:size(data,1)-1),median(data,2,'omitnan'),mad(data,1,2));hold on
%plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),data,'Color',c2);hold on
%plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data,2),'Color',c3);hold on
plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data,2),'Color',c1,'Marker','.','MarkerSize',10);hold on

tau = downsampling_f*median(RES_thr)*(0:size(data,1)-1);
dum=data;
g1 = mean(dum,2)+std(dum,[],2);
g2 = mean(dum,2)-std(dum,[],2);
p = patch([tau fliplr(tau)], [g1' fliplr(g2')], str2rgb(c1));hold on
p.FaceAlpha = 0.2;
p.EdgeAlpha = 0;

%p = plot(t,[s1],'k-','LineWidth',1.5);
%p = plot(t,[s2],'k--','LineWidth',1.5);

%' ','','data: SSH(x)','Location',
legend('sinc^2(x/L_x)','data','Location','northeast')

xlabel('along-track distance (m)')
title({'SSH' '(0.75 threshold)'})
grid on



ax1 = subplot(1,4,1)
plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',5);hold on

% t = DD_x(8001:end);
% s1 = DD_envelope(8001:end);
% s2 = DD_cut(8001:end); 
% p = patch([t fliplr(t)], [s1 fliplr(s2)], [0.8 0.8 0.8]);hold on
% p.FaceAlpha = 0.5;
% p.EdgeAlpha = 0;

data = ACF_r(1:downsampling_f:end,:);
%plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),data,'Color',c3);hold on
%plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data,2),'Color',c1);hold on
plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data,2),'Color',c2,'Marker','.','MarkerSize',10);hold on
%p = plot(t,[s1],'k-','LineWidth',1.5);
p = plot(t,[s2],'k--','LineWidth',1);hold on;

tau = downsampling_f*median(RES_thr)*(0:size(data,1)-1);
dum=ACF_r(1:downsampling_f:end,:);
g1 = mean(dum,2)+std(dum,[],2);
g2 = mean(dum,2)-std(dum,[],2);
p = patch([tau fliplr(tau)], [g1' fliplr(g2')], str2rgb(c2));hold on
p.FaceAlpha = 0.2;
p.EdgeAlpha = 0;

%'','','data: P(r,x)','prediction: C_P(0,x)',
legend('sinc^2(x/L_x)','data','model','Location','northeast')
ylabel('Noise autocorrelation function')

title({'waveform power' '(single rangebins)'})
grid on


ax2 = subplot(1,4,2)
plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',5);hold on

% t = DD_x(8001:end);
% s1 = DD_envelope(8001:end);
% s2 = DD_cut(8001:end); 
% p = patch([t fliplr(t)], [s1 fliplr(s2)], [0.8 0.8 0.8]);hold on
% p.FaceAlpha = 0.5;
% p.EdgeAlpha = 0;

data = ACF_sig(1:downsampling_f:end,:);
%errorbar(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data,2,'omitnan'),1.96*1/sqrt(numel(RES))*std(data,[],2,'omitnan'));hold on
%errorbar(downsampling_f*median(RES_thr)*(0:size(data,1)-1),median(data,2,'omitnan'),mad(data,1,2));hold on
%plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),data,'Color',c3);hold on
%plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data,2),'Color',c1);hold on
plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data,2),'Color',c2,'Marker','.','MarkerSize',10);hold on
p = plot(t,[s1],'k--','LineWidth',1.5);hold on;
%p = plot(t,[s2],'k--','LineWidth',1.5);

tau = downsampling_f*median(RES_thr)*(0:size(data,1)-1);
dum=data;
g1 = mean(dum,2)+std(dum,[],2);
g2 = mean(dum,2)-std(dum,[],2);
p = patch([tau fliplr(tau)], [g1' fliplr(g2')], str2rgb(c2));hold on
p.FaceAlpha = 0.2;
p.EdgeAlpha = 0;
%'','','data: \Sigma_rP(r,x)','prediction: \Sigma_r C_P(r,x)',
legend('sinc^2(x/L_x)','data','model','Location','northeast')


title({'cumulated waveform power' '(over range)'})

%plot(DD_x(8001:end),DD_envelope(8001:end),'k--');hold on
%plot(DD_x(8001:end),DD_cut(8001:end),'k-.');hold on
%plot(DD_x(8001:end),0.5*DD_cut(8001:end)+0.5*DD_envelope(8001:end),'k.');hold on
grid on
linkaxes([ax1 ax2 ax3 ax4],'xy')
xlim([0,320]);
ylim([-0.08,1])

%legend('','detrended SSH (SAMOSA2)','detrended SSH (0.75 threshold)','cumulated waveform power','sinc^2(x/L_x)','along-track PTR cut')

%exportgraphics(gcf,['/home/fehelers/PhD_Delft/2nd_paper/ACF_ssh_estimates_SAMOSA2_threshold_power_signal_seperate.pdf'],'Resolution',300)

%% then plot the data but sorted for sea state

str2rgb=@(x)get(line('color',x),'color');
t = DD_x(8001:end);
s1 = DD_envelope(8001:end);
s2 = DD_cut(8001:end); 
cs = {c1,c2,c3};
ms = {'v','o','^'};

figure('units','inch','position',[0,0,15,2.5]);

ax4 = subplot(1,5,4)
plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',5);hold on
grid
bin_borders = [0,1,2,3]
plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',5);hold on
for k =1:numel(bin_borders)-1
    bin_mask = (bin_borders(k)<SWH)&(SWH < bin_borders(k+1));
    plot(x,mean(ACF(:,bin_mask),2),'color',cs{k},'Marker',ms{k},'MarkerSize',5);hold on
end

grid on
%legend('sinc^2(x/L_x)','data','','','','','model','Location','northeast')
title({'SSH' '(SAMOSA2)'})

ax5 = subplot(1,5,5)
plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',5);hold on
grid
bin_borders = [0,1,2,3]
plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',5);hold on
for k =1:numel(bin_borders)-1
    bin_mask = (bin_borders(k)<SWH)&(SWH < bin_borders(k+1));
    plot(x,mean(ACF_swh(:,bin_mask),2),'color',cs{k},'Marker',ms{k},'MarkerSize',5);hold on
end

grid on
%legend('sinc^2(x/L_x)','data','','','','','model','Location','northeast')
title({'SWH' '(SAMOSA2)'})



ax3 = subplot(1,5,3)
plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',5);hold on

% downsample ACF
downsampling_f = 1;
data = ACF_thr(1:downsampling_f:end,:);
%plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data,2),'Color',c1,'Marker','.','MarkerSize',10);hold on
for k =1:numel(bin_borders)-1
    bin_mask = (bin_borders(k)<SWH)&(SWH < bin_borders(k+1));
    plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data(:,bin_mask),2),'color',cs{k},'Marker',ms{k},'MarkerSize',5);hold on
end
%legend('sinc^2(x/L_x)','data','Location','northeast')

xlabel('along-track distance (m)')
title({'SSH' '(0.75 threshold)'})
grid on

ax1 = subplot(1,5,1)
plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',5);hold on

data = ACF_r(1:downsampling_f:end,:);
for k =1:numel(bin_borders)-1
    bin_mask = (bin_borders(k)<SWH)&(SWH < bin_borders(k+1));
    plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data(:,bin_mask),2),'color',cs{k},'Marker',ms{k},'MarkerSize',5);hold on
end

%plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data,2),'Color',c2,'Marker','.','MarkerSize',10);hold on

p = plot(t,[s2],'k--','LineWidth',1);hold on;

%'','','data: P(r,x)','prediction: C_P(0,x)',
%legend('sinc^2(x/L_x)','data','model','Location','northeast')
ylabel('Noise autocorrelation function')

title({'waveform power' '(single rangebins)'})
grid on


ax2 = subplot(1,5,2)
plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',5);hold on

data = ACF_sig(1:downsampling_f:end,:);
%plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data,2),'Color',c2,'Marker','.','MarkerSize',10);hold on

for k =1:numel(bin_borders)-1
    bin_mask = (bin_borders(k)<SWH)&(SWH < bin_borders(k+1));
    plot(downsampling_f*median(RES_thr)*(0:size(data,1)-1),mean(data(:,bin_mask),2),'color',cs{k},'Marker',ms{k},'MarkerSize',5);hold on
end

p = plot(t,[s1],'k--','LineWidth',1.5);hold on;

legend('sinc^2(x/L_x)','SWH 0-1 m','SWH 1-2 m','SWH 2-3 m','model','Location','northeast')


title({'cumulated waveform power' '(over range)'})

%plot(DD_x(8001:end),DD_envelope(8001:end),'k--');hold on
%plot(DD_x(8001:end),DD_cut(8001:end),'k-.');hold on
%plot(DD_x(8001:end),0.5*DD_cut(8001:end)+0.5*DD_envelope(8001:end),'k.');hold on
grid on
linkaxes([ax1 ax2 ax3 ax4 ax5],'xy')
xlim([0,330]);
ylim([-0.08,1])

%legend('','detrended SSH (SAMOSA2)','detrended SSH (0.75 threshold)','cumulated waveform power','sinc^2(x/L_x)','along-track PTR cut')

exportgraphics(gcf,['/home/fehelers/PhD_Delft/2nd_paper/measured_ACFs_sea_state.pdf'],'Resolution',300)



%% now plot the data in the along-track distance 0:320 m and compare to modeled data

tau = mean(RES)*(0:size(ACF,1)-1)

ACF_dyn_model = zeros(size(ACF))
ACF_mean_model = zeros(size(ACF))
ACF_sinc_model = zeros(size(ACF))

for i = 1:numel(SWH)
    ACF_dyn_model(:,i) = C_epoch(tau,SWH(i)*ones(size(tau)));
    ACF_mean_model(:,i) = C_epoch(tau,mean(SWH)*ones(size(tau)));
    %ACF_mean_model(:,i) = C_epoch(tau,2*ones(size(tau)));
    ACF_sinc_model(:,i) = sinc(tau*DD_scale);
end

% figure;imagesc(ACF_dyn_model(1:30,:))
% figure;imagesc(ACF(1:30,:))

%find maximum lag index to consider:
maxind = max(find(tau<320));
maxind = 12;

% make scatterplots of models

figure('units','inch','position',[0,0,5,5]);
mod = ACF_sinc_model(1:maxind,:);
dat = ACF(1:maxind,:);
scatter(mod(:),dat(:),'MarkerEdgeColor',c3,'Marker','o');hold on
rmse_sinc = sqrt(mean((mod(:)-dat(:)).^2))
% R_sinc = corrcoef(mod(:),dat(:));
% R_sinc = R_sinc(1,2)

mod = ACF_mean_model(1:maxind,:);
scatter(mod(:),dat(:),'MarkerEdgeColor',c2,'Marker','^');hold on
rmse_mean = sqrt(mean((mod(:)-dat(:)).^2))
% R_mean = corrcoef(mod(:),dat(:));
% R_mean = R_mean(1,2)

mod = ACF_dyn_model(1:maxind,:);
scatter(mod(:),dat(:),60,'MarkerEdgeColor',c1,'Marker','.');hold on
rmse_dyn = sqrt(mean((mod(:)-dat(:)).^2))
% R_dyn = corrcoef(mod(:),dat(:));
% R_dyn = R_dyn(1,2)

pbaspect([1 1 1])
xlabel('predicted SSH noise correlation')
ylabel('measured SSH noise correlation')

plot([0,1],[0,1],'k-','LineWidth',2);hold on
grid on

legend(...
['sinc^2 model, RMSE=',num2str(rmse_sinc,2)],...
['model (SWH = 1.4 m), RMSE=',num2str(rmse_mean,2)],...
['model (dynamic SWH), RMSE=',num2str(rmse_dyn,2)],...
'',...
'Location','northoutside')

%exportgraphics(gcf,['/home/fehelers/PhD_Delft/2nd_paper/ACFs_model_sea_state_dependence.pdf'],'Resolution',300)

%% test coloring the data with SWH for plotting
x = median(RES)*(0:size(ACF,1)-1);

cmin = 0;
cmax = 4;

subplot(1,2,1)

[SWH_new,idx] = sort(SWH);
ACF_new = ACF_swh(:,idx);

str2rgb=@(x)get(line('color',x),'color');


plotcolors = pink(100)%custommap;

bin_borders = [0,1,2,3]

plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',5);hold on
for i = 1:numel(SWH)
    %plot(x,ACF_new(:,i),'color',plotcolors(round(SWH_new(i)/4*100),:),'Linewidth',2);hold on
    for k =1:numel(bin_borders)-1
        bin_mask = (bin_borders(k)<SWH_new)&(SWH_new < bin_borders(k+1));
        plot(x,mean(ACF_new(:,bin_mask),2),'color',cs{k},'Linewidth',1,'Marker','.');hold on
    end
end
%plot(x,C_SWH(x,0.2*ones(size(x))),'k--','Linewidth',2);hold on
%plot(x,C_SWH(x,4.*ones(size(x))),'k--','Linewidth',2);hold on



xlim([0,330])
ylim([-0.1,1])
grid on

% colorbar
% colormap('pink')
caxis([cmin,cmax])


subplot(1,2,2)

[SWH_new,idx] = sort(SWH);
ACF_new = ACF(:,idx);

plotcolors = pink(100);
cs = {c1,c2,c3};

plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',5);hold on
for i = 1:numel(SWH)
    %plot(x,ACF_new(:,i),'color',plotcolors(round(SWH_new(i)/4*100),:),'Linewidth',2);hold on
    for k =1:numel(bin_borders)-1
        bin_mask = (bin_borders(k)<SWH_new)&(SWH_new < bin_borders(k+1));
        plot(x,mean(ACF_new(:,bin_mask),2),'color',cs{k},'Linewidth',1);hold on
        mean(SWH_new(bin_mask))
    end
end
%plot(x,C_epoch(x,0.2*ones(size(x))),'k--','Linewidth',2);hold on
%plot(x,C_epoch(x,4.*ones(size(x))),'k--','Linewidth',2);hold on


xlim([0,330])
ylim([-0.1,1])
grid on

% colorbar
% colormap('pink')
caxis([cmin,cmax])


%% plot Pu, epoch and SWH dependence on sea state according to full SAMOSA 2 model:

SWH_vals = [0.2 1.5 3 5 10];


figure('units','inch','position',[0,0,12,3.5]);

M = 15;
M2 = 1500;


subplot(1,3,1)
plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',3);hold on
for h = SWH_vals
    plot(t(1:M:M2),C_Pu(t(1:M:M2),h*ones(size(t(1:M:M2)))),'k'); hold on
end
title('Pu / sigma0')
set(gca,'xtick',[0:30:330])
set(gca,'ytick',[0:0.1:1])
legend('sinc^2','model')
xlim([0,330])

ylabel('modelled noise autocorrelation')
grid on

subplot(1,3,2)
plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',3);hold on

for h = SWH_vals
    plot(t(1:M:M2),C_epoch(t(1:M:M2),h*ones(size(t(1:M:M2)))),'k'); hold on
end
set(gca,'xtick',[0:30:330])
set(gca,'ytick',[0:0.1:1])
legend('sinc^2','model')

title('epoch / range / SSH')
xlim([0,330])
xlabel('along-track distance (m)')
grid on


subplot(1,3,3)
plot(t,s1,'Color',[0.8 0.8 0.8],'Linewidth',3);hold on
for h = SWH_vals
    plot(t(1:M:M2),C_SWH(t(1:M:M2),h*ones(size(t(1:M:M2)))),'k'); hold on
end
set(gca,'xtick',[0:30:330])
set(gca,'ytick',[0:0.1:1])
legend('sinc^2','model')

title('SWH')
xlim([0,330])
grid on

%exportgraphics(gcf,['/home/fehelers/PhD_Delft/2nd_paper/ACFs_model_sea_state_dependence_nodata.pdf'],'Resolution',300)

%% calculate stds and acf for reduced 20 Hz products

L2_files(1:end)
n = numel(L2_files);
var = 'HEI'
lat_threshold
noise_threshold

CORR20 = zeros(n,6);
STD20 = zeros(n,6);

% order is uf80,ff80,uf40,ff40,uf20,ff20, uf20smooth

for i = 1:n
    load(fullfile(L2_files(i).folder,L2_files(i).name));

    CS2.f_080Hz.UF.HEI_d = detrend_and_filter(CS2.f_080Hz.UF,var,lat_threshold, noise_threshold, 80);
    CS2.f_080Hz.FF.HEI_d = detrend_and_filter(CS2.f_080Hz.FF,var,lat_threshold, noise_threshold, 80);
    CS2.f_040Hz.UF.HEI_d = detrend_and_filter(CS2.f_040Hz.UF,var,lat_threshold, noise_threshold, 40);
    CS2.f_040Hz.FF.HEI_d = detrend_and_filter(CS2.f_040Hz.FF,var,lat_threshold, noise_threshold, 40);
    CS2.f_020Hz.UF.HEI_d = detrend_and_filter(CS2.f_020Hz.UF,var,lat_threshold, noise_threshold, 20);
    CS2.f_020Hz.FF.HEI_d = detrend_and_filter(CS2.f_020Hz.FF,var,lat_threshold, noise_threshold, 20);

%     figure;
%     plot(CS2.f_080Hz.UF.LAT,CS2.f_080Hz.UF.HEI_d);hold on
%     plot(CS2.f_080Hz.FF.LAT,CS2.f_080Hz.FF.HEI_d);hold on
%     plot(CS2.f_020Hz.UF.LAT,CS2.f_020Hz.UF.HEI_d);hold on
%     plot(CS2.f_020Hz.FF.LAT,CS2.f_020Hz.FF.HEI_d);hold on

%     % compute acfs on original sampling
%     CS2.f_080Hz.UF.acf = autocorr(CS2.f_080Hz.UF.HEI_d,80); hold on
%     CS2.f_080Hz.FF.acf = autocorr(CS2.f_080Hz.FF.HEI_d,80); hold on
%     CS2.f_020Hz.UF.acf = autocorr(CS2.f_020Hz.UF.HEI_d,20); hold on
%     CS2.f_020Hz.FF.acf = autocorr(CS2.f_020Hz.FF.HEI_d,20); hold on

    % compute stds (no surprises, the FF-SAR 20 is the best because of highest correlation)
%     std(CS2.f_080Hz.UF.HEI_d,'omitnan')
%     std(CS2.f_080Hz.FF.HEI_d,'omitnan')
%     std(CS2.f_020Hz.UF.HEI_d,'omitnan')
%     std(CS2.f_020Hz.FF.HEI_d,'omitnan')

    % perform reduction on 20 Hz
    CS2.f_080Hz.UF.HEI_d20 = movmean(CS2.f_080Hz.UF.HEI_d,4);
    CS2.f_080Hz.FF.HEI_d20 = movmean(CS2.f_080Hz.FF.HEI_d,4);
    CS2.f_040Hz.UF.HEI_d20 = movmean(CS2.f_040Hz.UF.HEI_d,2);
    CS2.f_040Hz.FF.HEI_d20 = movmean(CS2.f_040Hz.FF.HEI_d,2);
    CS2.f_020Hz.UF.HEI_d20 = CS2.f_020Hz.UF.HEI_d;
    CS2.f_020Hz.FF.HEI_d20 = CS2.f_020Hz.FF.HEI_d;
    CS2.f_020Hz.UF.HEI_d20_smooth = conv(CS2.f_020Hz.UF.HEI_d20,[0.88,0.12]);

%     figure;
%     plot(CS2.f_080Hz.UF.LAT,CS2.f_080Hz.UF.HEI_d20);hold on
%     plot(CS2.f_080Hz.FF.LAT,CS2.f_080Hz.FF.HEI_d20);hold on
%     plot(CS2.f_020Hz.UF.LAT,CS2.f_020Hz.UF.HEI_d20);hold on
%     plot(CS2.f_020Hz.FF.LAT,CS2.f_020Hz.FF.HEI_d20);hold on

    % compute acfs on original sampling
    [CS2.f_080Hz.UF.acf,CS2.f_080Hz.UF.lags] = autocorr(CS2.f_080Hz.UF.HEI_d20,80); hold on
    [CS2.f_080Hz.FF.acf,CS2.f_080Hz.FF.lags] = autocorr(CS2.f_080Hz.FF.HEI_d20,80); hold on
    [CS2.f_040Hz.UF.acf,CS2.f_040Hz.UF.lags] = autocorr(CS2.f_040Hz.UF.HEI_d20,40); hold on
    [CS2.f_040Hz.FF.acf,CS2.f_040Hz.FF.lags] = autocorr(CS2.f_040Hz.FF.HEI_d20,40); hold on
    [CS2.f_020Hz.UF.acf,CS2.f_020Hz.UF.lags] = autocorr(CS2.f_020Hz.UF.HEI_d20,20); hold on
    [CS2.f_020Hz.FF.acf,CS2.f_020Hz.FF.lags] = autocorr(CS2.f_020Hz.FF.HEI_d20,20); hold on
    [CS2.f_020Hz.UF.acf_smooth,CS2.f_020Hz.UF.lags_smooth] = autocorr(CS2.f_020Hz.UF.HEI_d20_smooth,20); hold on

    CORR20(i,1) = CS2.f_080Hz.UF.acf(5) % 0 0.25 0.5 0.75 1
    CORR20(i,2) = CS2.f_080Hz.FF.acf(5)
    CORR20(i,3) = CS2.f_040Hz.UF.acf(3) % 0 0.5 1
    CORR20(i,4) = CS2.f_040Hz.FF.acf(3)
    CORR20(i,5) = CS2.f_020Hz.UF.acf(2) % 0 1
    CORR20(i,6) = CS2.f_020Hz.FF.acf(2)
    %CORR20(i,7) = CS2.f_020Hz.UF.acf_smooth(2)
    
    
%     figure;
%     plot(0.25*CS2.f_080Hz.UF.lags,CS2.f_080Hz.UF.acf,'o'); hold on
%     plot(0.25*CS2.f_080Hz.FF.lags,CS2.f_080Hz.FF.acf,'o'); hold on
%     plot(CS2.f_020Hz.UF.lags,CS2.f_020Hz.UF.acf); hold on
%     plot(CS2.f_020Hz.FF.lags,CS2.f_020Hz.FF.acf); hold on
%     plot(CS2.f_020Hz.UF.lags_smooth,CS2.f_020Hz.UF.acf_smooth); hold on
% 
%     legend('80to20 UF','80to20 FF','20UF','20FF','20UF smoothed')

    % compute stds
    STD20(i,1) = std(CS2.f_080Hz.UF.HEI_d20,'omitnan')
    STD20(i,2) = std(CS2.f_080Hz.FF.HEI_d20,'omitnan')
    STD20(i,3) = std(CS2.f_040Hz.UF.HEI_d20,'omitnan')
    STD20(i,4) = std(CS2.f_040Hz.FF.HEI_d20,'omitnan')
    STD20(i,5) = std(CS2.f_020Hz.UF.HEI_d20,'omitnan')
    STD20(i,6) = std(CS2.f_020Hz.FF.HEI_d20,'omitnan')
    %STD20(i,7) = std(CS2.f_020Hz.UF.HEI_d20_smooth,'omitnan')
end

%% illustrate the results

figure('units','inch','position',[0,0,10.5,3.5]);

% sort for decreasing std
[val,ind] = sort(median(STD20),'descend')
names = {'80Hz UF'; '80Hz FF'; '40Hz UF'; '40Hz FF'; '20Hz UF';'20Hz FF'};

%figure;

subplot(1,2,1)
data = median(STD20);
err = mad(STD20,1);

%errorbar(1:numel(data),data(ind),err(ind),'.');
boxplot(100*STD20(:,ind),'Color','k');
set(gca,'xtick',[1:6],'xticklabel',names(ind))
grid on
ylabel({'noise estimate (cm):';'20 Hz SSH std after 1 Hz detrending'})
ylim([2,5.5])

subplot(1,2,2)
data = median(CORR20);
err = mad(CORR20,1);
%errorbar(1:numel(data),data(ind),err(ind),'.');
boxplot(CORR20(:,ind),'Color','k');
set(gca,'xtick',[1:6],'xticklabel',names(ind))
grid on
ylabel({'correlation between';'consecutive 20 Hz SSH noise estimates'})
ylim([-0.2,0.35])
exportgraphics(gcf,['/home/fehelers/PhD_Delft/2nd_paper/20Hz_noise_and_correlation.pdf'],'Resolution',300)

%% ridiculed example
% generate some white noise and find kernel that establishes 0.1
% correlation for trailing edge
z = randn(1e6,1);
z_s = conv(z,[0.9,0.1]);

figure;
plot(z(1:100));hold on
plot(z_s(1:100));hold on

figure;
autocorr(z_s)

std(z)
std(z_s)

%%
function detrended_est = detrend_and_filter(data,var,lat_threshold, noise_threshold, movmedian_win_size)

    % apply different masks
    % latitude
    mask_lat = (data.LAT < lat_threshold)';
    data.(var)(mask_lat) = NaN;
    data.mask = mask_lat;
    
    % outlier with respect to moving median of k samples
    mask_outlier = (abs(data.(var)-movmedian(data.(var),movmedian_win_size,'omitnan')) > noise_threshold);
    data.(var)(mask_outlier) = NaN;
    data.mask = mask_outlier|data.mask|isnan(data.(var));
   
%     % substract second order polynomial function from data
%     P = polyfit(data.LAT(~data.mask),data.(var)(~data.mask),2);
%     data.(var)_slope = P(1)*data.LAT.^2 + P(2)*data.LAT+P(3);
%     data.(var)_d = data.(var) - data.(var)_slope'; 

    % substract moving median with the given window size
    detrended_est = data.(var) - movmedian(data.(var),movmedian_win_size,'omitnan');
end

function [ACF,RES,SWH] = compute_acfs(L2_files,var,num_lags,lat_threshold,noise_threshold,movmedian_win_size)
    n = numel(L2_files);

    %num_lags = 100;
    ACF = zeros(num_lags+1,n);
    RES = zeros(n,1);
    SWH = zeros(n,1);

    for i=1:numel(L2_files)
        load(fullfile(L2_files(i).folder,L2_files(i).name));
        data = CS2.f_240Hz.UF;
        data.(var) = double(data.(var));
        data.SWH = double(data.SWH);

        % ############ for HEI ###############
        data.detrended = detrend_and_filter(data,var,lat_threshold, noise_threshold, movmedian_win_size);
        
        along_track_distance = cumsum(deg2km(distance(data.LAT(1:end-1),data.LON(1:end-1),data.LAT(2:end),data.LON(2:end))));
%         figure;
%         plot(along_track_distance,data.detrended(1:end-1),'k-');hold on % detrended SSHi
%         title(i)
%         xlabel('km')

        [acf,lags] = autocorr(data.detrended,'NumLags',num_lags);
        azimuth_res = 1e3*median(deg2km(distance(data.LAT(1:end-1),data.LON(1:end-1),data.LAT(2:end),data.LON(2:end))));

    %     figure;
    %     plot( azimuth_res*lags, acf,'k-');hold on
    %     grid on
    %     if i==1
    %         xlabel('lag in m')
    %     end
        ACF(:,i) = acf;
        RES(i) = azimuth_res;
        SWH(i) = median(data.SWH(data.LAT>lat_threshold));
    end
end

function [ACF,ACF_s,ACF_r,RES] = compute_acfs_thr(L1_files,num_lags,lat_threshold,noise_threshold,movmedian_win_size) 
    n = numel(L1_files);
    movmedian_win_size = 1.4572e+04;
    num_lags = 6000/50;

    ACF = zeros(num_lags+1,n);
    ACF_s = zeros(num_lags+1,n); % sum of signal
    ACF_r = zeros(num_lags+1,n); % signal along a range bin
    RES = zeros(n,1);

    for i=1:numel(L1_files)
        i
        CS1b = load(fullfile(L1_files(i).folder,L1_files(i).name));
        CS1b = CS1b.CS1b;
        
        azimuth_res = 1e3*median(deg2km(distance(CS1b.GEO.LAT(1:end-1),CS1b.GEO.LON(1:end-1),CS1b.GEO.LAT(2:end),CS1b.GEO.LON(2:end))));
        DDAcf = DDA_ConfigFile('S3B','SAR');
        range_res = (DDAcf.c/2/DDAcf.B)/2;
        
        r = threshold_retracker(CS1b.SAR.data_pseudoDD,87,range_res);
        r = CS1b.GEO.H(:) - (CS1b.MEA.tracker_range(:) + r(:)); % to mimic HEI
        
        sig_sum = sum(CS1b.SAR.data_pseudoDD(105:160,:));
        sig_r = CS1b.SAR.data_pseudoDD(105:160,:);
        
        % lat filter
        mask_lat = (CS1b.GEO.LAT < lat_threshold);
        r(mask_lat) = NaN;
        sig_sum(mask_lat) = NaN;
        sig_r(:,mask_lat) = NaN;
        mask = mask_lat;
        mask = mask(:);
        
        clear CS1b
        
        % outlier filters
        mask_outlier = (abs(r-movmedian(r,movmedian_win_size,'omitnan')) > noise_threshold);
        r(mask_outlier) = NaN;
        
        s_movmed = movmedian(sig_sum,movmedian_win_size,'omitnan');
        s_movmad = movmad(sig_sum,movmedian_win_size,'omitnan');
        r_movmed = movmedian(sig_r,movmedian_win_size,2,'omitnan');
        r_movmad = movmad(sig_r,movmedian_win_size,2,'omitnan');
        
        mask_outlier_s = (sig_sum > s_movmed+4*s_movmad)|(sig_sum < s_movmed-4*s_movmad);
        mask_outlier_r = (sig_r > r_movmed+4*r_movmad)|(sig_r < r_movmed-4*r_movmad);
        
        sig_sum(mask_outlier_s) = NaN;
        sig_r(mask_outlier_r) = NaN;
        mask = mask_outlier|mask|isnan(r);
        
        % remove moving median to detrend
        r = r-movmedian(r,movmedian_win_size,'omitnan');
        sig_sum = sig_sum-movmedian(sig_sum,movmedian_win_size,'omitnan');
        sig_r = sig_r-movmedian(sig_r,movmedian_win_size,2,'omitnan');
        
        % remove additional outliers after the detrending
        outlier_mask_s2 = (sig_sum > 4*mad(sig_sum,1,2))|(sig_sum < -4*mad(sig_sum,1,2));
        sig_sum(outlier_mask_s2) = NaN;
        
        
        % calculate acf
        [acf,lags] = autocorr(r(25:50:end),'NumLags',num_lags);
        [acf_s,lags_s] = autocorr(sig_sum(25:50:end),'NumLags',num_lags);
        
        for n = 1:size(sig_r,1)
            if n==1
                [acf_r,lags_r] = autocorr(sig_r(n,25:50:end),'NumLags',num_lags);
            else
                [acf_r_new,lags_r] = autocorr(sig_r(n,25:50:end),'NumLags',num_lags);
                acf_r = acf_r + acf_r_new;
            end
        end
       acf_r = acf_r./size(sig_r,1);
        %[acf_r,lags_r] = 

    %     figure;
    %     plot( azimuth_res*lags, acf,'k-');hold on
    %     grid on
    %     if i==1
    %         xlabel('lag in m')
    %     end
        ACF(:,i) = acf;
        ACF_s(:,i) = acf_s;
        ACF_r(:,i) = acf_r;
        RES(i) = azimuth_res*50;
    end

end