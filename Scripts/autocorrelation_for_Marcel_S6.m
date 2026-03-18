%% ############### Sentinel-6! #################
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);

% choose color palette
cFF = [207, 0, 110]/255;
cDD = [0, 185, 227]/255;
cT = [214, 214, 214]/255;
%%

load('/home/fehelers/PhD_Delft/pseudoDDprocessing/S6test_FF_with_zD.mat')
%load('/home/fehelers/PhD_Delft/pseudoDDprocessing/swh1.5m_FF.mat')
CS1b = ffproc.CS1b
clear ffproc

%%
idxmax = 40000;
snippet = CS1b.SAR.data(200:400,1:idxmax);
%snippet = CS1b.SAR.data(100:180,1:idxmax);
figure;
imagesc(CS1b.SAR.data(:,1:idxmax));
title("FFSAR radargram R")
colormap('pink')
colorbar()

snippet = snippet./mean(snippet,1);

snippet_m = movmean(snippet,6000,2);
figure;
imagesc(snippet_m);
colormap('pink')
title("moving mean M")
colorbar()

snippet_norm = snippet-snippet_m;
figure;
imagesc(snippet_norm);
colormap('pink')
title("residual R - M")
colorbar()

%% calculate autocorrelation along these bins
lag = 7000;
Cptr = zeros(200,2*lag+1);
for i = 1:200
    Cptr(i,:) = xcorr(snippet_norm(i,:),snippet_norm(i,:),lag);
    Cptr(i,:) = Cptr(i,:)./max(Cptr(i,:));
end


%%
figure;
plot(mean(Cptr,1,'omitnan'))

%%
azimuth_res = 1e3*deg2km(distance(CS1b.GEO.LAT(1),CS1b.GEO.LON(1),CS1b.GEO.LAT(end),CS1b.GEO.LON(end)))/numel(CS1b.GEO.LON);

% ax1 = subplot(2,1,1)
% plot(azimuth_res*(-lag:lag),mean(Cptr(1:200,:)));hold on
% %plot(azimuth_res*(-lag:lag),sinc((-lag:lag)/0.7).^2);hold on
% %plot(azimuth_res*(-lag:lag),sinc((-lag:lag)./azimuth_res).^2);hold on
% plot(azimuth_res*(-lag:0.1:lag),sinc((-lag:0.1:lag)./azimuth_res).^2);hold on
% legend('autocorrelation','expected PTR shape: sinc^2(x/L_x), L_x = 1.07 m')
% xlabel('along-track lag (m)')
% title({['Sentinel-6 MF: FF-SAR radargram'] ['along-track autocorrelation function'] ['SWH ~ 1.8 m']})
% grid on
% xlim([-1000 1000])

%ax2 = subplot(2,1,2)
plot(azimuth_res*(-lag:lag),mean(Cptr(1:200,:)));hold on
%plot(azimuth_res*(-lag:lag),sinc((-lag:lag)/0.7).^2);hold on
%plot(azimuth_res*(-lag:lag),sinc((-lag:lag)./azimuth_res).^2);hold on
plot(azimuth_res*(-lag:0.1:lag),sinc((-lag:0.1:lag)./azimuth_res).^2);hold on
legend('autocorrelation','expected PTR shape: sinc^2(x/L_x), L_x = 1.07 m')
% title({['Sentinel-6 MF: FF-SAR radargram'] ['along-track autocorrelation function'] ['SWH ~ 1.8 m']})
xlabel('along-track lag (m)')
title({['Sentinel-6 MF: FF-SAR radargram'] ['along-track autocorrelation function']})
grid on

xlim([-1000 1000])
ylim([-0.005 0.01])
% ax2 = subplot(6,1,2)
% plot((-lag:lag),mean(Cptr(50:100,:)));hold on
% plot((-lag:lag),sinc((-lag:lag)/0.7).^2);hold on
% plot((-lag:lag),sinc((-lag:lag)/1.3).^2);hold on
% 
% ax3 = subplot(6,1,3)
% plot((-lag:lag),mean(Cptr(100:150,:)));hold on
% plot((-lag:lag),sinc((-lag:lag)/0.7).^2);hold on
% plot((-lag:lag),sinc((-lag:lag)/1.3).^2);hold on
% 
% ax4 = subplot(6,1,4)
% plot((-lag:lag),mean(Cptr(150:200,:)));hold on
% plot((-lag:lag),sinc((-lag:lag)/0.7).^2);hold on
% plot((-lag:lag),sinc((-lag:lag)/1.3).^2);hold on
% 
% ax5 = subplot(6,1,5)
% plot((-lag:lag),mean(Cptr(200:250,:)));hold on
% plot((-lag:lag),sinc((-lag:lag)/0.7).^2);hold on
% plot((-lag:lag),sinc((-lag:lag)/1.3).^2);hold on
% 
% ax6 = subplot(6,1,6)
% plot((-lag:lag),mean(Cptr(250:300,:)));hold on
% plot((-lag:lag),sinc((-lag:lag)/0.7).^2);hold on
% plot((-lag:lag),sinc((-lag:lag)/1.3).^2);hold on
% 
% linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'xy')
xlim([-1000 1000])
ylim([-0.005 0.01])