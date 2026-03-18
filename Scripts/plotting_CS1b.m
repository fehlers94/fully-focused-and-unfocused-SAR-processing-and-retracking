%%
% make default plot backgroudn white
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);

% choose color palette
cFF = [207, 0, 110]/255;
cDD = [0, 185, 227]/255;
cT = [214, 214, 214]/255;
%%
figure;
ax1 = subplot(5,1,1)
plot(CS2.f_020Hz.FF.LAT, CS2.f_020Hz.FF.SSHi)

ax2 = subplot(5,1,2)
plot(CS2.f_040Hz.FF.LAT, CS2.f_040Hz.FF.SSHi)

ax3 = subplot(5,1,3)
plot(CS2.f_080Hz.FF.LAT, CS2.f_080Hz.FF.SSHi)
ylabel('SSHi (m)')

ax4 = subplot(5,1,4)
plot(CS2.f_120Hz.FF.LAT, CS2.f_120Hz.FF.SSHi)

ax5 = subplot(5,1,5)
plot(CS2.f_240Hz.FF.LAT, CS2.f_240Hz.FF.SSHi)
xlabel('latitude (deg)')

linkaxes([ax1 ax2 ax3 ax4 ax5],'xy')

%% for conference, example of altimetry data
figure;

t_multilook = 1/(2*20):1/20:CS1b.GEO.Elapsed_Time(end); % 20 Hz timing
closestIndex = interp1(CS1b.GEO.Elapsed_Time, 1:numel(CS1b.GEO.Elapsed_Time), t_multilook, 'nearest', 'extrap')

ax1 = subplot(1,3,1)
imagesc(CS1b.GEO.LAT(closestIndex),1:256,log10(CS1b.SAR.data_pseudoDD(:,closestIndex)))
title('20 Hz UF-SAR (as from EUMETSAT)')
colormap('pink')
colorbar()
caxis([9,13])

ax2 = subplot(1,3,2)
imagesc(CS1b.GEO.LAT,1:256,log10(CS1b.SAR.data_pseudoDD))
title('own UF-SAR implementation')
colormap('pink')
colorbar()
caxis([9,13])

ax3 = subplot(1,3,3)
imagesc(CS1b.GEO.LAT,1:256,log10(movmean(CS1b.SAR.data,20,2)))
title('FF-SAR multilooked on 10 m')
colormap('pink')
colorbar()
caxis([9,13])
linkaxes([ax1 ax2 ax3],'xy')
%% for conference, example of altimetry data
figure;
rb=82
t_multilook = 1/(2*20):1/20:CS1b.GEO.Elapsed_Time(end); % 20 Hz timing
closestIndex = interp1(CS1b.GEO.Elapsed_Time, 1:numel(CS1b.GEO.Elapsed_Time), t_multilook, 'nearest', 'extrap')

ax1 = subplot(1,3,1)
plot(CS1b.GEO.LAT(closestIndex),log10(CS1b.SAR.data_pseudoDD(rb,closestIndex)))
title('20 Hz UF-SAR (as from EUMETSAT)')
colormap('pink')
colorbar()
caxis([8,14])

ax2 = subplot(1,3,2)
plot(CS1b.GEO.LAT,log10(CS1b.SAR.data_pseudoDD(rb,:)))
title('own UF-SAR implementation')
colormap('pink')
colorbar()
caxis([8,14])

ax3 = subplot(1,3,3)
plot(CS1b.GEO.LAT,log10(movmean(CS1b.SAR.data(rb,:),50,2)))
title('FF-SAR multilooked on 10 m')
colormap('pink')
colorbar()
caxis([8,14])
linkaxes([ax1 ax2 ax3],'xy')