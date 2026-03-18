




%%

datacase = 'f_240Hz'

figure;

plot(CS2.(datacase).FF.LAT,CS2.(datacase).FF.SSHi);hold on
plot(CS2.(datacase).UF.LAT,CS2.(datacase).UF.SSHi);hold on
%plot(CS2.(datacase).EUM.lat_20_ku,CS2.(datacase).EUM.SSHi)

legend('FF SAR','UF SAR','Eumetsat')

xlabel('latitude (deg)')
ylabel('instantaneous Sea Suface Height over ref. ell. WGS84 (m)')

%% analysis planned as a start:
%
% - comparing data to the nearest tide gauges
% - checking for biases versus EUMETSAT data: Is it constant in time?
% - 


%%
%load korea5c
%worldmap(korea5c,korea5cR);
%geoshow(korea5c,korea5cR,'DisplayType','texturemap')
%demcmap(korea5c)

latlim=[51,53];
lonlim=[0,6];
worldmap(latlim,lonlim)
load coastlines
%plotm(coastlat,coastlon);hold on
geoshow(coastlat,coastlon);hold on
plotm(CS2.f_020Hz_matched.UF.LAT,CS2.f_020Hz_matched.UF.LON)
