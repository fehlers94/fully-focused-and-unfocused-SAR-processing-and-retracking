% this script is supposed as a test-bed for the flat-earth approximation
% validity in case of off-nadir altimetry
% for plot of the geometry: see manuscript

clear all;
RE = 6371000;   % Earth radius
Hs = 1336000;   % satellite altitude

% parameters:
% input river height over round earth
Hr = -50:10:50;
% maximum cross track distance
cross_track_dist = (0:1000:10000)'; 

% calculate range R to satellite assuming input parameters

phi_r = cross_track_dist/RE; % angle towards river in radians
xr = (Hr+RE).*sin(phi_r);
yr = (Hr+RE).*cos(phi_r) - RE;

R = sqrt( (Hs-yr).^2 + (xr).^2 );


% ########## reverse the range measurement to a river height using (xr0,
% yr0) = (xr(Hr=0), yr(Hr=0)) and enu coordinates
xr0 = (RE)*sin(phi_r);
yr0 = (RE)*cos(phi_r) - RE;

d = xr0;
yr0;

yr_est = (Hs + xr0.^2/(2*Hs) - R);
Hr_est = yr_est - yr0;
%Hr_est = ((yr_est+RE)./cos(xr0/RE)-RE);

leg = {}
sym = {'.','o','x','+','*','s','d','v','^','<','>','p','h'}
for i = 1:numel(cross_track_dist)
    ax = plot(Hr,1000*(Hr_est(i,:)-Hr),['-' sym{i}]);hold on
    leg{i} = ['x_R=' num2str(cross_track_dist(i)/1000) ' km']
end

legend(leg,'NumColumns',2)
xlabel('true WSE (m)')
ylabel({'systematic difference between obtained','WSE and true WSE (mm)'})
grid on
set(gcf,'color','w');

%exportgraphics(gcf,'systematic_bias.pdf') 
