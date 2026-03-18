3e9
13.575e9*0.022

2.355*800000./(2*7000)*2*1




%% ######################################
% random waves with total amplitude and what velocity variance
% deep water:
% eta = a sin(wt-kx)
% ux = wa cos(wt-kx)
% w^2 = gk

% consider a snapshot at a given time t = 0

x = linspace(0,100,4000);

a = 0.1;
k = 2*pi/0.4;
w = sqrt(9.81*k);

eta = a*sin(-k*x);
ux = w*a*cos(-k*x);

a = 0.2;
k = 2*pi/1;
w = sqrt(9.81*k);

eta = eta + a*sin(-k*x);
ux = ux + w*a*cos(-k*x);

a = 0.3;
k = 2*pi/1.5;
w = sqrt(9.81*k);

eta = eta + a*sin(-k*x);
ux = ux + w*a*cos(-k*x);




4*std(eta)
std(ux)





