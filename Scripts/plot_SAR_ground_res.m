clc

R = 700000;
dr = 0.47;
n = 0:128-40;

x = sqrt( (R+n*dr).^2 - R^2 );

% plot concentric circles
th = linspace(0,2*pi); 

figure;
for r=x
    xx = r*cos(th);
    yy = r*sin(th);
    plot(xx,yy,'k'); hold on
end

for i = 0:28
    plot([300*i,300*i],[-8000,8000],'k'); hold on
end

% the total width until the last bin at distance d from the center is
d = linspace(1,max(x),200); % distance from center
for i = 1:numel(d)
    n_d(i) = sum(x>d(i));
end
h = sqrt(max(x)^2 - d.^2)

% figure;
% plot(d,h)
% 
% figure;
% plot(d,n_d)
% 
% figure;
% plot(d,n_d./h)
