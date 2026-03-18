clc
clear all
%%

res = 10;
N = 7000./res;
x = res.*(-N:N);
y = res.*(-N:N);

x_p = repmat(x,2*N+1,1);
y_p = x_p';

z_s = 800e3;
V = 7000;

% plot x and y planes
figure;
subplot(1,2,1)
imagesc(x,y,x_p)
colorbar()

subplot(1,2,2)
imagesc(x,y,y_p)
colorbar()

R_xy = sqrt( x_p.^2 + y_p.^2 + z_s.^2);

figure;
subplot(1,2,1)
imagesc(x,y,R_xy)
colorbar()
subplot(1,2,2)
imagesc(x,y,R_xy.^4)
colorbar()

f_D = (x_p.*V)./sqrt( x_p.^2 + y_p.^2 + z_s.^2);

figure;
subplot(1,2,1)
imagesc(x,y,f_D)
colorbar()

figure;
subplot(1,2,1)
imagesc(x,y,x_p.^2/(2*z_s))
colorbar()

%% look at the spotlight caused by the range sinc:

DDAcf = DDA_ConfigFile('S3B','SAR');
T = 2.1;

range = z_s + 5;
term = sinc(2*DDAcf.B./DDAcf.c * (range - z_s - x_p.^2/(2*z_s) - y_p.^2/(2*z_s))  ).^2.*sinc(64*DDAcf.PRI*2*DDAcf.fc.*(x_p-3000).*V./DDAcf.c./z_s).^2;

figure; 
imagesc(x,y,term);
colorbar()



