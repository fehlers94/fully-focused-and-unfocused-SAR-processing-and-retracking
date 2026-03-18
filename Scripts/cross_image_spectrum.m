%% cross spectrum

kx = 0.5;
ky = 0.25;
[x,y] = meshgrid(1:64,1:64)

subplot(1,2,1);imagesc(x);subplot(1,2,2);imagesc(y);

image = cos(kx.*x+ky.*y).^2.*rand(size(x)).^3;
image2 = cos(kx.*x+ky.*y+1).^2.*rand(size(x)).^3;

image = image - mean(image(:));
image2 = image2 - mean(image2(:));

figure;
imagesc(image)

figure;
imagesc(image2)

%% calculate cross image spectrum (because we have different noise realisations, it is better than autocorrelation hich is strongly peaked around the center due to the PTR!)
%% cross-covariance function is C(R) = mean( I1(r)*I2(r + R)  )
C = xcorr2_fft(image,image2);

imagesc(C)

%% cross image spectrum:
Y = fft2(C)

subplot(1,2,1)
imagesc(abs(Y))

subplot(1,2,2)
imagesc(angle(Y))








