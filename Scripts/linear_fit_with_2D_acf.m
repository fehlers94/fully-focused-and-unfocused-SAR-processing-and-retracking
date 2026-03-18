% initialise the kernel and predict resulting acfs (projected and cut)
k2 = zeros(3,3)
k2(2,2) = 0.8;
k2(3,2) = 0;
k2(1,2) = 0;
k2(:,1) = 0.2;
k2(:,3) = 0.2;
k = zeros(5,5);
k(2:4,2:4) = k2;

% plot the averaging kernel
figure;
imagesc(k);
colormap('pink')
colorbar()
%caxis([0 1])

% since kernel is symmetric, the expected acf is the 2d convolution of the kernel with itself
acf_2d = conv2(k,k,'full');
acf_2d = acf_2d./max(acf_2d(:))


% now, calculate both the projected acf (vertical sum): acf_proj
% and the cut through the middle layer of the 2d ACF: acf_x
acf_proj = sum(acf_2d);
acf_proj = acf_proj./max(acf_proj);
acf_x = acf_2d(5,:);
acf_x = acf_x./max(acf_x);


figure;
imagesc(acf_2d);
colormap('pink')
%colorbar()
%caxis([0 1])

ax1 = subplot(2,1,1)
plot(acf_proj);hold on
plot(acf_x);hold on
title('1D noise ACFs')
legend('projected','along-track cut')
ylim([-0.1,1.2])
ax2 = subplot(2,1,2)
imagesc(acf_2d)
colorbar()
xlabel('along-track')
ylabel('range')
title('2D noise ACF')
linkaxes([ax1,ax2],'x')


%% image parameters
m = 1;
b = 0;

x=1:64;
x = x-mean(x);
y = m*x+b;

%% build noisy "image" of those functions
n=10000;
y_im = repmat(y,n,1)';
y_im_true = y_im; 
imagesc(y_im_true)
colorbar()

%% random noise
noise_im = randn(size(y_im));
imagesc(noise_im)
colormap('pink')

%% convolute noise with kernel of acf
noise_im = conv2(noise_im,k,'same');
imagesc(noise_im);
colormap('pink')

%% check 2d acf of resulting noise aginst prediction
acf_test = xcorr2_fft(noise_im); % fft optimized but no standard function
%acf_test = xcorr2(noise_im); %slow
acf_test = acf_test./max(acf_test(:));

figure;
subplot(1,2,1)
imagesc(acf_2d)
colorbar()
caxis([0,1])

c = (size(acf_test)-1)/2+1;

subplot(1,2,2)
imagesc(acf_test(c(1)-4:c(1)+4,c(2)-4:c(2)+4))
colorbar()
%caxis([0,1])
colormap('pink')

%% add correlated noise to image
y_im = y_im + noise_im;
imagesc(y_im);
colorbar()
colormap('pink')

%% make linear fit on all columns independently and evaluate misfit
b_est = zeros(1,n);
m_est = zeros(1,n);
misfit = zeros(1,n);
integrated_noise = zeros(1,n);

for i=1:n
    P = polyfit(x,y_im(:,i),1);
    m_est(i) = P(1);
    b_est(i) = P(2);
    
    residual = (y_im(:,i)' - (P(1)*x+P(2))).^2;
    %integrated_noise(i) = sum((residual));
end

% plot autocorrelation of m_est
[m_acf,lags] = autocorr(m_est);
[b_acf,lags] = autocorr(b_est);
%[noise_acf,lags] = autocorr(integrated_noise);
[noise_along_track_acf,lags] = autocorr(y_im(30,:));
[summed_noise_along_track_acf,lags] = autocorr(sum(noise_im));


figure;
plot(1:10,m_acf(1:10));hold on
plot(1:10,b_acf(1:10));hold on
plot(1:10,noise_along_track_acf(1:10));hold on
plot(1:10,summed_noise_along_track_acf(1:10));hold on
plot(acf_x(5:end),'b.');hold on
plot(acf_proj(5:end),'bo');hold on


legend('ACF of m','ACF of b','ACF of image noise in along-track','ACF of summed column noise in along-track','along-track ACF','projected ACF')

% %% ############### second part with multi-column fit to y_im #####################
% 
% 
% % we have y = m*x+b;
% % and we want to solve A*params = y in a least squares sense. That means
% % the column-wise design matrix A is given by:
% A = [x;ones(size(x))]'
% % because then multiplication with params column-vector = (m b)' results in
% % the y= m * x + b;
% 
% % Now, acf_2d has a length of 5, so let's choose 5 columns at once (natural choice)
% cov = acf_2d(3:end-2,3:end-2)
% 
% A5 = repmat(A,5,1);
% 
% %% however, we want to solve this with a given correlation matrix (covariance is divided out, doesn't affect the mutal weighting)
% nn = size(y_im,1);
% e = ones(nn,1);
% e0 = zeros(nn,1);
% % build sparse sub-matrices for the problem
% S0 = spdiags(e0,0,nn,nn);
% S1 = spdiags(e*cov(:,3)',-2:2,nn,nn);
% S2 = spdiags(e*cov(:,4)',-2:2,nn,nn);
% S3 = spdiags(e*cov(:,5)',-2:2,nn,nn);
% 
% %display full matrix
% %full(S3)
% 
% % compose correlation matrix R for all elements in y5 saved as a sparse
% % matrix
% R5 = [S1 S2 S3 S0 S0;
%      S2 S1 S2 S3 S0;
%      S3 S2 S1 S2 S3;
%      S0 S3 S2 S1 S2;
%      S0 S0 S3 S2 S1];
% 
% 
% %% now compute the column-spanning estimates along the n columns
% 
% b5_est = zeros(1,n-5);
% m5_est = zeros(1,n-5);
% 
% for i=1:n-5
%     y5 = y_im(:,i+1:i+5);
%     y5 = y5(:);
% 
%     % lscov(A,y) solves A*params = y:
%     params5 = lscov(A5,y5,R);
% 
%     m5_est(i) = params5(1);
%     b5_est(i) = params5(2);
% end
% 
% % plot autocorrelation of m_est
% [m5_acf,lags] = autocorr(m5_est);
% [b5_acf,lags] = autocorr(b5_est);
% %%
% figure;
% plot(1:10,m5_acf(1:10));hold on
% plot(1:10,b5_acf(1:10));hold on
% plot(acf_x(5:end),'b.');hold on
% plot(acf_proj(5:end),'bo');hold on
% 
% legend('m acf','b acf','cut acf','proj acf')





