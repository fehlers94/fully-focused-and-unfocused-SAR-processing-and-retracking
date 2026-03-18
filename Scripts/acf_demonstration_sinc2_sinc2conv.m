clear all

%% preliminaries

% time spacing
dt = 0.1;
% width of sinc transfer function sinc(x/L)
L = 4;

% ### generate time series with sinc^2 acf and approx. gaussian distribution ####
%L = 1;
N=16000000;
%x = randn(1,N);
%y = conv(x,ones(1,4),'same');
y = zeros(1,N);
for i = 1%250
    %x = rand(1,N).^1.*exp(1j*2*pi*rand(1,N));%.^31;%.^200;
    %x = rand(1,N).^31;%.^31;%.^200;
    %x = randn(1,N) + 1j*randn(1,N);
    %x_1d = rand(1,N).^31;
    %x = rand(1,N).^31 + 1j*rand(1,N).^31;
    x_1d = 2*randn(1,N).^3;
    x = 2*randn(1,N).^3 + 2*1j*randn(1,N).^3;
    
    %x = exprnd(10,1,N).^2;
    x = x-mean(x);
    x_1d = x_1d-mean(x_1d);
    %y = y + sign(rand(1,1)-0.5)*(conv(x,sinc(-100:0.25:100),'same').^2);
    y = abs(conv(x,sinc((-100:dt:100)/L),'same')).^2;
    %y = y + sign(rand(1,1)-0.5)*(conv(x,sinc(-100:0.25:100).^2,'same')); %
    %convolution with a sinc square leads to autocorrelation that does not
    %match! ALso doesn't resemble reality! In fact, the sinc is the
    %response and the power detection comes after the rest
    %y = y + sign(rand(1,1)-0.5)*(conv(x,sinc(-2:0.25:2),'same').^2);
end

y = (y-mean(y));
% this summation makes y approx. gaussian

[acf,lags] = autocorr(y,round(15/dt*L/4));

% ### plot acf and histogram for a quick test ###
figure;
subplot(1,3,1)
plot(y(1:1000))
ylabel('noise value')
title('synthetic noise')

subplot(1,3,2)
hist(y,100);
xlabel('noise value')
title('histogram')

subplot(1,3,3)
plot(dt*lags,acf,'.'); hold on
plot(dt*lags,sinc(lags*dt/L).^2,'k')
%plot(dt*lags,sinc(lags*dt/L))
sinc2conv = conv(sinc((-100:dt:100)/L).^2,sinc((-100:dt:100)/L).^2,'same');
sinc2conv = sinc2conv./max(sinc2conv);
plot(dt*lags,sinc2conv((numel(sinc2conv)-1)/2+1:(numel(sinc2conv)-1)/2+1+numel(lags)-1),'k');

% residual = acf-sinc2conv((numel(sinc2conv)-1)/2+1:(numel(sinc2conv)-1)/2+1+numel(lags)-1);
% plot(residual);

title('Autocorrelation')
%legend('noise','sinc^2','sinc')


% calculate criteria that determines the acf shape
%L = 1;
K = 2*pi/L;
weight1 = 2*(mean(x_1d.^4)-3*mean(x_1d.^2).^2)*(4/6)*L/dt; % extra factor of 1/4 seemed to solve the issue!!! Where does it come from? L and L^2?
weight2 = 4*(mean(x_1d.^2)^2*(L/dt)^2);

% model the acf with the obtained formula
t = lags*dt;
formula1 = (K*t - sin(K*t))./(K*t).^3;
formula1(1) = 1/6;
formula1 = formula1*6; %normalised to 1 at 0


formula2 = sinc(t/L).^2; %normalised to 1 at 0

% figure;
% plot(t,6*formula1);hold on
% plot(t,formula2);hold on

%proper addition of terms for normalization
acf_pred = (weight1*formula1 + weight2*formula2)./(weight1+weight2);
%acf_pred = acf_pred./max(abs(acf_pred));
plot(dt*lags,acf_pred,'o');hold on


% evaluate the difference between acf and acf_pred depending on the weights w1/(w1+w2) and w2/(w1+w2):

w1 = -1:0.001:1;
w2 = 1-w1;

acf_pred_w1 = w1'*formula1 + w2'*formula2;

res = sum((acf - acf_pred_w1).^2,2);

[m,arg] = min(res);
w1_fit = w1(arg)
w1_calc = weight1./(weight1+weight2)

%weight_ratio_1_2 = 1/(1-w1_fit)
%weight1/weight2

plot(dt*lags,w1(arg)*formula1+w2(arg)*formula2)% best acf function fit

legend('noise','','','calc','fit')

