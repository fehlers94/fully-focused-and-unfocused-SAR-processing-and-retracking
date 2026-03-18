% this script uses a 100 km snippte to illustrate the problem with
% correlations in the 20 Hz along-track averages of real data

clc
clear all
close all
%%
dataset = 'S6test_FF'
mission = 'S6A'
LoadCommonSettings
DDAcf   = DDA_ConfigFile(mission,'SAR');
%[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
DOM = [55 57; -180 180];

%%
load(['/home/fehelers/PhD_Delft/pseudoDDprocessing/' dataset '.mat'],'ffproc')
CS1b = ffproc.CS1b;
clear ffproc

%% average onto ~80 Hz
% average on 80 Hz
N_avg = round((1/80)/median(diff(CS1b.GEO.Elapsed_Time)))-1;
CS1b_80 = FF_SAR_Average_Waveforms(CS1b,N_avg,mission);

% %% average on ~20 Hz
% N_avg = round((1/20)/median(diff(CS1b.GEO.Elapsed_Time)));
% CS1b_20 = FF_SAR_Average_Waveforms(CS1b,N_avg,mission);

%% get some things straight before retracking
CS1b_80.SAR.data_FF = CS1b_80.SAR.data;
CS1b_80.MEA.ref_range = CS1b_80.MEA.tracker_range;
CS1b_80.SAR.scale_factor_ku = CS1b_80.SAR.echo_scale_power;

% CS1b_20.SAR.data_FF = CS1b_20.SAR.data;
% CS1b_20.MEA.ref_range = CS1b_20.MEA.tracker_range;
% CS1b_20.SAR.scale_factor_ku = CS1b_20.SAR.echo_scale_power;

%% do the retracking

% FFSAR with SAMOSA zero Doppler beam
CS1b_80.SAR.data = CS1b_80.SAR.data_FF
%[L2_FF,~] = SAR_L1b_to_L2(mission,CS1b_avg,DOM,'SAMOSA2',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'})
[L2_FF,~] = SAR_L1b_to_L2(mission,CS1b_80,DOM,'SAMOSA2FF',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'});
%%
% UFSAR
CS1b_80.SAR.data = CS1b_80.SAR.data_pseudoDD
[L2_UF,~] = SAR_L1b_to_L2(mission,CS1b_80,DOM,'SAMOSA2',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'})

%% UFSAR threshold retracker as well
% CS1b_80.SAR.data = CS1b_80.SAR.data_pseudoDD
% [L2_UF,~] = SAR_L1b_to_L2(mission,CS1b_80,DOM,'SAMOSA2',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'})


%% check for outliers
figure;
plot(L2_FF.SSHi);hold on;plot(L2_UF.SSHi)

figure;
plot(L2_FF.SWH);hold on;plot(L2_UF.SWH)

L2_FF.SSHi(817) = L2_UF.SSHi(817);
L2_FF.SWH(817) = L2_UF.SWH(817);

%% detrend SSHi and SWH using 1 Hz moving average of uncorrelated FF-SAR estimates
k = 80;

L2_FF.SSHi_d = L2_FF.SSHi - movmean(L2_FF.SSHi,k);
%L2_FF.SSHi_d(1:k/2) = [];
%L2_FF.SSHi_d(end-k/2:end) = [];

L2_UF.SSHi_d = L2_UF.SSHi - movmean(L2_FF.SSHi,k);
%L2_UF.SSHi_d(1:k/2) = [];
%L2_UF.SSHi_d(end-k/2:end) = [];

L2_FF.SWH_d = L2_FF.SWH - movmean(L2_FF.SWH,k);
%L2_FF.SWH_d(1:k/2) = [];
%L2_FF.SWH_d(end-k/2:end) = [];

L2_UF.SWH_d = L2_UF.SWH - movmean(L2_FF.SWH,k);
%L2_UF.SWH_d(1:k/2) = [];
%L2_UF.SWH_d(end-k/2:end) = [];

%% plot the time series
figure;
plot(L2_FF.SSHi_d);hold on
plot(L2_UF.SSHi_d)

figure;
plot(L2_FF.SWH_d);hold on
plot(L2_UF.SWH_d)

%% compute and plot acfs for SSH and SWH
figure;
subplot(1,2,1)
[acf,lags] = autocorr(L2_FF.SSHi_d,100);
plot(lags,acf,'b.');hold on

[acf,lags] = autocorr(L2_UF.SSHi_d,100);
plot(lags,acf,'ro');hold on
grid on

% plot analytic sinc (DO MORE EXACT SOON!)
lags = 0:0.025:2;
plot(4*lags, sinc(lags).^2,'k-')

%plot 1/sqrt(N) error bars
plot([0,20],1.96/sqrt(numel(L2_FF.SSHi_d))*ones(2,1));hold on
plot([0,20],-1.96/sqrt(numel(L2_FF.SSHi_d))*ones(2,1));hold on


legend('FF-SAR','UF-SAR','sinc^2')
title('uncorrected SSH')

subplot(1,2,2)
[acf,lags] = autocorr(L2_FF.SWH_d,100);
plot(lags,acf,'b.');hold on

[acf,lags] = autocorr(L2_UF.SWH_d,100);
plot(lags,acf,'ro');hold on

% plot analytic sinc (DO MORE EXACT SOON!)
lags = 0:0.025:2;
plot(4*lags, sinc(lags).^2,'k-')

%plot 1/sqrt(N) error bars
plot([0,20],1.96/sqrt(numel(L2_FF.SWH_d))*ones(2,1));hold on
plot([0,20],-1.96/sqrt(numel(L2_FF.SWH_d))*ones(2,1));hold on


grid on
legend('FF-SAR','UF-SAR','sinc^2')
title('SWH')

%% compute some 20 Hz estimates and evaluate acf "precision"
% UF
L2_UF_20.SSHi_d_a = movmean(L2_UF.SSHi_d,4);
L2_UF_20.SSHi_d_a = L2_UF_20.SSHi_d_a(1:4:end); %averaged

L2_UF_20.SSHi_d_s = L2_UF.SSHi_d(1:4:end); %subsampled

L2_UF_20.SWH_d_a = movmean(L2_UF.SWH_d,4);
L2_UF_20.SWH_d_a = L2_UF_20.SWH_d_a(1:4:end);

L2_UF_20.SWH_d_s = L2_UF.SWH_d(1:4:end);

% FF
L2_FF_20.SSHi_d_a = movmean(L2_FF.SSHi_d,4);
L2_FF_20.SSHi_d_a = L2_FF_20.SSHi_d_a(1:4:end); %averaged

L2_FF_20.SSHi_d_s = L2_FF.SSHi_d(1:4:end); %subsampled

L2_FF_20.SWH_d_a = movmean(L2_FF.SWH_d,4);
L2_FF_20.SWH_d_a = L2_FF_20.SWH_d_a(1:4:end);

L2_FF_20.SWH_d_s = L2_FF.SWH_d(1:4:end);

%% compute fair comparison of precision

std(L2_FF_20.SSHi_d_a)
std(L2_UF_20.SSHi_d_s)
std(L2_UF_20.SSHi_d_a)
%%
std(L2_FF_20.SWH_d_a)
std(L2_UF_20.SWH_d_s)
std(L2_UF_20.SWH_d_a)

%% 20 Hz acfs

figure;
subplot(1,2,1)
[acf,lags] = autocorr(L2_FF_20.SSHi_d_a,10);
plot(lags,acf,'b.');hold on

[acf,lags] = autocorr(L2_UF_20.SSHi_d_a,10);
plot(lags,acf,'ro');hold on

[acf,lags] = autocorr(L2_UF_20.SSHi_d_s,10);
plot(lags,acf,'go');hold on

%plot 1/sqrt(N) error bars
plot([0,10],1.96/sqrt(numel(L2_FF_20.SSHi_d_a))*ones(2,1));hold on
plot([0,10],-1.96/sqrt(numel(L2_FF_20.SSHi_d_a))*ones(2,1));hold on

grid on

legend('av. FF-SAR','av. UF-SAR','subs. UF-SAR')

title('uncorrected SSH')

subplot(1,2,2)
[acf,lags] = autocorr(L2_FF_20.SWH_d_a,10);
plot(lags,acf,'b.');hold on

[acf,lags] = autocorr(L2_UF_20.SWH_d_a,10);
plot(lags,acf,'ro');hold on

[acf,lags] = autocorr(L2_UF_20.SWH_d_s,10);
plot(lags,acf,'go');hold on

%plot 1/sqrt(N) error bars
plot([0,10],1.96/sqrt(numel(L2_FF_20.SWH_d_a))*ones(2,1));hold on
plot([0,10],-1.96/sqrt(numel(L2_FF_20.SWH_d_a))*ones(2,1));hold on

grid on

legend('av. FF-SAR','av. UF-SAR','subs. UF-SAR')

title('SWH')


%% white noise acf errors

z = (2*rand(2500,1)-1).^3;
hist(z.^3,20)

figure;
[acf,lags] = autocorr(z,20);
plot(lags,acf);hold on
plot([0,20],1.96/sqrt(numel(z))*ones(2,1));hold on
plot([0,20],-1.96/sqrt(numel(z))*ones(2,1));hold on

%% try to bootstrap some uncertainty assuming acf is sinc^2? Doesn't work, dependent on real signal inside and averaging beforehand
k
% given sample size and knowing the true acf, compute the median and
% percentiles from 1000 tries
figure;
for j=1:200

    N=2500;
    z = zeros(1,N);
    for i = 1:250
        x = randn(1,N);
        z = z + sign(rand(1,1)-0.5)*(conv(x,sinc(-20:0.25:20),'same').^2);
    end

    z = (z-mean(z))./std(z);
    % add some smooth 'signal' and remove moving mean
    z_signal = conv(movmean(x,1000),sinc(-20:0.01:20),'same');
    z = z + z_signal;
    z = z - movmean(z,120);
    
    

    [acf,lags] = autocorr(z,100);
    plot(lags,acf,'r');hold on
end

plot(lags,sinc(0.25*lags).^2,'k')

%% use small snippets to evaluate acf and average them to get an error?

y = L2_UF.SSHi_d;
k=80
n = round(numel(y)/80)-1
acf_a = zeros(21,n);

for j = 0:n-1
    j
    [acf,lags] = autocorr(y(j*k+1:(j+1)*k),20);
    acf_a(:,j+1) = acf;
end

figure;
errorbar(mean(acf_a,2),1.96*std(acf_a,[],2)./sqrt(n));hold on;plot(sinc(0:0.25:2).^2)