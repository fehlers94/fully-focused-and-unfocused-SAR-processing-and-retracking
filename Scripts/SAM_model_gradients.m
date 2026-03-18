clc; clear all;
close all;
%%
% make default plot backgroudn white
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);

CS1b = load('/home/fehelers/PhD Delft/Projects/FFSAR4ROFI/Data/ROFI/L1b/S3A_SR_1_SRA_A__20180109T101705_20180109T110735_20180204T003249_3029_026_279______MAR_O_NT_003.SEN3.mat')
CS1b=CS1b.CS1b;

%% evaluate along-track correlation function of the resulting SAMOSA SSH estimate:

% calculate covariance model:

mission = 'S3B'
mode = 'all'
T = CS1b.GEO.integration_time
DDAcf = DDA_ConfigFile('S3B','SAR');

%% acf model:

H = median(CS1b.MEA.tracker_range);
V = median(CS1b.GEO.V.V);
dHdt = median(CS1b.GEO.H_rate);

% DDAcf = DDA_ConfigFile('S3A','SAR');
% H = 800000;         % m % or better take tracker range instead, because we want the range history of the surface.
% V = 7000;           % m/s
% dHdt = -10;         % m/s
% T = 2;

multi_target = false;
num_targets = 21; % number of point targets in the center range bin, if multi_target=true; this is important to assume, because given a signal value in a range bin, we cannot know where the target is; take something odd for symmetrical filling of rangebin;

%os_zp = 2; %oversampling of the waveform
n_bins = DDAcf.os_ZP*DDAcf.Np;

% calculate exact range history in flat earth approximation and compare it to its Taylor approximation
R0 = H;
%t = -T/2:0.02:T/2;
if strcmp(mode,'all')
    %t = (1:401)/401*T - 401/400*T/2;
    t = (1:181)/181*T - 181/180*T/2;
elseif strcmp(mode,'left')
    t = (1:401)/401*T - 401/400*T/2;
    t = t(1 : (size(t,2)+1)./2);
elseif strcmp(mode,'right')
    t = (1:401)/401*T - 401/400*T/2;
    t = t((size(t,2)+1)./2 : end);
end
alpha = asin(dHdt/V);

Rt = sqrt(  (H + V*t*sin(alpha)).^2 + (V*t*cos(alpha)).^2 );
Rt2 = sqrt(  H.^2 + 2*H*V*t*sin(alpha) + (V*t).^2 ); 
Rt3 = H.*sqrt(  1 + 2/H*V*t*sin(alpha) + (V*t/H).^2 ); 
Rt4 = H.*(1 + 0.5*(2/H*V*t*sin(alpha) + (V*t/H).^2));
Rt4 = H.*(1 + 0.5*(2/H*V*t*sin(alpha) + (V*t/H).^2));
Rt5 = H + V*t*sin(alpha) + 0.5*(V*t).^2/H; % take this as model for range history
% magnitude analysis
% H.^2 = 6.4000e+11
% 2*H*V*t*sin(alpha) = 1.6e+7
% (V*t).^2 = 4.9e+7


% % plotting Rt, to assess goodness of approximation
% figure;
% subplot(1,2,1)
% plot(t,Rt);hold on
% plot(t,Rt2,'ro');hold on
% plot(t,Rt3,'r.');hold on
% plot(t,Rt4,'bx');hold on
% plot(t,Rt5,'bo');hold on
% subplot(1,2,2)
% plot(t,1000*(Rt - Rt4),'bx');hold on % exact to a milimeter
% plot(t,1000*(Rt5 - Rt4),'b.');hold on % exact to a milimeter
% xlabel('time')
% ylabel('disagrement of approximation with "truth" in mm')

% %% now, Rt5 = H + V*t*sin(alpha) + 0.5*(V*t).^2/H is range history to target in good approximation and the difference between two adjacent range histories is then given by a simple linear function
% x = cos(alpha)*V*t;             % ground projected distance
% theta_p = atan(x./(H-dHdt*t));  % pitch angles to focus point in flat earth approximation, used for antenna pattern
% %% check error due to small angle approximation:
% x_s = V*t;             % ground projected distance
% theta_p_s = x./(H-dHdt*t);  % pitch angles to focus point in flat earth approximation, used for antenna pattern
% figure;
% plot(theta_p); hold on
% plot(theta_p_s)
% figure;
% plot(100*(theta_p-theta_p_s)./theta_p);hold on % relative error of 3e-3% very acceptable
% %% illustrate calculated vr (radial velocity from target at time t)
% vr = V^2*t./(H+dHdt*t) - dHdt;
% vr_s = V^2*t./H - dHdt;
% figure;
% plot(vr); hold on
% plot(vr_s)
%
gamma_x = 8*log(2)/DDAcf.theta_x^2; 
gamma_y = 8*log(2)/DDAcf.theta_y^2; %are the same

%antenna_pattern = exp(-gamma_x*theta_p.^2); define antenna pattern in the
%loop depending on changing theta_p

center = n_bins/2; %center bin
range_res_sample = 1/DDAcf.B*DDAcf.c/2/DDAcf.os_ZP; %oversampled range resolution
if strcmp(mission,'S6A')
    range_res_transmit = 1/DDAcf.Bt*DDAcf.c/2; %actual range resolution and scale of the sinc function
else
    range_res_transmit = 1/DDAcf.B*DDAcf.c/2; %actual range resolution and scale of the sinc function
end

rangebin = (1:n_bins); % floating point range bin number
d_res = 1;
N = 1200;
% so that the sidelobes from +- 1 km distance will be included
IRF = zeros(n_bins,2*N+1);
IRF_new = zeros(n_bins,2*N+1);
M = size(t,2);
M_c = round(size(t,2)/2);

rel_target_pos(1,1,:) = 1/num_targets*((1:num_targets)-mean(1:num_targets)); % relative position of several point targets within the center range bin

% calculate IRFs new
x = d_res*(-N:N); % ground-projected along-track distance of the scatterer to nadir
t_p = t;  % determined/given as some inform time vector above already
r = range_res_sample*((1:n_bins)-round(n_bins/2)); % positions the scatterer exactly in the middle of a bin

% make 3D meshgrid of required coordinates (screws memory)
% [x,r,t_p] = meshgrid(   d_res*(-N:N),...
%                         range_res_sample*((1:n_bins)-round(n_bins/2)),...
%                         t)
% calculate IRF (easy but slower)
% 
% tic
% for i = 1:n_bins
%     IRF_new(i,:) =      sum(...
%                     exp(-gamma_x*((V*t_p'+x)./H).^2).*...
%                     sinc((r(i) - V*t_p'.*x/H - x.^2/(2*H) + DDAcf.fc*V*x./(H*DDAcf.s)) ./ range_res_transmit).^2,...
%                     1);
%     i
% end            
% toc
% %% plot alignment errors
% 
% % make default plot backgroudn white
% get(0,'Factory');
% set(0,'defaultfigurecolor',[1 1 1]);
% 
% % choose color palette
% cFF = [207, 0, 110]/255;
% cDD = [0, 185, 227]/255;
% cT = [214, 214, 214]/255;
% figure('units','inch','position',[0,0,4,5]);
% 
% 
% ax4 = subplot(2,4,8);
% 
% index = [0,372,2*372];
% for j = [1,2,3];
%     i = index(j);
%     d = d_res*i % distance from focal point [m]
%     dt = d / (V*cos(alpha));
%     %Rt = H + V*t*sin(alpha) + 0.5*(V*t).^2/H;
%     %Rt_off = H + V*(t+dt)*sin(alpha) + 0.5*(V*(t+dt)).^2/H;
%     %D_rcmc = (Rt-Rt(M_c))-(Rt_off-Rt_off(M_c)); % relative error in the rcmc at the ground distance d with respect to the grating lobe target from distance d
%                         % R_rcmc is approximately linear
%     %D_rcmc = -V^2*t*dt/H - V^2*dt^2/(2*H) - dt*dHdt + DDAcf.fc*V^2*dt./(H*DDAcf.s); % replace with analytic formulation from calculations
%     %D_rcmc = -V^2*t*dt/H - V^2*dt^2/(2*H) + DDAcf.fc*V^2*dt./(H*DDAcf.s); % replace with analytic formulation from calculations and dropping one term
%     D_rcmc = V^2*t*dt/H + V^2*dt^2/(2*H) - DDAcf.fc*V^2*dt./(H*DDAcf.s); % replace with analytic formulation from calculations and dropping one term
% %     figure;
% %     plot(t,(D_rcmc)/0.47); hold on
% %     plot(t,(D_rcmc_analytic)/0.47,'r.'); hold on
% %     figure;
% %     plot(t,antenna_pattern)
% 
%     % get sum of sinc^2's blurred at distance d from the focal point (approximate convolution coarsly by just adding up shifted functions)
%     D_rangebin = D_rcmc/range_res_sample;
%     kernel = zeros(size(rangebin));
% 
%     
%     %sinc_scale = DDAcf.Bt*2/DDAcf.c*range_res_sample;
%     sinc_scale = range_res_sample/range_res_transmit;
%     %sinc_scale = DDAcf.Bt*DDAcf.c/2*range_res;
% 
%     
%     % rewrite this into a matrix operation!
%     theta_p = -V*(t+dt)./H; %modulo a plus/minus error, which should not matter due to squaring
%     antenna_pattern = exp(-gamma_x*theta_p.^2);
%     
%     kernel = kernel + antenna_pattern'.*sinc(sinc_scale*(rangebin - D_rangebin' - center)).^2;
%     ax = subplot(2,4,4+j)
%     img_axj(j) = ax;
%     kernel = kernel'./max(kernel(:));
%     imagesc([-T/2,T/2],1:256,10*log10(kernel));
%     colormap('pink')
%     %colorbar()
%     caxis([-30 5])
%     
%     axes(ax4)
%     plot(sum(kernel',1)./sum(sum(kernel',1)),0.25*(1:128*8),'k'); hold on
%     grid on
%     
%     ax = subplot(4,4,4+j)
%     ant_axj(j) = ax;
%     plot(T*(-200:200)./400,sum(kernel)./max(sum(kernel)),'k');
%     grid on
% end
% 
% linkaxes([ant_axj img_axj],'x')
% linkaxes([ant_axj],'y')
% linkaxes([img_axj ax4],'y')
% img_axj(1).YLim = [128-30,128+30]
% ax4.XLim = [-0.01,0.13]
% 
% img_axj(2).XLabel.String = 'relative pulse timing (s)'
% img_axj(1).YLabel.String = 'range bin index'
% ant_axj(2).Title.String = 'along-track antenna pattern'
% 
% img_axj(1).Title.String = ['x=' num2str(d_res*index(1)) ' m']
% img_axj(2).Title.String = ['x=' num2str(d_res*index(2)) ' m']
% img_axj(3).Title.String = ['x=' num2str(d_res*index(3)) ' m']
% 
% ax4.YTickLabel = []
% img_axj(2).YTickLabel = []
% img_axj(3).YTickLabel = []
% 
% %exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/paper_range_alignment_' mission '_model.pdf'],'Resolution',300)
% 

% calculate IRFs (looks complicated but faster)
tic
for i = -N:N
    d = d_res*i; % distance from focal point [m]
    dt = d / (V*cos(alpha));
    %Rt = H + V*t*sin(alpha) + 0.5*(V*t).^2/H;
    %Rt_off = H + V*(t+dt)*sin(alpha) + 0.5*(V*(t+dt)).^2/H;
    %D_rcmc = (Rt-Rt(M_c))-(Rt_off-Rt_off(M_c)); % relative error in the rcmc at the ground distance d with respect to the grating lobe target from distance d
                        % R_rcmc is approximately linear
    %D_rcmc = -V^2*t*dt/H - V^2*dt^2/(2*H) - dt*dHdt + DDAcf.fc*V^2*dt./(H*DDAcf.s); % replace with analytic formulation from calculations
    %D_rcmc = -V^2*t*dt/H - V^2*dt^2/(2*H) + DDAcf.fc*V^2*dt./(H*DDAcf.s); % replace with analytic formulation from calculations and dropping one term
    D_rcmc = V^2*t*dt/H + V^2*dt^2/(2*H) - DDAcf.fc*V^2*dt./(H*DDAcf.s); % replace with analytic formulation from calculations and dropping one term
%     figure;
%     plot(t,(D_rcmc)/0.47); hold on
%     plot(t,(D_rcmc_analytic)/0.47,'r.'); hold on
%     figure;
%     plot(t,antenna_pattern)

    % get sum of sinc^2's blurred at distance d from the focal point (approximate convolution coarsly by just adding up shifted functions)
    D_rangebin = D_rcmc/range_res_sample;
    kernel = zeros(size(rangebin));


    %sinc_scale = DDAcf.Bt*2/DDAcf.c*range_res_sample;
    sinc_scale = range_res_sample/range_res_transmit;
    %sinc_scale = DDAcf.Bt*DDAcf.c/2*range_res;


    % rewrite this into a matrix operation!
    theta_p = -V*(t+dt)./H; %modulo a plus/minus error, which should not matter due to squaring
    %theta_p = -0.8899*V*(t+dt)./H; % when comparing to data to correct for look angle inaccuracy
    %antenna_pattern = exp(-gamma_x*theta_p.^2);
    antenna_pattern = exp(-gamma_x*theta_p.^2).^2; % this needs to be squared in case of the ACF

    if ~multi_target %place single point scatterer in the center bin
        kernel = kernel + antenna_pattern*sinc(sinc_scale*(rangebin - D_rangebin' - center)).^2;
    else
%         for j = 1:num_targets %place multiple point scatterers in the center bin
%             kernel = kernel + antenna_pattern*sinc(sinc_scale*(rangebin - D_rangebin' - center - rel_target_pos(j))).^2;
%         end
        kernel = kernel + antenna_pattern*sum(sinc(sinc_scale*(rangebin - D_rangebin' - center - rel_target_pos)).^2,3);
    end

    IRF(:,i+N+1) = kernel';%/sum(kernel);
    %IRF(:,N+i+1) = kernel'/sum(kernel);
    %IRF(:,N-i+1) = kernel'/sum(kernel);
    %figure;plot(log10(kernel/max(kernel)));hold on
    %plot(log10(IRF(:,i)))
    %clear kernel
    %plot(log10(IRF(:,i)))
    i
end
toc

% DD_scale = (DDAcf.Nb*DDAcf.PRI)*2*DDAcf.fc*V/(DDAcf.c*H); %now, multiplied with d (y_a in Egido), we get along track the sinc argument
% IRF_DD = IRF.*sinc(DD_scale*d_res*(-N:N)).^2;
% IRF_DD = IRF_DD./sum(IRF_DD(:));

% FF_scale = (t(end)-t(1))*2*DDAcf.fc*V/(DDAcf.c*H); %now, multiplied with d (y_a in Egido), we get along track the sinc argument
% IRF_FF = IRF.*sinc(FF_scale*d_res*(-N:N)).^2;
% IRF_FF = IRF_FF./sum(IRF_FF(:));

% calculate grating lobe positions
f_D_rate = 2*DDAcf.fc*(V^2/H)/DDAcf.c % doppler rate [Hz/s]
BP = f_D_rate*DDAcf.BRI % units 1/s
ghost_spacing_m = V/BP; % units in m

DD_scale = (DDAcf.Nb*DDAcf.PRI)*2*DDAcf.fc*V/(DDAcf.c*H); %now, multiplied with d (y_a in Egido), we get along track the sinc argument
%DD_scale = 1/250;
envelope_DD = sinc(DD_scale*d_res*(-N:N)).^2;
IRF_DD = IRF.*envelope_DD;
IRF_DD = IRF_DD./sum(IRF_DD(:));


FF_scale = (t(end)-t(1))*2*DDAcf.fc*V/(DDAcf.c*H); %now, multiplied with d (y_a in Egido), we get along track the sinc argument

%
grating_lobes_FF = sum(sinc(FF_scale*(d_res*(-N:N) - ghost_spacing_m*(-40:40)')).^2,1);
envelope_FF = envelope_DD.*grating_lobes_FF;
IRF_FF = IRF.*envelope_FF;
IRF_FF = IRF_FF./sum(IRF_FF(:));

% %% florian resolution check
% 
% L=T*V/2;
% res_a=(((L^2+H^2)^0.5+DDAcf.lambda0/4)^2-H^2)^0.5-L

%%
Pu = 1;
epoch = -0;

%HSs = 0.2:0.05:10 % for paper plot
HSs = 1:0.5:7;%1%:0.05:1

ACF_mod = zeros(2*N+1,numel(HSs),3);

for m = 1:numel(HSs)
    Hs = HSs(m)
    waveform = SAMOSA(Pu,epoch,Hs,CS1b);
    %plot(waveform);hold on

    % evaluate partial derivatives
    dp = 1e-8;
    dwdHs = (SAMOSA(Pu,epoch,Hs+dp,CS1b)-SAMOSA(Pu,epoch,Hs,CS1b))./dp;
    dwdepoch = (SAMOSA(Pu,epoch+dp,Hs,CS1b)-SAMOSA(Pu,epoch,Hs,CS1b))./dp;
    dwdPu = (SAMOSA(Pu+dp,epoch,Hs,CS1b)-SAMOSA(Pu,epoch,Hs,CS1b))./dp;
    % figure;
    % plot(dwdHs);hold on
    % plot(dwdepoch);hold on;
    % plot(dwdPu,'ro');hold on;

    % build jacobian matrix with derivatives in columns
    J = [dwdPu,dwdepoch,dwdHs];

    % now, weighting matrix params = W * measured_signal + const. is given by
    W = (J'*J)^(-1)*J';

    % % plot all the different weighting vectors for each parameter
%     figure('units','inch','position',[0,0,12,3.2]);
%     ax1 = subplot(1,3,1)
% 
%     plot(waveform);hold on
%     plot(W(1,:)./max(abs(W(1,:))))
%     title('Pu')
%     ylabel('normalized power')
%     grid on
% 
%     ax2 = subplot(1,3,2)
% 
%     plot(waveform);hold on
%     plot(W(2,:)./max(abs(W(2,:))))
%     title('epoch')
%     xlabel('range gate')
%     grid on
% 
% 
%     ax3 = subplot(1,3,3)
% 
%     plot(waveform);hold on
%     plot(W(3,:)./max(abs(W(3,:))))
%     title('Hs')
%     grid on
%     legend('waveform power','weighting vector')
% 
% 
%     linkaxes([ax1,ax2,ax3]);
%     ax1.XLim = [60,120];
%     ax1.YLim = [-1.1,1.1];

    % test the whole procedure, with a 0.1 m increaeed waveheight, SAMOSA model
    % the linear change should be predictable: Dparams = W* (S1-S2) should
    % diagnose e.g. 0.1 m more waveheight, if correct

    dparams = W*(SAMOSA(Pu-0.01,epoch+0.1,Hs+0.1,CS1b)-SAMOSA(Pu,epoch,Hs,CS1b));
    % check! this is exactly what we want! Plus some errors from the linear
    % approximation...

    % now use the weights and waveform to convolute with the calculated ACF
    %figure;
    %imagesc(IRF_DD);
    for n = 1:3
        w = W(n,:)'; % weight from before
        sigma = waveform; % noise std proportional to waveform power

        Q = w.*sigma;

        % now convolve the IRF in range direction with Q and then take the weighted
        % sum with Q to get the ACF of the SSH estimate:

        dum = conv2(IRF_DD,Q,'same');
        %figure; imagesc(dum)

        ACF_mod(:,m,n) = sum(dum.*Q)./max(sum(dum.*Q));
    end    
%     % now plot the resulting acf together with the estimated ones:
%     %figure;
%     plot(d_res*(-N:N),sum(IRF_DD)./max(sum(IRF_DD)),'k');hold on
%     plot(d_res*(-N:N),IRF_DD(128,:)./max(IRF_DD(128,:)),'k');hold on
%     plot(d_res*(-N:N),Cepoch./max(Cepoch),'r');hold on
%     xlim([0,600])
end
%%
figure;
imagesc(HSs,d_res*(-N:N),(ACF_mod(:,:,1)));colormap('pink');xlabel('significant wave height (m)');ylabel('along-track distance (m)');title('Pu')
figure;
imagesc(HSs,d_res*(-N:N),(ACF_mod(:,:,2)));colormap('pink');xlabel('significant wave height (m)');ylabel('along-track distance (m)');title('SSH')
figure;
imagesc(HSs,d_res*(-N:N),(ACF_mod(:,:,3)));colormap('pink');xlabel('significant wave height (m)');ylabel('along-track distance (m)');title('SWH')

%% plot the functions with black lines
figure;
subplot(1,3,1)
plot(d_res*(-N:N),ACF_mod(:,:,1),'k')
xlim([0,320])
title('Pu')
grid on
ylabel('noise auto-correlation function')

subplot(1,3,2)
plot(d_res*(-N:N),ACF_mod(:,:,2),'k')
xlim([0,320])
grid on
title('SSH')

subplot(1,3,3)
plot(d_res*(-N:N),ACF_mod(:,:,3),'k')
xlim([0,320])
title('SWH')
grid on
xlabel('along-track distance')

%% estimate the frequencies within. ACF and PSD (square of fourier amplitudes) are related by fourier transform.

figure('units','inch','position',[0,0,4.5,3.5]);

x = d_res*(-N:N);
Lx = 1/DD_scale;
Rpu = ACF_mod(:,:,1);%./sum(ACF_mod(:,:,1));
Rssh = ACF_mod(:,:,2);%./sum(ACF_mod(:,:,1));
Rswh = ACF_mod(:,:,3);%./sum(ACF_mod(:,:,1));
Ropt = sinc(x/Lx).^2;%./sum(sinc(x/Lx).^2);

% periodogram plots
Nfft = 20000;
[Sopt,f] = periodogram(Ropt,[],Nfft,1/d_res);
[Spu,f] = periodogram(Rpu,[],Nfft,1/d_res); %hamming(numel(x))
[Sssh,f] = periodogram(Rssh,[],Nfft,1/d_res); %hamming(numel(x))
[Sswh,f] = periodogram(Rswh,[],Nfft,1/d_res); %hamming(numel(x))
% sampling is 1 sample per d_res meters
% that means psd is given with cycles/meter
% 1/Lx cycles per meter = 20 Hz, hence
% f*1/Lx*20 = frequency in hertz?
fnew = 20*f*Lx;
start_ind = 2;
Sopt_ = movmean(sqrt(Sopt),1)./sqrt(Sopt(start_ind));
Spu_ = movmean(sqrt(Spu),1)./sqrt(Spu(start_ind));
Sssh_ = movmean(sqrt(Sssh),1)./sqrt(Sssh(start_ind));
Sswh_ = movmean(sqrt(Sswh),1)./sqrt(Sswh(start_ind));

% semilogy(fnew,Sopt_/mean(Sopt_(end-100:end)-Sswh(end-100:end)));hold on
% semilogy(fnew,Spu_/mean(Spu_(end-100:end)-Sswh(end-100:end)));hold on
% semilogy(fnew,Sssh_/mean(Sssh_(end-100:end)-Sswh(end-100:end)));hold on
% semilogy(fnew,Sswh_/mean(Sswh_(end-100:end)-Sswh(end-100:end)));hold on

plot(fnew(start_ind:end),Sopt_(start_ind:end),'LineWidth',2,'LineStyle','-','Color','k');hold on
plot(fnew(start_ind:end),Spu_(start_ind:end),'LineWidth',2,'LineStyle','--','Color','k');hold on
plot(fnew(start_ind:end),Sssh_(start_ind:end),'LineWidth',2,'LineStyle','-.','Color','k');hold on
plot(fnew(start_ind:end),Sswh_(start_ind:end),'LineWidth',2,'LineStyle',':','Color','k');hold on
grid on

xlim([0,80])
xlabel('posting rate equivalent (Hz)')
ylabel('noise PSD (normalized to maximum)')

% methodology confirmed by observing that R = sinc^2 yields triangular
% function. Since F[R] is F[sinc*sinc] = F[sinc] \ast F[sinc] = box \ast
% box = triangle function

% calculate where the PSDs fall below 1%:
iopt1 = find(Sopt_ < 0.01,1,'first');
iopt2 = find(Sopt_ > 0.01,1,'last');
fopt = 0.5*(fnew(iopt1) + fnew(iopt2))

iopt1 = find(Spu_ < 0.01,1,'first');
iopt2 = find(Spu_ > 0.01,1,'last');
fpu = 0.5*(fnew(iopt1) + fnew(iopt2))

iopt1 = find(Sssh_ < 0.01,1,'first');
iopt2 = find(Sssh_ > 0.01,1,'last');
fssh = 0.5*(fnew(iopt1) + fnew(iopt2))

iopt1 = find(Sswh_ < 0.01,1,'first');
iopt2 = find(Sswh_ > 0.01,1,'last');
fswh = 0.5*(fnew(iopt1) + fnew(iopt2))

legend(['R=sinc^2(x/L_x)' newline 'f(PSD = 0.01) \approx ',num2str(fopt,'%.1f'), ' Hz'],...
       ['Pu, sigma0' newline 'f(PSD \approx 0.01) \approx ',num2str(fpu,'%.1f'), ' Hz'],...
       ['epoch, range, SSH' newline 'f(PSD \approx 0.01) \approx ',num2str(fssh,'%.1f'), ' Hz'],...
        ['SWH' newline 'f(PSD \approx 0.01) \approx ',num2str(fswh,'%.1f'), ' Hz'])

%legend(['a' newline 'b'],['a' newline 'c'])

exportgraphics(gcf,['/home/fehelers/PhD_Delft/2nd_paper/PSD_estimates_SWH=1.pdf'],'Resolution',300)

%% estimate the frequencies within. ACF and PSD (square of fourier amplitudes) are related by fourier transform.
S = periodogram(R);
plot(abs(S));hold on


R = ACF_mod(:,:,3);

S = periodogram(R);
plot(abs(S));hold on


S = periodogram(Ropt);
plot(abs(S));hold on



%% Samosa function at end of script
function [waveform] = SAMOSA(Pu,epoch,Hs,CS1b)
    % read data from TUM processor
    data_dir = [getenv('HOME') '/TUDTUM/'];
    %tum_data = load([data_dir '/sam_model_validation/samosa_model_validation_tu_delft_swh0-20.mat']);
    
    DDAcf = DDA_ConfigFile('S3A','SAR')
    % constants
    CONST_A = 6378137; %Equatorial radius [m]
    CONST_B = 6356752.3142; % Polar radius [m]
    CONST_C = 299792458; % speed of light [m/s]

    % prepare input for TUD model
    %dt = 1./tum_data.sensor_sets.B_r_Hz;

    x(1) = Pu;%tum_data.Pu; %Pu
    x(2) = epoch;%double(tum_data.epoch_ns); %epoch
    x(3) = Hs;%2; %Hs

    % some preliminaries
    % calc Re
    Re = sqrt(CONST_A.^2 * cos(CS1b.GEO.LAT(1)/360*2*pi).^2 + CONST_B.^2 .* sin(CS1b.GEO.LAT(1)/360*2*pi).^2);

    % calc orbit_slope
%     if tum_data.model_params.ascending
%         track_sign = -1;
%     else
%         track_sign = 1;
%     end
    orbit_slope = 0;%track_sign*((CONST_A^2 - CONST_B^2)./(2*Re.^2)).*sin(2*tum_data.model_params.lat_rad) - (-double(tum_data.model_params.h_rate_m_per_s)./double(tum_data.model_params.Vs_m_per_s));

    % fill n data vector with parameters
    % n = [n;TN;NaN;Re;h;vt;OrbSlp;theta_p;theta_r;SmodID;EstPar];
    n = zeros(138,1);
    n(1:266) = (1:266);

    n(end) = 3;
    n(end-1) = 1; %SmodID;
    n(end-2) = 0;%tum_data.model_params.ksiy_rad; %roll mispointing angle (across-track)
    n(end-3) = 0;%tum_data.model_params.ksix_rad; %pitch mispointing angle (along-track)
    n(end-4) = orbit_slope;
    n(end-5) = CS1b.GEO.V.V(1);%tum_data.model_params.Vs_m_per_s;
    n(end-6) = CS1b.GEO.H(1);%tum_data.model_params.alt_m;
    n(end-7) = Re;
    n(end-8) = NaN;
    n(end-9) = 0.0; %TN;

%     DDAcf.Np = double(tum_data.wf_sets.np);
%     DDAcf.RefBin = double(tum_data.model_params.epoch_ref_gate) + 1.0;
%     DDAcf.fp = tum_data.sensor_sets.prf_Hz; % 1.782531194295900e+04
%     DDAcf.Nb = double(tum_data.sensor_sets.n_b); %64
%     DDAcf.B = tum_data.sensor_sets.B_r_Hz; %320e6
%     DDAcf.shx = 1; %1
%     DDAcf.shy = 1; %1
%     DDAcf.theta_x = tum_data.sensor_sets.theta_x_rad; %0.023352505391684
%     DDAcf.theta_y = tum_data.sensor_sets.theta_y_rad; %0.023352505391684
%     DDAcf.c = CONST_C;
%     DDAcf.fc = tum_data.sensor_sets.f_c_Hz; %1.357500000000000e+10
%     DDAcf.alpha_PF = tum_data.model_sets.alpha_p_mean; %0.5
%     DDAcf.lambda0 = CONST_C .* 1/tum_data.sensor_sets.f_c_Hz; %0.022084158968692

    %DDAcf = DDA_ConfigFile('S3A','SAR');

    fid = fopen('alphap_table_SEN3_09_Nov_2017.txt'); LUT_AP = cell2mat(textscan(fid,'%f %f','HeaderLines',1)); fclose(fid);
    LUT_AP = griddedInterpolant(LUT_AP(:,1),LUT_AP(:,2),'linear','nearest');
    fun_x = (-19:.001:41.999)';
    LUT_F0 = (pi/4) * abs(fun_x).^0.5 .* (besseli(-1/4,(.5*fun_x).^2,1) + sign(fun_x).*besseli(1/4,(.5*fun_x).^2,1));
    LUT_F0(fun_x == 0) = (pi*2^(3/4)) / (4 * gamma(3/4));
    LUT_F0 = griddedInterpolant(fun_x,LUT_F0,'spline','nearest');
    LUT_F1 = (pi/8) * abs(fun_x).^1.5 .* ((besseli(1/4,(.5*fun_x).^2,1)-besseli(-3/4,(.5*fun_x).^2,1)) + sign(fun_x).*(besseli(-1/4,(.5*fun_x).^2,1)-besseli(3/4,(.5*fun_x).^2,1)));
    LUT_F1(fun_x == 0) = -(2^(3/4) * gamma(3/4)) / 4;
    LUT_F1 = griddedInterpolant(fun_x,LUT_F1,'spline','nearest');

    vt = CS1b.GEO.V.V(1);%double(tum_data.model_params.Vs_m_per_s);
    BAstck = CS1b.SAR.BeamAngle(:,1);%tum_data.model_params.beam_ang_stack_rad - pi/2;
    DFstck  = (2*vt/DDAcf.lambda0) * sin(BAstck);     %Eq. 90
    dfa     = DDAcf.fp/DDAcf.Nb;                     %Eq. 61
    %BeamIDX = double(unique(int32(DFstck/dfa)));
    BeamIDX = double(DFstck/dfa);

    stack_mask_start_stop = CS1b.SAR.stack_mask_start_stop(:,1);

    waveform = RETRACKER.SAMOSA_DPM_v2_5_2_fun(x,n,DDAcf,LUT_AP,LUT_F0,LUT_F1,BeamIDX',stack_mask_start_stop',false);
end