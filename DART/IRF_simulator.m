% The goal of this script is to plot the IRF of S6 unfocused SAR and S6
% focused SAR impulse response functions (IRF) and to derive a convolution
% kernel, that transforms S6 FF-SAR waveforms into S6 UF-SAR waveforms via:
% UF-SAR-wf = conv(FF-SAR-wf,kernel,'same');
%% get some orbit parameters from CS1b and choose some settings

%mission = 'S6A'
mission = 'S3B'
%mission = 'CS'

if strcmp(mission,'S6A')
    DDAcf = DDA_ConfigFile('S6A','SAR');
    DDAcf.BRI = median(CS1b.GEO.BRI);
    DDAcf.PRI = median(CS1b.MEA.PRI);
    %T = 3.6;
    %T = 3.413;
    T = 2;
end

if strcmp(mission,'S3B')
    DDAcf = DDA_ConfigFile('S3B','SAR');
    %DDAcf.BRI = CS1b.GEO.BRI;
    DDAcf.PRI = DDAcf.PRI
    %T = 2.2922;%2.1;
    T = CS1b.GEO.integration_time
end

if strcmp(mission, 'CS')
    DDAcf = DDA_ConfigFile('CS','SAR');
    %DDAcf.BRI = CS1b.GEO.BRI;
    DDAcf.PRI = DDAcf.PRI
    T = 1.9;
end

%exclusively for paper plots:
%DDAcf.os_ZP = 8;


%% get different integration times from stack mask:
% binary_mask = zeros(size(CS1b.SAR.stack_mask_start_stop(:,1),1),DDAcf.os_ZP*DDAcf.Np);
% 
% for i = 1:size(CS1b.SAR.stack_mask_start_stop(:,1),1)
%     binary_mask(i,1:CS1b.SAR.stack_mask_start_stop(i,1)) = 1;
% end
% 
% figure;imagesc(binary_mask)
% num_valid_entries = sum(binary_mask,1);
% figure; plot(num_valid_entries./max(num_valid_entries));
% 
% T_fraction = num_valid_entries./max(num_valid_entries);

%%

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

%% calculate exact range history in flat earth approximation and compare it to its Taylor approximation
R0 = H;
%t = -T/2:0.02:T/2;
t = (1:401)/401*T - 401/400*T/2
%t = -T/2;
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
%%
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
d_res = 0.25;
N = 8000;
% so that the sidelobes from +- 1 km distance will be included
IRF = zeros(n_bins,2*N+1);
IRF_new = zeros(n_bins,2*N+1);
M = size(t,2);
M_c = round(size(t,2)/2);

rel_target_pos(1,1,:) = 1/num_targets*((1:num_targets)-mean(1:num_targets)); % relative position of several point targets within the center range bin

%% calculate IRFs new
x = d_res*(-N:N); % ground-projected along-track distance of the scatterer to nadir
t_p = t;  % determined/given as some inform time vector above already
r = range_res_sample*((1:n_bins)-round(n_bins/2)); % positions the scatterer exactly in the middle of a bin

% make 3D meshgrid of required coordinates (screws memory)
% [x,r,t_p] = meshgrid(   d_res*(-N:N),...
%                         range_res_sample*((1:n_bins)-round(n_bins/2)),...
%                         t)
%% calculate IRF (easy but slower)
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

%% calculate IRFs (looks complicated but faster)
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
    antenna_pattern = exp(-gamma_x*theta_p.^2).^2; % squared in case of ACF instead of PTR
    
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

%%
 
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

%%
grating_lobes_FF = sum(sinc(FF_scale*(d_res*(-N:N) - ghost_spacing_m*(-40:40)')).^2,1);
envelope_FF = envelope_DD.*grating_lobes_FF;
IRF_FF = IRF.*envelope_FF;
IRF_FF = IRF_FF./sum(IRF_FF(:));
%% get convolution kernel
% for the calculation of the kernel, it does not matter the x-scale or any
% distortion in this direction, because the kernel is integrated over that
% direction anyway

DDkernel = sum(IRF_DD,2);
FFkernel = sum(IRF_FF,2);

FTInput = fft(FFkernel);
FtOutput = fft(DDkernel);
kernel = fftshift(ifft(FtOutput./FTInput));
kernel = kernel(2:end); %to have the kernel centered

%% plot IRFs for UF-SAR and FF-SAR and the corresponding convolution kernel
figure;
ax1 = subplot(1,3,1)
imagesc(log10(IRF_DD))
caxis([max(log10(IRF_DD(:)))-5,max(log10(IRF_DD(:)))]);
colorbar();
colormap('pink');
title({'UF-SAR IRF', 'multiple targets within range bin (0/1) = ', num2str(multi_target)})
ax2 = subplot(1,3,2)
imagesc(log10(IRF_FF))
caxis([max(log10(IRF_FF(:)))-5,max(log10(IRF_FF(:)))]);
colorbar();
colormap('pink');
title({'FF-SAR IRF', 'multiple targets within range bin (0/1) = ', num2str(multi_target)})
ax3 = subplot(1,3,3)
plot(DDkernel,1:n_bins); hold on
plot(FFkernel,1:n_bins);hold on
plot(kernel,1:n_bins-1);hold on
plot(conv2(FFkernel,kernel,'same'),1:n_bins,'ro');hold on
legend('Delay/Doppler widening', 'S6 FFSAR widening', 'transformation kernel')

linkaxes([ax1,ax2],'xy')
%linkaxes([ax1,ax2,ax3],'y')

%% comparison of a slice through IRF DD to argue for the along track resolution issue
figure;
DD_theory = sum(IRF_DD)./max(sum(IRF_DD));
DD_practice = IRF_DD(256,:)./max(IRF_DD(256,:));

plot(x,DD_theory); hold on
plot(x,DD_practice); hold on
plot([-1000,1000],[1/exp(1),1/exp(1)],'k'); hold on
xlabel('along-track distance (m)')
ylabel('normalized power')
grid on;
xlim([-150,150])

legend('UF-SAR theoretical','UF-SAR range smeared')



%% P(r) comparison
figure;
plot(DDkernel(1:4:end)); hold on
plot(FFkernel(1:4:end));hold on

figure;
plot(DDkernel(2:4:end)); hold on
plot(FFkernel(2:4:end));hold on

figure;
plot(DDkernel(3:4:end)); hold on
plot(FFkernel(3:4:end));hold on


%% theoretically, each bin of the waveform has a different effective integration time due to masking. Hence, the convolution with a constant kernel is not strictly correct!

% figure;
% %[val,ind] = max(irfDD,[],2);
% ind = 4001;
% plot(log10(IRF_DD(:,ind)));hold on
% plot(log10(IRF_DD(:,ind+1000)));hold on
% plot(log10(IRF_DD(:,ind-1000)));hold on
% plot(log10(IRF_DD(:,ind+3000)));hold on
% plot(log10(IRF_DD(:,ind-3000)));hold on
% 
% grid on

%%
% 
% figure;
% [val,ind] = max(IRF,[],2);
% ind = 4001;
% plot(log10(IRF(:,ind)));hold on
% plot(log10(IRF(:,ind+1000)));hold on
% plot(log10(IRF(:,ind-1000)));hold on
% plot(log10(IRF(:,ind+3000)));hold on
% plot(log10(IRF(:,ind-3000)));hold on
% 
% grid on

%%
% make default plot backgroudn white
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);

% choose color palette
cFF = [207, 0, 110]/255;
cDD = [0, 185, 227]/255;
cT = [214, 214, 214]/255;
%% paper plot 

%mission = 'S6A'
irfDD = IRF_DD;
irfFF = IRF_FF;

PtotDD = sum(sum(irfDD));
PtotFF = sum(sum(irfFF));

irfDD = irfDD/PtotDD;
irfFF = irfFF/PtotFF;

% calculate scale of DD along track PTR sinc function
DD_scale = (DDAcf.Nb*DDAcf.PRI)*2*DDAcf.fc*median(CS1b.GEO.V.V)/(DDAcf.c*median(CS1b.MEA.tracker_range)); %now, multiplied with d (y_a in Egido), we get along track the sinc argument
m_per_look = d_res
N_looks = 1:size(CS1b.SAR.data,2);

% range sinc parameter
Range_scale = 1/DDAcf.os_ZP; % correct for S3 but also for S6!?
if strcmp(mission,'S6A')
    Range_scale = DDAcf.Bt/DDAcf.B/DDAcf.os_ZP
end


% calculate distance to transponder projected in along track direction
dist2Tr = abs(x); %1e3*deg2km(distance(CS1b.GEO.LAT(:),CS1b.GEO.LON(:),TR.crete.lat,TR.crete.lon));
% this distance is an hyperpbolic, because of the remaining cross-track
% distance of the transponder, correct for this:
% dist2Tr = sqrt(dist2Tr.^2 - min(dist2Tr.^2));
[val,ind_x] = min(dist2Tr);
dist2Tr(1:ind_x) = -dist2Tr(1:ind_x);

% position of transponder in range:
[val,ind_r] = max(irfDD(:,ind_x))

DD_theo = sinc(DD_scale*dist2Tr).^2/sum(sinc(DD_scale*dist2Tr).^2);
dr = 0.1*(1:DDAcf.Np*DDAcf.os_ZP*10);
if strcmp(mission,'S6A')
    offset = 0.;
elseif strcmp(mission,'S3A')|strcmp(mission,'S3B')|strcmp(mission,'CS')
    offset = 0;
end
Range_theo = 10*sinc(0.1*Range_scale*(((1:DDAcf.Np*DDAcf.os_ZP*10) - 10*(ind_r+offset)))).^2/ sum(sinc(0.1*Range_scale*(((1:DDAcf.Np*DDAcf.os_ZP*10) - 10*ind_r))).^2);

figure('units','inch','position',[0,0,10,5]);
% first plot the integrated and cumulated power for DD:
ax1 = subplot(6,2,1);
plot(dist2Tr,cumsum(DD_theo),'LineWidth',4,'Color',cT); hold on
plot(dist2Tr,cumsum(sum(irfDD)),'Color',cDD); hold on
plot(dist2Tr,cumsum(sum(irfFF)),'Color',cFF); hold on
%legend('UF-SAR theoretical','UF-SAR','FF-SAR')
ylabel('\int P(x)/P_{tot}')%,'Interpreter','latex')
grid()
colorbar() % placeholder

ax2 = subplot(6,2,3);
plot(dist2Tr,10*log10(DD_theo),'LineWidth',4,'Color',cT); hold on
plot(dist2Tr,10*log10(sum(irfDD)),'Color',cDD); hold on
plot(dist2Tr,10*log10(sum(irfFF)),'Color',cFF); hold on
%int_irfFF = movsum(sum(irfFF),300);
int_irfFF = movsum(sum(irfFF),20);
[peak_val,peak_ind] = findpeaks(int_irfFF,'MinPeakProminence',0.00001);
%plot(dist2Tr,10*log10(int_irfFF),'Color',cFF); hold on
s = scatter(dist2Tr(peak_ind),10*log10(int_irfFF(peak_ind)/max(int_irfFF(peak_ind))*max(sum(irfDD))),50,'filled','MarkerEdgeColor',cFF,'MarkerFaceColor',cFF); hold on
s.MarkerEdgeAlpha = 0.5;
s.MarkerFaceAlpha = .5;

legend('UF-SAR theoretical','UF-SAR','FF-SAR','FF-SAR integrated peak energy')
ylabel('P(x)/P_{tot} (dB)')
ylim([-80,0])
grid()
colorbar() % placeholder

% now along range axis
axr1 = subplot(3,6,10);
plot(10*log10(Range_theo),dr,'k','LineWidth',2,'Color',cT); hold on
plot(10*log10(sum(irfDD,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cDD); hold on
xlim([-45,0])
set(gca, 'YDir','reverse');
%set(gca, 'XDir','reverse');
%legend('UF-SAR theoretical','UF-SAR')
xlabel('P(r)/P_{tot} (dB)')
grid()

axr2 = subplot(3,6,16);
plot(10*log10(Range_theo),dr,'k','LineWidth',2,'Color',cT); hold on
plot(10*log10(sum(irfFF,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cFF); hold on
%plot(10*log10(irfFF(:,11650)/max(irfFF(:,11650))),1:DDAcf.Np*DDAcf.os_ZP,'Color',cFF); hold on
set(gca, 'YDir','reverse');
xlim([-45,0])
xlabel('P(r)/P_{tot} (dB)')
%set(gca, 'XDir','reverse');
%legend('FF-SAR theoretical','FF-SAR')
grid()

% now along range axis
axr3 = subplot(3,6,11);
plot(cumsum(Range_theo/10),dr,'k','LineWidth',2,'Color',cT); hold on
plot(cumsum(sum(irfDD,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cDD); hold on
set(gca, 'YDir','reverse');
%set(gca, 'XDir','reverse');
%legend('UF-SAR')
xlabel('\int P(r)/P_{tot}')
grid()

axr4 = subplot(3,6,17);
plot(cumsum(Range_theo/10),dr,'k','LineWidth',2,'Color',cT); hold on
plot(cumsum(sum(irfFF,2)),1:DDAcf.Np*DDAcf.os_ZP,'Color',cFF); hold on
set(gca, 'YDir','reverse');
%set(gca, 'XDir','reverse');
%legend('FF-SAR')
xlabel('\int P(r)/P_{tot}')
grid()

% now plot irf's

ax3 = subplot(3,2,3);
imagesc(dist2Tr,1:DDAcf.Np*DDAcf.os_ZP,10*log10(irfDD))
cax = colorbar()
ylabel(cax,'P(r,x)/P_{tot} (dB)')
caxis([-80,0])
colormap('pink')
ylabel('range bin index')

ax4 = subplot(3,2,5);
imagesc(dist2Tr,1:DDAcf.Np*DDAcf.os_ZP,10*log10(irfFF))
cax = colorbar()
ylabel(cax,'P(r,x)/P_{tot} (dB)')
caxis([-80,0])
xlabel('transponder distance (m)')
ylabel('range bin index')

linkaxes([ax1,ax2,ax3,ax4],'x')
linkaxes([ax3,ax4,axr1,axr2,axr3,axr4],'y')
linkaxes([axr1,axr2],'x')

if strcmp(mission,'S6A')
    ax3.XLim = [-600,600]
    ax3.YLim = [256-30,256+30]%60*DDAcf.os_ZP]
end

if strcmp(mission,'S3B')
    ax3.XLim = [-650,650]
    ax3.YLim = [128-30,128+30]%60*DDAcf.os_ZP]
end

ax2.YTick = [-80 -60 -40 -20 0]
% 
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);

%exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/paperIRF_' mission '_model.pdf'],'Resolution',300)


%exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/paperIRF_' mission '_model.svg'],'Resolution',300)
%saveas(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/paperIRF_' 'S6_model.pdf'])

%% plot simple plot of P(r)/P_tot for S3 and S6


