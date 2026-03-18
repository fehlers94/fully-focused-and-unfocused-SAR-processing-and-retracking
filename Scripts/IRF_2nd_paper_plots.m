% The goal of this script is to plot the IRF of S6 unfocused SAR and S6
% focused SAR impulse response functions (IRF) and to derive a convolution
% kernel, that transforms S6 FF-SAR waveforms into S6 UF-SAR waveforms via:
% UF-SAR-wf = conv(FF-SAR-wf,kernel,'same');
clear all
% make default plot backgroudn white
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);

% choose color palette
cFF = [207, 0, 110]/255;
cDD = [0, 185, 227]/255;
cT = [214, 214, 214]/255;

%% get some orbit parameters from CS1b and choose some settings

IRFs = {};
%sampling = 90;%5; 
sampling = 0.025;

%mission = 'S6A'
%mission = 'S3B'
%mission = 'CS'

for mission = [{'S6A'},{'CS'},{'S3B'}]
    mission = mission{1}
    
    if strcmp(mission,'S6A')
        DDAcf = DDA_ConfigFile('S6A','SAR');
        DDAcf.BRI = 718.1024*1e-5;%median(CS1b.GEO.BRI);
        DDAcf.PRI = 1.0876e-04;%median(CS1b.MEA.PRI);
        H = 1.3441e+06;%median(CS1b.MEA.tracker_range);
        V = 6.9703e+03;%median(CS1b.GEO.V.V);
        dHdt = 0;%median(CS1b.GEO.H_rate);
        %T = 3.6;
        %T = 3.413;
        %T = 2;
        %T_list = [(0:sampling:180)*3.413/180]
        T_illu = 1.928*2;
        T_list = [(0:sampling:1.5)*T_illu]
    end

    if strcmp(mission,'S3B')
        H = 8.0553e+05;%median(CS1b.MEA.tracker_range);
        V = 7.5441e+03;%median(CS1b.GEO.V.V);
        dHdt = 0;%median(CS1b.GEO.H_rate);
        DDAcf = DDA_ConfigFile('S3B','SAR');
        %DDAcf.BRI = CS1b.GEO.BRI;
        DDAcf.PRI = DDAcf.PRI
        %T = 2.1;
        %T_list = [(0:sampling:180)*DDAcf.BRI];
        T_illu = 0.9931*2;
        T_list = [(0:sampling:1.5)*T_illu]
    end

    if strcmp(mission,'CS')
        H = 717242;
        V = 7498;
        dHdt = 0;%median(CS1b.GEO.H_rate);
        DDAcf = DDA_ConfigFile('CS','SAR');
        %DDAcf.BRI = CS1b.GEO.BRI;
        DDAcf.PRI = DDAcf.PRI;
        %T = 2.1;
        %T_list = [(0:sampling:180)*DDAcf.BRI]
        T_illu = 0.879*2;
        T_list = [(0:sampling:1.5)*T_illu]
    end

    DDAcf.os_ZP = 2

    % for plotting of 3 different cases:
    %T_list = [64*DDAcf.PRI,  90*DDAcf.BRI, 180*DDAcf.BRI]

    for T = T_list

        multi_target = false;
        num_targets = 21; % number of point targets in the center range bin, if multi_target=true; this is important to assume, because given a signal value in a range bin, we cannot know where the target is; take something odd for symmetrical filling of rangebin;

        n_bins = DDAcf.os_ZP*DDAcf.Np;

        R0 = H;
        t = (1:401)/401*T - 401/400*T/2
        alpha = asin(dHdt/V);

        gamma_x = 8*log(2)/DDAcf.theta_x^2; 
        gamma_y = 8*log(2)/DDAcf.theta_y^2; %are the same

        %antenna_pattern = exp(-gamma_x*theta_p.^2); define antenna pattern in the
        %loop depending on changing theta_p

        center = n_bins/2; %center bin
        range_res_sample = 1/DDAcf.B*DDAcf.c/2/DDAcf.os_ZP; %oversampled range resolution
        if strcmp(mission,'S6A')
            range_res_transmit = 1/DDAcf.Bt*DDAcf.c/2; %actual range resolution and scale of the sinc function
        elseif strcmp(mission(1:2),'S3')
            range_res_transmit = 1/DDAcf.B*DDAcf.c/2; %actual range resolution and scale of the sinc function
        elseif strcmp(mission(1:2),'CS')
            range_res_transmit = 1/DDAcf.B*DDAcf.c/2; %actual range resolution and scale of the sinc function
        end

        rangebin = (1:n_bins); % floating point range bin number
        d_res = 8*0.25;
        N = 4000/8;
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
            antenna_pattern = exp(-gamma_x*theta_p.^2);

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

        % calculate grating lobe positions
        f_D_rate = 2*DDAcf.fc*(V^2/H)/DDAcf.c % doppler rate [Hz/s]
        BP = f_D_rate*DDAcf.BRI % units 1/s
        ghost_spacing_m = V/BP; % units in m

        DD_scale = (DDAcf.Nb*DDAcf.PRI)*2*DDAcf.fc*V/(DDAcf.c*H); %now, multiplied with d (y_a in Egido), we get along track the sinc argument
        envelope_DD = sinc(DD_scale*d_res*(-N:N)).^2;
        IRF_DD = IRF.*envelope_DD;
        IRF_DD = IRF_DD./sum(IRF_DD(:));

        FF_scale = (t(end)-t(1))*2*DDAcf.fc*V/(DDAcf.c*H); %now, multiplied with d (y_a in Egido), we get along track the sinc argument

        grating_lobes_FF = sum(sinc(FF_scale*(d_res*(-N:N) - ghost_spacing_m*(-40:40)')).^2,1);
        envelope_FF = envelope_DD.*grating_lobes_FF;
        IRF_FF = IRF.*envelope_FF;
        IRF_FF = IRF_FF./sum(IRF_FF(:));

        IRFs.(mission).(['IRF_DD_' num2str(round(100*T))]) = IRF_DD
    end
end
%% plot IRFs for UF-SAR
figure('units','inch','position',[0,0,4.5,3.5]);

clim = 35;
IRF_DD = IRFs.('S3B').IRF_DD_0/max(IRFs.('S3B').IRF_DD_0(:));

ax1 = subplot(3,2,1)
%imagesc(x,1:256,10*log10(IRF_DD)); hold on
imagesc(x,-127:128,IRF_DD); hold on
plot([-1000,1000],[0,0],'r--')
%caxis([max(log10(IRF_DD(:)))-clim,max(log10(IRF_DD(:)))]);
caxis([0 1]);
%colorbar();
colormap('pink');
title({'PTR, 1 burst'})
xticklabels([])

IRF_DD = IRFs.('S3B').IRF_DD_114/max(IRFs.('S3B').IRF_DD_114(:));

ax2 = subplot(3,2,3)
%imagesc(x,1:256,10*log10(IRF_DD));hold on
imagesc(x,-127:128,IRF_DD);hold on
plot([-1000,1000],[0,0],'r-.')
%caxis([max(log10(IRF_DD(:)))-clim,max(log10(IRF_DD(:)))]);
caxis([0 1]);
%h = colorbar();
%ylabel(h,'power (dB)')
colormap('pink');
title({'       multilooked PTR, 90 bursts'})
xticklabels([])
ylabel('range gate number')



IRF_DD = IRFs.('S3B').IRF_DD_228/max(IRFs.('S3B').IRF_DD_228(:));

ax3 = subplot(3,2,5)
%imagesc(x,1:256,10*log10(IRF_DD)); hold on
imagesc(x,-127:128,IRF_DD); hold on
plot([-1000,1000],[0,0],'r-')
%caxis([max(log10(IRF_DD(:)))-clim,max(log10(IRF_DD(:)))]);
caxis([0 1]);
%colorbar();
colormap('pink');
title({'       multilooked PTR, 180 bursts'})
xlabel('along-track distance (m)')


linkaxes([ax1,ax2,ax3],'xy')
ax1.YLim = [0-18,0+18]
ax1.XLim = [-490,490]

ax4 = subplot(1,2,2)

plot(x,envelope_DD,'Color',cT,'LineWidth',3); hold on
plot(x,IRFs.('S3B').IRF_DD_0(128*DDAcf.os_ZP/2,:)/max(IRFs.('S3B').IRF_DD_0(:)),'k--'); hold on
plot(x,IRFs.('S3B').IRF_DD_114(128*DDAcf.os_ZP/2,:)/max(IRFs.('S3B').IRF_DD_114(:)),'k-.'); hold on
plot(x,IRFs.('S3B').IRF_DD_228(128*DDAcf.os_ZP/2,:)/max(IRFs.('S3B').IRF_DD_228(:)),'k-'); hold on
ylim([-0.05, 1.25])
xlabel('along-track distance (m)')
ylabel('normalized power (no unit)')
grid on
legend('sinc^2(x/L_x)','1 burst','90 bursts','180 bursts')


%linkaxes([ax1,ax2,ax3,ax4],'x')
ax4.XLim = [-350,350]

% linkaxes([ax1,ax2,ax3],'y')

exportgraphics(gcf,['/home/fehelers/PhD_Delft/2nd_paper/paperIRF_resolution_' mission '_v2_noscale.pdf'],'Resolution',300)

%% determine full width half maximum of the plots
figure('units','inch','position',[0,0,4.5,3.5]);

ls.S6A = 'k-'
ls.S3B = 'k-.'
ls.CS = 'k--'

colororder({'k','k'})

for mission = [{'S6A'},{'CS'},{'S3B'}];
    mission = mission{1};

    if strcmp(mission,'S6A')
        DDAcf = DDA_ConfigFile('S6A','SAR');
        DDAcf.BRI = 718.1024*1e-5;%median(CS1b.GEO.BRI);
        DDAcf.PRI = 1.0876e-04;%median(CS1b.MEA.PRI);
        H = 1.3441e+06;%median(CS1b.MEA.tracker_range);
        V = 6.9703e+03;%median(CS1b.GEO.V.V);
        dHdt = 0;%median(CS1b.GEO.H_rate);
        %T = 3.6;
        %T = 3.413;
        %T = 2;
        %T_list = [(0:sampling:180)*3.413/180]
        T_illu = 1.928*2;
        T_list = [(0:sampling:1.5)*T_illu]
    end

    if strcmp(mission,'S3B')
        H = 8.0553e+05;%median(CS1b.MEA.tracker_range);
        V = 7.5441e+03;%median(CS1b.GEO.V.V);
        dHdt = 0;%median(CS1b.GEO.H_rate);
        DDAcf = DDA_ConfigFile('S3B','SAR');
        %DDAcf.BRI = CS1b.GEO.BRI;
        DDAcf.PRI = DDAcf.PRI
        %T = 2.1;
        %T_list = [(0:sampling:180)*DDAcf.BRI];
        T_illu = 0.9931*2;
        T_list = [(0:sampling:1.5)*T_illu]
    end

    if strcmp(mission,'CS')
        H = 717242;
        V = 7498;
        dHdt = 0;%median(CS1b.GEO.H_rate);
        DDAcf = DDA_ConfigFile('CS','SAR');
        %DDAcf.BRI = CS1b.GEO.BRI;
        DDAcf.PRI = DDAcf.PRI;
        %T = 2.1;
        %T_list = [(0:sampling:180)*DDAcf.BRI]
        T_illu = 0.879*2;
        T_list = [(0:sampling:1.5)*T_illu]
    end
    
    name = fieldnames(IRFs.(mission));
    FWHMs = zeros(numel(T_list),1);

    for i = 1:numel(T_list)
        f = (IRFs.(mission).(name{i})(DDAcf.Np*DDAcf.os_ZP/2,:)/max(IRFs.(mission).IRF_DD_0(:)));
        
        %determine FWHM from a linear interpolation
        center_ind = (numel(f)+1)/2
        f_interp = griddedInterpolant(d_res*(0:center_ind-1),f(center_ind:end)-0.5);
        x_hr = 0:0.1:300;
        f_hr = f_interp(0:0.1:300);
        FWHMs(i) = 2*x_hr(find(f_hr<0,1,'first'));%d_res*(find(f >= 0.5, 1, 'last') - find(f >= 0.5, 1, 'first'));
    end
    
    ax = gca;
    
    %FWHMs = FWHMs./sqrt((6371e3 + H)/6371e3)
    %yyaxis left
    plot(T_list./T_illu,20./(FWHMs./max(FWHMs)),ls.(mission)); hold on
    lgd = legend('Sentinel-6 MF','CryoSat-2','Sentinel-3')
    lgd.Location = 'northwest'
    grid on
    %ylabel('FWHM(T)/FWHM(0)')
    ylabel('required posting rate (Hz)')
    xlabel('used fraction of illumination time T_{m}/T_{illu}')
    ylim([10 70])
    xlim([0 1.2])
end

% %plot some lines as fwhm frequencies
% freq = 20:10:80
% fwhm_freq = 20./(freq);
% for i = 1:numel(freq)
%     plot([0 3],[fwhm_freq(i) fwhm_freq(i)],'k');
% end

exportgraphics(gcf,['/home/fehelers/PhD_Delft/2nd_paper/estimated_posting_rate_v1.pdf'],'Resolution',300)

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

linkaxes([ax1,ax2],'x')
linkaxes([ax1,ax2,ax3],'y')

%% comparison of a slice through IRF DD to argue for the along track resolution issue
figure;
DD_theory = sum(IRF_DD)./max(sum(IRF_DD));
DD_practice = IRF_DD(255,:)./max(IRF_DD(255,:));

plot(x,DD_theory); hold on
plot(x,DD_practice); hold on
plot([-1000,1000],[1/exp(1),1/exp(1)],'k'); hold on
xlabel('along-track distance (m)')
ylabel('normalized power')
grid on;
xlim([-450,450])

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
elseif strcmp(mission,'S3A')|strcmp(mission,'S3B')
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

exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/paperIRF_' mission '_model.pdf'],'Resolution',300)
%exportgraphics(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/paperIRF_' mission '_model.svg'],'Resolution',300)
%saveas(gcf,['/home/fehelers/PhD_Delft/pseudoDDprocessing/paperIRF_' 'S6_model.pdf'])

%% plot simple plot of P(r)/P_tot for S3 and S6


