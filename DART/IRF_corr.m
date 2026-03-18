function kernel = IRF_corr(CS1b, FFSAR_processing_settings)

% The goal of this script is to plot the IRF of S6 unfocused SAR and S6
% focused SAR impulse response functions (IRF) and to derive a convolution
% kernel, that transforms S6 FF-SAR waveforms into S6 UF-SAR waveforms via:
% UF-SAR-wf = conv(FF-SAR-wf,kernel,'same');
%% get some orbit parameters from CS1b and choose some settings
DDAcf = DDA_ConfigFile('S6A','SAR');
%DDAcf.BRI = CS1b.GEO.BRI;
DDAcf.PRI = median(CS1b.MEA.PRI);

T = FFSAR_processing_settings.integration_time;

H = median(CS1b.MEA.tracker_range);
V = median(CS1b.GEO.V.V);
dHdt = median(CS1b.GEO.H_rate);

multi_target = false;
num_targets = 21; % number of point targets in the center range bin, if multi_target=true; this is important to assume, because given a signal value in a range bin, we cannot know where the target is; take something odd for symmetrical filling of rangebin;

%os_zp = 2; %oversampling of the waveform
n_bins = DDAcf.os_ZP*DDAcf.Np;

%% calculate exact range history in flat earth approximation and compare it to its Taylor approximation
% R0 = H;
%t = -T/2:0.02:T/2;
t = (1:401)/401*T - 401/400*T/2;
alpha = asin(dHdt/V);

% Rt = sqrt(  (H + V*t*sin(alpha)).^2 + (V*t*cos(alpha)).^2 );
% Rt2 = sqrt(  H.^2 + 2*H*V*t*sin(alpha) + (V*t).^2 ); 
% Rt3 = H.*sqrt(  1 + 2/H*V*t*sin(alpha) + (V*t/H).^2 ); 
% Rt4 = H.*(1 + 0.5*(2/H*V*t*sin(alpha) + (V*t/H).^2));
% Rt4 = H.*(1 + 0.5*(2/H*V*t*sin(alpha) + (V*t/H).^2));
% Rt5 = H + V*t*sin(alpha) + 0.5*(V*t).^2/H; % take this as model for range history
% % magnitude analysis
% % H.^2 = 6.4000e+11
% % 2*H*V*t*sin(alpha) = 1.6e+7
% % (V*t).^2 = 4.9e+7
% 
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
% gamma_y = 8*log(2)/DDAcf.theta_y^2; %are the same

%antenna_pattern = exp(-gamma_x*theta_p.^2); define antenna pattern in the
%loop depending on changing theta_p

center = n_bins/2; %center bin
range_res_sample = 1/DDAcf.B*DDAcf.c/2/DDAcf.os_ZP; %oversampled range resolution
range_res_transmit = 1/DDAcf.Bt*DDAcf.c/2; %actual range resolution and scale of the sinc function
rangebin = (1:n_bins); % floating point range bin number
d_res = 0.25;
N = 4000;
% so that the sidelobes from +- 1 km distance will be included
IRF = zeros(n_bins,2*N+1);
% IRF_new = zeros(n_bins,2*N+1);
% M = size(t,2);
% M_c = round(size(t,2)/2);

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
%     i
end
toc

%%

DD_scale = (DDAcf.Nb*DDAcf.PRI)*2*DDAcf.fc*V/(DDAcf.c*H); %now, multiplied with d (y_a in Egido), we get along track the sinc argument
IRF_DD = IRF.*sinc(DD_scale*d_res*(-N:N)).^2;
IRF_DD = IRF_DD./sum(IRF_DD(:));

FF_scale = (t(end)-t(1))*2*DDAcf.fc*V/(DDAcf.c*H); %now, multiplied with d (y_a in Egido), we get along track the since argument
IRF_FF = IRF.*sinc(FF_scale*d_res*(-N:N)).^2;
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

end
