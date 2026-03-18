clc; clear all;
close all;
%%
% make default plot backgroudn white
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);


%%
Pu = 1;
epoch = -0;
Hs = 5;

waveform = SAMOSA(Pu,epoch,Hs);
plot(waveform,'bo');hold on


%%

% evaluate partial derivatives
dp = 1e-8;
dwdHs = (SAMOSA(Pu,epoch,Hs+dp)-SAMOSA(Pu,epoch,Hs))./dp;
dwdepoch = (SAMOSA(Pu,epoch+dp,Hs)-SAMOSA(Pu,epoch,Hs))./dp;
dwdPu = (SAMOSA(Pu+dp,epoch,Hs)-SAMOSA(Pu,epoch,Hs))./dp;

plot(dwdHs);hold on
plot(dwdepoch);hold on;
plot(dwdPu,'ro');hold on;

% build jacobian matrix with derivatives in columns
J = [dwdPu,dwdepoch,dwdHs]

% now, weighting matrix params = W * measured_signal + const. is given by
W = (J'*J)^(-1)*J';

% plot all the different weighting vectors for each parameter
figure;
ax1 = subplot(1,3,1)

plot(waveform);hold on
plot(W(1,:)./max(abs(W(1,:))))
title('Pu')
ylabel('normalized power')
grid on

ax2 = subplot(1,3,2)

plot(waveform);hold on
plot(W(2,:)./max(abs(W(2,:))))
title('epoch')
xlabel('range gate')
grid on


ax3 = subplot(1,3,3)

plot(waveform);hold on
plot(W(3,:)./max(abs(W(3,:))))
title('Hs')
grid on
legend('waveform power','weighting vector')


linkaxes([ax1,ax2,ax3]);
ax1.XLim = [60,120];
ax1.YLim = [-1.1,1.1];

% test the whole procedure, with a 0.1 m increaeed waveheight, SAMOSA model
% the linear change should be predictable: Dparams = W* (S1-S2) should
% diagnose e.g. 0.1 m more waveheight, if correct

dparams = W*(SAMOSA(Pu-0.01,epoch+0.1,Hs+0.1)-SAMOSA(Pu,epoch,Hs))
% check! this is exactly what we want! Plus some errors from the linear
% approximation...

%% Samosa function at end of script
function [waveform] = SAMOSA(Pu,epoch,Hs)
    % read data from TUM processor
    data_dir = [getenv('HOME') '/TUDTUM/'];
    tum_data = load([data_dir '/sam_model_validation/samosa_model_validation_tu_delft_swh0-20.mat']);

    % constants
    CONST_A = 6378137; %Equatorial radius [m]
    CONST_B = 6356752.3142; % Polar radius [m]
    CONST_C = 299792458; % speed of light [m/s]

    % prepare input for TUD model
    dt = 1./tum_data.sensor_sets.B_r_Hz;

    x(1) = Pu;%tum_data.Pu; %Pu
    x(2) = epoch;%double(tum_data.epoch_ns); %epoch
    x(3) = Hs;%2; %Hs

    % some preliminaries
    % calc Re
    Re = sqrt(CONST_A.^2 * cos(tum_data.model_params.lat_rad).^2 + CONST_B.^2 .* sin(tum_data.model_params.lat_rad).^2);

    % calc orbit_slope
    if tum_data.model_params.ascending
        track_sign = -1;
    else
        track_sign = 1;
    end
    orbit_slope = track_sign*((CONST_A^2 - CONST_B^2)./(2*Re.^2)).*sin(2*tum_data.model_params.lat_rad) - (-double(tum_data.model_params.h_rate_m_per_s)./double(tum_data.model_params.Vs_m_per_s));

    % fill n data vector with parameters
    % n = [n;TN;NaN;Re;h;vt;OrbSlp;theta_p;theta_r;SmodID;EstPar];
    n = zeros(138,1);
    n(1:255) = (1:255);

    n(end) = 3;
    n(end-1) = 1; %SmodID;
    n(end-2) = tum_data.model_params.ksiy_rad; %roll mispointing angle (across-track)
    n(end-3) = tum_data.model_params.ksix_rad; %pitch mispointing angle (along-track)
    n(end-4) = orbit_slope;
    n(end-5) = tum_data.model_params.Vs_m_per_s;
    n(end-6) = tum_data.model_params.alt_m;
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

    DDAcf = DDA_ConfigFile('S3A','SAR');

    fid = fopen('alphap_table_SEN3_09_Nov_2017.txt'); LUT_AP = cell2mat(textscan(fid,'%f %f','HeaderLines',1)); fclose(fid);
    LUT_AP = griddedInterpolant(LUT_AP(:,1),LUT_AP(:,2),'linear','nearest');
    fun_x = (-19:.001:41.999)';
    LUT_F0 = (pi/4) * abs(fun_x).^0.5 .* (besseli(-1/4,(.5*fun_x).^2,1) + sign(fun_x).*besseli(1/4,(.5*fun_x).^2,1));
    LUT_F0(fun_x == 0) = (pi*2^(3/4)) / (4 * gamma(3/4));
    LUT_F0 = griddedInterpolant(fun_x,LUT_F0,'spline','nearest');
    LUT_F1 = (pi/8) * abs(fun_x).^1.5 .* ((besseli(1/4,(.5*fun_x).^2,1)-besseli(-3/4,(.5*fun_x).^2,1)) + sign(fun_x).*(besseli(-1/4,(.5*fun_x).^2,1)-besseli(3/4,(.5*fun_x).^2,1)));
    LUT_F1(fun_x == 0) = -(2^(3/4) * gamma(3/4)) / 4;
    LUT_F1 = griddedInterpolant(fun_x,LUT_F1,'spline','nearest');

    vt = double(tum_data.model_params.Vs_m_per_s);
    BAstck = tum_data.model_params.beam_ang_stack_rad - pi/2;
    DFstck  = (2*vt/DDAcf.lambda0) * sin(BAstck);     %Eq. 90
    dfa     = DDAcf.fp/DDAcf.Nb;                     %Eq. 61
    BeamIDX = double(unique(int32(DFstck/dfa)));

    stack_mask_start_stop = NaN(size(BeamIDX));

    waveform = RETRACKER.SAMOSA_DPM_v2_5_2_fun(x,n,DDAcf,LUT_AP,LUT_F0,LUT_F1,BeamIDX,stack_mask_start_stop,true);
end