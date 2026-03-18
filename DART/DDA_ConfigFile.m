function DDAcf = DDA_ConfigFile(SAT,MODE,Baseline_ID)

%DDA_CONFIGFILE includes universal constants, sensor parameters, and
%control data values

defval('Baseline_ID','D')

%% Definition of Universal Constants
DDAcf.c         = 299792458;                       %Light velocity [m/s]

%% Definition of reference ellipsoid, sensor parameters, and control data values
switch SAT
    case {'CS','CryoSat'}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Reference ellipsoid
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DDAcf.RefEll    = referenceEllipsoid('wgs84','m');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Sensor Parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DDAcf.fc        = 13.575E9;                %Radio frequency [Hz]
        DDAcf.B         = 320E6;                   %Receiver bandwidth [Hz]
        DDAcf.theta_x   = deg2rad(1.06);           %Along-track beamwidth (Scagliola 2013, CryoSat footprints - Aresys Technical Note, Table 1)
        DDAcf.theta_y   = deg2rad(1.1992);         %Cross-track beamwidth (Scagliola 2013, CryoSat footprints - Aresys Technical Note, Table 1)
        % DDAcf.theta_x   = deg2rad(1.10);           %Along-track beamwidth (value proposed by R. Scharroo used by Salvatore Dinardo)
        % DDAcf.theta_y   = deg2rad(1.22);           %Cross-track beamwidth (value proposed by R. Scharroo used by Salvatore Dinardo)
        switch MODE
            case {'LRM','SAR'}
                DDAcf.Np        = 128;             %Nr of bins/samples in any waveform (NO zero-padding)
            case 'SIN'
                DDAcf.Np        = 512;             %Nr of bins/samples in any waveform (NO zero-padding)
            otherwise
                error('MODE: %s not recognized',MODE)
        end
        DDAcf.fp        = 1/(4400*(1/80E6));       %Pulse repetition frequency [Hz]
        DDAcf.tau_u     = 44.8E-6;                 %Reception Chirp Duration / Usable pulse length (Scagliola, M. and M. Fornari (2016), CryoSat Characterization for FBR users)
        DDAcf.Nb        = 64;                      %Number of pulses per burst
        DDAcf.BRI       = (943437)*(1/80e6);       %Burst Repetition Interval [s]
        DDAcf.shx       = 1;                       %Antenna gain along-track shape factor
        DDAcf.shy       = 1;                       %Antenna gain across-track shape factor
        DDAcf.s         = DDAcf.B/DDAcf.tau_u;     %Chirp slope [Hz/s]
        % DDAcf.alpha_PF  = 0.47356;                 %Alpha power factor. Required to scale Pu
        DDAcf.alpha_PF  = 1/(0.886*sqrt(2*pi));    %Alpha power factor. Required to scale Pu (Gommenginger et al. 2017)
        DDAcf.PRI       = 1/DDAcf.fp;              %Pulse repitition interval
        DDAcf.lambda0   = DDAcf.c/DDAcf.fc;        %Radar wavelength [m]
        DDAcf.G0        = 10^(4.2814);             %Antenna Gain at Boresight (taken from CryoSat Characterization for FBR users v2.0 (see Eq. 8))
        DDAcf.PTR_width = 2.819E-9;                %3dB Range Point Target Response Temporal Width [s]
        DDAcf.tau_B     = DDAcf.Nb/DDAcf.fp;       %Burst Length [s]
        DDAcf.ABL       = 1.1676;                  %Antenna baseline length [m]
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %L1b processing Control Parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calibration (source: Scagliola, M. and M. Fornari (2016), CryoSat
        %Characterization for FBR users):
        DDAcf.ADC_MULT      = 1000;                %ADC multiplier
        DDAcf.p2p_amp       = ncread('cs_users_characterization_C002.nc','cal1_p2p_amplitude_sar');
        DDAcf.p2p_phs       = ncread('cs_users_characterization_C002.nc','cal1_p2p_phase_sar');
        DDAcf.lpf_sar       = ncread('cs_users_characterization_C002.nc','cal2_lpf_sar');
        DDAcf.lpf_sarin_rx1 = ncread('cs_users_characterization_C002.nc','cal2_lpf_sarin_rx1');
        DDAcf.p2p_amp_rx1   = ncread('cs_users_characterization_C002.nc','cal1_p2p_amplitude_sarin_rx1');
        DDAcf.p2p_phs_rx1   = ncread('cs_users_characterization_C002.nc','cal1_p2p_phase_sarin_rx1');
        DDAcf.lpf_sarin_rx2 = ncread('cs_users_characterization_C002.nc','cal2_lpf_sarin_rx2');
        DDAcf.p2p_amp_rx2   = ncread('cs_users_characterization_C002.nc','cal1_p2p_amplitude_sarin_rx2');
        DDAcf.p2p_phs_rx2   = ncread('cs_users_characterization_C002.nc','cal1_p2p_phase_sarin_rx2');
        
        %Apply Azimuth Processing:
        DDAcf.ApplyWindow   = true;                %Apply weighting before beam forming takes place to minimize the impact of side-lobe effects in the Doppler/azimuth PTR (true/false)
        DDAcf.WindowMode   = 'burst'; % or 'aperture'
        DDAcf.ApplyWindowRange = false;          
        DDAcf.Window        = 'hamming';           %Window funtion ('Hamming'/'Hanning')
        DDAcf.ThSurfVar     = 40;                  %Threshold surface variablity. If surface variability > threshold exact beam steering is applied, else approximated beam steering
        
        % %Generate stack
        % theta                    = 0.0005;
        % n                        = 20;
        % DDAcf.Uplt           = pi/2+(n*theta);
        % DDAcf.Lwlt           = pi/2-(n*theta);
        % DDAcf.beam_selection = true;
        
        %Geometry corrections
        DDAcf.os_ZP         = 2;                   %Interpolation/oversampling factor. If == 2, zero-padding is applied prior to range FFT
        
        %Along-track antenna pattern compensation (source: Scagliola, M. Dinardo,
        %S. and M. Fornari (2015), An extended analysis of along-track antenna
        %pattern compensation for SAR altimetry):
        DDAcf.ACP           = false;
        DDAcf.gamma1        = 0.0116;              %gamma1 is a function of the -3 dB width of the along-track antenna beam.

        %Set reference bin for NON-zero padded waveform
        DDAcf.RefBin        = (DDAcf.Np/2)+1;      %+1 in Matlab because first bin has index 1

        DDAcf.nTs=2; %time offset nTs
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Known biases in CS data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch Baseline_ID
            case {'C'}
                %Known biases in CryoSat Baseline C Level1b products
                %https://wiki.services.eoportal.org/tiki-index.php?page=CryoSat+Technical+Notes
                % ARESYS - Known biases in CryoSat Level1b products
                % ARESYS - Main evolutions and expected quality improvements in Baseline C Level 1b products
                switch MODE
                    case 'LRM'
                        DDAcf.Biases.Datation = 0;         %[s] TO BE ASSESSED
                        DDAcf.Biases.Range    = -22E-3;    %[m]
                        DDAcf.Biases.sigma0   = 0;         %[dB] UNKNOWN
                    case 'SAR'
                        DDAcf.Biases.Datation = -0.029E-3; %[s]
                        DDAcf.Biases.Range    = 47E-3;     %[m]
                        DDAcf.Biases.sigma0   = -3.04;     %[dB] (Salvatore Dinardo)
                    case 'SIN'
                        DDAcf.Biases.Datation = -0.022E-3; %[s]
                        DDAcf.Biases.Range    = 17E-3;     %[m]
                        DDAcf.Biases.Roll     = 0.0069;    %[deg] (Garcia-Mondejar et al., 2017)
                        DDAcf.Biases.sigma0   = -3.04;     %[dB] (Salvatore Dinardo)
                    otherwise
                        error('MODE: %s not recognized',MODE)
                end
            case {'D','OWN'}
                DDAcf.Biases.Datation = 0;                 %[s]
                DDAcf.Biases.Range    = 0;                 %[s]
                DDAcf.Biases.Roll     = 0;                 %[s]
                DDAcf.Biases.sigma0   = 0;                 %[s]
            otherwise
                error('Baseline_ID: %s not recognized',Baseline_ID)
        end
        
         % file handling
        DDAcf.ncgrp_l1a = '/';
        DDAcf.ncgrp_l1b_l2 = '/';
        DDAcf.var_mapping_l1b = struct('range_cor_com','range_cor_com');
        DDAcf.var_mapping_l2 = struct('lat','lat_20_ku','lon','lon_20_ku','alt','alt_20_ku','swh','swh_ocean_20_ku','range','range_ocean_20_ku','range_cor','','dist2coast', '');
        
    case {'S3A','Sentinel-3A','S3B','Sentinel-3B'}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Reference ellipsoid
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DDAcf.RefEll    = referenceEllipsoid('wgs84','m');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Sensor Parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DDAcf.fc        = 13.575E9;                %Radio frequency [Hz]
        DDAcf.B         = 320E6;                   %Receiver bandwidth [Hz]
        % DDAcf.theta_x   = deg2rad(1.35);           %Along-track beamwidth [rad] (./RadAlt/Sentinel3A/SR_1_SRA/2016/12/S3A_SR___CHDNAX_20000101T000000_20991231T235959_20160603T120000___________________MPC_O_AL_002.SEN3/S3A_SR_CCDB_CHAR_NOM.20160526140700_1.nc)
        % DDAcf.theta_y   = deg2rad(1.35);           %Cross-track beamwidth [rad] (./RadAlt/Sentinel3A/SR_1_SRA/2016/12/S3A_SR___CHDNAX_20000101T000000_20991231T235959_20160603T120000___________________MPC_O_AL_002.SEN3/S3A_SR_CCDB_CHAR_NOM.20160526140700_1.nc)
        DDAcf.theta_x   = deg2rad(1.338);          %3dB beamwidth along track [rad], from S.Dinardo (mail 29.05.20)
        DDAcf.theta_y   = deg2rad(1.338);          %3dB beamwidth across track [rad], from S.Dinardo (mail 29.05.20)
        switch MODE
            case 'SAR'
                DDAcf.Np        = 128;             %Nr of bins/samples in any waveform (NO zero-padding)
            otherwise
                error('MODE: %s not recognized',MODE)
        end
        DDAcf.fp        = 1/(4488/80E6);           %Pulse repetition frequency [Hz] (Determined by reverse engineering, confirmed by W. Smith)
        DDAcf.tau_u     = 44.8E-6;                 %Reception Chirp Duration / Usable pulse length (Gommenginger et al. 2017)
        DDAcf.Nb        = 64;                      %Number of pulses per burst (Gommenginger et al. 2017)
        DDAcf.BRI       = 1018710*1/80E6;          %Burst Repetition Interval [s] (Determined by reverse engineering, confirmed by W. Smith)
        DDAcf.shx       = 1;                       %Antenna gain along-track shape factor (Gommenginger et al. 2017)
        DDAcf.shy       = 1;                       %Antenna gain across-track shape factor (Gommenginger et al. 2017)
        DDAcf.s         = DDAcf.B/DDAcf.tau_u;     %Chirp slope [Hz/s]
        % DDAcf.alpha_PF  = 1/(0.886*sqrt(2*pi));    %Alpha power factor. Required to scale Pu (Gommenginger et al. 2017)
        DDAcf.alpha_PF  = 0.5;                     %As in Dinardo (2020), DPM 2.5.2: 1 / (0.886 * np.sqrt(2 * np.pi)) # ~0.45, PTR Gaussian approximation coefficient value 
        DDAcf.PRI       = 1/DDAcf.fp;              %Pulse repitition interval
        DDAcf.lambda0   = DDAcf.c/DDAcf.fc;        %Radar wavelength [m]
        DDAcf.PTR_width = 0.4162/DDAcf.c;          %3dB Range Point Target Response Temporal Width [s] (S3 - A SRAL Cyclic Performance Report Cycle No. 038)
        DDAcf.tau_B     = DDAcf.Nb/DDAcf.fp;       %Burst Length [s]
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %L1b processing Control Parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calibration:
        
        %Apply Azimuth Processing:
        DDAcf.ApplyWindow   = false;                %Apply weighting before beam forming takes place to minimize the impact of side-lobe effects in the Doppler/azimuth PTR (true/false)
        DDAcf.WindowMode   = 'burst'; % or 'aperture'
        DDAcf.ApplyWindowRange = false;          
        DDAcf.Window        = 'hanning';           %Window funtion ('Hamming'/'Hanning')
        DDAcf.ThSurfVar     = 40;                  %Threshold surface variablity. If surface variability > threshold exact beam steering is applied, else approximated beam steering
        
        %Geometry corrections
        DDAcf.os_ZP         = 2;                   %Interpolation/oversampling factor. If == 2, zero-padding is applied prior to range FFT
        
        %Set reference bin for NON-zero padded waveform
        DDAcf.RefBin        = 43+1;                %The tracking point is gate 44 i.e. bin index 43 for Ku-band; 43+1 in Matlab because first bin has index 1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Known biases in S3 data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DDAcf.Biases.Datation   = 0; %[s]
        DDAcf.Biases.Range      = 0; %[m]
        if any(strcmp(SAT,{'S3A','Sentinel-3A'}))
            DDAcf.Biases.sigma0 = -18.96;          %[dB] (Sentinel-3 Product Notice – STM L2 Marine, 6 February 2019, Notice #3)
        elseif any(strcmp(SAT,{'S3B','Sentinel-3B'}))
            DDAcf.Biases.sigma0 = -19.17;          %[dB] (Sentinel-3 Product Notice – STM L2 Marine, 6 February 2019, Notice #3)
        end
        
        DDAcf.nTs=0.3; %time offset nTs
        
        % file handling
        DDAcf.ncgrp_l1a = '/';
        DDAcf.ncgrp_l1b_l2 = '/';
        DDAcf.var_mapping_l1b = struct('range_cor_com','range_cor_com');
        DDAcf.var_mapping_l2 = struct('lat','lat_l1bs_echo_sar_ku','lon','lon_l1bs_echo_sar_ku','alt','alt_20_ku','swh','swh_ocean_20_ku','range','range_ocean_20_ku','range_cor','net_instr_cor_range_ocean','dist2coast', 'dist_coast_20_ku');
        
    case {'S6A','Sentinel-6'}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Reference ellipsoid
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DDAcf.RefEll    = referenceEllipsoid('wgs84','m');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Sensor Parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DDAcf.fc        = 13.575E9;                %Radio frequency [Hz]
        DDAcf.B         = 395E6;                   %Receiver bandwidth [Hz]
        DDAcf.Bt        = 320E6;                   %Transmitter bandwidth [Hz]
        DDAcf.theta_x   = deg2rad(1.33);           %3dB beamwidth along track [rad], 1.33 value: pers. comm. Remko Scharroo
        DDAcf.theta_y   = deg2rad(1.33);           %3dB beamwidth across track [rad], 1.33 value: pers. comm. Remko Scharroo
        DDAcf.alt_nom   = 1336e3;                  %nominal altitude [m]
        DDAcf.v_nom   = 7200;                      %nominal velocity [m/s]
        DDAcf.res_a=0.5; % approximate along-track resolution [m]
        switch MODE
            case 'SAR'
                DDAcf.Np        = 256;             %Nr of bins/samples in any waveform (NO zero-padding)
            otherwise
                error('MODE: %s not recognized',MODE)
        end
%         DDAcf.fp        = 9.1946E+3;           %Pulse repetition frequency [Hz], this is dynamically read from the S6 file
        DDAcf.tau_u     = 32E-6;                 %Reception Chirp Duration / Usable pulse length
        DDAcf.Nb        = 64;                      %Number of pulses per burst
%         DDAcf.BRI       = 0.0072;          %Burst Repetition Interval [s], likely not correct, check this
        DDAcf.shx       = 1;                       %Antenna gain along-track shape factor
        DDAcf.shy       = 1;                       %Antenna gain across-track shape factor
        DDAcf.s         = DDAcf.Bt/DDAcf.tau_u;     %Chirp slope [Hz/s]
        DDAcf.alpha_PF  = 0.55;                     % as documented in PDAP_AUX_CNF2_HR_S6A_010.xml (S6VT)
%         DDAcf.PRI       = 1/DDAcf.fp;              %Pulse repitition interval, this is dynamically read from the S6 file
        DDAcf.lambda0   = DDAcf.c/DDAcf.fc;        %Radar wavelength [m]
%         DDAcf.PTR_width = 0.4162/DDAcf.c;          %3dB Range Point Target Response Temporal Width [s]
%         DDAcf.tau_B     = DDAcf.Nb/DDAcf.fp;       %Burst Length [s]
        
        % Retracker Parameters
        DDAcf.Enable_Slope_Effect_Flag = false;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %L1b processing Control Parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calibration:
        
        %Apply Azimuth Processing:
        DDAcf.ApplyWindow   = true;                %Apply weighting before beam forming takes place to minimize the impact of side-lobe effects in the Doppler/azimuth PTR (true/false)
        DDAcf.WindowMode   = 'aperture'; % or 'aperture'
        DDAcf.ApplyWindowRange = true;          
        DDAcf.Window        = 'hamming';           %Window funtion ('Hamming'/'Hanning')
        DDAcf.ThSurfVar     = 40;                  %Threshold surface variablity. If surface variability > threshold exact beam steering is applied, else approximated beam steering
        
        %Geometry corrections
        DDAcf.os_ZP         = 2;                   %Interpolation/oversampling factor. If == 2, zero-padding is applied prior to range FFT
        
        %Set reference bin for NON-zero padded waveform
        DDAcf.RefBin        = 129;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Known biases in S6 data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DDAcf.Biases.Datation   = 0; %[s]
        DDAcf.Biases.Range      = 0; %[m]
        DDAcf.Biases.sigma0 = 0;     %[dB]
        
        DDAcf.nTs=-0.5; %time offset nTs
        
        % file handling
        DDAcf.ncgrp_l1a = 'data_140/ku/';
        DDAcf.ncgrp_l1b_l2 = 'data_20/ku/';
        DDAcf.var_mapping_l1b = struct('range_cor_com','range_cor_com');
        DDAcf.var_mapping_l2 = struct('lat','latitude','lon','longitude','alt','altitude','swh','swh_ocean','sigma0','sig0_ocean','time','time','range','range_ocean','range_cor','net_instr_cor_range_ocean','dist2coast', 'distance_to_coast');
        
    otherwise
        error('SAT: %s not recognized',SAT)
end

DDAcf.var_mapping_kuststrook = struct('lat','latitude','lon','longitude','swh','wave_height_hm0','range','waverlevel');

end