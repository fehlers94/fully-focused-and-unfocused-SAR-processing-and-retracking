function [DATA,CS] = SAR_L1b_to_L2(SAT,FName,DOM,Retracker,SetRetr,DEM,MSK,IDXpoi)

%SAR_L1B_TO_L2 processes level 2 data from CryoSat & Sentinel-3A/B level 1b
%SAR data.

%% Settings
%General settings
LoadCommonSettings

%Computation settings
defval('SAT','CS')                                                                %Satellite mission from which data are processed ('CS'/'S3A'/'S3B')
defval('FName','2016/03/CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C001.DBL') %*.DBL/*.nc file that contains level 1b data
defval('DOM',[59 84;-74 -10])                                                     %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','Threshold')                                                   %Retracker to be used
defval('SetRetr',{})
for i=1:size(SetRetr,1), eval(sprintf('%s = %s;',SetRetr{i,1},SetRetr{i,2})); end
defval('SAMOSAimpl','S+')                                                         %SAMOSA implementation to be used ('S+' = 'SAMOSA+' as described by Salvatore Dinardo et al. 2018; 'DPM_v2_5_2' = Gommenginger et al. 2017)
% defval('SAMOSAimpl','DPM_v2_5_2')                                                 %SAMOSA implementation to be used ('S+' = 'SAMOSA+' as described by Salvatore Dinardo et al. 2018; 'DPM_v2_5_2' = Gommenginger et al. 2017)
defval('SolverAlgo','TRR')                                                        %Solver Algorithm to be applied in case Retracker is an analytical retracker ('LM' = Levenberg-Marquardt; 'TRR' = trust-region-reflective)
defval('max_waveform_energy',150)                                                 %Maximum allowed energy of NORMALIZED waveform (i.e., area below the curve)
defval('ApplyWfClass',false)                                                      %Apply surface type classification based on waveform parameters [true/false].
defval('MAfromSTR',true)                                                          %Obtain mispointing angle from star tracker data (for CS baseline C data only).
defval('MINbin',10);                                                              %Minimum bin index of interval in which retracking point is accepted for NON-zero padded waveform
defval('MAXbin',125);                                                             %Maximum bin index of interval in which retracking point is accepted for NON-zero padded waveform
defval('Wf_Norm_Aw',1)                                                            %Width (in gates) of the sliding window for waveform normalization
defval('Wf_TN_First',3);                                                          %First gate of the NON-zero padded waveform to estimate the amplitude of the thermal noise floor
defval('Wf_TN_Last',8);                                                           %Last gate of the NON-zero padded waveform to estimate the amplitude of the thermal noise floor
defval('UseDEM2CompInitRP',false)                                                 %Use DEM to compute initial retracking point [true/false].
defval('DEM',[])                                                                  %DEM used to compute initial retracking point
defval('UseLSmskInWfC',false)                                                     %Use high-resolution land/sea mask in waveform type classification [true/false].
defval('MSK',[])                                                                  %Land/sea mask used in waveform classification
defval('UseFLAGS',true)                                                           %Determines whether FLAGS are used, for debugging
defval('UseCOR',true)                                                             %Determines whether COR is used for geophysical corrections from L2 file, unhandy for debugging/testing
defval('UseSurfTypeSeaOnly',false)    
defval('UseOwnPath',false) 
defval('UseEmptyBeamAngle',false) 
defval('do_subwaveform_fit',false)
defval('fit_zero_doppler',false)

%whether or not to set all surface types to ocean when ApplyWfClass is false
if UseDEM2CompInitRP && isempty(DEM), error('Provide DEM as input!'), end
if UseLSmskInWfC && isempty(MSK), error('Provide land/sea mask as input!'), end

%Create default optimization options for SolverAlgo with the named
%parameters altered with the specified values
if ~any(strcmp(Retracker,{'OCOG','Threshold'}))
    switch SolverAlgo
        case 'LM'
            %Levenberg-marquardt
            options = optimoptions('lsqcurvefit','Display','off','Algorithm','levenberg-marquardt','StepTolerance',1E-2,'FunctionTolerance',1E-2,'OptimalityTolerance',1E-2);
        case 'TRR'
            %Trust-region-reflective
            options = optimoptions('lsqcurvefit','Display','off','Algorithm','trust-region-reflective','StepTolerance',1E-5,'FunctionTolerance',1E-5,'OptimalityTolerance',1E-5);
        otherwise
            error('Solver not recognized!')
    end
end

%Physical constants & instrument characteristics
DDAcf   = DDA_ConfigFile(SAT,'SAR');

%Remaining settings
defval('MakePlots',exist('IDXpoi','var')) %Make plots
if MakePlots; if ~isscalar(IDXpoi); MakePlots = false; end; end
defval('FntSz',14)                        %Set fontsize to be used in figures

is_s6 = contains(SAT, 'S6A');

%% Read & Crop Level 1b SAR data
switch SAT
    case {'CS','CryoSat'}
        %Read data
        if isstruct(FName)
            CS     = FName;
        else
            [~,~,Fext] = fileparts(FName);
            if isequal(Fext,'.DBL')
                if ~UseOwnPath
                    [~,CS] = Cryo_L1b_read(fullfile(PathDATA,'RadAlt','CryoSat','SIR_SAR_L1',FName));
                else
                    [~,CS] = Cryo_L1b_read(FName);
                end
                
                %Apply scaling
                CS.SAR.echo_scaling = CS.SAR.echo_scaling.*1E-9;
            else
%                 [~,CS] = Cryo_L1b_read_nc(fullfile(PathDATA,'RadAlt','CryoSat','SIR_SAR_L1',FName));
                [~,CS] = Cryo_L1b_read_nc(FName);
            end
        end
        
        %Regenerate DDAcf - update Baseline_ID dependent biases
        DDAcf = DDA_ConfigFile(SAT,'SAR',CS.GEO.Baseline_ID);

        %Wrap angle in degrees to [-180 180]
        CS.GEO.LON = wrapTo180(CS.GEO.LON);
        %If segment crosses -180 meridian, wrap angle in degrees to [0 360]
        if max(abs(diff(CS.GEO.LON(CS.SAR.N_averaged_echoes(:) > 0)))) > 350
            CS.GEO.LON = wrapTo360(CS.GEO.LON);
        end
        
        %Remove field AVG to save RAM
        CS = rmfield(CS,'AVG');
        
        %If field BeamAngle does not exists, add field with beam angles.
        %Here, we define the beam angle as the angle at which the surface
        %sample is seen w.r.t. the normal to the satellite velocity vector
        %(sometimes it is defined as the angle between the satellite
        %velocity vector and the surface sample)
        if ~isfield(CS.SAR,'BeamAngle')
            CS.SAR.('BeamAngle') = cell(numel(CS.GEO.LAT),1);
            for i=1:numel(CS.GEO.LAT)
                CS.SAR.BeamAngle{i} = linspace(CS.SAR.beam_param(8,i),CS.SAR.beam_param(9,i),CS.SAR.beam_param(12,i))';
            end
        end
        
        %Correct window delay for DORIS USO drift (pp 20 CryoSat Product Handbook)
        if isequal(CS.GEO.Baseline_ID,'C')
            CS.MEA.win_delay = CS.MEA.win_delay.*(CS.GEO.USO+1);
        end
        
        %Compute reference range
        CS.MEA.('ref_range') = 0.5*DDAcf.c*CS.MEA.win_delay(:);
        
        %Transform acquisition time to datenum format
        CS.('TIME') = datenum('2000','yyyy') + CS.GEO.TAI.days(:) + CS.GEO.TAI.secs(:)./86400 + CS.GEO.TAI.microsecs(:)./1e6./86400 - DDAcf.Biases.Datation/86400;
        
        %Get Roll, Yaw, and Pitch angles from star tracker data
        if (isequal(CS.GEO.Baseline_ID,'C') || isequal(CS.GEO.Baseline_ID,'OWN')) && MAfromSTR
            try
                tSTR          = datetime(CS.GEO.Start_Time./(24.*60.*60) + datenum('01-Jan-2000 00:00:00') - 1,'ConvertFrom','datenum');
                mpSTR         = READ_STR_Mispointing_Angles(tSTR);
                DUM           = interp1(mpSTR.Time,[mpSTR.Pitch mpSTR.Roll mpSTR.Yaw],CS.TIME+DDAcf.Biases.Datation/86400,'spline');
                IDXnan        = isnan(interp1(mpSTR.Time,mpSTR.Roll,CS.TIME+DDAcf.Biases.Datation/86400));
                DUM(IDXnan,:) = NaN;
                CS.GEO.Antenna_Bench_Roll(~isnan(DUM(:,2)))  = DUM(~isnan(DUM(:,2)),2);
                %Note that yaw angle is defined in different sign convention!
                CS.GEO.Antenna_Bench_Yaw(~isnan(DUM(:,3)))   = -DUM(~isnan(DUM(:,3)),3);
                %Note that pitch angle is defined in different sign convention!
                CS.GEO.Antenna_Bench_Pitch(~isnan(DUM(:,1))) = -DUM(~isnan(DUM(:,1)),1);
            catch exception
                fprintf('%s\n',exception.message)
            end
        end
        
        % %Salvatore Dinardo applies the pitch and roll biases proposed by Remko
        % %Scharroo. In doing so, he undoes the biases applied in producing baseline
        % %C (Main evolutions and expected quality improvements in Baseline C Level
        % %1b products - ARESYS / ESA, V1.3, https://wiki.services.eoportal.org/tiki-
        % %index.php?page=CryoSat+Technical+Notes)
        % CS.GEO.Antenna_Bench_Roll  = CS.GEO.Antenna_Bench_Roll  + 0.1062 - 0.086;
        % CS.GEO.Antenna_Bench_Pitch = CS.GEO.Antenna_Bench_Pitch + 0.0550 - 0.096;
        
        %Convert Roll, Yaw & Pitch angles to radians
        CS.GEO.Antenna_Bench_Roll  = deg2rad(CS.GEO.Antenna_Bench_Roll);
        CS.GEO.Antenna_Bench_Yaw   = deg2rad(CS.GEO.Antenna_Bench_Yaw);
        CS.GEO.Antenna_Bench_Pitch = deg2rad(CS.GEO.Antenna_Bench_Pitch);
        
        %The power echo	sample values are all scaled to	fit	between	0 and 65535.
        %The scaling factors can change	for	each waveform. To convert these back to
        %values in Watts the following equation should be used (CryoSat Product
        %Handbook - Baseline D 1.1, section 2.4.4):
        %Power in Watts	= scaled value * (scale	factor) * 2^scale power
        CS.SAR.data = bsxfun(@times,CS.SAR.data,reshape(CS.SAR.echo_scaling .* 2.^(CS.SAR.echo_scale_power),1,size(CS.SAR.data,2),size(CS.SAR.data,3)));
        
        %Convert amplitude of Gaussian that fits the range integrated power
        %of the single look echoes within a stack from dB to Watt
        CS.SAR.beam_param(3,:,:) = 10.^(CS.SAR.beam_param(3,:,:)./10);
        if UseFLAGS
            %Identify flagged data
            IDXfd = CS.GEO.MCD_FLAG.Block_Degraded(:) == 1 | CS.GEO.MCD_FLAG.Blank_Block(:) == 1 | ...
                CS.GEO.MCD_FLAG.Datation_Degraded(:) == 1 | CS.GEO.MCD_FLAG.Orbit_Propag_Err(:) == 1 | ...
                CS.GEO.MCD_FLAG.Echo_Saturation(:) == 1 | CS.GEO.MCD_FLAG.Other_Echo_Err(:) == 1 | ...
                CS.GEO.MCD_FLAG.Rx1_Err_SARin(:) == 1 | CS.GEO.MCD_FLAG.Rx2_Err_SARin(:) == 1 | ...
                CS.GEO.MCD_FLAG.Wind_Delay_Incon(:) == 1 | CS.GEO.MCD_FLAG.AGC_Incon(:) == 1 | ...
                CS.GEO.MCD_FLAG.TRK_ECHO_Err(:) == 1 | CS.GEO.MCD_FLAG.RX1_ECHO_Err(:) == 1 | ...
                CS.GEO.MCD_FLAG.RX2_ECHO_Err(:) == 1 | CS.GEO.MCD_FLAG.NPM_Incon(:) == 1 | ...
                CS.GEO.MCD_FLAG.Power_Scaling_Err(:) == 1;
        else
            IDXfd = false(1,size(CS.SAR.data,2));
        end

        if UseCOR
            %Compute sum of instrumental, propagation and geophysical corrections
            %Corrections to be applied in case (i) surface == open oceans or
            %semi‐enclosed seas OR (ii) surface == enclosed seas or lakes (SSB corr. is
            %still lacking). Results in SSH as defined in Eq. 1 (Dinardo et al. 2018)
            CorrST01  = CS.COR.ocean_equilibrium_tide+CS.COR.ocean_longperiod_tide+CS.COR.ocean_loading_tide+...
                CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide+CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.gim_ion+...
                CS.COR.dac;
            %Corrections to be applied in case the instantaneous SSHs (SSHi, see Eq. 4
            %Dinardo et al. 2018) is needed. Note that the SSB correction still has to
            %be applied.
            CorrST01i = CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.ocean_loading_tide+CS.COR.gim_ion+...
                CS.COR.solidearth_tide+0.468*CS.COR.geocentric_polar_tide;
            %Corrections to be applied in case (i) surface == continental ice OR (ii)
            %surface == land
            CorrST23  = CS.COR.ocean_loading_tide+CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide+...
                CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.gim_ion;
            warning('In case of sea ice, different set of corrections should be applied!')

            %Interpolate corrections to 20 Hz
            if isequal(CS.GEO.Baseline_ID,'C') || isequal(CS.GEO.Baseline_ID,'OWN')
                CS.COR.time_cor = CS.GEO.TAI.days*86400 + CS.GEO.TAI.secs + CS.GEO.TAI.microsecs/1e6 - DDAcf.Biases.Datation;
                CS.COR.time_cor = CS.COR.time_cor(1,:);
            end
            TIME_COR         = datenum('2000','yyyy') + CS.COR.time_cor'./86400;
            CS.('surf_type') = interp1(TIME_COR,CS.COR.surf_type',CS.TIME,'nearest','extrap');
            CorrST01         = interp1(TIME_COR,CorrST01',CS.TIME,'linear','extrap');  CorrST01(IDXfd)  = NaN;
            CorrST01i        = interp1(TIME_COR,CorrST01i',CS.TIME,'linear','extrap'); CorrST01i(IDXfd) = NaN;
            CorrST23         = interp1(TIME_COR,CorrST23',CS.TIME,'linear','extrap');  CorrST23(IDXfd)  = NaN;
            FlgC01           = isnan(CorrST01); FlgC01i = isnan(CorrST01i); FlgC23 = isnan(CorrST23);
        end
        %Load look-up-table alpha_p parameter (provided by S. Dinardo)
        LUT_AP    = load('alphaPTR_table_05_Aug_2014.txt');
        LUT_AP    = griddedInterpolant(LUT_AP(:,1),LUT_AP(:,2),'linear','nearest');
        
    case {'S3A','Sentinel-3A','S3B','Sentinel-3B', 'S6A'}
        %Read data
        if isstruct(FName)
            CS     = FName;
        else
            if contains(FName,'S3A_SR_1_SRA___') || contains(FName,'S3B_SR_1_SRA___')
                [~,CS] = S3_L1b_read(FName);
            elseif contains (FName, 'S6A_P4_1B_')
                CS = S6_L1bs_read(FName);
            else
                [~,CS] = S3_L1bs_read(FName,DOM);
            end
            if isempty(CS), DATA.('LAT') = []; return; end
        end

        %Wrap angle in degrees to [-180 180]
        CS.GEO.LON = wrapTo180(CS.GEO.LON);
        %If segment crosses -180 meridian, wrap angle in degrees to [0 360]
        if max(abs(diff(CS.GEO.LON(CS.SAR.N_averaged_echoes(:) > 0)))) > 350
            CS.GEO.LON = wrapTo360(CS.GEO.LON);
        end
        
        %Compute tangential velocities
        CS.GEO.V.('V') = sqrt(CS.GEO.V.Vx.^2 + CS.GEO.V.Vy.^2 + CS.GEO.V.Vz.^2);

        %Transform acquisition time to datenum format
        CS.('TIME') = datenum('2000','yyyy') + CS.GEO.TAI.days(:) + CS.GEO.TAI.secs(:)./86400;

        %Convert Roll, Yaw & Pitch angles to radians
        CS.GEO.Antenna_Bench_Roll  = deg2rad(CS.GEO.Antenna_Bench_Roll);
        CS.GEO.Antenna_Bench_Yaw   = deg2rad(CS.GEO.Antenna_Bench_Yaw);
        CS.GEO.Antenna_Bench_Pitch = deg2rad(CS.GEO.Antenna_Bench_Pitch);

        
        
        %Though the netcdf field name refers to these angles as beam
        %angles, the comment in the netcdf file specifies: these are
        %LOOKING angles (angle at which the surface position is seen with
        %respect to the nadir direction of the altimeter) and not beam
        %angles (angle between the normal to the satellite velocity vector
        %and the satellite to surface direction). The latter is needed as
        %input to the SAMOSA retracker. To obtain the beam angles, we need
        %to correct for the pitch.
        % CS.SAR.BeamAngle = CS.SAR.BeamAngle + mean([CS.SAR.beam_param.start_beam_ang_stack_l1bs_echo_sar_ku-CS.SAR.beam_param.start_look_angle_stack_l1bs_echo_sar_ku,...
        %     CS.SAR.beam_param.stop_beam_ang_stack_l1bs_echo_sar_ku-CS.SAR.beam_param.stop_look_angle_stack_l1bs_echo_sar_ku],2)';
        if ~is_s6
            try
                CS.SAR.BeamAngle = CS.SAR.BeamAngle + CS.GEO.Antenna_Bench_Pitch';
            catch
                CS.SAR.BeamAngle = CS.SAR.BeamAngle + CS.GEO.Antenna_Bench_Pitch;
            end
        end
        
        %Convert matrix with beam angles to cell array
        CS.SAR.BeamAngle = mat2cell(CS.SAR.BeamAngle, size(CS.SAR.BeamAngle,1), ones(1,size(CS.SAR.BeamAngle,2)));

        if UseEmptyBeamAngle
           for jj = 1:numel(CS.SAR.BeamAngle)
               CS.SAR.BeamAngle{jj} = zeros(0,0);
           end
        end
        %The echo is corrected for Doppler range effect, phase/power burst
        %calibration and GPRW effect . the echo is scaled using the
        %corrected AGC (agc_ku). Note that ''it is not possible to convert
        %power waveform in watts because it requires some instrumental
        %information from the industry who designed the altimeter and this
        %information is not available to the users. Nevertheless, the S3A
        %waveforms can be compared to each other, by applying the agc
        %values provided in the SRAL L1 products.''

        
        if UseFLAGS
            %Identify flagged data
            IDXfd = CS.GEO.FLAGS.time_corr_val(:) == 2 | ...                   %time correlation info invalid % CS.GEO.FLAGS.man_pres(:) ~= 0 | ... %ongoing manoeuvre
                CS.GEO.FLAGS.man_thrust(:) == 1 | ...                          %ongoing thrust % CS.GEO.FLAGS.man_plane(:) == 1 | ...               %out of plane
                CS.GEO.FLAGS.gnss_status(:) == 1 | ...                         %navigation message GNSS receiver not valid/available
                CS.GEO.FLAGS.loss_track(:) == 1;                               %loss of track
        else
            IDXfd = false(1,size(CS.SAR.data,2));
        end
        
        if UseCOR
            %Compute sum of instrumental, propagation and geophysical corrections
            %Corrections to be applied in case (i) surface == open oceans or
            %semi‐enclosed seas OR (ii) surface == enclosed seas or lakes.
            %Results in SSH as defined in Eq. 1 (Dinardo et al. 2018). Note
            %that the SSB correction still has to be applied as the one that is
            %in the L2 file has not been tuned for S3A/S3B and contains Jason-2
            %SSB solution (Sentinel-3 Product Notice – STM L2 Marine, 6
            %February 2019, Notice #5)
            CorrST01           = CS.COR.ocean_equilibrium_tide+CS.COR.ocean_longperiod_tide+CS.COR.ocean_loading_tide+...
                CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide+CS.COR.dry_trop+CS.COR.wet_trop+...
                CS.COR.inv_bar+CS.COR.hf_fluct_cor+CS.COR.iono_cor_gim_01_ku+...
                CS.COR.cog_cor_01+CS.COR.mod_instr_cor_range_01_ku;
            IDXnn              = knnsearch(CS.COR.TIME1Hz,CS.TIME);
            FlgC01             = abs(CS.TIME-CS.COR.TIME1Hz(IDXnn)) > 1/86400 | isnan(CorrST01(IDXnn)) | isnan(CS.COR.LAT1Hz(IDXnn));
            if nnz(~isnan(CorrST01)) >= 2
                CorrST01       = interp1(CS.COR.TIME1Hz(~isnan(CorrST01)),CorrST01(~isnan(CorrST01)),CS.TIME,'pchip',NaN);
            else
                CorrST01       = zeros(size(CS.TIME));
            end
            CorrST01(FlgC01)   = 0;
            %Corrections to be applied in case the instantaneous SSHs (SSHi, see Eq. 4
            %Dinardo et al. 2018) are needed. Note that the SSB correction
            %still has to be applied as the one that is in the L2 file has not
            %been tuned for S3A/S3B and contains Jason-2 SSB solution
            %(Sentinel-3 Product Notice – STM L2 Marine, 6 February 2019,
            %Notice #5)
            CorrST01i          = CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.ocean_loading_tide+...
                CS.COR.solidearth_tide+0.468*CS.COR.geocentric_polar_tide+CS.COR.iono_cor_gim_01_ku+...
                CS.COR.cog_cor_01+CS.COR.mod_instr_cor_range_01_ku;
            FlgC01i            = abs(CS.TIME-CS.COR.TIME1Hz(IDXnn)) > 1/86400 | isnan(CorrST01i(IDXnn)) | isnan(CS.COR.LAT1Hz(IDXnn));
            if nnz(~isnan(CorrST01i)) >= 2
                CorrST01i      = interp1(CS.COR.TIME1Hz(~isnan(CorrST01i)),CorrST01i(~isnan(CorrST01i)),CS.TIME,'pchip',NaN);
            else
                CorrST01i      = zeros(size(CS.TIME));
            end
            CorrST01i(FlgC01i) = 0;
            %Corrections to be applied in case (i) surface == continental ice OR (ii)
            %surface == land. Since the surface type parameter is ver coarse,
            %we copy the values of CorrST01i to CorrST23 to avoid that in the
            %end we loose useful data in the coastal areas.
            CorrST23           = CorrST01i;
            FlgC23             = FlgC01i;
            warning('In case of sea ice, different set of corrections should be applied!')
        end

        %Load look-up-table alpha_p parameter (provided by S. Dinardo)
        if is_s6
            if ~isfield(CS.GEO, 'FName')
                disp('warning: CS.GEO.FName not found. Run again L1b processing!')
            end
            
            if contains(Retracker,'FF') || fit_zero_doppler
                s6_lut_file = 'S6A_TEST_AUX_FLUT___00000000T000000_99999999T999999_0001.NC';
            elseif contains(CS.GEO.FName, 'F06')
                s6_lut_file = 'AUX_RLUT_S6A_003.nc';
            else
                s6_lut_file = 'AUX_RLUT_S6A_002.nc';
            end
            alpha_p_x = ncread(s6_lut_file, 'AlphaP_X');
            alpha_p_y = ncread(s6_lut_file, 'AlphaP_Y');
            LUT_AP   = griddedInterpolant(alpha_p_x,alpha_p_y,'linear','nearest');
        else
            fid      = fopen('alphap_table_SEN3_09_Nov_2017.txt'); LUT_AP = cell2mat(textscan(fid,'%f %f','HeaderLines',1)); fclose(fid);
            LUT_AP   = griddedInterpolant(LUT_AP(:,1),LUT_AP(:,2),'linear','nearest');
        end

    otherwise
        error('SAT: %s not recognized',SAT)
end

if ~isfield(CS.SAR, 'stack_mask_start_stop')
    CS.SAR.stack_mask_start_stop = NaN(1, numel(CS.TIME));
end

%Determine Zero-Padding Oversampling Factor and Waveform sampling interval.
%If os_ZP ~= 1, adjust reference bin.
DDAcf.Ns = size(CS.SAR.data,1);        %Nr of bins/samples in any waveform
os_ZP    = DDAcf.Ns/DDAcf.Np;          %Zero-Padding Oversampling Factor
RefBin   = (DDAcf.RefBin-1)*os_ZP + 1;
dt       = 1/(os_ZP*DDAcf.B);          %Waveform sampling interval [s]

% %Compute mispointing angles
% CS.GEO.('MPA') = acos(cos(CS.GEO.Antenna_Bench_Pitch)./sqrt(sum(cat(3,sin(CS.GEO.Antenna_Bench_Roll),sin(CS.GEO.Antenna_Bench_Yaw),cos(CS.GEO.Antenna_Bench_Pitch)).^2,3)));

%Compute the local radii of curvature of the Earth's surface (Maulik Jain,
%Improved sea level determination in the Arctic regions through development
%of tolerant altimetry retracking, Eq. 5.2, pp. 47)
% CS.GEO.('Re')  = sqrt(DDAcf.RefEll.SemimajorAxis^2*cosd(geocentricLatitude(CS.GEO.LAT,DDAcf.RefEll.Flattening)).^2 + DDAcf.RefEll.SemiminorAxis^2*sind(geocentricLatitude(CS.GEO.LAT,DDAcf.RefEll.Flattening)).^2);

% use Re calculus from SAMOSA DPM 2.5.2
CONST_B = 6356752.3142;
% CS.GEO.Re = sqrt(DDAcf.RefEll.SemimajorAxis^2 * cos(deg2rad(CS.GEO.LAT)).^2 + DDAcf.RefEll.SemiminorAxis^2 * sin(deg2rad(CS.GEO.LAT)).^2);
CS.GEO.Re = sqrt(DDAcf.RefEll.SemimajorAxis^2 * cos(deg2rad(CS.GEO.LAT)).^2 + CONST_B^2 * sin(deg2rad(CS.GEO.LAT)).^2);

%Compute slope of orbit
if CS.GEO.LAT(find(~IDXfd,1,'first')) < CS.GEO.LAT(find(~IDXfd,1,'last'))
    track_sign = -1; %ascending track
else
    track_sign = 1;  %descending track
end
CS.GEO.('orbit_slope') = track_sign*((DDAcf.RefEll.SemimajorAxis^2 - DDAcf.RefEll.SemiminorAxis^2)./(2*CS.GEO.Re.^2)).*sin(2*deg2rad(CS.GEO.LAT)) - (-CS.GEO.H_rate./CS.GEO.V.V);

%Normalize each waveform by max value of waveform
NORMfactor = 1./max(movmean(CS.SAR.data,Wf_Norm_Aw,1),[],1);
NORMdata   = NORMfactor .* CS.SAR.data;

% %Normalize amplitude of Gaussian that fits the range integrated power
% CS.SAR.beam_param(3,:,:) = NORMfactor .* CS.SAR.beam_param(3,:,:);

NORMdata_mean = size(mean(NORMdata, 2));
is_rmc_mode_active_s6 = (sum(NORMdata_mean(end-(length(NORMdata_mean)/2-30):end)) < 1) & is_s6;

%Estimate the normalized thermal noise
if is_s6
    Wf_TN_First = 13; %0-based
    Wf_TN_Last = 17; %0-based
end

CS.MEA.('TN') = mean(NORMdata(Wf_TN_First*os_ZP+1:os_ZP*Wf_TN_Last+os_ZP,:,:));

%Crop data to area of interest
if ischar(DOM)
    switch DOM
        case {'JakobsHavn','jakobshavn'}
            %Load polygon that outlines the JakobsHavn glacier
            load(fullfile(PathDATA,'Geography','jakobshavn_basin.mat'))
            IDX = inpolygon(CS.GEO.LON(:),CS.GEO.LAT(:),polyLon,polyLat);
        otherwise
            error('DOMain of interest not reckognized')
    end
else
    eps_lat = 0.002;
%     IDX = ingeoquad(CS.GEO.LAT(:),CS.GEO.LON(:),DOM(1,:),DOM(2,:));
    IDX = ingeoquad(CS.GEO.LAT(:),CS.GEO.LON(:),DOM(1,:)+[-eps_lat,eps_lat],DOM(2,:));
end

%% Data editing
%Exclude flagged data
IDX(IDXfd) = false;

%Identify waveforms for which power == 0 for all entries
IDX(squeeze(all(CS.SAR.data == 0,1))) = false;

%Energy NORMALIZED waveform < max_waveform_energy && power at first bins of
%NORMALIZED waveform should be at noise level
% IDX(squeeze(trapz(1:DDAcf.Ns,NORMdata,1) >= max_waveform_energy) | squeeze(any(NORMdata(1:MINbin*os_ZP,:,:) > .1))) = false;

%Return if no data remain
if ~any(IDX)
    DATA = struct('TIME',[],'LAT',[],'LON',[],'HEI',[],'SurfT',[]);
    return
end

%Select points of interest (by default entire track is retracked)
defval('IDXpoi',find(IDX)')

%% Classify waveforms
if ApplyWfClass
    %Preliminaries
    [CS.('n_ret'),CS.('Pu'),CS.('WFc')] = deal(nan(numel(CS.GEO.LAT),1));
    
    %Apply threshold retracker to obtain retracking point and Pu
    [CS.n_ret(IDXpoi),CS.Pu(IDXpoi)]    = RETRACKER.Threshold_mat(NORMdata(:,IDXpoi));
    
    %Compute range and backscatter coefficient (sigma0 / sigma-naught)
    CS.('range')   = CS.MEA.ref_range + (0.5*DDAcf.c*(CS.n_ret*dt - RefBin*dt)) - DDAcf.Biases.Range;
    switch SAT
        case {'CS','CryoSat'}
            CS.Pu          = CS.Pu ./ NORMfactor(:);
            CS.('sigma0')  = Compute_Sigma0(CS,DDAcf);
            CS.Pu          = CS.Pu .* NORMfactor(:);
        case {'S3A','Sentinel-3A','S3B','Sentinel-3B'}    
            %Surface Topography Mission (STM) SRAL/MWR L2 Algorithms
            %Definition, Accuracy and Specification, Section 2.15.3.3
            CS.Pu          = CS.Pu ./ NORMfactor(:);
            CS.('sigma0')  = CS.SAR.scale_factor_ku + 10*log10(CS.Pu) + DDAcf.Biases.sigma0;
            CS.Pu          = CS.Pu .* NORMfactor(:);            
        otherwise
            error('SAT: %s not recognized',SAT)
    end
    
    %Apply land/sea mask to find points wrongly classified as land
    if UseLSmskInWfC
        classIN              = CS.surf_type(IDXpoi(:));
        LorS                 = MSK(CS.GEO.LON(IDXpoi(:)),CS.GEO.LAT(IDXpoi(:)));
        IDXs                 = LorS == 0 & ismember(classIN,[2,3]);
        classIN(IDXs)        = 0;
        CS.surf_type(IDXpoi) = classIN;
        clear('classIN','LorS','IDXs')
    end
    
    %Classify waveforms
    CS.WFc(IDXpoi) = Classify_Waveforms(SAT,NORMdata(:,IDXpoi),CS.sigma0(IDXpoi)',CS.surf_type(IDXpoi)');
elseif UseSurfTypeSeaOnly
    CS.WFc(IDXpoi) = 0;
else
    CS.WFc(IDXpoi) = CS.surf_type(IDXpoi);
end
    
%% Preliminaries
%Declare arrays
[CS.('n_ret'),CS.('Pu'),CS.('SWH'),CS.('qual_flag')] = deal(nan(numel(CS.GEO.LAT),1));
[CS.('nu'),CS.('ExitF'),CS.('MF')]  = deal(nan(numel(CS.GEO.LAT),1));
CS.('PCorr')                        = nan(numel(CS.GEO.LAT),1);
BinIDs                              = (1:DDAcf.Ns)';
CS.('WDrecon')                      = nan(size(CS.SAR.data));

%Generate look-up-tables for fast evaluation of f0 and f1 terms
if strncmp(Retracker,'SAMOSA',6)
    x              = (-19:.001:41.999)';
    LUT_F0         = (pi/4) * abs(x).^0.5 .* (besseli(-1/4,(.5*x).^2,1) + sign(x).*besseli(1/4,(.5*x).^2,1));
    LUT_F0(x == 0) = (pi*2^(3/4)) / (4 * gamma(3/4));
    LUT_F0         = griddedInterpolant(x,LUT_F0,'spline','nearest');
    LUT_F1         = (pi/8) * abs(x).^1.5 .* ((besseli(1/4,(.5*x).^2,1)-besseli(-3/4,(.5*x).^2,1)) + sign(x).*(besseli(-1/4,(.5*x).^2,1)-besseli(3/4,(.5*x).^2,1)));
    LUT_F1(x == 0) = -(2^(3/4) * gamma(3/4)) / 4;
    LUT_F1         = griddedInterpolant(x,LUT_F1,'spline','nearest');
end

%Use DEM to compute initial retracking point
Xpk_ALL = nan(size(CS.GEO.H(:)));
if UseDEM2CompInitRP
    Xpk_ALL = ((CS.GEO.H(:) - DEM(CS.GEO.LON(:),CS.GEO.LAT(:))) - CorrST01(:) + DDAcf.Biases.Range - CS.MEA.ref_range(:) + (0.5*DDAcf.c*RefBin*dt))/(0.5*DDAcf.c*dt);
    Xpk_ALL(Xpk_ALL < 1 | Xpk_ALL > DDAcf.Ns) = NaN;
end

% profile on -detail builtin -history

%% Retrack waveforms
for i = IDXpoi
    %Copy normalized waveform i to vector WD
    WD        = NORMdata(:,i);
    
    if ~IDX(i) || any(isnan(WD))
        continue;
    end

    %Find largest peak in waveform and return associated bin index
    [Ypk,Xpk] = max(WD(MINbin*os_ZP:MAXbin*os_ZP));
    Xpk       = Xpk+(MINbin*os_ZP)-1;
    if ~isnan(Xpk_ALL(i)), Xpk = Xpk_ALL(i); end
    Wpk       = 10; %Just a number

    %Set initial values for SAMOSA retracker
    % IDXrm     = max([1 i-10]):min([numel(CS.GEO.LAT) i+9]);
    % [~,Xpk]   = max(prod(NORMdata(MINbin*os_ZP:MAXbin*os_ZP,IDXrm),2));
    % Xpk       = Xpk+(MINbin*os_ZP)-1;
    t0_0      = (Xpk - RefBin)*dt*1E9;
%     IDXrm     = max([1 i-20]):max([1 i-1]);
%     Pu0       = nanmean(CS.Pu(IDXrm));
%     SWH0      = nanmean(CS.SWH(IDXrm));
%     nu0       = nanmean(CS.nu(IDXrm));
    Pu0       = 1;
    SWH0      = 2;
    nu0       = 10;
    NrBins    = DDAcf.Ns;
    
    if is_s6
        DDAcf.fp = 1 / CS.MEA.PRI(i);
        DDAcf.BRI = CS.GEO.BRI;
    end
    
    %Apply retracking
    switch Retracker
        case 'BetaX'
            %X-parameter Beta-retracker with a linear trailing edge (Martin
            %et al., 1983)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.Beta5(WD,BinIDs,Ypk,Xpk,'linear');
        case 'BetaX_expTE'
            %X-parameter Beta-retracker with an exponential trailing edge
            %(Deng & Featherstone, 2006)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.Beta5(WD,BinIDs,Ypk,Xpk,'exponential');
        case 'Brown'
            %Brown Theoretical Ocean Model (Passaro et al., 2014)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.Brown(WD,BinIDs,Ypk,Xpk,CS.GEO.H(i),dt*1E9);
        case 'BrownHayne'
            %Brown-Hayne Theoretical Ocean Model (Gommenginger et al., 2011)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.BrownHayne(WD,BinIDs,Ypk,Xpk,CS.GEO.H(i),dt*1E9);
        case 'D2P'
            %D2P (Delay/Doppler Phase-monopulse Radar Altimeter) retracker
            %(Giles et al., 2007)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.D2P(WD,BinIDs,Ypk,Xpk,Wpk);
        case 'FunctionFit'
            %"Function Fit" retracker (Surface Topography Mission (STM)
            %SRAL/MWR L2 Algorithms Definition, Accuracy and Specification
            %[SD-03] [SD-07])
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.FunctionFit(WD,BinIDs,Ypk,Xpk,Wpk);
        case 'OCOG'
            %Offset Centre Of Gravity retracker (Wingham et al., 1986)
            [CS.n_ret(i),CS.Pu(i)] = RETRACKER.OCOG(WD,BinIDs');
        case 'Threshold'
            %Threshold retracker (Davis, 1997).
            [CS.n_ret(i),CS.Pu(i)] = RETRACKER.Threshold(WD,BinIDs',WD);
        case 'SAMOSA2'
            %SAMOSA2 retracker (Ray et al. 2015; Dinardo et al., 2018)
            [x0,lb,ub,XDATA,IDXnr,BeamIDX,stack_mask_start_stop] = RETRACKER.SAMOSA(WD,BinIDs,DDAcf,CS.MEA.TN(i),Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),CS.GEO.H(i),CS.GEO.V.V(i),CS.GEO.orbit_slope(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),1,0,CS.SAR.BeamAngle{i}',CS.SAR.stack_mask_start_stop(:,i));
        case 'SAMOSA2FF'
            %SAMOSA2 retracker (Ray et al. 2015; Dinardo et al., 2018)
            [x0,lb,ub,XDATA,IDXnr,BeamIDX,stack_mask_start_stop] = RETRACKER.SAMOSA(WD,BinIDs,DDAcf,CS.MEA.TN(i),Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),CS.GEO.H(i),CS.GEO.V.V(i),CS.GEO.orbit_slope(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),1,0,0,CS.SAR.stack_mask_start_stop(:,i));
        case 'SAMOSA3'
            %SAMOSA3 retracker (Ray et al. 2015; Dinardo et al., 2018)
            [x0,lb,ub,XDATA,IDXnr,BeamIDX,stack_mask_start_stop] = RETRACKER.SAMOSA(WD,BinIDs,DDAcf,CS.MEA.TN(i),Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),CS.GEO.H(i),CS.GEO.V.V(i),CS.GEO.orbit_slope(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),0,0,CS.SAR.BeamAngle{i}',CS.SAR.stack_mask_start_stop(:,i));
        case 'SAMOSA3FF'
            %SAMOSA3 retracker (Ray et al. 2015; Dinardo et al., 2018)
            [x0,lb,ub,XDATA,IDXnr,BeamIDX,stack_mask_start_stop] = RETRACKER.SAMOSA(WD,BinIDs,DDAcf,CS.MEA.TN(i),Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),CS.GEO.H(i),CS.GEO.V.V(i),CS.GEO.orbit_slope(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),0,0,0,CS.SAR.stack_mask_start_stop(:,i));
        otherwise
            error('Retracker not implemented!')
    end

    %Set function handle
    if strncmp(Retracker,'SAMOSA',6) && strcmp(SAMOSAimpl,'S+')
        fun = @(x,n) RETRACKER.SAMOSAfun(x,n,DDAcf,LUT_AP,LUT_F0,LUT_F1,BeamIDX, stack_mask_start_stop);
    elseif strncmp(Retracker,'SAMOSA',6) && strcmp(SAMOSAimpl,'DPM_v2_5_2')
        fun = @(x,n) RETRACKER.SAMOSA_DPM_v2_5_2_fun(x,n,DDAcf,LUT_AP,LUT_F0,LUT_F1,BeamIDX,stack_mask_start_stop,fit_zero_doppler);
    end
    
    %In case of analytical retrackers, solve the non-linear obs. eq.
    if ~any(strcmp(Retracker,{'OCOG','Threshold'}))
        %In case of the SAMOSA retracker, the 4 unknown parameters are not
        %solved simultaneously. For ocean waveforms, nu (the inverse of the
        %mean-square slope of the sea surface) is set to 0 and not
        %estimated. In case of lead waveforms, the SWH is set to 0 and not
        %estimated. For land contaminated waveforms, Dinardo et al. (2018)
        %apply a dual step retracking where first Pu, t0, and SWH are
        %estimated. Thereafter, the SWH is set 0 and Pu, t0, and nu are
        %(re-)estimated.
        
        %Solve for parameters in case of ocean waveform
        if any(CS.WFc(i) == [0 1 5])
            %Estimate Pu, t0, and SWH
            XDATA(end)           = 3;
            x0_tmp               = x0(1:3); lb_tmp = lb(1:3); ub_tmp = ub(1:3);
            if strcmp(SolverAlgo,'LM'), lb_tmp = []; ub_tmp = []; end
            
%             do_subwaveform_fit = true;
            % ignore first few and second halve of the range gates (RMC mode for Sentinel-6)
            if is_rmc_mode_active_s6     
                w = zeros(length(WD),1);
                w(11 * os_ZP + 1:132 * os_ZP + 1) = 1;
                fun = @(x,n) fun(x,n) .* w;
                WD = WD .* w;
            elseif do_subwaveform_fit
                [~,k_max] = max(WD);
                stopgates = 20;
                last_gate = k_max + stopgates;

                w = zeros(length(WD),1);
                w(1:min(length(WD),last_gate)) = 1;
                fun = @(x,n) fun(x,n) .* w;

                WD = WD .* w;
            end
            
            try
                [x,Resnorm,exitflag] = DARTutils.SolveNLLS(fun,x0_tmp,XDATA,WD,lb_tmp,ub_tmp,options);
            catch
                x = [NaN, NaN, NaN];
                Resnorm = NaN;
                exitflag = -1;
            end
            
            x(4)                 = 0;

            %Verify reason why algorithm is terminated
            CS.ExitF(i) = exitflag;
            if exitflag <= 0
                continue;
            end
            
            %Assess whether waveform is land contaminated based on quality
            %of fit
            CS.PCorr(i) = corr(WD,fun(x,XDATA),'type','Pearson');
            if CS.PCorr(i) <= 0.9 && ApplyWfClass
                CS.WFc(i) = 6;
                XDATA(end-8) = x(3);
            end
            %if 100*sqrt(Resnorm/DDAcf.Ns) > 4, CS.WFc(i) = 6; XDATA(end-8) = x(3); end
%             if CS.PCorr(i) <= 0.9
%                 CS.WFc(i)            = 6;
%                 [~,Xpk_wf]           = findpeaks(WD(MINbin*os_ZP:MAXbin*os_ZP),'MinPeakProminence',.1,'MINPEAKDISTANCE',5,'SortStr','none');
%                 Xpk_wf               = Xpk_wf+(MINbin*os_ZP)-1;
%                 [~,IDpk]             = min(abs(Xpk_wf-Xpk));
%                 %Set number of bins beyond Xpk
%                 BBP                  = 10;
%                 if IDpk < numel(Xpk_wf), [~,BBP] = min(WD(ceil(Xpk):Xpk_wf(IDpk+1))); end
%                 if floor(Xpk)+BBP > DDAcf.Ns, BBP = DDAcf.Ns - floor(Xpk); end
% %                 [~,Xpk_max]          = max(WD);
% %                 %Set number of bins beyond Xpk
% %                 BBP                  = 10;
% %                 if Xpk < Xpk_max, [~,BBP] = min(WD(ceil(Xpk):Xpk_max)); end
% %                 if floor(Xpk)+BBP > DDAcf.Ns, BBP = DDAcf.Ns - floor(Xpk); end
%                 IDnp                 = 1:floor(Xpk)+BBP; NrBins = numel(IDnp);
%                 XDATAsub             = [XDATA(IDnp);XDATA(DDAcf.Ns+1:end)];
%                 [x,Resnorm,exitflag] = DARTutils.SolveNLLS(fun,x0_tmp,XDATAsub,WD(IDnp),lb_tmp,ub_tmp,options);
%                 x(4)                 = 0;
% 
%                 %Verify reason why algorithm is terminated
%                 CS.ExitF(i) = exitflag;
%                 if exitflag <= 0, continue, end
%                 
%                 CS.PCorr(i) = corr(WD(IDnp),fun(x,XDATAsub),'type','Pearson');
%                 if CS.PCorr(i) <= 0.9, XDATA(end-8) = x(3); end
%             end
        end

        %Solve for parameters in case of lead or land contaminated waveform
        if any(CS.WFc(i) == [4 6])
%         if CS.WFc(i) == 4 || (CS.WFc(i) == 6 && CS.PCorr(i) <= 0.9)
            if CS.WFc(i) == 4, XDATA(end-8) = 0; end
            
            %(Re-)estimate Pu, t0, and nu (inverse of the mean-square slope of
            %the sea surface).
            XDATA(end)           = 4;
            x0_tmp               = x0([1 2 4]); lb_tmp = lb([1 2 4]); ub_tmp = ub([1 2 4]);
            if strcmp(SolverAlgo,'LM'), lb_tmp = []; ub_tmp = []; end
            [x,Resnorm,exitflag] = DARTutils.SolveNLLS(fun,x0_tmp,XDATA,WD,lb_tmp,ub_tmp,options);
            x                    = [x(1:2) XDATA(end-8) x(3)];

            %Verify again the reason why algorithm is terminated
            CS.ExitF(i) = exitflag;
            if exitflag <= 0
                continue;
            end

            %Assess whether waveform has indeed a specular shape
            CS.PCorr(i) = corr(WD,fun(x,XDATA),'type','Pearson');
            if CS.PCorr(i) <= 0.95, CS.WFc(i) = 99; CS.PCorr(i) = NaN; end
        end

        %Apply threshold retracker in case waveform belongs to classes 2,
        %3, or 99
        if any(CS.WFc(i) == [2 3 99])
            x           = nan(1,4);
            [x(2),x(1)] = RETRACKER.Threshold(WD,BinIDs',WD);
            x(2)        = (x(2) - RefBin)*(dt*1E9);
            Resnorm     = NaN;
        end
        
        %Copy retracking location [bins] and other estimated parameters to CS
        if any(strcmp(Retracker,{'Brown','BrownHayne'}))
            CS.n_ret(i) = x(IDXnr)/(dt*1E9);
            CS.Pu(i)    = x(IDXnr-1);
        elseif any(strcmp(Retracker,{'SAMOSA2','SAMOSA2FF','SAMOSA3','SAMOSA3FF'}))
            CS.n_ret(i) = (x(IDXnr)/(dt*1E9)) + RefBin;
            CS.Pu(i)    = x(IDXnr-1);
            CS.SWH(i)   = x(IDXnr+1);
            CS.nu(i)    = x(IDXnr+2);
            CS.MF(i)    = 100*sqrt(Resnorm/NrBins);
            CS.qual_flag(i)    = CS.MF(i) > 4;  % flag sample as 'bad' if MF is > 4 (qual_flag = true)
        else
            CS.n_ret(i) = x(IDXnr);
            CS.Pu(i)    = x(IDXnr-1);
        end
        
        %Reconstruct waveform
        CS.WDrecon(:,i) = fun(x,XDATA);
        
        %plot fitted waveform
%         figure; plot(WD); hold on; plot(CS.WDrecon(:,i)); legend('w_r', 'w_{opt}'); grid on; title(['index i=' num2str(i)]);xlabel('range gates k'); ylabel('normalised power');
    end
end
clear('WD','sub_n','sub_WD','Ypk','Xpk','IDXpks','i','j')

% profile viewer
% profile off

%% Data editing
CS.n_ret(CS.n_ret < MINbin*os_ZP | CS.n_ret > MAXbin*os_ZP) = NaN;

%% Compute corrected range (corrections include instrumental, propagation and geophysical corrections)
%Compute range Eqn 2.8‐1 CryoSat Product Handbook
CS.('range')      = CS.MEA.ref_range(:) + (0.5*DDAcf.c*(CS.n_ret*dt - RefBin*dt)) - DDAcf.Biases.Range;
CS.('range_SSHi') = CS.range;

if UseCOR
    %Apply propagation/geophysical corrections
    ST                     = CS.surf_type;
    CS.('SumCorrST')       = nan(size(CorrST01)); CS.('FlgC')       = true(size(CorrST01));
    CS.SumCorrST(ST <= 1)  = CorrST01(ST <= 1);   CS.FlgC(ST <= 1)  = FlgC01(ST <= 1);
    CS.SumCorrST(ST >= 2)  = CorrST23(ST >= 2);   CS.FlgC(ST >= 2)  = FlgC23(ST >= 2);
    CS.('SumCorrSTi')      = CS.SumCorrST;        CS.('FlgCi')      = CS.FlgC;
    CS.SumCorrSTi(ST <= 1) = CorrST01i(ST <= 1);  CS.FlgCi(ST <= 1) = FlgC01i(ST <= 1);
    CS.range               = CS.range + CS.SumCorrST;
    CS.range_SSHi          = CS.range_SSHi + CS.SumCorrSTi;
    clear('CorrST01','CorrST01i','CorrST23','ST','FlgC01','FlgC01i','FlgC23')
end

%Compute (instantaneous) (sea) surface heights relative to the reference ellipsoid
CS.('HEI')  = CS.GEO.H(:) - CS.range;
CS.('SSHi') = CS.GEO.H(:) - CS.range_SSHi;

%% Compute backscatter coefficient (sigma0 / sigma-naught)
switch SAT
    case {'CS','CryoSat'}
        CS.Pu         = CS.Pu ./ NORMfactor(:);
        CS.('sigma0') = Compute_Sigma0(CS,DDAcf);
        CS.Pu         = CS.Pu .* NORMfactor(:);
    case {'S3A','Sentinel-3A','S3B','Sentinel-3B'}
        %Surface Topography Mission (STM) SRAL/MWR L2 Algorithms
        %Definition, Accuracy and Specification, Section 2.15.3.3
        CS.Pu         = CS.Pu ./ NORMfactor(:);
        CS.('sigma0') = CS.SAR.scale_factor_ku + 10*log10(CS.Pu) + DDAcf.Biases.sigma0;
        CS.Pu         = CS.Pu .* NORMfactor(:);
    case {'S6A'}
        CS.Pu         = CS.Pu ./ NORMfactor(:);
        CS.('sigma0') = CS.SAR.scale_factor_ku + 10*log10(CS.Pu) + DDAcf.Biases.sigma0;
    otherwise
        error('SAT: %s not recognized',SAT)
end

%% Analyze output
if MakePlots
    for i = IDXpoi
        if all(isnan(CS.n_ret(i,:))), continue, end
    
        figure('Position',get(0,'Screensize'));

        %Select valid entries
        IDXvalid = ~isnan(CS.n_ret(i,:));

        %Plot observed/reconstructed waveform
        WD        = CS.SAR.data(:,IDXpoi);
        WDnorm    = NORMdata(:,i);
        StrYlabel = 'Power (Watts)';
        if ~any(strcmp(Retracker,{'OCOG','Threshold'})); WD = WDnorm; StrYlabel = 'Normalized power'; end
        subplot(2,3,1),plot(1:DDAcf.Ns,WD,'.-','LineWidth',2),hold on
        if ~all(isnan(CS.WDrecon(:,IDXpoi)))
            plot(1:DDAcf.Ns,CS.WDrecon(:,IDXpoi),'g-','LineWidth',2)
            plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:DDAcf.Ns,CS.WDrecon(:,IDXpoi),CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
            LegendEntries = {'Observed waveform','Reconstructed waveform','Retracking points'};
        else
            plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:DDAcf.Ns,WD,CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
            LegendEntries = {'Observed waveform','Retracking points'};
        end
        axis([0 DDAcf.Ns 0 max(WD)])
        title('Observed waveform','fontsize',FntSz)
        legend(LegendEntries,'box','off')
        xlabel('Bin','fontsize',FntSz)
        ylabel(StrYlabel,'fontsize',FntSz)
        set(gca,'fontsize',FntSz)

        WD        = CS.SAR.data(:,IDXpoi);
        subplot(2,3,2),plot(1:DDAcf.Ns,WD,'.-','LineWidth',2),hold on
        if ~all(isnan(CS.WDrecon(:,IDXpoi)))
            plot(1:DDAcf.Ns,fun(x,XDATA)/NORMfactor(IDXpoi),'g-','LineWidth',2)
            plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:DDAcf.Ns,fun(x,XDATA)/NORMfactor(IDXpoi),CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
            LegendEntries = {'Observed waveform','Reconstructed waveform','Retracking points'};
        else
            plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:DDAcf.Ns,WD,CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
            LegendEntries = {'Observed waveform','Retracking points'};
        end
        plot(1:DDAcf.Ns,ones(1,DDAcf.Ns)*CS.Pu(IDXpoi)/NORMfactor(IDXpoi),'k--');
        LegendEntries = horzcat(LegendEntries,{'Pu'});
        axis([0 DDAcf.Ns 0 max(WD)])
        title('Observed waveform','fontsize',FntSz)
        legend(LegendEntries,'box','off')
        xlabel('Bin','fontsize',FntSz)
        ylabel('Power (Watts)','fontsize',FntSz)
        set(gca,'fontsize',FntSz)

        %Plot coastline
        LOLA = Get_River_Border_data('GSHHS','coastlines','f',DOM(1,:),DOM(2,:));
        subplot(2,3,4),plot(LOLA(:,1),LOLA(:,2),'k-'),hold on
        plot(CS.GEO.LON(CS.GEO.LON~=0),CS.GEO.LAT(CS.GEO.LON~=0),'.')
        plot(CS.GEO.LON(i),CS.GEO.LAT(i),'ro','LineWidth',2)
        axis([minmax(LOLA(:,1)) minmax(LOLA(:,2))])
        xlabel('Longitude (deg)','fontsize',FntSz)
        ylabel('Latitude (deg)','fontsize',FntSz)
        set(gca,'fontsize',FntSz)

        %Plot DEM around point of interest including ground track
        if ~isempty(DEM)
        AZIM   = azimuth(CS.GEO.LAT(1:end-1),CS.GEO.LON(1:end-1),CS.GEO.LAT(2:end),CS.GEO.LON(2:end),DDAcf.RefEll);
        MaxDoN = km2deg(20,6371); %Maximum distance-of-nadir [km]
        DEMx   = DEM.GridVectors{1}; DEMy = DEM.GridVectors{2};
        [~,Ix] = min(abs(DEMx-CS.GEO.LON(i))); [~,Iy] = min(abs(DEMy-CS.GEO.LAT(i))); 
        MDidx  = ceil(MaxDoN/mode(diff(DEMx)));
        subplot(2,3,5),imagesc(DEMx(Ix+(-MDidx:MDidx)),DEMy(Iy+(-MDidx:MDidx)),DEM.Values((Iy+(-MDidx:MDidx)),(Ix+(-MDidx:MDidx)))),hold on
        axis square, colorbar, set(gca,'YDir','normal')
        set(get(findobj(gcf,'tag','Colorbar'),'ylabel'),'String','meters','FontSize',FntSz);
        %Plot coastline
        plot(LOLA(:,1),LOLA(:,2),'-','Color',[.5 .5 .5]),hold on
        %Plot direction of flight
        quiver(CS.GEO.LON(i),CS.GEO.LAT(i),(CS.GEO.LON(i+1)-CS.GEO.LON(i)),(CS.GEO.LAT(i+1)-CS.GEO.LAT(i)),20,'Color','k','LineWidth',2,'MaxHeadSize',0.8);
        %Plot across-track points
        [lat_ac,lon_ac] = reckon(CS.GEO.LAT(i),CS.GEO.LON(i),deg2km(-MaxDoN:mode(diff(DEMx))/2:MaxDoN,6371)*1000,AZIM(i)-90,DDAcf.RefEll);
        plot(lon_ac,lat_ac,'k.-')
        %Plot identified scatterers
        plot(CS.GEO.LON(i),CS.GEO.LAT(i),'o','Color',[.5 .5 .5],'LineWidth',2,'MarkerSize',12)
        title('DEM','fontsize',FntSz)
        xlabel('Longitude (deg)','fontsize',FntSz)
        ylabel('Latitude (deg)','fontsize',FntSz)
        set(gca,'fontsize',FntSz)
        end

        % subplot(2,3,6),plot(CS.GEO.LON(CS.GEO.LON~=0),CS.GEO.LAT(CS.GEO.LON~=0),'.'),hold on
        % plot(CS.GEO.LON(i),CS.GEO.LAT(i),'ro','LineWidth',2)
        % axis([minmax(DEMx(Ix+(-MDidx:MDidx))) minmax(DEMy(Iy+(-MDidx:MDidx)))])
        % plot_google_map
        % xlabel('Longitude (deg)','fontsize',FntSz)
        % ylabel('Latitude (deg)','fontsize',FntSz)
        % set(gca,'fontsize',FntSz)

        %Set background color to white
        set(gcf, 'Color',[1 1 1])
    end
end

%% Save output
DATA              = struct;
% IDX               = ~isnan(CS.HEI);
DATA.('TIME')     = CS.TIME(IDX);
DATA.('LAT')      = CS.GEO.LAT(IDX);
DATA.('LON')      = CS.GEO.LON(IDX);
DATA.('HEI')      = single(CS.HEI(IDX));
DATA.('SSHi')     = single(CS.SSHi(IDX));
DATA.('ALT')     = CS.GEO.H(IDX);
DATA.('range') = CS.range(IDX);
if UseCOR
    DATA.('sumCORR')  = single(CS.SumCorrST(IDX));
    DATA.('FlgC')     = CS.FlgC(IDX);
    DATA.('sumCORRi') = single(CS.SumCorrSTi(IDX));
    DATA.('FlgCi')    = CS.FlgCi(IDX);
end
DATA.('sigma0')   = single(CS.sigma0(IDX));
DATA.('SWH')      = single(CS.SWH(IDX));
if ~UseSurfTypeSeaOnly
DATA.('SurfT')    = int8(CS.surf_type(IDX));
end
DATA.('WFc')      = int8(CS.WFc(IDX));
DATA.('nu')       = single(CS.nu(IDX));
DATA.('ExitF')    = int8(CS.ExitF(IDX));
DATA.('MF')       = single(CS.MF(IDX));
if isfield(CS,'qual_flag')
    DATA.('qual_flag')= single(CS.qual_flag(IDX));
end
DATA.('PCorr')    = single(CS.PCorr(IDX));
DATA.('IDXpoi')   = IDXpoi;

end
