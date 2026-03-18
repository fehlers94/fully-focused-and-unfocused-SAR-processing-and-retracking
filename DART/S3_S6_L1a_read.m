function CS1a = S3_S6_L1a_read(FName,startLoc,count)

nc_grps = {'/'};

LoadCommonSettings

if ~exist('FName','var')
    FName='S3A_SR_1_SRA____20180227T200232_20180227T205254_20180325T110545_3021_028_213______MAR_O_NT_003.SEN3';
    
    %Set paths
    if contains(FName,'S3A')
        s3_str = 'Sentinel3A';
    elseif contains(FName,'S3B')
        s3_str = 'Sentinel3B';
    end
    
    PathL1a = fullfile(PathDATA,'RadAlt',s3_str,'SR_1_SRA');
    FName = fullfile(PathL1a,FName,'measurement.nc');
end

mission = mission_from_fname(FName);
is_s3 = contains(mission,'S3');
is_s6 = contains(mission,'S6');

if is_s6
    nc_grps = {'data_140/ku/', 'global/ku/'};
end


for i_grp=1:length(nc_grps)
    nc_grp = nc_grps{i_grp};

    %Read netCDF file info
    finfo = ncinfo(FName, nc_grp);
    % ncid = netcdf.open(FName,'NC_NOWRITE');

    %Get number of variables in the netCDF file.
    nvars = length(finfo.Variables);

    %Create list of variable names associated to SAR ku echoes
    varname = cell(nvars,1);
    varsize = cell(nvars,1);
    
    for i=1:nvars
        varname{i} = finfo.Variables(i).Name;
        varsize{i} = finfo.Variables(i).Size;

        if contains(varname{i}, 'lat')
            latvar = varname{i};
        end
        if is_s3
            if contains(varname{i}, 'i_meas_ku')
                ivar = varname{i};
            end       
        elseif is_s6
            if contains(varname{i}, 'i_')
                ivar = varname{i};
            end
        end
    end
    qvar = replace(ivar,'i_','q_');

    %Read all variables. replace fill values by NaNs
    for i=1:numel(varname)
        if is_s6 || (is_s3 && contains(varname{i},'echo_sar_ku'))
            try FV = ncreadatt(FName,[nc_grp varname{i}],'_FillValue'); catch, FV = NaN; end

            try
                if length(varsize{i}) > 1
                    vs = length(varsize{i})-1;
                    DUM = ncread(FName,[nc_grp varname{i}],[ones(1,vs) startLoc],[inf*ones(1,vs) count]);
                else
                    DUM = ncread(FName,[nc_grp varname{i}],startLoc,count);    
                end
            catch
                DUM = ncread(FName,[nc_grp varname{i}]); % bug fixing due to some new S3B file variables from coda.eumetsat.int
                %error(sprintf('Error reading var %s',varname{i}));
            end

            DUM(DUM == FV) = NaN;
            eval(sprintf('ncdata.(''%s'') = DUM;',varname{i}));
        end
    end   
end

vecsize = size(ncdata.(latvar));
    
%Copy fields of the struct CS1a to a struct that has the same format as the
CS1a.FBR.data                                   = complex(ncdata.(ivar),ncdata.(qvar));

if is_s3
    CS1a.GEO.Start_Time = datetime(2000,1,1,0,0,ncdata.time_l1a_echo_sar_ku(1), 'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
    
    CS1a.GEO.TAI.days                               = ncdata.UTC_day_l1a_echo_sar_ku;
    CS1a.GEO.TAI.secs                               = ncdata.UTC_sec_l1a_echo_sar_ku;
    CS1a.GEO.LAT                                    = ncdata.lat_l1a_echo_sar_ku;
    CS1a.GEO.LON                                    = ncdata.lon_l1a_echo_sar_ku;
    CS1a.GEO.H                                      = ncdata.alt_l1a_echo_sar_ku;
    CS1a.GEO.H_rate                                 = ncdata.orb_alt_rate_l1a_echo_sar_ku;
    
    CS1a.GEO.FLAGS.time_status                      = ncdata.flag_time_status_l1a_echo_sar_ku;
    % CS1a.GEO.FLAGS.time_corr_val                    = ncdata.flag_time_corr_val_l1a_echo_sar_ku;
    % CS1a.GEO.FLAGS.man_pres                         = ncdata.flag_man_pres_l1a_echo_sar_ku;
    % CS1a.GEO.FLAGS.man_thrust                       = ncdata.flag_man_thrust_l1a_echo_sar_ku;
    % CS1a.GEO.FLAGS.man_plane                        = ncdata.flag_man_plane_l1a_echo_sar_ku;
    % CS1a.GEO.FLAGS.gnss_status                      = ncdata.flag_gnss_status_l1a_echo_sar_ku;
    CS1a.GEO.FLAGS.nav_bul_status                   = ncdata.nav_bul_status_l1a_echo_sar_ku;
    CS1a.GEO.FLAGS.nav_bul_source                   = ncdata.nav_bul_source_l1a_echo_sar_ku;
    % CS1a.GEO.FLAGS.isp_time_status                  = ncdata.isp_time_status_l1a_echo_sar_ku;
    CS1a.GEO.FLAGS.oper_instr                       = ncdata.oper_instr_l1a_echo_sar_ku;
    CS1a.GEO.FLAGS.SAR_mode                         = ncdata.SAR_mode_l1a_echo_sar_ku;
    CS1a.GEO.FLAGS.cl_gain                          = ncdata.cl_gain_l1a_echo_sar_ku;
    CS1a.GEO.FLAGS.acq_stat                         = ncdata.acq_stat_l1a_echo_sar_ku;
    CS1a.GEO.FLAGS.dem_eeprom                       = ncdata.dem_eeprom_l1a_echo_sar_ku;
    CS1a.GEO.FLAGS.weighting                        = ncdata.weighting_l1a_echo_sar_ku;
    CS1a.GEO.FLAGS.loss_track                       = ncdata.loss_track_l1a_echo_sar_ku;
    
    CS1a.GEO.V.Vx                                   = ncdata.x_vel_l1a_echo_sar_ku;
    CS1a.GEO.V.Vy                                   = ncdata.y_vel_l1a_echo_sar_ku;
    CS1a.GEO.V.Vz                                   = ncdata.z_vel_l1a_echo_sar_ku;
    CS1a.GEO.V.V                                    = sqrt(CS1a.GEO.V.Vx.^2+CS1a.GEO.V.Vy.^2+CS1a.GEO.V.Vz.^2);
    CS1a.GEO.Antenna_Bench_Roll                     = ncdata.roll_sral_mispointing_l1a_echo_sar_ku;
    CS1a.GEO.Antenna_Bench_Pitch                    = ncdata.pitch_sral_mispointing_l1a_echo_sar_ku;
    CS1a.GEO.Antenna_Bench_Yaw                      = zeros(vecsize);

    CS1a.MEA.ref_range                              = ncdata.range_ku_l1a_echo_sar_ku;
    CS1a.MEA.win_delay                              = ncdata.range_ku_l1a_echo_sar_ku ./ (CONST.c./2);
    CS1a.MEA.LAI                                    = zeros(vecsize);
    CS1a.MEA.burst_phase_cor                        = ncdata.burst_phase_cor_ku_l1a_echo_sar_ku;
    CS1a.MEA.burst_power_cor                        = ncdata.burst_power_cor_ku_l1a_echo_sar_ku;
    
    CS1a.SAR.agc_ku                                 = ncdata.agc_ku_l1a_echo_sar_ku;
    CS1a.SAR.scale_factor_ku                        = ncdata.scale_factor_ku_l1a_echo_sar_ku;
    % CS1a.SAR.N_averaged_echoes                      = ncdata.nb_stack_l1a_echo_sar_ku;
    % CS1a.SAR.beam_param.max_stack                   = ncdata.max_stack_l1a_echo_sar_ku;
    % CS1a.SAR.beam_param.stdev_stack                 = ncdata.stdev_stack_l1a_echo_sar_ku;
    % CS1a.SAR.beam_param.skew_stack                  = ncdata.skew_stack_l1a_echo_sar_ku;
    % CS1a.SAR.beam_param.kurt_stack                  = ncdata.kurt_stack_l1a_echo_sar_ku;
    % CS1a.SAR.BeamAngle                              = ncdata.beam_ang_stack_l1a_echo_sar_ku; %These are looking angles!
    % CS1a.SAR.data                                   = ncdata.i2q2_meas_ku_l1a_echo_sar_ku;

    CS1a.COR.surf_type                              = ncdata.surf_type_l1a_echo_sar_ku;
    CS1a.surf_type                                  = ncdata.surf_type_l1a_echo_sar_ku;

    CS1a.TIME                                       = ncdata.time_l1a_echo_sar_ku;
elseif is_s6
    CS1a.GEO.Start_Time = datetime(2000,1,1,0,0,ncdata.time(1), 'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');

    CS1a.GEO.TAI.days                               = ncdata.time_tai / (60*60*24);
    CS1a.GEO.TAI.secs                               = ncdata.time_tai;
    CS1a.GEO.LAT                                    = ncdata.latitude;
    CS1a.GEO.LON                                    = ncdata.longitude;
    CS1a.GEO.H                                      = ncdata.altitude;
    CS1a.GEO.H_rate                                 = ncdata.altitude_rate;
    
    CS1a.GEO.MCD_FLAG                               = ncdata.mcd_flags;
    if contains(FName, 'GPP')
        CS1a.GEO.MCD_FLAG                               = ~ncdata.mcd_flags;
    end
    
%     CS1a.GEO.FLAGS.time_status                      = ncdata.flag_time_status_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.time_corr_val                    = ncdata.flag_time_corr_val_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.man_pres                         = ncdata.flag_man_pres_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.man_thrust                       = ncdata.flag_man_thrust_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.man_plane                        = ncdata.flag_man_plane_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.gnss_status                      = ncdata.flag_gnss_status_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.nav_bul_status                   = ncdata.nav_bul_status_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.nav_bul_source                   = ncdata.nav_bul_source_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.isp_time_status                  = ncdata.isp_time_status_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.oper_instr                       = ncdata.oper_instr_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.SAR_mode                         = ncdata.SAR_mode_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.cl_gain                          = ncdata.cl_gain_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.acq_stat                         = ncdata.acq_stat_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.dem_eeprom                       = ncdata.dem_eeprom_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.weighting                        = ncdata.weighting_l1a_echo_sar_ku;
%     CS1a.GEO.FLAGS.loss_track                       = ncdata.loss_track_l1a_echo_sar_ku;
    
    CS1a.GEO.V.Vx                                   = ncdata.velocity_vector(1,:)';
    CS1a.GEO.V.Vy                                   = ncdata.velocity_vector(2,:)';
    CS1a.GEO.V.Vz                                   = ncdata.velocity_vector(3,:)';
    CS1a.GEO.V.V                                    = sqrt(CS1a.GEO.V.Vx.^2+CS1a.GEO.V.Vy.^2+CS1a.GEO.V.Vz.^2);
    CS1a.GEO.Antenna_Bench_Roll                     = ncdata.off_nadir_roll_angle_pf;
    CS1a.GEO.Antenna_Bench_Pitch                    = ncdata.off_nadir_pitch_angle_pf;
    CS1a.GEO.Antenna_Bench_Yaw                      = zeros(vecsize);

    CS1a.MEA.altimeter_clock                        = ncdata.altimeter_clock;
    CS1a.MEA.PRI                                    = ncdata.tm_pri .* 4 .* 1./ncdata.altimeter_clock;
    CS1a.MEA.ref_range                              = ncdata.tracker_range_calibrated;
    CS1a.MEA.win_delay                              = ncdata.tracker_range_calibrated ./ (CONST.c./2);
    CS1a.MEA.LAI                                    = zeros(vecsize);
    CS1a.MEA.p4_mode_flag                           = ncdata.p4_mode_flag;
    CS1a.MEA.tm_burst_num                           = ncdata.tm_burst_num;
    if ~contains(FName, 'GPP')
        CS1a.MEA.burst_phase_cor                        = ncdata.burst_phase_array_cor;
        CS1a.MEA.burst_power_cor                        = ncdata.burst_power_array_cor;
        CS1a.MEA.cal2_correction                        = ncdata.cal2_correction;
    end
    
    CS1a.SAR.agc_ku                                 = ncdata.variable_digital_gain;
    CS1a.SAR.scale_factor_ku                        = ncdata.power_scaling_to_antenna;
    % CS1a.SAR.N_averaged_echoes                      = ncdata.nb_stack_l1a_echo_sar_ku;
    % CS1a.SAR.beam_param.max_stack                   = ncdata.max_stack_l1a_echo_sar_ku;
    % CS1a.SAR.beam_param.stdev_stack                 = ncdata.stdev_stack_l1a_echo_sar_ku;
    % CS1a.SAR.beam_param.skew_stack                  = ncdata.skew_stack_l1a_echo_sar_ku;
    % CS1a.SAR.beam_param.kurt_stack                  = ncdata.kurt_stack_l1a_echo_sar_ku;
    % CS1a.SAR.BeamAngle                              = ncdata.beam_ang_stack_l1a_echo_sar_ku; %These are looking angles!
    % CS1a.SAR.data                                   = ncdata.i2q2_meas_ku_l1a_echo_sar_ku;

    % [0 1 2 3 4 5 6]=[open_ocean land continental_water aquatic_vegetation continental_ice_snow floating_ice salted_basin] (S6 baseline F03)
    CS1a.COR.surf_type                              = ncdata.surface_classification_flag;
    CS1a.surf_type                                  = ncdata.surface_classification_flag;

    CS1a.TIME                                       = ncdata.time;   
end

% some dummy stuff
CS1a.GEO.FName = FName;
CS1a.GEO.BaseLine.X = zeros(vecsize);
CS1a.GEO.BaseLine.Y = zeros(vecsize);
CS1a.GEO.BaseLine.Z = zeros(vecsize);
CS1a.GEO.Beam.X = zeros(vecsize);
CS1a.GEO.Beam.Y = zeros(vecsize);
CS1a.GEO.Beam.Z = zeros(vecsize);

end