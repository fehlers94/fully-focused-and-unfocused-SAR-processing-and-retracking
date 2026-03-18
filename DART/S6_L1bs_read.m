function CS = S6_L1bs_read(FName,DOM,ReadStack,FNameL2)

LoadCommonSettings

% defval('FName','S6A_SR_1_SRA_BS_20180211T201725_20180211T210749_20180308T213530_3024_027_370______LN3_O_NT_003.SEN3')
defval('ReadStack',false);
defval('FNameL2',[]);

%Read netCDF file info
nc_grp = 'data_20/ku/';
finfo = ncinfo(FName, nc_grp);
% ncid = netcdf.open(FName,'NC_NOWRITE');

%Get number of variables in the netCDF file.
nvars = length(finfo.Variables);

%Create list of variable names associated to SAR ku echoes
varname = cell(nvars,1);
varsize = cell(nvars,1);
latvar = '';
for i=1:nvars
    varname{i} = finfo.Variables(i).Name;
    varsize{i} = finfo.Variables(i).Size;
    
    if contains(varname{i}, 'lat')
        latvar = varname{i};
    end
%     if contains(varname{i}, 'i_')
%         ivar = varname{i};
%     end
end
% qvar = replace(ivar,'i_','q_');

%Read all variables. replace fill values by NaNs
for i=1:numel(varname)
    try FV = ncreadatt(FName,[nc_grp varname{i}],'_FillValue'); catch, FV = NaN; end

    try
        if length(varsize{i}) > 1
            vs = length(varsize{i})-1;
%             DUM = ncread(FName,[nc_grp varname{i}],[ones(1,vs) startLoc],[inf*ones(1,vs) count]);
            DUM = ncread(FName,[nc_grp varname{i}]);
        else
%             DUM = ncread(FName,[nc_grp varname{i}],startLoc,count);
            DUM = ncread(FName,[nc_grp varname{i}]);
        end
    catch
        error(sprintf('Error reading var %s',varname{i}));
    end

    DUM(DUM == FV) = NaN;
    eval(sprintf('ncdata.(''%s'') = DUM;',varname{i}));
end

%Get start time
var_first_time = 'first_measurement_time';
StrtTime = ncreadatt(FName,'/',var_first_time);

vecsize = size(ncdata.(latvar));

CS.GEO.TAI.days                               = ncdata.time_tai / (60*60*24);
CS.GEO.TAI.secs                               = ncdata.time_tai;
CS.GEO.LAT                                    = ncdata.latitude;
CS.GEO.LON                                    = ncdata.longitude;
CS.GEO.H                                      = ncdata.altitude;
CS.GEO.H_rate                                 = ncdata.altitude_rate;
CS.GEO.BRI                                    = mean(diff(ncdata.time));

CS.GEO.MCD_FLAG                               = ncdata.mcd_flags;
% CS.GEO.FLAGS.time_status                      = ncdata.flag_time_status_l1a_echo_sar_ku;
% CS.GEO.FLAGS.time_corr_val                    = ncdata.flag_time_corr_val_l1a_echo_sar_ku;
% CS.GEO.FLAGS.man_pres                         = ncdata.flag_man_pres_l1a_echo_sar_ku;
% CS.GEO.FLAGS.man_thrust                       = ncdata.flag_man_thrust_l1a_echo_sar_ku;
% CS.GEO.FLAGS.man_plane                        = ncdata.flag_man_plane_l1a_echo_sar_ku;
% CS.GEO.FLAGS.gnss_status                      = ncdata.flag_gnss_status_l1a_echo_sar_ku;
% CS.GEO.FLAGS.nav_bul_status                   = ncdata.nav_bul_status_l1a_echo_sar_ku;
% CS.GEO.FLAGS.nav_bul_source                   = ncdata.nav_bul_source_l1a_echo_sar_ku;
% CS.GEO.FLAGS.isp_time_status                  = ncdata.isp_time_status_l1a_echo_sar_ku;
% CS.GEO.FLAGS.oper_instr                       = ncdata.oper_instr_l1a_echo_sar_ku;
% CS.GEO.FLAGS.SAR_mode                         = ncdata.SAR_mode_l1a_echo_sar_ku;
% CS.GEO.FLAGS.cl_gain                          = ncdata.cl_gain_l1a_echo_sar_ku;
% CS.GEO.FLAGS.acq_stat                         = ncdata.acq_stat_l1a_echo_sar_ku;
% CS.GEO.FLAGS.dem_eeprom                       = ncdata.dem_eeprom_l1a_echo_sar_ku;
% CS.GEO.FLAGS.weighting                        = ncdata.weighting_l1a_echo_sar_ku;
% CS.GEO.FLAGS.loss_track                       = ncdata.loss_track_l1a_echo_sar_ku;

CS.GEO.V.Vx                                   = ncdata.velocity_vector(1,:)';
CS.GEO.V.Vy                                   = ncdata.velocity_vector(2,:)';
CS.GEO.V.Vz                                   = ncdata.velocity_vector(3,:)';
CS.GEO.V.V                                    = sqrt(CS.GEO.V.Vx.^2+CS.GEO.V.Vy.^2+CS.GEO.V.Vz.^2);
CS.GEO.Antenna_Bench_Roll                     = ncdata.off_nadir_roll_angle_pf;
CS.GEO.Antenna_Bench_Pitch                    = ncdata.off_nadir_pitch_angle_pf;
CS.GEO.Antenna_Bench_Yaw                      = ncdata.off_nadir_yaw_angle_pf;

CS.GEO.FName = FName;

CS.MEA.ALTIMETER_CLOCK                        = mean(ncdata.altimeter_clock); % TODO might cause problems
CS.MEA.PRI                                    = ncdata.pulse_repetition_interval;
CS.MEA.ref_range                              = ncdata.tracker_range_calibrated;
CS.MEA.win_delay                              = ncdata.tracker_range_calibrated ./ (CONST.c./2);
CS.MEA.LAI                                    = zeros(vecsize);
CS.MEA.p4_mode_flag                               = ncdata.p4_mode_flag;

CS.SAR.agc_ku                                 = ones(vecsize);
% CS.SAR.agc_ku                                 = ncdata.variable_digital_gain;
CS.SAR.scale_factor_ku                        = ones(vecsize);
% CS.SAR.scale_factor_ku                        = ncdata.power_scaling_to_antenna;
%CS.SAR.StackMask                      = ncdata.stack_mask_start_stop;
CS.SAR.N_averaged_echoes                      = min(ncdata.num_looks_start_stop, ncdata.num_looks_multilooking);
% CS.SAR.beam_param.max_stack                   = ncdata.max_stack_l1a_echo_sar_ku;
% CS.SAR.beam_param.stdev_stack                 = ncdata.stdev_stack_l1a_echo_sar_ku;
% CS.SAR.beam_param.skew_stack                  = ncdata.skew_stack_l1a_echo_sar_ku;
% CS.SAR.beam_param.kurt_stack                  = ncdata.kurt_stack_l1a_echo_sar_ku;
% CS.SAR.BeamAngle                              = ncdata.beam_ang_stack_l1a_echo_sar_ku; %These are looking angles!
max_looks = ncinfo(FName).Dimensions(5).Length;
angle_per_look_rad = (ncdata.look_angle_stop - ncdata.look_angle_start) ./ (ncdata.num_looks_start_stop - 1);
CS.SAR.BeamAngle = NaN(max_looks, vecsize(1));
CS.SAR.stack_mask_start_stop                  = ncdata.stack_mask_start_stop;
% CS.SAR.stack_mask_start_stop                  = ones(size(ncdata.stack_mask_start_stop)) * NaN;

for i=1:vecsize(1)
    angs = ncdata.look_angle_start(i) + (0:ncdata.num_looks_start_stop(i)-1) * angle_per_look_rad(i);
    num_eff_looks = min(length(angs), ncdata.num_looks_multilooking(i));
    CS.SAR.BeamAngle(1:num_eff_looks, i) = angs(1:num_eff_looks);
end
CS.SAR.data                                   = ncdata.power_waveform .* transpose(ncdata.waveform_scale_factor);

% TODO: check how definition of this differs from S3's surf_type_l1a_echo_sar_ku
CS.COR.surf_type                              = ncdata.surface_classification_flag;
CS.surf_type                                  = ncdata.surface_classification_flag;

CS.TIME                                       = ncdata.time;

%Load corrections from level-2 file
if ~isempty(FNameL2)
    CS.COR.TIME1Hz                   = datenum('2000','yyyy') + double(ncread(FNameL2,'UTC_day_01')) + ncread(FNameL2,'UTC_sec_01')/86400;
    CS.COR.TIME20Hz                  = datenum('2000','yyyy') + double(ncread(FNameL2,'UTC_day_20_ku')) + ncread(FNameL2,'UTC_sec_20_ku')/86400;
    CS.COR.LAT1Hz                    = ncread(FNameL2,'lat_01');
    CS.COR.LON1Hz                    = ncread(FNameL2,'lon_01');
    CS.COR.LAT20Hz                   = ncread(FNameL2,'lat_20_ku');
    CS.COR.LON20Hz                   = ncread(FNameL2,'lon_20_ku');
    CS.COR.iono_cor_alt_20_ku        = ncread(FNameL2,'iono_cor_alt_20_ku');
    CS.COR.iono_cor_gim_01_ku        = ncread(FNameL2,'iono_cor_gim_01_ku');
    CS.COR.dry_trop                  = ncread(FNameL2,'mod_dry_tropo_cor_zero_altitude_01');
    CS.COR.wet_trop                  = ncread(FNameL2,'rad_wet_tropo_cor_01_ku');
    CS.COR.SSB                       = ncread(FNameL2,'sea_state_bias_01_ku');
    CS.COR.solidearth_tide           = ncread(FNameL2,'solid_earth_tide_01');
    CS.COR.ocean_equilibrium_tide    = ncread(FNameL2,'ocean_tide_sol1_01');
    CS.COR.ocean_longperiod_tide     = zeros(size(CS.COR.ocean_equilibrium_tide)); %This value has already been added to the two geocentric ocean tide height values recorded in the product (ocean_tide_sol1_01 and ocean_tide_sol2_01).
    CS.COR.ocean_loading_tide        = zeros(size(CS.COR.ocean_equilibrium_tide)); %This value has already been added to the corresponding ocean tide height value recorded in the product (ocean_tide_sol2_01)."
    CS.COR.geocentric_polar_tide     = ncread(FNameL2,'pole_tide_01');
    CS.COR.inv_bar                   = ncread(FNameL2,'inv_bar_cor_01');
    CS.COR.hf_fluct_cor              = ncread(FNameL2,'hf_fluct_cor_01');
    CS.COR.cog_cor_01                = ncread(FNameL2,'cog_cor_01');
    CS.COR.mod_instr_cor_range_01_ku = ncread(FNameL2,'mod_instr_cor_range_01_ku');
else
    CS.COR.TIME1Hz                   = datenum('2000','yyyy') + CS.GEO.TAI.days + CS.GEO.TAI.secs/86400;
    CS.COR.TIME20Hz                  = CS.COR.TIME1Hz;
    CS.COR.LAT1Hz                    = nan(numel(CS.COR.TIME1Hz),1);
    CS.COR.LON1Hz                    = nan(numel(CS.COR.TIME1Hz),1);
    CS.COR.LAT20Hz                   = nan(numel(CS.COR.TIME1Hz),1);
    CS.COR.LON20Hz                   = nan(numel(CS.COR.TIME1Hz),1);
    CS.COR.iono_cor_alt_20_ku        = zeros(numel(CS.COR.TIME1Hz),1);
    CS.COR.iono_cor_gim_01_ku        = zeros(numel(CS.COR.TIME1Hz),1);
    CS.COR.dry_trop                  = zeros(numel(CS.COR.TIME1Hz),1);
    CS.COR.wet_trop                  = zeros(numel(CS.COR.TIME1Hz),1);
    CS.COR.SSB                       = zeros(numel(CS.COR.TIME1Hz),1);
    CS.COR.solidearth_tide           = zeros(numel(CS.COR.TIME1Hz),1);
    CS.COR.ocean_equilibrium_tide    = zeros(numel(CS.COR.TIME1Hz),1);
    CS.COR.ocean_longperiod_tide     = zeros(numel(CS.COR.TIME1Hz),1);
    CS.COR.ocean_loading_tide        = zeros(numel(CS.COR.TIME1Hz),1);
    CS.COR.geocentric_polar_tide     = zeros(numel(CS.COR.TIME1Hz),1);
    CS.COR.inv_bar                   = zeros(numel(CS.COR.TIME1Hz),1);
    CS.COR.hf_fluct_cor              = zeros(numel(CS.COR.TIME1Hz),1);
    CS.COR.cog_cor_01                = zeros(numel(CS.COR.TIME1Hz),1);
    CS.COR.mod_instr_cor_range_01_ku = zeros(numel(CS.COR.TIME1Hz),1);
end

end