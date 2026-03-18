function [S3,CS] = S3_L1bs_read(FName,DOM,ReadStack)

LoadCommonSettings

% defval('FName','S3A_SR_1_SRA_BS_20180211T201725_20180211T210749_20180308T213530_3024_027_370______LN3_O_NT_003.SEN3')
% defval('ReadStack',false);
% 
% %Set paths
% if contains(FName,'S3A')
%     PathL1b = fullfile(PathDATA,'RadAlt','Sentinel3A','SR_1_SRA_BS');
%     PathL2  = fullfile(PathDATA,'RadAlt','Sentinel3A','SR_2_WAT');
% elseif contains(FName,'S3B')
%     PathL1b = fullfile(PathDATA,'RadAlt','Sentinel3B','SR_1_SRA_BS');
%     PathL2  = fullfile(PathDATA,'RadAlt','Sentinel3B','SR_2_WAT');
% end

%Open netCDF file
ncid = netcdf.open(FName,'NC_NOWRITE');

%Return the number of variables in a netCDF file.
[~,nvars,~,~] = netcdf.inq(ncid);

%Create list of variable names associated to SAR ku echoes
varname = cell(nvars,1);
for i=0:nvars-1
    [varname{i+1},~,~,~] = netcdf.inqVar(ncid,i);
end
varids  = find(contains(varname,'echo_sar_ku'));
varname = varname(varids);

%Exclude fields i_echoes_ku_l1bs_echo_sar_ku and q_echoes_ku_l1bs_echo_sar_ku
if ~ReadStack
    varids  = varids(~ismember(varname,{'i_echoes_ku_l1bs_echo_sar_ku','q_echoes_ku_l1bs_echo_sar_ku'}));
    varname = varname(~ismember(varname,{'i_echoes_ku_l1bs_echo_sar_ku','q_echoes_ku_l1bs_echo_sar_ku'}));
end

%Identify indices of observations inside the area of interest
LAT  = ncread(FName,'lat_l1bs_echo_sar_ku');
LON  = ncread(FName,'lon_l1bs_echo_sar_ku');
IDX  = ingeoquad(LAT,LON,DOM(1,:),DOM(2,:));
IDXs = find(IDX,1,'first')-1;
IDXc = sum(IDX);

%Return if file does not contain 20Hz sar data
if IDXc == 0, [S3,CS] = deal([]); netcdf.close(ncid); return, end

%Read all variables. Apply scaling/add offset and replace fill values by NaNs
for i=1:numel(varids)
    try SF = netcdf.getAtt(ncid,varids(i)-1,'scale_factor'); catch, SF = 1; end
    try AO = netcdf.getAtt(ncid,varids(i)-1,'add_offset'); catch, AO = 0; end
    try FV = netcdf.getAtt(ncid,varids(i)-1,'_FillValue'); catch, FV = NaN; end
    try 
        DUM  = netcdf.getVar(ncid,varids(i)-1,IDXs,IDXc,'double');
    catch
        [~,~,dimids,~] = netcdf.inqVar(ncid,varids(i)-1);
        [~,dimlen1]    = netcdf.inqDim(ncid,dimids(1));
        if numel(dimids) == 2
            DUM            = netcdf.getVar(ncid,varids(i)-1,[0,IDXs],[dimlen1,IDXc],'double');
        else
            [~,dimlen2]    = netcdf.inqDim(ncid,dimids(2));
            DUM            = netcdf.getVar(ncid,varids(i)-1,[0,0,IDXs],[dimlen1,dimlen2,IDXc],'double');
        end            
    end
    DUM(DUM == FV) = NaN;
    eval(sprintf('S3.(''%s'') = (SF*DUM) + AO;',varname{i}));
end

%Return if file does not contain 20Hz sar data
if isempty(S3.time_l1bs_echo_sar_ku), [S3,CS] = deal([]); netcdf.close(ncid); return, end

%Get start time
StrtTime = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'first_meas_time');

%Close netCDF file
netcdf.close(ncid);

%If number of output arguments == 1, return. Else return data in format
%used to store CryoSat data
if nargout == 1, return, end

%Copy fields of the struct S3 to a struct that has the same format as the
%one used for CryoSat data
CS.GEO.('Start_Time')                                         = StrtTime;
CS.GEO.TAI.('days')                                           = S3.UTC_day_l1bs_echo_sar_ku;
CS.GEO.TAI.('secs')                                           = S3.UTC_sec_l1bs_echo_sar_ku;
CS.GEO.('LAT')                                                = S3.lat_l1bs_echo_sar_ku;
CS.GEO.('LON')                                                = S3.lon_l1bs_echo_sar_ku;
CS.GEO.('H')                                                  = S3.alt_l1bs_echo_sar_ku;
CS.GEO.('H_rate')                                             = S3.orb_alt_rate_l1bs_echo_sar_ku;
CS.GEO.FLAGS.('time_corr_val')                                = zeros(size(S3.lat_l1bs_echo_sar_ku)); %Info not available
CS.GEO.FLAGS.('man_thrust')                                   = zeros(size(S3.lat_l1bs_echo_sar_ku)); %Info not available
CS.GEO.FLAGS.('gnss_status')                                  = zeros(size(S3.lat_l1bs_echo_sar_ku)); %Info not available
CS.GEO.FLAGS.('loss_track')                                   = zeros(size(S3.lat_l1bs_echo_sar_ku)); %Info not available
CS.GEO.V.('Vx')                                               = S3.x_vel_l1bs_echo_sar_ku;
CS.GEO.V.('Vy')                                               = S3.y_vel_l1bs_echo_sar_ku;
CS.GEO.V.('Vz')                                               = S3.z_vel_l1bs_echo_sar_ku;
CS.GEO.('Antenna_Bench_Roll')                                 = S3.roll_sral_mispointing_l1bs_echo_sar_ku;
CS.GEO.('Antenna_Bench_Pitch')                                = S3.pitch_sral_mispointing_l1bs_echo_sar_ku;
CS.GEO.('Antenna_Bench_Yaw')                                  = S3.yaw_sral_mispointing_l1bs_echo_sar_ku;
CS.MEA.('ref_range')                                          = S3.range_ku_l1bs_echo_sar_ku;
CS.SAR.('agc_ku_l1bs_echo_sar_ku')                            = S3.agc_ku_l1bs_echo_sar_ku;
CS.SAR.('scale_factor_ku')                                    = S3.scale_factor_ku_l1bs_echo_sar_ku;
CS.SAR.('N_averaged_echoes')                                  = S3.nb_stack_l1bs_echo_sar_ku;
CS.SAR.beam_param.('max_stack')                               = S3.max_stack_l1bs_echo_sar_ku;
CS.SAR.beam_param.('max_loc_stack')                           = S3.max_loc_stack_l1bs_echo_sar_ku;
CS.SAR.beam_param.('stdev_stack')                             = S3.stdev_stack_l1bs_echo_sar_ku;
CS.SAR.beam_param.('stdev_stack')                             = S3.stdev_stack_l1bs_echo_sar_ku;
CS.SAR.beam_param.('skew_stack')                              = S3.skew_stack_l1bs_echo_sar_ku;
CS.SAR.beam_param.('kurt_stack')                              = S3.kurt_stack_l1bs_echo_sar_ku;
CS.SAR.beam_param.('start_beam_ang_stack')                    = S3.start_beam_ang_stack_l1bs_echo_sar_ku;
CS.SAR.beam_param.('stop_beam_ang_stack')                     = S3.stop_beam_ang_stack_l1bs_echo_sar_ku;
CS.SAR.beam_param.('start_look_angle_stack')                  = S3.start_look_angle_stack_l1bs_echo_sar_ku;
CS.SAR.beam_param.('stop_look_angle_stack')                   = S3.stop_look_angle_stack_l1bs_echo_sar_ku;
CS.SAR.('BeamAngle')                                          = S3.beam_ang_stack_l1bs_echo_sar_ku; %These are looking angles!
CS.SAR.('data')                                               = S3.i2q2_meas_ku_l1bs_echo_sar_ku;
CS.COR.('surf_type')                                          = S3.surf_type_l1bs_echo_sar_ku;
CS.('surf_type')                                              = S3.surf_type_l1bs_echo_sar_ku;
% 
% %Load corrections from level-2 file
% FNameL2 = regexprep(regexprep(FName,'SR_1_SRA_BS','SR_2_WAT___'),'LN3','MAR');
% TMP     = dir(fullfile(PathL2,sprintf('%s*%s',FNameL2(1:55),FNameL2(end-29:end))));
% i       = 0;
% while isempty(TMP)
%     i   = i+1;
%     TMP = dir(fullfile(PathL2,sprintf('%s*%s',FNameL2(1:55-i),FNameL2(end-29:end))));
%     if i == 54, break, end
% end
% FNameL2 = TMP; clear i TMP
FNameL2 = [];
if ~isempty(FNameL2)
    try
        FNameL2                      = dir(fullfile(FNameL2.folder,FNameL2.name,'enhanced_measurement.nc'));
    catch
        FNameL2                      = dir(fullfile(FNameL2.folder,FNameL2.name,'standard_measurement.nc'));
    end
    CS.COR.TIME1Hz                   = datenum('2000','yyyy') + double(ncread(fullfile(FNameL2.folder,FNameL2.name),'UTC_day_01')) + ncread(fullfile(FNameL2.folder,FNameL2.name),'UTC_sec_01')/86400;
    CS.COR.TIME20Hz                  = datenum('2000','yyyy') + double(ncread(fullfile(FNameL2.folder,FNameL2.name),'UTC_day_20_ku')) + ncread(fullfile(FNameL2.folder,FNameL2.name),'UTC_sec_20_ku')/86400;
    CS.COR.LAT1Hz                    = ncread(fullfile(FNameL2.folder,FNameL2.name),'lat_01');
    CS.COR.LON1Hz                    = ncread(fullfile(FNameL2.folder,FNameL2.name),'lon_01');
    CS.COR.LAT20Hz                   = ncread(fullfile(FNameL2.folder,FNameL2.name),'lat_20_ku');
    CS.COR.LON20Hz                   = ncread(fullfile(FNameL2.folder,FNameL2.name),'lon_20_ku');
    CS.COR.iono_cor_alt_20_ku        = ncread(fullfile(FNameL2.folder,FNameL2.name),'iono_cor_alt_20_ku');
    CS.COR.iono_cor_gim_01_ku        = ncread(fullfile(FNameL2.folder,FNameL2.name),'iono_cor_gim_01_ku');
    CS.COR.dry_trop                  = ncread(fullfile(FNameL2.folder,FNameL2.name),'mod_dry_tropo_cor_zero_altitude_01');
    CS.COR.wet_trop                  = ncread(fullfile(FNameL2.folder,FNameL2.name),'rad_wet_tropo_cor_01_ku');
    CS.COR.SSB                       = ncread(fullfile(FNameL2.folder,FNameL2.name),'sea_state_bias_01_ku');
    CS.COR.solidearth_tide           = ncread(fullfile(FNameL2.folder,FNameL2.name),'solid_earth_tide_01');
    CS.COR.ocean_equilibrium_tide    = ncread(fullfile(FNameL2.folder,FNameL2.name),'ocean_tide_sol1_01');
    CS.COR.ocean_longperiod_tide     = zeros(size(CS.COR.ocean_equilibrium_tide)); %This value has already been added to the two geocentric ocean tide height values recorded in the product (ocean_tide_sol1_01 and ocean_tide_sol2_01).
    CS.COR.ocean_loading_tide        = zeros(size(CS.COR.ocean_equilibrium_tide)); %This value has already been added to the corresponding ocean tide height value recorded in the product (ocean_tide_sol2_01)."
    CS.COR.geocentric_polar_tide     = ncread(fullfile(FNameL2.folder,FNameL2.name),'pole_tide_01');
    CS.COR.inv_bar                   = ncread(fullfile(FNameL2.folder,FNameL2.name),'inv_bar_cor_01');
    CS.COR.hf_fluct_cor              = ncread(fullfile(FNameL2.folder,FNameL2.name),'hf_fluct_cor_01');
    CS.COR.cog_cor_01                = ncread(fullfile(FNameL2.folder,FNameL2.name),'cog_cor_01');
    CS.COR.mod_instr_cor_range_01_ku = ncread(fullfile(FNameL2.folder,FNameL2.name),'mod_instr_cor_range_01_ku');
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