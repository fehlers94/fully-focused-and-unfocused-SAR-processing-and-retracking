function [HDR,CS,CS_nc] = Cryo_L2_read_nc(full_filename)

% Scope: Matlab Function to ingest CryoSat L2 netcdf data products to matlab workspace
% Data Level: 1a & 1b
% Supported Modes: FBR, LRM, SAR, SARin, FDM
%
% Input Argument: <full_filename> is the full pathname (path + file) where the CryoSat .nc
% file is stored in your local drive
%
% Output Arguments: the structures <HDR>, containing the header of the read file, and the structure <CS>,
% containing the read data fields
%
% Author: Cornelis Slobbe
% Date: 05/01/2021
% Version: 1.0
% Debugging: for any issues, please write to D.C.Slobbe@tudelft.nl
% Track Change Log:
%
%    - version 1.0: read baseline D

%Return information about a NetCDF source
ginfo = ncinfo(full_filename);

%%%%%%%%%%%%%%%%%%%%%%%%  DATA HEADER READING %%%%%%%%%%%%%%%%%%%

HDR   = struct;
for i=1:numel(ginfo.Attributes)
    HDR.(ginfo.Attributes(i).Name) = ginfo.Attributes(i).Value;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Open netCDF file
ncid = netcdf.open(full_filename,'NOWRITE');

%Return the number of variables in a netCDF file.
[~,nvars,~,~] = netcdf.inq(ncid);

%Create list of variable names associated to SAR ku echoes
varname = cell(nvars,1);
for i=0:nvars-1
    [varname{i+1},~,~,~] = netcdf.inqVar(ncid,i);
end
varids  = 1:nvars;

%Read all variables. Apply scaling/add offset and replace fill values by NaNs
for i=1:numel(varids)
    try SF = netcdf.getAtt(ncid,varids(i)-1,'scale_factor','double'); catch, SF = 1; end
    try AO = netcdf.getAtt(ncid,varids(i)-1,'add_offset','double'); catch, AO = 0; end
    try FV = netcdf.getAtt(ncid,varids(i)-1,'_FillValue','double'); catch, FV = NaN; end
    DUM            = netcdf.getVar(ncid,varids(i)-1,'double');
    DUM(DUM == FV) = NaN;
    eval(sprintf('CS_nc.(''%s'') = (SF*DUM) + AO;',regexprep(varname{i},{'20_ku','85_ku'},'NN_ku')));
end

%Return if file does not contain 20Hz/85Hz data
if isempty(CS_nc.time_NN_ku), [HDR,CS] = deal([]); netcdf.close(ncid); return, end
NrDATA = numel(CS_nc.time_NN_ku);

%Close netCDF file
netcdf.close(ncid);

%%%%%%%%%%%%%%%%%%%%%%%%% OPERATIVE MODE IDENTIFICATION  %%%%%%%%%%%%%%%%%%

CS.GEO.OPERATION_MODE = HDR.product_name(9:18);
DUM                   = textscan(CS.GEO.OPERATION_MODE,'%s','Delimiter','_');
CS.GEO.OPERATION_MODE = sprintf('%s_%s_L%s',DUM{1}{1},deblank(HDR.sir_op_mode),DUM{1}{3});
CS.GEO.Start_Time     = floor((datenum(regexprep(HDR.first_record_time,{'TAI=','T'},' '))-datenum('01-Jan-2000 00:00:00')).*24.*60.*60);
DUM                   = strsplit(HDR.first_record_time,'.');
CS.GEO.Start_Time     = CS.GEO.Start_Time + str2double(sprintf('0.%s',DUM{2}));
tmp                   = deblank(HDR.product_name);
CS.GEO.Baseline_ID    = tmp(end-3);

%%%%%%%%%%%%%%%%%%%%%%%%  DATA STRUCTURE INFORMATION %%%%%%%%%%%%%%%%%%%%%%

switch CS.GEO.OPERATION_MODE
    
    case {'SIR_LRM_L2','SIR_SAR_L2A','SIR_SAR_L2B','SIR_SIN_L2','SIR_SID_L2','SIR_GDR_2A','SIR_GDR_2B','SIR_SAR_L2','SIR_GDR_2_','SIR_SARIN_L2'}
        
        N_block=20;
        
    case {'SIR_FDM_L2'}
        
        error('Mode not implemented');
        
    otherwise
        
        disp('Mode not Supported or File Product not recognized');
        CS=[];HDR=[];return;
        
end

n_recs = ceil(numel(CS_nc.time_NN_ku)/20);

if ~~(mod(n_recs,1))
    
    disp('File Product corrupt');
    HDR=[];CS=[];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%  DATA STRUCTURE INITIALIZATION %%%%%%%%%%%%%%%%%%%

switch CS.GEO.OPERATION_MODE
    
    case  {'SIR_LRM_L2','SIR_SAR_L2A','SIR_SAR_L2B','SIR_SIN_L2','SIR_SID_L2','SIR_GDR_2A','SIR_GDR_2B','SIR_SAR_L2','SIR_GDR_2_','SIR_SARIN_L2'}
        
        CS.GEO.MEAS_MODE=zeros(N_block,n_recs);
        CS.COR.SURF_TYPE=zeros(N_block,n_recs);
        
        CS.MEA.time_20Hz=zeros(N_block,n_recs);
        CS.MEA.delta_time_20Hz=zeros(N_block,n_recs);
        CS.MEA.LAT_20Hz=zeros(N_block,n_recs);
        CS.MEA.LON_20Hz=zeros(N_block,n_recs);
        CS.MEA.surf_height_r1_20Hz=zeros(N_block,n_recs);
        CS.MEA.surf_height_r2_20Hz=zeros(N_block,n_recs);
        CS.MEA.surf_height_r3_20Hz=zeros(N_block,n_recs);
        CS.MEA.height_sea_ice_floe=zeros(N_block,n_recs);
        CS.MEA.height_sea_ice_lead=zeros(N_block,n_recs);
        CS.MEA.backsc_sig_r1_20Hz=zeros(N_block,n_recs);
        CS.MEA.backsc_sig_r2_20Hz=zeros(N_block,n_recs);
        CS.MEA.backsc_sig_r3_20Hz=zeros(N_block,n_recs);
        CS.MEA.freeboard_20Hz=zeros(N_block,n_recs);
        CS.MEA.SLA_interp_20Hz=zeros(N_block,n_recs);
        CS.MEA.SLA_N_rec_20Hz=zeros(N_block,n_recs);
        CS.MEA.SLA_interp_qual_20Hz=zeros(N_block,n_recs);
        CS.MEA.peakiness_20Hz=zeros(N_block,n_recs);
        CS.MEA.beam_avg_N_20Hz=zeros(N_block,n_recs);
        CS.MEA.Qual_r1_Value_20Hz=zeros(N_block,n_recs);
        CS.MEA.Qual_r2_Value_20Hz=zeros(N_block,n_recs);
        CS.MEA.Qual_r3_Value_20Hz=zeros(N_block,n_recs);
        
        CS.MEA.Qual_flag_20Hz.height_sea_ice_error=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.record_degraded=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.orbit_err=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.orbit_discont=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.height_err_r1=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.height_err_r2=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.height_err_r3=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.backsc_err_r1=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.backsc_err_r2=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.backsc_err_r3=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.interp_SSHA_err=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.peakiness_err=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.freeboard_err=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.SAR_discriminator_ocean=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.SAR_discriminator_lead=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.SAR_discriminator_ice=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.SAR_discriminator_unknown=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.SARin_Xtrack_angle_err=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.SARin_ch1_err=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.SARin_ch2_err=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.SIRAL_id=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.surf_mod_available=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.misp_err=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.delta_time_err=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.LRM_Slope_Model_valid=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.SIN_Baseline_bad=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.SIN_Out_of_Range=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.SIN_Velocity_bad=zeros(N_block,n_recs);
        CS.MEA.Qual_flag_20Hz.cal_wrng=zeros(N_block,n_recs);
        
        CS.MEA.Corr_Applic_flag_20Hz.corr_int_cal=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_rad_doppler=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_dry_tropo=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_wet_tropo=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_inv_bar=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_hf_atmo=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_iono_gim=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_iono_model=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_ocean_tide=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_lpe_tide=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_ocean_loading=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_se_tide=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_geo_polar_tide=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_doppler_slope_correction=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.mode_window_offset_applied=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.SAR_retrack_applied=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.SIN_retrack_applied=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.LRM_retrack_applied=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.LRM_ocean_bias_applied=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.LRM_ice_bias_applied=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.SAR_ocean_bias_applied=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.SAR_ice_bias_applied=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.SIN_ocean_bias_applied=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.SIN_ice_bias_applied=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_slope=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.sarin_bad_baseline=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.sarin_out_of_range=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.sarin_bad_velocity=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.corr_ssb=zeros(N_block,n_recs);
        CS.MEA.Corr_Applic_flag_20Hz.master_failure=zeros(N_block,n_recs);
            
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading the product

switch CS.GEO.OPERATION_MODE
    
    case  {'SIR_LRM_L2','SIR_SAR_L2A','SIR_SAR_L2B','SIR_SIN_L2','SIR_SID_L2','SIR_GDR_2A','SIR_GDR_2B','SIR_SAR_L2','SIR_GDR_2_','SIR_SARIN_L2'}
        
        % Time and Orbit Group
        CS.GEO.TAI.days                                    = floor(CS_nc.time_cor_01/86400);
        CS.GEO.TAI.secs                                    = floor(CS_nc.time_cor_01 - (CS.GEO.TAI.days*86400));
        CS.GEO.TAI.microsecs                               = (CS_nc.time_cor_01 - floor(CS_nc.time_cor_01)).*1E6;
        CS.GEO.MEAS_MODE(1:NrDATA)                         = CS_nc.flag_instr_mode_op_NN_ku;
        CS.GEO.LAT                                         = CS_nc.lat_01;
        CS.GEO.LON                                         = CS_nc.lon_01;
        CS.GEO.H                                           = CS_nc.alt_01;
        CS.GEO.Spacecraft_Roll_Angle                       = CS_nc.off_nadir_roll_angle_str_01;
        CS.GEO.Spacecraft_Pitch_Angle                      = CS_nc.off_nadir_pitch_angle_str_01;
        CS.GEO.Spacecraft_Yaw_Angle                        = CS_nc.off_nadir_yaw_angle_str_01;
        CS.GEO.valid_meas                                  = CS_nc.num_valid_01;
        
        %External Correction Group
        CS.COR.dry_tropo                                   = CS_nc.mod_dry_tropo_cor_01;
        CS.COR.wet_tropo                                   = CS_nc.mod_wet_tropo_cor_01;
        CS.COR.inv_baro                                    = CS_nc.inv_bar_cor_01;
        CS.COR.dac                                         = CS_nc.hf_fluct_total_cor_01;
        CS.COR.iono                                        = CS_nc.iono_cor_01;
        CS.COR.gim_ion                                     = CS_nc.iono_cor_gim_01;
        CS.COR.ssb                                         = CS_nc.sea_state_bias_01_ku;
        CS.COR.ocean_tide                                  = CS_nc.ocean_tide_01;
        CS.COR.lpe_ocean                                   = CS_nc.ocean_tide_eq_01;
        CS.COR.ocean_loading                               = CS_nc.load_tide_01;
        CS.COR.solid_earth                                 = CS_nc.solid_earth_tide_01;
        CS.COR.geoc_polar                                  = CS_nc.pole_tide_01;
        CS.COR.SURF_TYPE(1:NrDATA)                         = CS_nc.surf_type_NN_ku;
        CS.COR.mean_sea_surf_sea_ice                       = CS_nc.mean_sea_surf_sea_ice_01;
        CS.COR.geoid                                       = CS_nc.geoid_01;
        CS.COR.odle                                        = CS_nc.odle_01;
        CS.COR.ice_concentration                           = CS_nc.sea_ice_concentration_01;
        CS.COR.snow_depth                                  = CS_nc.snow_depth_01;
        CS.COR.snow_density                                = CS_nc.snow_density_01;
        dummy                                              = dec2bin(CS_nc.flag_cor_err_01,24);
        CS.COR.corr_status.dry_tropo                       = bin2dec(dummy(:,13));  % 0=OK 1=invalid
        CS.COR.corr_status.wet_tropo                       = bin2dec(dummy(:,14));  % 0=OK 1=invalid
        CS.COR.corr_status.inv_baro                        = bin2dec(dummy(:,15));  % 0=OK 1=invalid
        CS.COR.corr_status.dac                             = bin2dec(dummy(:,16));  % 0=OK 1=invalid
        CS.COR.corr_status.iono                            = bin2dec(dummy(:,18));  % 0=OK 1=invalid
        CS.COR.corr_status.gim_ion                         = bin2dec(dummy(:,17));  % 0=OK 1=invalid
        CS.COR.corr_status.ssb                             = bin2dec(dummy(:,12));  % 0=OK 1=invalid
        CS.COR.corr_status.ocean_tide                      = bin2dec(dummy(:,19));  % 0=OK 1=invalid
        CS.COR.corr_status.lpe_ocean                       = bin2dec(dummy(:,20));  % 0=OK 1=invalid
        CS.COR.corr_status.ocean_loading                   = bin2dec(dummy(:,21));  % 0=OK 1=invalid
        CS.COR.corr_status.solid_earth                     = bin2dec(dummy(:,22));  % 0=OK 1=invalid
        CS.COR.corr_status.geoc_polar                      = bin2dec(dummy(:,23));  % 0=OK 1=invalid
        CS.COR.corr_status.surf_type                       = bin2dec(dummy(:,24));  % 0=OK 1=invalid
        CS.COR.corr_status.mss                             = bin2dec(dummy(:,7));   % 0=OK 1=invalid
        CS.COR.corr_status.geoid                           = bin2dec(dummy(:,8));   % 0=OK 1=invalid
        CS.COR.corr_status.odle_model                      = bin2dec(dummy(:,9));   % 0=OK 1=invalid
        CS.COR.corr_status.ice_conc                        = bin2dec(dummy(:,4));   % 0=OK 1=invalid
        CS.COR.corr_status.snow_depth                      = bin2dec(dummy(:,5));   % 0=OK 1=invalid
        CS.COR.corr_status.snow_dens                       = bin2dec(dummy(:,6));   % 0=OK 1=invalid
        CS.COR.corr_status.swh                             = bin2dec(dummy(:,3));   % 0=OK 1=invalid
        CS.COR.corr_status.wind_spd                        = bin2dec(dummy(:,2));   % 0=OK 1=invalid
        CS.COR.corr_status.slope_model_error               = bin2dec(dummy(:,11));  % 0=OK 1=invalid
        CS.COR.corr_status.dem_error                       = bin2dec(dummy(:,10));  % 0=OK 1=invalid
        CS.COR.SWH                                         = CS_nc.swh_ocean_01_ku;
        CS.COR.wind_speed                                  = CS_nc.wind_speed_alt_01_ku;
        
        CS.COR.total_ocean=CS.COR.dry_tropo+CS.COR.wet_tropo+CS.COR.dac+CS.COR.iono+...
            CS.COR.ssb+CS.COR.ocean_tide+CS.COR.lpe_ocean+CS.COR.ocean_loading+CS.COR.solid_earth+CS.COR.geoc_polar;
        
        CS.COR.total_ocean_gim=CS.COR.dry_tropo+CS.COR.wet_tropo+CS.COR.dac+CS.COR.gim_ion+...
            CS.COR.ssb+CS.COR.ocean_tide+CS.COR.lpe_ocean+CS.COR.ocean_loading+CS.COR.solid_earth+CS.COR.geoc_polar;
        
        %Measurements Group
        CS.MEA.time_20Hz(1:NrDATA)                                           = CS_nc.time_NN_ku;
        CS.MEA.delta_time_20Hz                                               = CS.MEA.time_20Hz-repmat(CS.MEA.time_20Hz(1,:),N_block,1);
        CS.MEA.delta_time_20Hz(NrDATA+1:end)                                 = 0;
        CS.MEA.LAT_20Hz(1:NrDATA)                                            = CS_nc.lat_poca_NN_ku;
        CS.MEA.LON_20Hz(1:NrDATA)                                            = CS_nc.lon_poca_NN_ku;
        CS.MEA.surf_height_r1_20Hz(1:NrDATA)                                 = CS_nc.height_1_NN_ku;
        CS.MEA.surf_height_r2_20Hz(1:NrDATA)                                 = CS_nc.height_2_NN_ku;
        CS.MEA.surf_height_r3_20Hz(1:NrDATA)                                 = CS_nc.height_3_NN_ku;
        CS.MEA.height_sea_ice_floe(1:NrDATA)                                 = CS_nc.height_sea_ice_floe_NN_ku;
        CS.MEA.height_sea_ice_lead(1:NrDATA)                                 = CS_nc.height_sea_ice_lead_NN_ku;
        CS.MEA.backsc_sig_r1_20Hz(1:NrDATA)                                  = CS_nc.sig0_1_NN_ku;
        CS.MEA.backsc_sig_r2_20Hz(1:NrDATA)                                  = CS_nc.sig0_2_NN_ku;
        CS.MEA.backsc_sig_r3_20Hz(1:NrDATA)                                  = CS_nc.sig0_3_NN_ku;
        CS.MEA.freeboard_20Hz(1:NrDATA)                                      = CS_nc.freeboard_NN_ku;
        CS.MEA.SLA_interp_20Hz(1:NrDATA)                                     = CS_nc.ssha_interp_NN_ku;
        CS.MEA.SLA_N_rec_20Hz(1:NrDATA)                                      = CS_nc.ssha_interp_numval_NN_ku;
        CS.MEA.SLA_interp_qual_20Hz(1:NrDATA)                                = CS_nc.ssha_interp_rms_NN_ku;
        CS.MEA.peakiness_20Hz(1:NrDATA)                                      = CS_nc.peakiness_NN_ku;
        CS.MEA.beam_avg_N_20Hz(1:NrDATA)                                     = CS_nc.echo_avg_numval_NN_ku;
        dummy                                                                = dec2bin(CS_nc.flag_prod_status_NN_ku,32);
        CS.MEA.Qual_flag_20Hz.height_sea_ice_error(1:NrDATA)                 = bin2dec(dummy(:,4));
        CS.MEA.Qual_flag_20Hz.record_degraded(1:NrDATA)                      = bin2dec(dummy(:,5));  % 0=OK 1=block should not be processed
        CS.MEA.Qual_flag_20Hz.orbit_err(1:NrDATA)                            = bin2dec(dummy(:,6));  % 0=OK 1=error detected
        CS.MEA.Qual_flag_20Hz.orbit_discont(1:NrDATA)                        = bin2dec(dummy(:,7));  % 0=OK 1=orbit discontinuity occureg (e.g. gap)
        CS.MEA.Qual_flag_20Hz.height_err_r1(1:NrDATA)                        = bin2dec(dummy(:,8));  % 0=no 1=error in height derivation
        CS.MEA.Qual_flag_20Hz.height_err_r2(1:NrDATA)                        = bin2dec(dummy(:,9));  % 0=no 1=error in height derivation
        CS.MEA.Qual_flag_20Hz.height_err_r2(1:NrDATA)                        = bin2dec(dummy(:,10)); % 0=no 1=error in height derivation
        CS.MEA.Qual_flag_20Hz.backsc_err_r1(1:NrDATA)                        = bin2dec(dummy(:,11)); % 0=no 1=error
        CS.MEA.Qual_flag_20Hz.backsc_err_r2(1:NrDATA)                        = bin2dec(dummy(:,12)); % 0=no 1=error
        CS.MEA.Qual_flag_20Hz.backsc_err_r3(1:NrDATA)                        = bin2dec(dummy(:,13)); % 0=no 1=error
        CS.MEA.Qual_flag_20Hz.interp_SSHA_err(1:NrDATA)                      = bin2dec(dummy(:,14)); % 0=no 1=error
        CS.MEA.Qual_flag_20Hz.peakiness_err(1:NrDATA)                        = bin2dec(dummy(:,15)); % 0=no 1=error detected
        CS.MEA.Qual_flag_20Hz.freeboard_err(1:NrDATA)                        = bin2dec(dummy(:,16)); % 0=OK 1=invalid
        CS.MEA.Qual_flag_20Hz.SAR_discriminator_ocean(1:NrDATA)              = bin2dec(dummy(:,17)); % 0=no 1=yes
        CS.MEA.Qual_flag_20Hz.SAR_discriminator_lead(1:NrDATA)               = bin2dec(dummy(:,18)); % 0=no 1=yes
        CS.MEA.Qual_flag_20Hz.SAR_discriminator_ice(1:NrDATA)                = bin2dec(dummy(:,19)); % 0=no 1=yes
        CS.MEA.Qual_flag_20Hz.SAR_discriminator_unknown(1:NrDATA)            = bin2dec(dummy(:,20)); % 0=no 1=yes
        CS.MEA.Qual_flag_20Hz.SARin_Xtrack_angle_err(1:NrDATA)               = bin2dec(dummy(:,21)); % 0=no 1=ambiguous angle
        CS.MEA.Qual_flag_20Hz.SARin_ch1_err(1:NrDATA)                        = bin2dec(dummy(:,22)); % 0=OK 1=degraded or missing
        CS.MEA.Qual_flag_20Hz.SARin_ch2_err(1:NrDATA)                        = bin2dec(dummy(:,23)); % 0=OK 1=degraded or missing
        CS.MEA.Qual_flag_20Hz.SIRAL_id(1:NrDATA)                             = bin2dec(dummy(:,24)); % 0=nominal 1=redundant
        CS.GEO.MEAS_MODE_Siral_id                                            = CS.MEA.Qual_flag_20Hz.SIRAL_id(1,:)';
        CS.MEA.Qual_flag_20Hz.surf_mod_available(1:NrDATA)                   = bin2dec(dummy(:,25)); % 0=OK 1=no DEM/SLOPE model for location
        CS.MEA.Qual_flag_20Hz.misp_err(1:NrDATA)                             = bin2dec(dummy(:,26)); % 0=OK 1=error during calculation
        CS.MEA.Qual_flag_20Hz.delta_time_err(1:NrDATA)                       = bin2dec(dummy(:,27)); % 0=OK 1=error during calculation
        CS.MEA.Qual_flag_20Hz.LRM_Slope_Model_valid(1:NrDATA)                = bin2dec(dummy(:,28)); % 0=OK 1=error during calculation
        CS.MEA.Qual_flag_20Hz.SIN_Baseline_bad(1:NrDATA)                     = bin2dec(dummy(:,29)); % 0=OK 1=error during calculation
        CS.MEA.Qual_flag_20Hz.SIN_Out_of_Range(1:NrDATA)                     = bin2dec(dummy(:,30)); % 0=OK 1=error during calculation
        CS.MEA.Qual_flag_20Hz.SIN_Velocity_bad(1:NrDATA)                     = bin2dec(dummy(:,31)); % 0=OK 1=error during calculation
        CS.MEA.Qual_flag_20Hz.cal_wrng(1:NrDATA)                             = bin2dec(dummy(:,32)); % 0=no 1=non-nominal calibration correction
        dummy                                                                = dec2bin(CS_nc.flag_cor_applied_NN_ku,32);
        CS.MEA.Corr_Applic_flag_20Hz.corr_int_cal(1:NrDATA)                  = bin2dec(dummy(:,3));
        CS.MEA.Corr_Applic_flag_20Hz.corr_rad_doppler(1:NrDATA)              = bin2dec(dummy(:,4));
        CS.MEA.Corr_Applic_flag_20Hz.corr_dry_tropo(1:NrDATA)                = bin2dec(dummy(:,5));
        CS.MEA.Corr_Applic_flag_20Hz.corr_wet_tropo(1:NrDATA)                = bin2dec(dummy(:,6));
        CS.MEA.Corr_Applic_flag_20Hz.corr_inv_bar(1:NrDATA)                  = bin2dec(dummy(:,7));
        CS.MEA.Corr_Applic_flag_20Hz.corr_hf_atmo(1:NrDATA)                  = bin2dec(dummy(:,8));
        CS.MEA.Corr_Applic_flag_20Hz.corr_iono_gim(1:NrDATA)                 = bin2dec(dummy(:,9));
        CS.MEA.Corr_Applic_flag_20Hz.corr_iono_model(1:NrDATA)               = bin2dec(dummy(:,10));
        CS.MEA.Corr_Applic_flag_20Hz.corr_ocean_tide(1:NrDATA)               = bin2dec(dummy(:,11));
        CS.MEA.Corr_Applic_flag_20Hz.corr_lpe_tide(1:NrDATA)                 = bin2dec(dummy(:,12));
        CS.MEA.Corr_Applic_flag_20Hz.corr_ocean_loading(1:NrDATA)            = bin2dec(dummy(:,13));
        CS.MEA.Corr_Applic_flag_20Hz.corr_se_tide(1:NrDATA)                  = bin2dec(dummy(:,14));
        CS.MEA.Corr_Applic_flag_20Hz.corr_geo_polar_tide(1:NrDATA)           = bin2dec(dummy(:,15));
        CS.MEA.Corr_Applic_flag_20Hz.corr_doppler_slope_correction(1:NrDATA) = bin2dec(dummy(:,16));
        CS.MEA.Corr_Applic_flag_20Hz.mode_window_offset_applied(1:NrDATA)    = bin2dec(dummy(:,17));
        CS.MEA.Corr_Applic_flag_20Hz.SAR_retrack_applied(1:NrDATA)           = bin2dec(dummy(:,18));
        CS.MEA.Corr_Applic_flag_20Hz.SIN_retrack_applied(1:NrDATA)           = bin2dec(dummy(:,19));
        CS.MEA.Corr_Applic_flag_20Hz.LRM_retrack_applied(1:NrDATA)           = bin2dec(dummy(:,20));
        CS.MEA.Corr_Applic_flag_20Hz.LRM_ocean_bias_applied(1:NrDATA)        = bin2dec(dummy(:,21));
        CS.MEA.Corr_Applic_flag_20Hz.LRM_ice_bias_applied(1:NrDATA)          = bin2dec(dummy(:,22));
        CS.MEA.Corr_Applic_flag_20Hz.SAR_ocean_bias_applied(1:NrDATA)        = bin2dec(dummy(:,23));
        CS.MEA.Corr_Applic_flag_20Hz.SAR_ice_bias_applied(1:NrDATA)          = bin2dec(dummy(:,24));
        CS.MEA.Corr_Applic_flag_20Hz.SIN_ocean_bias_applied(1:NrDATA)        = bin2dec(dummy(:,25));    
        CS.MEA.Corr_Applic_flag_20Hz.SIN_ice_bias_applied(1:NrDATA)          = bin2dec(dummy(:,26));
        CS.MEA.Corr_Applic_flag_20Hz.corr_slope(1:NrDATA)                    = bin2dec(dummy(:,27));
        CS.MEA.Corr_Applic_flag_20Hz.sarin_bad_baseline(1:NrDATA)            = bin2dec(dummy(:,28));
        CS.MEA.Corr_Applic_flag_20Hz.sarin_out_of_range(1:NrDATA)            = bin2dec(dummy(:,29));
        CS.MEA.Corr_Applic_flag_20Hz.sarin_bad_velocity(1:NrDATA)            = bin2dec(dummy(:,30));
        CS.MEA.Corr_Applic_flag_20Hz.corr_ssb(1:NrDATA)                      = bin2dec(dummy(:,31));
        CS.MEA.Corr_Applic_flag_20Hz.master_failure(1:NrDATA)                = bin2dec(dummy(:,32));
        CS.MEA.Qual_r1_Value_20Hz(1:NrDATA)                                  = CS_nc.retracker_1_quality_NN_ku;
        CS.MEA.Qual_r2_Value_20Hz(1:NrDATA)                                  = CS_nc.retracker_2_quality_NN_ku;
        CS.MEA.Qual_r3_Value_20Hz(1:NrDATA)                                  = CS_nc.retracker_3_quality_NN_ku;
        
    case {'SIR_FDM_L2'}
        
        %%%%%%%%% Time and Orbit Group Reading %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%% Range Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%% Corrections Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%% SWH & Sigma0 Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%% Geophysical Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
end

CS.GEO.Serial_Sec_Num=CS.GEO.TAI.days.*24.*60.*60+CS.GEO.TAI.secs+CS.GEO.TAI.microsecs./1e6;
CS.GEO.Elapsed_Time=zeros(size(CS.GEO.Serial_Sec_Num));
CS.GEO.Elapsed_Time(CS.GEO.Serial_Sec_Num~=0)=CS.GEO.Serial_Sec_Num(CS.GEO.Serial_Sec_Num~=0)-CS.GEO.Start_Time;
