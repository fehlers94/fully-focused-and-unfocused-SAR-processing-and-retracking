function [HDR,CS,CS_nc] = Cryo_L1b_read_nc(full_filename)

% Scope: Matlab Function to ingest CryoSat L1b netcdf data products to matlab workspace
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
% Date: 07/04/2020
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
if isfield(CS_nc,'time_avg_01_ku'), NrAVG = numel(CS_nc.time_avg_01_ku); end
NrCOR  = numel(CS_nc.time_cor_01);

%Close netCDF file
netcdf.close(ncid);

%%%%%%%%%%%%%%%%%%%%%%%%% OPERATIVE MODE IDENTIFICATION  %%%%%%%%%%%%%%%%%%

CS.GEO.OPERATION_MODE = HDR.product_name(9:18);
DUM                   = textscan(CS.GEO.OPERATION_MODE,'%s','Delimiter','_');
CS.GEO.OPERATION_MODE = sprintf('%s_L%s_%s',DUM{1}{1},DUM{1}{3},deblank(HDR.sir_op_mode));
CS.GEO.Start_Time     = floor((datenum(regexprep(HDR.first_record_time,{'TAI=','T'},' '))-datenum('01-Jan-2000 00:00:00')).*24.*60.*60);
DUM                   = strsplit(HDR.first_record_time,'.');
CS.GEO.Start_Time     = CS.GEO.Start_Time + str2double(sprintf('0.%s',DUM{2}));
tmp                   = deblank(HDR.product_name);
CS.GEO.Baseline_ID    = tmp(end-3);

%%%%%%%%%%%%%%%%%%%%%%%%  DATA STRUCTURE INFORMATION %%%%%%%%%%%%%%%%%%%%%%

N_block=20;
SAR_pulses_burst=64;

switch  CS.GEO.Baseline_ID
    
    case  'D'
        
        N_samples=128;
        N_samples_SAR=256;
        
    otherwise
        
        error('Baseline_ID not recognized!')
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch CS.GEO.OPERATION_MODE
    
    case {'SIR_L1B_LRM','SIR_L1B_FDM'}
        
        n_points=N_samples;
        
    case 'SIR_FBR_SAR'
        
        n_points=2*N_samples*SAR_pulses_burst;
        
    case 'SIR_FBR_SARIN'
        
        n_points=2*(4.*N_samples)*SAR_pulses_burst;
        
    case 'SIR_L1B_SAR'
        
        n_points=N_samples_SAR;
        
    case 'SIR_L1B_SARIN'
        
        n_points=N_samples_SAR.*4; % 512 bins instead of 128 bins
        
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

CS.GEO.TAI.days=zeros(N_block,n_recs);
CS.GEO.TAI.secs=zeros(N_block,n_recs);
CS.GEO.TAI.microsecs=zeros(N_block,n_recs);
CS.GEO.USO=zeros(N_block,n_recs);
CS.GEO.SRC_CNT=zeros(N_block,n_recs);
CS.GEO.BURST_CNT=zeros(N_block,n_recs);
CS.GEO.LAT=zeros(N_block,n_recs);
CS.GEO.LON=zeros(N_block,n_recs);
CS.GEO.H=zeros(N_block,n_recs);
CS.GEO.H_rate=zeros(N_block,n_recs);
CS.GEO.V.Vx=zeros(N_block,n_recs);
CS.GEO.V.Vy=zeros(N_block,n_recs);
CS.GEO.V.Vz=zeros(N_block,n_recs);
CS.GEO.Beam.X=zeros(N_block,n_recs);
CS.GEO.Beam.Y=zeros(N_block,n_recs);
CS.GEO.Beam.Z=zeros(N_block,n_recs);
CS.GEO.BaseLine.X=zeros(N_block,n_recs);
CS.GEO.BaseLine.Y=zeros(N_block,n_recs);
CS.GEO.BaseLine.Z=zeros(N_block,n_recs);
CS.GEO.ind_first_meas_20hz_01=zeros(1,n_recs);
CS.GEO.ind_meas_1hz_20_ku=zeros(N_block,n_recs);

switch CS.GEO.OPERATION_MODE
    
    case {'SIR_L1B_LRM','SIR_L1B_FDM','SIR_L1B_SAR','SIR_L1B_SARIN'}
        
        CS.GEO.Star_Tracker_Usage=zeros(N_block,n_recs);
        CS.GEO.Antenna_Bench_Roll=zeros(N_block,n_recs);
        CS.GEO.Antenna_Bench_Pitch=zeros(N_block,n_recs);
        CS.GEO.Antenna_Bench_Yaw=zeros(N_block,n_recs);
        
end

CS.GEO.MODE_ID.Ins_Oper_Mode=zeros(N_block,n_recs);
CS.GEO.MODE_ID.Sarin_Degr_Case=zeros(N_block,n_recs);
CS.GEO.MODE_ID.Cal4_Flag=zeros(N_block,n_recs);
CS.GEO.MODE_ID.Plat_Att_Ctrl=zeros(N_block,n_recs);

CS.GEO.INS_CFG.Rx_Chain_Use=zeros(N_block,n_recs);
CS.GEO.INS_CFG.SIRAL_ID=zeros(N_block,n_recs);
CS.GEO.INS_CFG.Band_FLAG=zeros(N_block,n_recs);
CS.GEO.INS_CFG.Tracking_Mode=zeros(N_block,n_recs);
CS.GEO.INS_CFG.Ext_Cal=zeros(N_block,n_recs);
CS.GEO.INS_CFG.Loop_Status=zeros(N_block,n_recs);
CS.GEO.INS_CFG.Echo_Loss=zeros(N_block,n_recs);
CS.GEO.INS_CFG.Real_Time_Err=zeros(N_block,n_recs);
CS.GEO.INS_CFG.Echo_Satur_Err=zeros(N_block,n_recs);
CS.GEO.INS_CFG.Rx_Band_Atten=zeros(N_block,n_recs);
CS.GEO.INS_CFG.Cycle_Report_Err=zeros(N_block,n_recs);
CS.GEO.INS_CFG.STR_ATTREF=zeros(N_block,n_recs);

CS.GEO.MCD_FLAG.Block_Degraded=ones(N_block,n_recs);
CS.GEO.MCD_FLAG.Blank_Block=ones(N_block,n_recs);
CS.GEO.MCD_FLAG.Datation_Degraded=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.Orbit_Propag_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.Orbit_File_Change=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.Orbit_Discontinuity=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.Echo_Saturation=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.Other_Echo_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.Rx1_Err_SARin=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.Rx2_Err_SARin=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.Wind_Delay_Incon=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.AGC_Incon=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.CAL1_Corr_Miss=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.CAL1_Corr_IPF=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.DORIS_USO_Corr=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.Complex_CAL1_Corr_IPF=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.TRK_ECHO_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.RX1_ECHO_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.RX2_ECHO_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG.NPM_Incon=zeros(N_block,n_recs);

switch  CS.GEO.OPERATION_MODE
    
    case   {'SIR_L1B_LRM','SIR_L1B_FDM','SIR_L1B_SAR','SIR_L1B_SARIN'}
        
        CS.GEO.MCD_FLAG.CAL2_Corr_Missing=zeros(N_block,n_recs);
        CS.GEO.MCD_FLAG.Phase_Pertubation_Corr=zeros(N_block,n_recs);
        CS.GEO.MCD_FLAG.CAL2_Corr_Miss=zeros(N_block,n_recs);
        CS.GEO.MCD_FLAG.CAL2_Corr_IPF=zeros(N_block,n_recs);
        CS.GEO.MCD_FLAG.Power_Scaling_Err=zeros(N_block,n_recs);
        CS.GEO.MCD_FLAG.Att_Corr_Miss=zeros(N_block,n_recs);
        CS.GEO.MCD_FLAG.Phase_Pertubation_Mode=zeros(N_block,n_recs);
        
    case  {'SIR_FBR_SAR', 'SIR_FBR_SARIN'}
        
        CS.GEO.MCD_FLAG.Attitude_Corr_Missing=zeros(N_block,n_recs);
        CS.GEO.MCD_FLAG.CAL1_Corr_Type=zeros(N_block,n_recs);
        
end

CS.MEA.win_delay=zeros(N_block,n_recs);
CS.MEA.Ho=zeros(N_block,n_recs);
CS.MEA.Trk_H_rate=zeros(N_block,n_recs);
CS.MEA.LAI=zeros(N_block,n_recs);
CS.MEA.FAI=zeros(N_block,n_recs);
CS.MEA.AGC_1=zeros(N_block,n_recs);
CS.MEA.AGC_2=zeros(N_block,n_recs);
CS.MEA.Gain_Rx1=zeros(N_block,n_recs);
CS.MEA.Gain_Rx2=zeros(N_block,n_recs);
CS.MEA.Tx_Pwr=zeros(N_block,n_recs);
CS.MEA.dpl_range_corr=zeros(N_block,n_recs);
CS.MEA.ins_range_corr_rx_tx=zeros(N_block,n_recs);
CS.MEA.ins_range_corr_rx=zeros(N_block,n_recs);
CS.MEA.ins_gain_corr_rx_tx=zeros(N_block,n_recs);
CS.MEA.ins_gain_corr_rx=zeros(N_block,n_recs);
CS.MEA.int_phase_corr=zeros(N_block,n_recs);
CS.MEA.ext_phase_corr=zeros(N_block,n_recs);
CS.MEA.noise_power=zeros(N_block,n_recs);
CS.MEA.phase_slope_corr=zeros(N_block,n_recs);

CS.COR.dry_trop=zeros(1,n_recs);
CS.COR.wet_trop=zeros(1,n_recs);
CS.COR.inv_bar=zeros(1,n_recs);
CS.COR.dac=zeros(1,n_recs);
CS.COR.gim_ion=zeros(1,n_recs);
CS.COR.model_ion=zeros(1,n_recs);
CS.COR.ocean_equilibrium_tide=zeros(1,n_recs);
CS.COR.ocean_longperiod_tide=zeros(1,n_recs);
CS.COR.ocean_loading_tide=zeros(1,n_recs);
CS.COR.solidearth_tide=zeros(1,n_recs);
CS.COR.geocentric_polar_tide=zeros(1,n_recs);
CS.COR.surf_type=zeros(1,n_recs);

CS.COR.corr_status.dry_trop=zeros(1,n_recs);
CS.COR.corr_status.wet_trop=zeros(1,n_recs);
CS.COR.corr_status.inv_bar=zeros(1,n_recs);
CS.COR.corr_status.dac=zeros(1,n_recs);
CS.COR.corr_status.gim_iono=zeros(1,n_recs);
CS.COR.corr_status.model_iono=zeros(1,n_recs);
CS.COR.corr_status.ocean_equilibrium_tide=zeros(1,n_recs);
CS.COR.corr_status.ocean_longperiod_tide=zeros(1,n_recs);
CS.COR.corr_status.ocean_loading_tide=zeros(1,n_recs);
CS.COR.corr_status.solidearth_tide=zeros(1,n_recs);
CS.COR.corr_status.geocentric_polar_tide=zeros(1,n_recs);
CS.COR.corr_status.surface_type=zeros(1,n_recs);

CS.COR.corr_error.dry_trop=zeros(1,n_recs);
CS.COR.corr_error.wet_trop=zeros(1,n_recs);
CS.COR.corr_error.inv_bar=zeros(1,n_recs);
CS.COR.corr_error.dac=zeros(1,n_recs);
CS.COR.corr_error.gim_iono=zeros(1,n_recs);
CS.COR.corr_error.model_iono=zeros(1,n_recs);
CS.COR.corr_error.ocean_equilibrium_tide=zeros(1,n_recs);
CS.COR.corr_error.ocean_longperiod_tide=zeros(1,n_recs);
CS.COR.corr_error.ocean_loading_tide=zeros(1,n_recs);
CS.COR.corr_error.solidearth_tide=zeros(1,n_recs);
CS.COR.corr_error.geocentric_polar_tide=zeros(1,n_recs);
CS.COR.corr_error.surface_type=zeros(1,n_recs);

switch CS.GEO.OPERATION_MODE
    
    case {'SIR_L1B_LRM','SIR_L1B_FDM'}
        
        CS.AVG.TAI.days=zeros(1,NrAVG);
        CS.AVG.TAI.secs=zeros(1,NrAVG);
        CS.AVG.TAI.microsecs=zeros(1,NrAVG);
        CS.AVG.lat=zeros(1,NrAVG);
        CS.AVG.lon=zeros(1,NrAVG);
        CS.AVG.H=zeros(1,NrAVG);
        CS.AVG.win_delay=zeros(1,NrAVG);
        CS.AVG.data=zeros(n_points,NrAVG);
        CS.AVG.echo_scaling=zeros(1,NrAVG);
        CS.AVG.echo_scale_power=zeros(1,NrAVG);
        CS.AVG.N_averaged_echoes=zeros(1,NrAVG);
        CS.AVG.OneHz_Echo_Err=zeros(1,NrAVG);
        CS.LRM.data=zeros(n_points,N_block,n_recs);
        CS.LRM.echo_scaling=zeros(N_block,n_recs);
        CS.LRM.echo_scale_power=zeros(N_block,n_recs);
        CS.LRM.N_averaged_echoes=zeros(N_block,n_recs);
        CS.LRM.FLAG=zeros(N_block,n_recs);
        
    case {'SIR_L1B_SAR'}
        
        CS.AVG.TAI.days=zeros(1,NrAVG);
        CS.AVG.TAI.secs=zeros(1,NrAVG);
        CS.AVG.TAI.microsecs=zeros(1,NrAVG);
        CS.AVG.lat=zeros(1,NrAVG);
        CS.AVG.lon=zeros(1,NrAVG);
        CS.AVG.H=zeros(1,NrAVG);
        CS.AVG.win_delay=zeros(1,NrAVG);
        CS.AVG.data=zeros(N_samples,NrAVG);
        CS.AVG.echo_scaling=zeros(1,NrAVG);
        CS.AVG.echo_scale_power=zeros(1,NrAVG);
        CS.AVG.N_averaged_echoes=zeros(1,NrAVG);
        CS.AVG.OneHz_Echo_Err=zeros(1,NrAVG);
        CS.AVG.Mispointing_Err=zeros(1,NrAVG);
        
        CS.SAR.data=zeros(n_points,N_block,n_recs);
        CS.SAR.echo_scaling=zeros(N_block,n_recs);
        CS.SAR.echo_scale_power=zeros(N_block,n_recs);
        CS.SAR.N_averaged_echoes=zeros(N_block,n_recs);
        CS.SAR.FLAG.Approximate_Beam_Steering=zeros(N_block,n_recs);
        CS.SAR.FLAG.Exact_Beam_Steering=zeros(N_block,n_recs);
        CS.SAR.FLAG.Doppler_Weighting_Computed=zeros(N_block,n_recs);
        CS.SAR.FLAG.Doppler_Weighting_Applied_Before_Stack=zeros(N_block,n_recs);
        CS.SAR.FLAG.Multilook_Incomplete=zeros(N_block,n_recs);
        CS.SAR.FLAG.Beam_Angle_Steering_Err=zeros(N_block,n_recs);
        CS.SAR.FLAG.AntiAliased_Power_Echo=zeros(N_block,n_recs);
        CS.SAR.FLAG.Auto_Beam_Steering=zeros(N_block,n_recs);
        CS.SAR.beam_param=zeros(16,N_block,n_recs);
        
    case {'SIR_FBR_SAR'}
        
        CS.FBR.data=zeros(N_samples,SAR_pulses_burst,N_block,n_recs);
        CS.FBR.N_pulses=zeros(N_block,n_recs);
        CS.FBR.FLAG=char(zeros(16,N_block,n_recs));
        
    case {'SIR_FBR_SARIN'}
        
        CS.FBR.data_ch1=zeros(4.*N_samples,SAR_pulses_burst,N_block,n_recs);
        CS.FBR.data_ch2=zeros(4.*N_samples,SAR_pulses_burst,N_block,n_recs);
        CS.FBR.N_pulses=zeros(N_block,n_recs);
        CS.FBR.FLAG=char(zeros(16,N_block,n_recs));
        
    case 'SIR_L1B_SARIN'
        
        CS.AVG.TAI.days=zeros(1,NrAVG);
        CS.AVG.TAI.secs=zeros(1,NrAVG);
        CS.AVG.TAI.microsecs=zeros(1,NrAVG);
        CS.AVG.lat=zeros(1,NrAVG);
        CS.AVG.lon=zeros(1,NrAVG);
        CS.AVG.H=zeros(1,NrAVG);
        CS.AVG.win_delay=zeros(1,NrAVG);
        CS.AVG.data=zeros(N_samples.*4,NrAVG);
        CS.AVG.echo_scaling=zeros(1,NrAVG);
        CS.AVG.echo_scale_power=zeros(1,NrAVG);
        CS.AVG.N_averaged_echoes=zeros(1,NrAVG);
        CS.AVG.OneHz_Echo_Err=zeros(1,NrAVG);
        CS.AVG.Mispointing_Err=zeros(1,NrAVG);
        
        CS.SIN.data=zeros(n_points,N_block,n_recs);
        CS.SIN.echo_scaling=zeros(N_block,n_recs);
        CS.SIN.echo_scale_power=zeros(N_block,n_recs);
        CS.SIN.N_averaged_echoes=zeros(N_block,n_recs);
        CS.SIN.FLAG.Approximate_Beam_Steering=zeros(N_block,n_recs);
        CS.SIN.FLAG.Exact_Beam_Steering=zeros(N_block,n_recs);
        CS.SIN.FLAG.Doppler_Weighting_Computed=zeros(N_block,n_recs);
        CS.SIN.FLAG.Doppler_Weighting_Applied_Before_Stack=zeros(N_block,n_recs);
        CS.SIN.FLAG.Multilook_Incomplete=zeros(N_block,n_recs);
        CS.SIN.FLAG.Beam_Angle_Steering_Err=zeros(N_block,n_recs);
        CS.SIN.FLAG.AntiAliased_Power_Echo=zeros(N_block,n_recs);
        CS.SIN.FLAG.Auto_Beam_Steering=zeros(N_block,n_recs);
        CS.SIN.beam_param=zeros(16,N_block,n_recs);
        CS.SIN.coherence=zeros(n_points,N_block,n_recs);
        CS.SIN.phase_difference=zeros(n_points,N_block,n_recs);
        
end

%%%%%%%%% Time and Orbit Group Reading %%%%%%%%%%%%%%%%%%%%%%%%

CS.GEO.TAI.days(1:NrDATA)                          = floor(CS_nc.time_NN_ku/86400);
CS.GEO.TAI.secs(1:NrDATA)                          = floor(CS_nc.time_NN_ku - (CS.GEO.TAI.days(1:NrDATA)'*86400));
CS.GEO.TAI.microsecs(1:NrDATA)                     = (CS_nc.time_NN_ku - floor(CS_nc.time_NN_ku)).*1E6;
CS.GEO.USO(1:NrDATA)                               = CS_nc.uso_cor_NN_ku;
CS.GEO.SRC_CNT(1:NrDATA)                           = CS_nc.seq_count_NN_ku;
CS.GEO.BURST_CNT(1:NrDATA)                         = CS_nc.rec_count_NN_ku;
CS.GEO.LAT(1:NrDATA)                               = CS_nc.lat_NN_ku;
CS.GEO.LON(1:NrDATA)                               = CS_nc.lon_NN_ku;
CS.GEO.H(1:NrDATA)                                 = CS_nc.alt_NN_ku;
CS.GEO.H_rate(1:NrDATA)                            = CS_nc.orb_alt_rate_NN_ku;
CS.GEO.V.Vx(1:NrDATA)                              = CS_nc.sat_vel_vec_NN_ku(1,:);
CS.GEO.V.Vy(1:NrDATA)                              = CS_nc.sat_vel_vec_NN_ku(2,:);
CS.GEO.V.Vz(1:NrDATA)                              = CS_nc.sat_vel_vec_NN_ku(3,:);
CS.GEO.Beam.X(1:NrDATA)                            = CS_nc.beam_dir_vec_NN_ku(1,:);
CS.GEO.Beam.Y(1:NrDATA)                            = CS_nc.beam_dir_vec_NN_ku(2,:);
CS.GEO.Beam.Z(1:NrDATA)                            = CS_nc.beam_dir_vec_NN_ku(3,:);
CS.GEO.BaseLine.X(1:NrDATA)                        = CS_nc.inter_base_vec_NN_ku(1,:);
CS.GEO.BaseLine.Y(1:NrDATA)                        = CS_nc.inter_base_vec_NN_ku(2,:);
CS.GEO.BaseLine.Z(1:NrDATA)                        = CS_nc.inter_base_vec_NN_ku(3,:);
CS.GEO.ind_first_meas_20hz_01(:)                   = CS_nc.ind_first_meas_20hz_01;
CS.GEO.ind_meas_1hz_20_ku(1:NrDATA)                = CS_nc.ind_meas_1hz_NN_ku;

switch CS.GEO.OPERATION_MODE
    
    case {'SIR_L1B_LRM','SIR_L1B_FDM','SIR_L1B_SAR','SIR_L1B_SARIN'}
        
        CS.GEO.Star_Tracker_Usage(1:NrDATA)              = CS_nc.flag_instr_conf_rx_str_in_use_NN_ku;
        CS.GEO.Antenna_Bench_Roll(1:NrDATA)              = CS_nc.off_nadir_roll_angle_str_NN_ku;
        CS.GEO.Antenna_Bench_Pitch(1:NrDATA)             = CS_nc.off_nadir_pitch_angle_str_NN_ku;
        CS.GEO.Antenna_Bench_Yaw(1:NrDATA)               = CS_nc.off_nadir_yaw_angle_str_NN_ku;
        
end

CS.GEO.MODE_ID.Ins_Oper_Mode(1:NrDATA)             = CS_nc.flag_instr_mode_op_NN_ku;
dummy                                              = dec2bin(CS_nc.flag_instr_mode_flags_NN_ku,2);
CS.GEO.MODE_ID.Sarin_Degr_Case(1:NrDATA)           = str2num(dummy(:,1));
CS.GEO.MODE_ID.Cal4_Flag(1:NrDATA)                 = str2num(dummy(:,2));
CS.GEO.MODE_ID.Plat_Att_Ctrl(1:NrDATA)             = CS_nc.flag_instr_mode_att_ctrl_NN_ku;
CS.GEO.INS_CFG.Rx_Chain_Use(1:NrDATA)              = CS_nc.flag_instr_conf_rx_in_use_NN_ku;
CS.GEO.INS_CFG.Band_FLAG(1:NrDATA)                 = CS_nc.flag_instr_conf_rx_bwdt_NN_ku;
CS.GEO.INS_CFG.Tracking_Mode(1:NrDATA)             = CS_nc.flag_instr_conf_rx_trk_mode_NN_ku;
dummy                                              = dec2bin(abs(CS_nc.flag_instr_conf_rx_flags_NN_ku),8);
CS.GEO.INS_CFG.SIRAL_ID(1:NrDATA)                  = bin2dec(dummy(:,1));
CS.GEO.INS_CFG.Ext_Cal(1:NrDATA)                   = bin2dec(dummy(:,2));
CS.GEO.INS_CFG.Loop_Status(1:NrDATA)               = bin2dec(dummy(:,3));
CS.GEO.INS_CFG.Echo_Loss(1:NrDATA)                 = bin2dec(dummy(:,4));
CS.GEO.INS_CFG.Real_Time_Err(1:NrDATA)             = bin2dec(dummy(:,5));
CS.GEO.INS_CFG.Echo_Satur_Err(1:NrDATA)            = bin2dec(dummy(:,6));
CS.GEO.INS_CFG.Rx_Band_Atten(1:NrDATA)             = bin2dec(dummy(:,7));
CS.GEO.INS_CFG.Cycle_Report_Err(1:NrDATA)          = bin2dec(dummy(:,8));
CS.GEO.INS_CFG.STR_ATTREF(1:NrDATA)                = double(CS_nc.flag_instr_conf_rx_str_in_use_NN_ku > 0);

dummy                                              = dec2bin(abs(CS_nc.flag_mcd_NN_ku),32);
CS.GEO.MCD_FLAG.Block_Degraded(1:NrDATA)           = bin2dec(dummy(:,1));
CS.GEO.MCD_FLAG.Blank_Block(1:NrDATA)              = bin2dec(dummy(:,2));
CS.GEO.MCD_FLAG.Datation_Degraded(1:NrDATA)        = bin2dec(dummy(:,3));
CS.GEO.MCD_FLAG.Orbit_Propag_Err(1:NrDATA)         = bin2dec(dummy(:,4));
CS.GEO.MCD_FLAG.Orbit_File_Change(1:NrDATA)        = bin2dec(dummy(:,5));
CS.GEO.MCD_FLAG.Orbit_Discontinuity(1:NrDATA)      = bin2dec(dummy(:,6));
CS.GEO.MCD_FLAG.Echo_Saturation(1:NrDATA)          = bin2dec(dummy(:,7));
CS.GEO.MCD_FLAG.Other_Echo_Err(1:NrDATA)           = bin2dec(dummy(:,8));
CS.GEO.MCD_FLAG.Rx1_Err_SARin(1:NrDATA)            = bin2dec(dummy(:,9));
CS.GEO.MCD_FLAG.Rx2_Err_SARin(1:NrDATA)            = bin2dec(dummy(:,10));
CS.GEO.MCD_FLAG.Wind_Delay_Incon(1:NrDATA)         = bin2dec(dummy(:,11));
CS.GEO.MCD_FLAG.AGC_Incon(1:NrDATA)                = bin2dec(dummy(:,12));
CS.GEO.MCD_FLAG.CAL1_Corr_Miss(1:NrDATA)           = bin2dec(dummy(:,13));
CS.GEO.MCD_FLAG.CAL1_Corr_IPF(1:NrDATA)            = bin2dec(dummy(:,14));
CS.GEO.MCD_FLAG.DORIS_USO_Corr(1:NrDATA)           = bin2dec(dummy(:,15));
CS.GEO.MCD_FLAG.Complex_CAL1_Corr_IPF(1:NrDATA)    = bin2dec(dummy(:,16));
CS.GEO.MCD_FLAG.TRK_ECHO_Err(1:NrDATA)             = bin2dec(dummy(:,17));
CS.GEO.MCD_FLAG.RX1_ECHO_Err(1:NrDATA)             = bin2dec(dummy(:,18));
CS.GEO.MCD_FLAG.RX2_ECHO_Err(1:NrDATA)             = bin2dec(dummy(:,19));
CS.GEO.MCD_FLAG.NPM_Incon(1:NrDATA)                = bin2dec(dummy(:,20));

switch  CS.GEO.OPERATION_MODE
    
    case   {'SIR_L1B_LRM','SIR_L1B_FDM','SIR_L1B_SAR','SIR_L1B_SARIN'}
        
        CS.GEO.MCD_FLAG.CAL1_Corr_Type(1:NrDATA)         = bin2dec(dummy(:,21));
        CS.GEO.MCD_FLAG.Phase_Pertubation_Corr(1:NrDATA) = bin2dec(dummy(:,25));
        CS.GEO.MCD_FLAG.CAL2_Corr_Miss(1:NrDATA)         = bin2dec(dummy(:,26));
        CS.GEO.MCD_FLAG.CAL2_Corr_IPF(1:NrDATA)          = bin2dec(dummy(:,27));
        CS.GEO.MCD_FLAG.Power_Scaling_Err(1:NrDATA)      = bin2dec(dummy(:,28));
        CS.GEO.MCD_FLAG.Att_Corr_Miss(1:NrDATA)          = bin2dec(dummy(:,29));
        CS.GEO.MCD_FLAG.Phase_Pertubation_Mode(1:NrDATA) = bin2dec(dummy(:,32));
        
    case  {'SIR_FBR_SAR', 'SIR_FBR_SARIN'}
        
        CS.GEO.MCD_FLAG.CAL1_Corr_Type(1:NrDATA)         = bin2dec(dummy(:,21));
        CS.GEO.MCD_FLAG.Attitude_Corr_Missing(1:NrDATA)  = bin2dec(dummy(:,29));
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Measurements Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%

CS.MEA.win_delay(1:NrDATA)                         = CS_nc.window_del_NN_ku;
CS.MEA.Ho(1:NrDATA)                                = CS_nc.h0_applied_NN_ku;
CS.MEA.Trk_H_rate(1:NrDATA)                        = CS_nc.cor2_applied_NN_ku; %Note that in Cryo_L1b_read no scaling is applied
CS.MEA.LAI(1:NrDATA)                               = CS_nc.h0_lai_word_NN_ku;
CS.MEA.FAI(1:NrDATA)                               = CS_nc.h0_fai_word_NN_ku;
switch  CS.GEO.OPERATION_MODE
    case   {'SIR_L1B_LRM','SIR_L1B_FDM','SIR_L1B_SAR','SIR_L1B_SARIN'}
        CS.MEA.AGC_1(1:NrDATA)                             = CS_nc.agc_ch1_NN_ku;
        CS.MEA.AGC_2(1:NrDATA)                             = CS_nc.agc_ch2_NN_ku;
    case  {'SIR_FBR_SAR', 'SIR_FBR_SARIN'}
        CS.MEA.AGC_1(1:NrDATA)                             = CS_nc.agc_1_NN_ku;
        CS.MEA.AGC_2(1:NrDATA)                             = CS_nc.agc_2_NN_ku;
end
CS.MEA.Gain_Rx1(1:NrDATA)                          = CS_nc.tot_gain_ch1_NN_ku;
CS.MEA.Gain_Rx2(1:NrDATA)                          = CS_nc.tot_gain_ch2_NN_ku;
CS.MEA.Tx_Pwr(1:NrDATA)                            = CS_nc.transmit_pwr_NN_ku;
CS.MEA.dpl_range_corr(1:NrDATA)                    = CS_nc.dop_cor_NN_ku;
CS.MEA.ins_range_corr_rx_tx(1:NrDATA)              = CS_nc.instr_cor_range_tx_rx_NN_ku;
CS.MEA.ins_range_corr_rx(1:NrDATA)                 = CS_nc.instr_cor_range_rx_NN_ku;
CS.MEA.ins_gain_corr_rx_tx(1:NrDATA)               = CS_nc.instr_cor_gain_tx_rx_NN_ku;
CS.MEA.ins_gain_corr_rx(1:NrDATA)                  = CS_nc.instr_cor_gain_rx_NN_ku;
switch  CS.GEO.OPERATION_MODE
    case   {'SIR_L1B_SARIN','SIR_FBR_SARIN'}
        CS.MEA.int_phase_corr(1:NrDATA)                    = CS_nc.instr_int_ph_cor_NN_ku;
        CS.MEA.ext_phase_corr(1:NrDATA)                    = CS_nc.instr_ext_ph_cor_NN_ku;
        CS.MEA.phase_slope_corr(1:NrDATA)                  = CS_nc.ph_slope_cor_NN_ku;
end
CS.MEA.noise_power(1:NrDATA)                       = CS_nc.noise_power_NN_ku;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Corrections Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%

CS.COR.time_cor(1:NrCOR)                           = CS_nc.time_cor_01;
CS.COR.dry_trop(1:NrCOR)                           = CS_nc.mod_dry_tropo_cor_01;
CS.COR.wet_trop(1:NrCOR)                           = CS_nc.mod_wet_tropo_cor_01;
CS.COR.inv_bar(1:NrCOR)                            = CS_nc.inv_bar_cor_01;
CS.COR.dac(1:NrCOR)                                = CS_nc.hf_fluct_total_cor_01;
CS.COR.gim_ion(1:NrCOR)                            = CS_nc.iono_cor_gim_01;
CS.COR.model_ion(1:NrCOR)                          = CS_nc.iono_cor_01;               % some values are different wrt original file
CS.COR.ocean_equilibrium_tide(1:NrCOR)             = CS_nc.ocean_tide_01;
CS.COR.ocean_longperiod_tide(1:NrCOR)              = CS_nc.ocean_tide_eq_01;
CS.COR.ocean_loading_tide(1:NrCOR)                 = CS_nc.load_tide_01;
CS.COR.solidearth_tide(1:NrCOR)                    = CS_nc.solid_earth_tide_01;
CS.COR.geocentric_polar_tide(1:NrCOR)              = CS_nc.pole_tide_01;
CS.COR.surf_type(1:NrCOR)                          = CS_nc.surf_type_01;

dummy                                              = dec2bin(CS_nc.flag_cor_status_01,12);
CS.COR.corr_status.dry_trop(1:NrCOR)               = bin2dec(dummy(:,1));
CS.COR.corr_status.wet_trop(1:NrCOR)               = bin2dec(dummy(:,2));
CS.COR.corr_status.inv_bar(1:NrCOR)                = bin2dec(dummy(:,3));
CS.COR.corr_status.dac(1:NrCOR)                    = bin2dec(dummy(:,4));
CS.COR.corr_status.gim_iono(1:NrCOR)               = bin2dec(dummy(:,5));
CS.COR.corr_status.model_iono(1:NrCOR)             = bin2dec(dummy(:,6));
CS.COR.corr_status.ocean_equilibrium_tide(1:NrCOR) = bin2dec(dummy(:,7));
CS.COR.corr_status.ocean_longperiod_tide(1:NrCOR)  = bin2dec(dummy(:,8));
CS.COR.corr_status.ocean_loading_tide(1:NrCOR)     = bin2dec(dummy(:,9));
CS.COR.corr_status.solidearth_tide(1:NrCOR)        = bin2dec(dummy(:,10));
CS.COR.corr_status.geocentric_polar_tide(1:NrCOR)  = bin2dec(dummy(:,11));
CS.COR.corr_status.surface_type(1:NrCOR)           = bin2dec(dummy(:,12));

dummy                                              = dec2bin(CS_nc.flag_cor_err_01,12);
CS.COR.corr_error.dry_trop(1:NrCOR)                = bin2dec(dummy(:,1));
CS.COR.corr_error.wet_trop(1:NrCOR)                = bin2dec(dummy(:,2));
CS.COR.corr_error.inv_bar(1:NrCOR)                 = bin2dec(dummy(:,3));
CS.COR.corr_error.dac(1:NrCOR)                     = bin2dec(dummy(:,4));
CS.COR.corr_error.gim_iono(1:NrCOR)                = bin2dec(dummy(:,5));
CS.COR.corr_error.model_iono(1:NrCOR)              = bin2dec(dummy(:,6));
CS.COR.corr_error.ocean_equilibrium_tide(1:NrCOR)  = bin2dec(dummy(:,7));
CS.COR.corr_error.ocean_longperiod_tide(1:NrCOR)   = bin2dec(dummy(:,8));
CS.COR.corr_error.ocean_loading_tide(1:NrCOR)      = bin2dec(dummy(:,9));
CS.COR.corr_error.solidearth_tide(1:NrCOR)         = bin2dec(dummy(:,10));
CS.COR.corr_error.geocentric_polar_tide(1:NrCOR)   = bin2dec(dummy(:,11));
CS.COR.corr_error.surface_type(1:NrCOR)            = bin2dec(dummy(:,12));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch CS.GEO.OPERATION_MODE
    
    case {'SIR_L1B_LRM','SIR_L1B_FDM'}
        
        %%%%%%%%% LRM Average Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        CS.AVG.TAI.days(1:NrAVG)                         = floor(CS_nc.time_avg_01_ku/86400);
        CS.AVG.TAI.secs(1:NrAVG)                         = floor(CS_nc.time_avg_01_ku - (CS.AVG.TAI.days(1:NrAVG)'*86400));
        CS.AVG.TAI.microsecs(1:NrAVG)                    = (CS_nc.time_avg_01_ku - floor(CS_nc.time_avg_01_ku)).*1E6;
        CS.AVG.USO(1:NrAVG)                              = CS_nc.uso_cor_avg_01_ku;
        CS.AVG.lat(1:NrAVG)                              = CS_nc.lat_avg_01_ku;
        CS.AVG.lon(1:NrAVG)                              = CS_nc.lon_avg_01_ku;
        CS.AVG.H(1:NrAVG)                                = CS_nc.alt_avg_01_ku;
        CS.AVG.win_delay(1:NrAVG)                        = CS_nc.window_del_avg_01_ku;
        CS.AVG.data(:,1:NrAVG)                           = CS_nc.pwr_waveform_avg_01_ku;           % different wrt original file
        CS.AVG.echo_scaling(1:NrAVG)                     = CS_nc.echo_scale_factor_avg_01_ku;      % SF of 1e-09 applied, but not in original file
        CS.AVG.echo_scale_power(1:NrAVG)                 = CS_nc.echo_scale_pwr_avg_01_ku;
        CS.AVG.N_averaged_echoes(1:NrAVG)                = CS_nc.echo_numval_avg_01_ku;
        
        dummy                                            = dec2bin(abs(CS_nc.flag_echo_avg_01_ku),16);
        CS.AVG.OneHz_Echo_Err(1:NrAVG)                   = bin2dec(dummy(:,1));
        CS.AVG.Mispointing_Err(1:NrAVG)                  = bin2dec(dummy(:,16));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%% LRM Full Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        CS.LRM.data(:,1:NrDATA)                          = CS_nc.pwr_waveform_NN_ku;
        CS.LRM.echo_scaling(1:NrDATA)                    = CS_nc.echo_scale_factor_NN_ku;
        CS.LRM.echo_scale_power(1:NrDATA)                = CS_nc.echo_scale_pwr_NN_ku;
        CS.LRM.N_averaged_echoes(1:NrDATA)               = CS_nc.echo_numval_NN_ku;
        dummy                                            = dec2bin(CS_nc.flag_trk_cycle_NN_ku,16);
        CS.LRM.FLAG(1:NrDATA)                            = bin2dec(dummy(:,14:16));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case {'SIR_L1B_SAR'}
        
        %%%%%%%%% SAR Average Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        CS.AVG.TAI.days(1:NrAVG)                         = floor(CS_nc.time_avg_01_ku/86400);
        CS.AVG.TAI.secs(1:NrAVG)                         = floor(CS_nc.time_avg_01_ku - (CS.AVG.TAI.days(1:NrAVG)'*86400));
        CS.AVG.TAI.microsecs(1:NrAVG)                    = (CS_nc.time_avg_01_ku - floor(CS_nc.time_avg_01_ku)).*1E6;
        CS.AVG.USO(1:NrAVG)                              = CS_nc.uso_cor_avg_01_ku;
        CS.AVG.lat(1:NrAVG)                              = CS_nc.lat_avg_01_ku;
        CS.AVG.lon(1:NrAVG)                              = CS_nc.lon_avg_01_ku;
        CS.AVG.H(1:NrAVG)                                = CS_nc.alt_avg_01_ku;
        CS.AVG.win_delay(1:NrAVG)                        = CS_nc.window_del_avg_01_ku;
        CS.AVG.data(:,1:NrAVG)                           = CS_nc.pwr_waveform_avg_01_ku;
        CS.AVG.echo_scaling(1:NrAVG)                     = CS_nc.echo_scale_factor_avg_01_ku;
        CS.AVG.echo_scale_power(1:NrAVG)                 = CS_nc.echo_scale_pwr_avg_01_ku;
        CS.AVG.N_averaged_echoes(1:NrAVG)                = CS_nc.echo_numval_avg_01_ku;
        
        dummy                                            = dec2bin(abs(CS_nc.flag_echo_avg_01_ku),16);
        CS.AVG.OneHz_Echo_Err(1:NrAVG)                   = bin2dec(dummy(:,1));
        CS.AVG.Mispointing_Err(1:NrAVG)                  = bin2dec(dummy(:,16));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%% SAR Full Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        CS.SAR.data(:,1:NrDATA)                          = CS_nc.pwr_waveform_NN_ku;
        CS.SAR.echo_scaling(1:NrDATA)                    = CS_nc.echo_scale_factor_NN_ku;
        CS.SAR.echo_scale_power(1:NrDATA)                = CS_nc.echo_scale_pwr_NN_ku;
        CS.SAR.N_averaged_echoes(1:NrDATA)               = CS_nc.echo_numval_NN_ku;
        dummy                                            = dec2bin(abs(CS_nc.flag_echo_NN_ku),16);
        CS.SAR.FLAG.Approximate_Beam_Steering(1:NrDATA)  = bin2dec(dummy(:,1));
        CS.SAR.FLAG.Exact_Beam_Steering(1:NrDATA)        = bin2dec(dummy(:,2));
        CS.SAR.FLAG.Doppler_Weighting_Computed(1:NrDATA) = bin2dec(dummy(:,3));
        CS.SAR.FLAG.Doppler_Weighting_Applied_Before_Stack(1:NrDATA) = bin2dec(dummy(:,4));
        CS.SAR.FLAG.Multilook_Incomplete(1:NrDATA)       = bin2dec(dummy(:,5));
        CS.SAR.FLAG.Beam_Angle_Steering_Err(1:NrDATA)    = bin2dec(dummy(:,6));
        CS.SAR.FLAG.AntiAliased_Power_Echo(1:NrDATA)     = bin2dec(dummy(:,7));
        CS.SAR.FLAG.Auto_Beam_Steering(1:NrDATA)         = bin2dec(dummy(:,8));
        
        CS.SAR.beam_param(1,1:NrDATA)                    = CS_nc.stack_std_NN_ku;
        CS.SAR.beam_param(2,1:NrDATA)                    = CS_nc.stack_centre_NN_ku;
        CS.SAR.beam_param(3,1:NrDATA)                    = CS_nc.stack_scaled_amplitude_NN_ku;
        CS.SAR.beam_param(4,1:NrDATA)                    = CS_nc.stack_skewness_NN_ku;
        CS.SAR.beam_param(5,1:NrDATA)                    = CS_nc.stack_kurtosis_NN_ku;
        CS.SAR.beam_param(6,1:NrDATA)                    = CS_nc.stack_std_angle_NN_ku;
        CS.SAR.beam_param(7,1:NrDATA)                    = CS_nc.stack_centre_angle_NN_ku;
        CS.SAR.beam_param(8,1:NrDATA)                    = CS_nc.dop_angle_start_NN_ku;
        CS.SAR.beam_param(9,1:NrDATA)                    = CS_nc.dop_angle_stop_NN_ku;
        CS.SAR.beam_param(10,1:NrDATA)                   = CS_nc.look_angle_start_NN_ku;
        CS.SAR.beam_param(11,1:NrDATA)                   = CS_nc.look_angle_stop_NN_ku;
        CS.SAR.beam_param(12,1:NrDATA)                   = CS_nc.stack_number_after_weighting_NN_ku;
        CS.SAR.beam_param(13,1:NrDATA)                   = CS_nc.stack_number_before_weighting_NN_ku;
        CS.SAR.beam_param(14,1:NrDATA)                   = CS_nc.stack_centre_look_angle_NN_ku;
        CS.SAR.beam_param(15,1:NrDATA)                   = CS_nc.stack_gaussian_fitting_residuals_NN_ku;
        CS.SAR.beam_param(16,1:NrDATA)                   = CS_nc.stack_peakiness_NN_ku;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case {'SIR_FBR_SAR'}
        
        %%%%%%%%% FBR Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CS.FBR.data(:,:,1:NrDATA) = complex(CS_nc.cplx_waveform_ch1_i_NN_ku,CS_nc.cplx_waveform_ch1_q_NN_ku);
        CS.FBR.N_pulses(1:NrDATA) = CS_nc.echo_numval_NN_ku;
        CS.FBR.FLAG               = dec2bin(CS_nc.flag_echo_NN_ku,16);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case {'SIR_FBR_SARIN'}
        
        %%%%%%%%% FBR Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        CS.FBR.data_ch1(:,:,1:NrDATA) = complex(CS_nc.cplx_waveform_ch1_i_NN_ku,CS_nc.cplx_waveform_ch1_q_NN_ku);
        CS.FBR.data_ch2(:,:,1:NrDATA) = complex(CS_nc.cplx_waveform_ch2_i_NN_ku,CS_nc.cplx_waveform_ch2_q_NN_ku);
        CS.FBR.N_pulses(1:NrDATA)     = CS_nc.echo_numval_NN_ku;
        CS.FBR.FLAG                   = dec2bin(CS_nc.flag_echo_NN_ku,16);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case {'SIR_L1B_SARIN'}
        
        %%%%%%%%% SIN Average Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        CS.AVG.TAI.days(1:NrAVG)                          = floor(CS_nc.time_avg_01_ku/86400);
        CS.AVG.TAI.secs(1:NrAVG)                          = floor(CS_nc.time_avg_01_ku - (CS.AVG.TAI.days(1:NrAVG)'*86400));
        CS.AVG.TAI.microsecs(1:NrAVG)                     = (CS_nc.time_avg_01_ku - floor(CS_nc.time_avg_01_ku)).*1E6;
        CS.AVG.USO(1:NrAVG)                               = CS_nc.uso_cor_avg_01_ku;
        CS.AVG.lat(1:NrAVG)                               = CS_nc.lat_avg_01_ku;
        CS.AVG.lon(1:NrAVG)                               = CS_nc.lon_avg_01_ku;
        CS.AVG.H(1:NrAVG)                                 = CS_nc.alt_avg_01_ku;
        CS.AVG.win_delay(1:NrAVG)                         = CS_nc.window_del_avg_01_ku;
        CS.AVG.data(:,1:NrAVG)                            = CS_nc.pwr_waveform_avg_01_ku;
        CS.AVG.echo_scaling(1:NrAVG)                      = CS_nc.echo_scale_factor_avg_01_ku;
        CS.AVG.echo_scale_power(1:NrAVG)                  = CS_nc.echo_scale_pwr_avg_01_ku;
        CS.AVG.N_averaged_echoes(1:NrAVG)                 = CS_nc.echo_numval_avg_01_ku;
        
        dummy                                             = dec2bin(abs(CS_nc.flag_echo_avg_01_ku),16);
        CS.AVG.OneHz_Echo_Err(1:NrAVG)                    = bin2dec(dummy(:,1));
        CS.AVG.Mispointing_Err(1:NrAVG)                   = bin2dec(dummy(:,16));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%% SIN Full Waveform Group Reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        CS.SIN.data(:,1:NrDATA)                          = CS_nc.pwr_waveform_NN_ku;
        CS.SIN.echo_scaling(1:NrDATA)                    = CS_nc.echo_scale_factor_NN_ku;
        CS.SIN.echo_scale_power(1:NrDATA)                = CS_nc.echo_scale_pwr_NN_ku;
        CS.SIN.N_averaged_echoes(1:NrDATA)               = CS_nc.echo_numval_NN_ku;
        dummy                                            = dec2bin(abs(CS_nc.flag_echo_NN_ku),16);
        CS.SIN.FLAG.Approximate_Beam_Steering(1:NrDATA)  = bin2dec(dummy(:,1));
        CS.SIN.FLAG.Exact_Beam_Steering(1:NrDATA)        = bin2dec(dummy(:,2));
        CS.SIN.FLAG.Doppler_Weighting_Computed(1:NrDATA) = bin2dec(dummy(:,3));
        CS.SIN.FLAG.Doppler_Weighting_Applied_Before_Stack(1:NrDATA) = bin2dec(dummy(:,4));
        CS.SIN.FLAG.Multilook_Incomplete(1:NrDATA)       = bin2dec(dummy(:,5));
        CS.SIN.FLAG.Beam_Angle_Steering_Err(1:NrDATA)    = bin2dec(dummy(:,6));
        CS.SIN.FLAG.AntiAliased_Power_Echo(1:NrDATA)     = bin2dec(dummy(:,7));
        CS.SIN.FLAG.Auto_Beam_Steering(1:NrDATA)         = bin2dec(dummy(:,8));
        
        CS.SIN.beam_param(1,1:NrDATA)                    = CS_nc.stack_std_NN_ku;
        CS.SIN.beam_param(2,1:NrDATA)                    = CS_nc.stack_centre_NN_ku;
        CS.SIN.beam_param(3,1:NrDATA)                    = CS_nc.stack_scaled_amplitude_NN_ku;
        CS.SIN.beam_param(4,1:NrDATA)                    = CS_nc.stack_skewness_NN_ku;
        CS.SIN.beam_param(5,1:NrDATA)                    = CS_nc.stack_kurtosis_NN_ku;
        CS.SIN.beam_param(6,1:NrDATA)                    = CS_nc.stack_std_angle_NN_ku;
        CS.SIN.beam_param(7,1:NrDATA)                    = CS_nc.stack_centre_angle_NN_ku;
        CS.SIN.beam_param(8,1:NrDATA)                    = CS_nc.dop_angle_start_NN_ku;
        CS.SIN.beam_param(9,1:NrDATA)                    = CS_nc.dop_angle_stop_NN_ku;
        CS.SIN.beam_param(10,1:NrDATA)                   = CS_nc.look_angle_start_NN_ku;
        CS.SIN.beam_param(11,1:NrDATA)                   = CS_nc.look_angle_stop_NN_ku;
        CS.SIN.beam_param(12,1:NrDATA)                   = CS_nc.stack_number_after_weighting_NN_ku;
        CS.SIN.beam_param(13,1:NrDATA)                   = CS_nc.stack_number_before_weighting_NN_ku;
        CS.SIN.beam_param(14,1:NrDATA)                   = CS_nc.stack_centre_look_angle_NN_ku;
        CS.SIN.beam_param(15,1:NrDATA)                   = CS_nc.stack_gaussian_fitting_residuals_NN_ku;
        CS.SIN.beam_param(16,1:NrDATA)                   = CS_nc.stack_peakiness_NN_ku;
        CS.SIN.coherence(:,1:NrDATA)                     = CS_nc.coherence_waveform_NN_ku;
        CS.SIN.phase_difference(:,1:NrDATA)              = CS_nc.ph_diff_waveform_NN_ku;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
end

CS.GEO.V.V=sqrt(CS.GEO.V.Vx.^2+CS.GEO.V.Vy.^2+CS.GEO.V.Vz.^2);
CS.GEO.Serial_Sec_Num=CS.GEO.TAI.days.*24.*60.*60+CS.GEO.TAI.secs+CS.GEO.TAI.microsecs./1e6;

CS.COR.TOTAL_gim=CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.dac+CS.COR.gim_ion+CS.COR.ocean_equilibrium_tide+CS.COR.ocean_longperiod_tide+...
    CS.COR.ocean_loading_tide+CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide;

CS.COR.TOTAL_model=CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.dac+CS.COR.model_ion+CS.COR.ocean_equilibrium_tide+CS.COR.ocean_longperiod_tide+...
    CS.COR.ocean_loading_tide+CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide;

CS.GEO.Elapsed_Time=zeros(size(CS.GEO.Serial_Sec_Num));
CS.GEO.Elapsed_Time(CS.GEO.Serial_Sec_Num~=0)=CS.GEO.Serial_Sec_Num(CS.GEO.Serial_Sec_Num~=0)-CS.GEO.Start_Time;

switch CS.GEO.OPERATION_MODE
    
    case {'SIR_L1B_SAR'}
        
        CS.AVG.Serial_Sec_Num=CS.AVG.TAI.days.*24.*60.*60+CS.AVG.TAI.secs+CS.AVG.TAI.microsecs./1e6;
        CS.AVG.Elapsed_Time=CS.AVG.Serial_Sec_Num-CS.GEO.Start_Time;
               
    case {'SIR_L1B_SARIN'}
        
        CS.AVG.Serial_Sec_Num=CS.AVG.TAI.days.*24.*60.*60+CS.AVG.TAI.secs+CS.AVG.TAI.microsecs./1e6;
        CS.AVG.Elapsed_Time=CS.AVG.Serial_Sec_Num-CS.GEO.Start_Time;
        
end

