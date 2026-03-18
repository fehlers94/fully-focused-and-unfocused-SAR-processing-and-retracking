function export_L1b_to_nc(CS1b, FName, l1a_srcfile)
% export_L1b_to_nc exports the CS1b struct to destination netCDF file FName
% 
% according the following doc https://de.mathworks.com/help/matlab/import_export/exporting-to-network-common-data-form-netcdf-files.html

% init file
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLOBBER'));

if exist(FName, 'file')==2
  delete(FName);
end
ncid = netcdf.create(FName, mode);

%% definition mode
% dimensions
time_dimid = netcdf.defDim(ncid,'time', numel(CS1b.GEO.Elapsed_Time));
samples_ov_dimid = netcdf.defDim(ncid,'samples_ov', size(CS1b.SAR.data,1));
looks_dimid = netcdf.defDim(ncid,'looks', size(CS1b.SAR.stack_mask_start_stop, 1));
space_3d_dimid = netcdf.defDim(ncid,'space_3d', 3);

% variables
time_id = netcdf.defVar(ncid,'time','NC_DOUBLE',time_dimid);
power_waveform_id = netcdf.defVar(ncid,'power_waveform','NC_DOUBLE',[samples_ov_dimid time_dimid]);
power_waveform_opt_id = netcdf.defVar(ncid,'power_waveform_opt','NC_DOUBLE',[samples_ov_dimid time_dimid]);
power_waveform_opt2_id = netcdf.defVar(ncid,'power_waveform_opt2','NC_DOUBLE',[samples_ov_dimid time_dimid]);
power_waveform_pseudo_dda_id = netcdf.defVar(ncid,'power_waveform_pseudo_dda','NC_DOUBLE',[samples_ov_dimid time_dimid]);
lat_id = netcdf.defVar(ncid,'latitude','NC_DOUBLE',time_dimid);
lon_id = netcdf.defVar(ncid,'longitude','NC_DOUBLE',time_dimid);
alt_id = netcdf.defVar(ncid,'altitude','NC_DOUBLE',time_dimid);
alt_rate_id = netcdf.defVar(ncid,'altitude_rate','NC_DOUBLE',time_dimid);
vel_id = netcdf.defVar(ncid,'velocity_vector','NC_DOUBLE',[space_3d_dimid time_dimid]);
pitch_id = netcdf.defVar(ncid,'off_nadir_pitch_angle_pf','NC_DOUBLE', time_dimid);
roll_id = netcdf.defVar(ncid,'off_nadir_roll_angle_pf','NC_DOUBLE', time_dimid);
tracker_range_id = netcdf.defVar(ncid,'tracker_range_calibrated','NC_DOUBLE', time_dimid);
stack_mask_start_stop_id = netcdf.defVar(ncid,'stack_mask_start_stop','NC_DOUBLE', [looks_dimid time_dimid]);
look_angles_id = netcdf.defVar(ncid,'look_angles','NC_DOUBLE', [looks_dimid time_dimid]);
epoch_ref_gate_id = netcdf.defVar(ncid,'epoch_ref_gate_var','NC_DOUBLE', []);

if isfield(CS1b.MEA, 'PRI')
    pri_id = netcdf.defVar(ncid,'pulse_repetition_interval','NC_DOUBLE', time_dimid);
end
if isfield(CS1b, 'IDXpoi')
    record_ind_id = netcdf.defVar(ncid,'record_ind','NC_INT', time_dimid);
end
    
if isfield(CS1b.MEA, 'p4_mode_flag')
    p4_mode_flag_id = netcdf.defVar(ncid,'p4_mode_flag','NC_DOUBLE', time_dimid);
    netcdf.putAtt(ncid, p4_mode_flag_id, 'flag_meanings', 'Acquisition LRM_CL LRM_OL LX_CL LX_OL LRMC_CL LRMC_OL LX2_CL TRANSPONDER Self_test CAL1_INSTR CAL1_LRM CAL1_SAR CAL1_RMC CAL2');
    netcdf.putAtt(ncid, p4_mode_flag_id, 'flag_values', '[ 5  6  7  8  9 10 11 12 13 14 15 16 17 18 19]');
end

if isfield(CS1b.SAR.FLAG, 'cleanse_pdm')
    cleanse_pdm_id = netcdf.defVar(ncid,'cleanse_pdm_flag','NC_DOUBLE', time_dimid);
    netcdf.putAtt(ncid, cleanse_pdm_id, 'flag_meanings', 'Indicates, whether individual pulses during single-look formation, (0=good, 1=bad) ');
    netcdf.putAtt(ncid, cleanse_pdm_id, 'flag_values', '[0 1]');
end

% global attributes
netcdf.putAtt(ncid,netcdf.getConstant('GLOBAL'),'title','FF-SAR DART processor from TU Delft');
netcdf.putAtt(ncid,netcdf.getConstant('GLOBAL'),'contact','cornelis.slobbe@tu-delft.com');
if exist(l1a_srcfile)
    netcdf.putAtt(ncid,netcdf.getConstant('GLOBAL'), 'L1A source file', l1a_srcfile);
end

netcdf.putAtt(ncid,netcdf.getConstant('GLOBAL'),'first_measurement_time',string(CS1b.GEO.Start_Time));

os_zp = 1;
if isfield(CS1b.GEO, 'os_zp_factor')
    os_zp = CS1b.GEO.os_zp_factor;
    netcdf.putAtt(ncid,netcdf.getConstant('GLOBAL'),'os_zp_factor', int8(os_zp));
end

if isfield(CS1b.GEO, 'integration_time')
    netcdf.putAtt(ncid,netcdf.getConstant('GLOBAL'),'integration_time',CS1b.GEO.integration_time);
end

netcdf.putAtt(ncid,netcdf.getConstant('GLOBAL'),'epoch_ref_gate', int16(128*os_zp));

netcdf.endDef(ncid)

%% writing mode
% start_time = datetime(CS1b.GEO.Start_Time(1:end-1),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSS');
% times = start_time + seconds(CS1b.GEO.Elapsed_Time);
netcdf.putVar(ncid,time_id,CS1b.GEO.Elapsed_Time)
netcdf.putVar(ncid,power_waveform_id,CS1b.SAR.data)
if isfield(CS1b.SAR, 'data_opt')
    netcdf.putVar(ncid,power_waveform_opt_id,CS1b.SAR.data_opt)
end
if isfield(CS1b.SAR, 'data_opt2')
    netcdf.putVar(ncid,power_waveform_opt2_id,CS1b.SAR.data_opt2)
end
if isfield(CS1b.SAR, 'data_pseudoDD')
    is_pseudoDD_on = ~all(all(isnan(CS1b.SAR.data_pseudoDD)));
    
    if is_pseudoDD_on
        netcdf.putVar(ncid,power_waveform_pseudo_dda_id,CS1b.SAR.data_pseudoDD)
    end
end
netcdf.putVar(ncid,lat_id,CS1b.GEO.LAT)
netcdf.putVar(ncid,lon_id,CS1b.GEO.LON)
netcdf.putVar(ncid,alt_id,CS1b.GEO.H)
netcdf.putVar(ncid,alt_rate_id,CS1b.GEO.H_rate)
netcdf.putVar(ncid,vel_id,[CS1b.GEO.V.Vx;CS1b.GEO.V.Vy;CS1b.GEO.V.Vz])
netcdf.putVar(ncid,pitch_id,CS1b.GEO.Antenna_Bench_Pitch)
netcdf.putVar(ncid,roll_id,CS1b.GEO.Antenna_Bench_Roll)
netcdf.putVar(ncid,tracker_range_id,CS1b.MEA.tracker_range)
netcdf.putVar(ncid,stack_mask_start_stop_id,round(CS1b.SAR.stack_mask_start_stop))
netcdf.putVar(ncid,epoch_ref_gate_id,int16(128*os_zp))
% epoch_ref_gate_id

if iscell(CS1b.SAR.BeamAngle)
    look_angles = zeros([numel(CS1b.SAR.BeamAngle{1}), numel(CS1b.SAR.BeamAngle)]);
    for i=1:size(look_angles,2)
        look_angles(:,i) = CS1b.SAR.BeamAngle{1};
    end
    netcdf.putVar(ncid,look_angles_id,look_angles)
else
    netcdf.putVar(ncid,look_angles_id,CS1b.SAR.BeamAngle)
end

if isfield(CS1b, 'IDXpoi') && (numel(CS1b.GEO.Elapsed_Time) == numel(CS1b.IDXpoi))
    netcdf.putVar(ncid,record_ind_id,CS1b.IDXpoi - 1) %make it 0-based
    netcdf.putAtt(ncid,record_ind_id,'comment', '0-based');
elseif isfield(CS1b, 'IDXpoi') && (rem(numel(CS1b.GEO.Elapsed_Time), numel(CS1b.IDXpoi)) == 0)
    n_rep = numel(CS1b.GEO.Elapsed_Time) / numel(CS1b.IDXpoi);
    netcdf.putVar(ncid,record_ind_id,repelem(CS1b.IDXpoi - 1, n_rep)) %make it 0-based
    netcdf.putAtt(ncid,record_ind_id,'comment', '0-based');
    
end
    
if isfield(CS1b.MEA, 'PRI')
    netcdf.putVar(ncid,pri_id,CS1b.MEA.PRI)
end

if isfield(CS1b.MEA, 'p4_mode_flag')
    netcdf.putVar(ncid,p4_mode_flag_id,CS1b.MEA.p4_mode_flag)
end

if isfield(CS1b.SAR.FLAG, 'cleanse_pdm')
    netcdf.putVar(ncid,cleanse_pdm_id,CS1b.SAR.FLAG.cleanse_pdm)
end

netcdf.close(ncid)

end