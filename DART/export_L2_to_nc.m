function export_L2_to_nc(CS1b, DATA_L2, FName, l1a_srcfile)
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
len_time = numel(DATA_L2.TIME);
time_dimid = netcdf.defDim(ncid,'time', len_time);

% variables
time_id = netcdf.defVar(ncid,'time','NC_DOUBLE',time_dimid);
lat_id = netcdf.defVar(ncid,'latitude','NC_DOUBLE',time_dimid);
lon_id = netcdf.defVar(ncid,'longitude','NC_DOUBLE',time_dimid);
alt_id = netcdf.defVar(ncid,'altitude','NC_DOUBLE',time_dimid);
% alt_rate_id = netcdf.defVar(ncid,'altitude_rate','NC_DOUBLE',time_dimid);
tracker_range_id = netcdf.defVar(ncid,'tracker_range_calibrated','NC_DOUBLE', time_dimid);
record_ind_id = netcdf.defVar(ncid,'record_ind','NC_INT', time_dimid);

swh_id = netcdf.defVar(ncid,'swh_ocean','NC_DOUBLE',time_dimid);
swh_qual_id = netcdf.defVar(ncid,'swh_ocean_qual','NC_DOUBLE',time_dimid);
% amplitude_id = netcdf.defVar(ncid,'amplitude_ocean','NC_DOUBLE',time_dimid);
% epoch_id = netcdf.defVar(ncid,'epoch_ocean','NC_DOUBLE',time_dimid);
range_id = netcdf.defVar(ncid,'range_ocean','NC_DOUBLE',time_dimid);
ssha_id = netcdf.defVar(ncid,'ssha_ocean','NC_DOUBLE',time_dimid);
% n_iter_id = netcdf.defVar(ncid,'num_iterations_ocean','NC_DOUBLE',time_dimid);
% misfit_id = netcdf.defVar(ncid,'misfit_ocean','NC_DOUBLE',time_dimid);
sig0_id = netcdf.defVar(ncid,'sig0_ocean','NC_DOUBLE',time_dimid);

if isfield(CS1b.MEA, 'p4_mode_flag')
    p4_mode_flag_id = netcdf.defVar(ncid,'p4_mode_flag','NC_DOUBLE', time_dimid);
    netcdf.putAtt(ncid, p4_mode_flag_id, 'flag_meanings', 'Acquisition LRM_CL LRM_OL LX_CL LX_OL LRMC_CL LRMC_OL LX2_CL TRANSPONDER Self_test CAL1_INSTR CAL1_LRM CAL1_SAR CAL1_RMC CAL2');
    netcdf.putAtt(ncid, p4_mode_flag_id, 'flag_values', '[ 5  6  7  8  9 10 11 12 13 14 15 16 17 18 19]');
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
netcdf.putVar(ncid,time_id,DATA_L2.TIME)
netcdf.putVar(ncid,lat_id,DATA_L2.LAT)
netcdf.putVar(ncid,lon_id,DATA_L2.LON)
netcdf.putVar(ncid,alt_id,DATA_L2.ALT)
% netcdf.putVar(ncid,alt_rate_id,CS1b.GEO.H_rate)
netcdf.putVar(ncid,tracker_range_id,CS1b.MEA.tracker_range(1:len_time))
if isfield(CS1b, 'IDXpoi') && (numel(CS1b.GEO.Elapsed_Time) == numel(CS1b.IDXpoi))
    netcdf.putVar(ncid,record_ind_id,CS1b.IDXpoi-1)  %make it 0-based
    netcdf.putAtt(ncid,record_ind_id,'comment', '0-based');
elseif isfield(CS1b, 'IDXpoi') && (rem(numel(CS1b.GEO.Elapsed_Time), numel(CS1b.IDXpoi)) == 0)
    n_rep = round(numel(DATA_L2.TIME) / numel(CS1b.IDXpoi));
    inds_rep = repelem(CS1b.IDXpoi - 1, n_rep);
    netcdf.putVar(ncid,record_ind_id,inds_rep(1:len_time))
    netcdf.putAtt(ncid,record_ind_id,'comment', '0-based');
end

netcdf.putVar(ncid,swh_id,DATA_L2.SWH)
if isfield(DATA_L2, 'swh_qual')
    netcdf.putVar(ncid,swh_qual_id,DATA_L2.qual_flag)
end
netcdf.putVar(ncid,range_id,DATA_L2.ALT' - DATA_L2.HEI)
netcdf.putVar(ncid,ssha_id, DATA_L2.HEI)
netcdf.putVar(ncid,sig0_id,DATA_L2.sigma0)

if isfield(CS1b.MEA, 'p4_mode_flag')
    netcdf.putVar(ncid,p4_mode_flag_id,CS1b.MEA.p4_mode_flag(1:len_time))
end

netcdf.close(ncid)

end