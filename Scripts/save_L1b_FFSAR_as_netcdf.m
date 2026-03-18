clear all
LoadCommonSettings

PathDATA = '/home/fehelers/PhD Delft/Projects/hydrology river/Data/L1b_Creuse/';
%PathDATA = '/home/fehelers/PhD Delft/Projects/hydrology river/Data/L1b_Garonne/';
fnames = dir([PathDATA '*.mat'])

for i=1:numel(fnames)
    fname = fnames(i).name

    load([PathDATA fname])
    
    tracker_range = CS1b.MEA.tracker_range;
    latitude = CS1b.GEO.LAT;
    longitude = CS1b.GEO.LON;
    altitude = CS1b.GEO.H;
    waveform_power = CS1b.SAR.data;
    waveform_index = 1:size(CS1b.GEO.LAT,2);
    range_gate_index = 1:size(waveform_power, 1);

    %CS1b.GEO.Start_Time;

    % Define the dimensions
    nx = size(waveform_power, 1);
    ny = size(waveform_power, 2);

    % Create a NetCDF file
    ncid = netcdf.create([PathDATA 'nc/' fname(1:end-4) '_TUDelft.nc'], 'NC_WRITE');
    
    % Define dimensions
    dimid_x = netcdf.defDim(ncid, 'range_gate_index', nx);
    dimid_y = netcdf.defDim(ncid, 'waveform_index', ny);
    dimid_time = netcdf.defDim(ncid, 'time', 26); % Define time dimension

    % Define variables
    varid_x = netcdf.defVar(ncid, 'range_gate_index', 'double', dimid_x);
    varid_y = netcdf.defVar(ncid, 'waveform_index', 'double', dimid_y);
    varid_time = netcdf.defVar(ncid, 'start_time', 'char', dimid_time);

    varid_tracker_range = netcdf.defVar(ncid, 'tracker_range', 'double', [dimid_y]);
    varid_lat = netcdf.defVar(ncid, 'latitude', 'double', [dimid_y]);
    varid_lon = netcdf.defVar(ncid, 'longitude', 'double', [dimid_y]);
    varid_alt = netcdf.defVar(ncid, 'altitude', 'double', [dimid_y]);
    varid_power = netcdf.defVar(ncid, 'waveform_power', 'double', [dimid_x, dimid_y]);


    netcdf.putAtt(ncid, varid_time, 'description', 'UTC Time of the overpass. Strictly speaking, this is the time corresponding to the waveform_index 0, but the satellite records the whole scene within a few seconds anyhow (7 km/s velocity).');
    netcdf.putAtt(ncid, varid_time, 'units', 'yyyy-MM-dd hh:mm:ss.sss');

    netcdf.putAtt(ncid, varid_x, 'description', 'Range gate index of the zero-padded waveforms');
    netcdf.putAtt(ncid, varid_x, 'units', 'integer');

    netcdf.putAtt(ncid, varid_y, 'description', 'This index represents the satellite positions / waveform number in along track direction. The posting rate of the initial FFSAR waveforms is roughly 1 m.');
    netcdf.putAtt(ncid, varid_y, 'units', 'integer');

    netcdf.putAtt(ncid, varid_tracker_range, 'description', 'The calibrated tracker range measured from the CoM of the satellite platform, as in the EUMETSAT L1a files, but interpolated.');
    netcdf.putAtt(ncid, varid_tracker_range, 'units', 'meter');

    netcdf.putAtt(ncid, varid_lat, 'description', 'The latitude of the satellite CoM.');
    netcdf.putAtt(ncid, varid_lat, 'units', 'degree North');

    netcdf.putAtt(ncid, varid_lon, 'description', 'The longitude of the satellite CoM.');
    netcdf.putAtt(ncid, varid_lon, 'units', 'degree East');

    netcdf.putAtt(ncid, varid_alt, 'description', 'The altitude of the satellite CoM above the WGS84 reference ellipsoid.');
    netcdf.putAtt(ncid, varid_alt, 'units', 'meter');

    netcdf.putAtt(ncid, varid_power, 'description', 'The FFSAR-processed power waveforms.');
    netcdf.putAtt(ncid, varid_power, 'units', 'intensity');

    % End definitions
    netcdf.endDef(ncid);

    % Write data to variables
    netcdf.putVar(ncid, varid_x, 1:nx);
    netcdf.putVar(ncid, varid_y, 1:ny);
    netcdf.putVar(ncid, varid_tracker_range, tracker_range);
    netcdf.putVar(ncid, varid_lat, latitude);
    netcdf.putVar(ncid, varid_lon, longitude);
    netcdf.putVar(ncid, varid_alt, altitude);
    netcdf.putVar(ncid, varid_power, waveform_power);
    netcdf.putVar(ncid, varid_time, char(CS1b.GEO.Start_Time)); % Write time value

    % Close the NetCDF file
    netcdf.close(ncid);
end

%% some test-reading
%clear all

power = ncread(['nc/' fname(1:end-4) '_TUDelft.nc'],"waveform_power");
lon = ncread(['nc/' fname(1:end-4) '_TUDelft.nc'],"longitude")
lat = ncread(['nc/' fname(1:end-4) '_TUDelft.nc'],"latitude")
%imagesc(log10(power))

%%
scatter(lon(:),lat(:),[],wse(:));hold on
