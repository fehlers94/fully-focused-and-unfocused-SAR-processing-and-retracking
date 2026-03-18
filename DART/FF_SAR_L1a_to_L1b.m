function [CS1b,run_check] = FF_SAR_L1a_to_L1b(FName,DOM)

%CRYOSAT_SAR_FF_L1A_TO_L1B produces L1b data from CryoSat baseline C L1a
%SAR data based on a Fully Focused SAR processing algorithm

% close all

data_dir = fullfile([getenv('HOME') '/TUDTUM/ffsar_data/']);

%TODO: consider the following for future implementation into class
% integrate global TUD paths
% [~,CS1a] = Cryo_L1b_read(fullfile(PathDATA,'RadAlt','CryoSat','SIR_SAR_FR',FName));

% Read a polygon to minimize the data
% poly=dlmread('/data/Projects/Versatile_Hydrodynamics/Data/Tide_Gauge/GriddedData/Lake_IJssel/IJsselmeer.txt');

%% Settings
%General settings
LoadCommonSettings
% PathDATA = '/data/Projects/Versatile_Hydrodynamics/Data';

%Computation settings
defval('FName',[data_dir 'cs/l1a/CS_LTA__SIR1SAR_FR_20120529T025446_20120529T025502_C001.DBL'])
% defval('FName',[data_dir 's3a/l1a/crete/S3A_SR_1_SRA_A__20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement_l1a.nc'])
% defval('FName',[data_dir 's3b/l1a/crete/S3B_SR_1_SRA_A__20190620T083403_20190620T092432_20191107T073220_3029_026_335______MR1_R_NT_004.SEN3/measurement_l1a.nc'])
% defval('FName',[data_dir 's6a/S6A_GPP__P4__FR__1A_20210116T201350_20210116T201413_0001.NC'])

%region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
% defval('DOM',[52.3 53.3; 4.5 7.5]) 
defval('DOM',[77.2 78.7; 14.9 15.9]) % CS
% defval('DOM',[35.15 35.4; 23 24]) %crete transponder, lat 35.3379, lon 23.7795
% defval('DOM',[])

%Remaining settings
defval('MakePlots',false)                 %Make plots
defval('FntSz',14)                        %Set fontsize to be used in figures

% %Load Parameter Settings (Reference ellipsoid, SIRAL Characteristics, etc.)
% Constants_CryoSat;
% Constants_CryoSat_FF;
% t=t'; % fast time samples [s]

%Read a polygon to minimize the data
% poly=dlmread('/data/Projects/Versatile_Hydrodynamics/Data/Tide_Gauge/GriddedData/Lake_IJssel/IJsselmeer.txt');

%% Determine mission based on FName
mission = mission_from_fname(FName);

is_cs = strcmp(mission,'CS');
is_s3 = contains(mission,'S3');
is_s6 = contains(mission,'S6');

nc_grp='/';
if is_s6
    nc_grp='data_140/ku/';
end

%% Read & Crop CryoSat Level 1a SAR data
%Get indices to read cropped data (for all missions except CS)
if ~is_cs
    finfo=ncinfo(FName,nc_grp);
    for i=1:numel(finfo.Variables)
        if contains(finfo.Variables(i).Name,'lat')
            latvar=finfo.Variables(i).Name;
        elseif contains(finfo.Variables(i).Name,'lon')
            lonvar=finfo.Variables(i).Name;
        end
    end

    nc_lat=ncread(FName,[nc_grp latvar]);
    nc_lon=ncread(FName,[nc_grp lonvar]);
    if isempty(DOM)
        IDX_mask = true(numel(nc_lat),1);
    else
        if ischar(DOM)
            switch DOM
                case {'JakobsHavn','jakobshavn'}
                    %Load polygon that outlines the JakobsHavn glacier
                    load(fullfile(PathDATA,'Geography','jakobshavn_basin.mat'));
                    IDX_mask = inpolygon(nc_lon(:),nc_lat(:),polyLon,polyLat);
                otherwise
                    error('DOMain of interest not reckognized')
            end
        else
            IDX_mask = ingeoquad(nc_lat(:),nc_lon(:),DOM(1,:),DOM(2,:));
        end
    end

    IDX = find(IDX_mask);
    startLoc=min(IDX);
    count=max(IDX)-startLoc+1;
    IDX_mask=IDX_mask(IDX_mask);
end
clear nc_lat nc_lon


%Read data
% [~,CS1a] = Cryo_L1b_read(fullfile(PathDATA,'RadAlt','CryoSat','SIR_SAR_FR',FName));
switch mission
    case 'CS'
        [~,CS1a] = Cryo_L1b_read(FName);
    case {'S3A', 'S3B', 'S6A'}
        CS1a = S3_S6_L1a_read(FName,startLoc,count);
    otherwise
        error('mission: %s not recognized',mission)
end

%Load Parameter Settings (Reference ellipsoid, SIRAL Characteristics, etc.)
if is_cs
    DDAcf = DDA_ConfigFile(mission,'SAR',CS1a.GEO.Baseline_ID);
else
    DDAcf = DDA_ConfigFile(mission,'SAR');
end

%If segment crosses -180 meridian, wrap angle in degrees to [0 360]
if is_cs
    if max(abs(diff(CS1a.GEO.LON(CS1a.FBR.N_pulses(:) > 0)))) > 350
        CS1a.GEO.LON = wrapTo360(CS1a.GEO.LON);
    end
else
    if max(abs(diff(CS1a.GEO.LON))) > 350
        CS1a.GEO.LON = wrapTo360(CS1a.GEO.LON);
    end
end

%Correct window delay for DORIS USO drift (pp 20 CryoSat Product Handbook)
if is_cs
    if isequal(CS1a.GEO.Baseline_ID,'C')
        CS1a.MEA.win_delay = CS1a.MEA.win_delay.*(CS1a.GEO.USO+1);
    end
end

%Correct window delay for internal path delay
if is_cs
    CS1a.MEA.win_delay = CS1a.MEA.win_delay + (2*CS1a.MEA.ins_range_corr_rx_tx/DDAcf.c);
end

%Transform acquisition time to datenum format
if is_cs
    CS1a.('TIME') = datenum('2000','yyyy') + CS1a.GEO.TAI.days(:) + CS1a.GEO.TAI.secs(:)./86400 + CS1a.GEO.TAI.microsecs(:)./1e6./86400;
% else %do nothing here, because we have the time vector explicitely given
% in the netCDF file
end

%Crop data
if is_cs
    if isempty(DOM)
        IDX_mask = true(numel(CS1a.GEO.LAT),1);
    else
        if ischar(DOM)
            switch DOM
                case {'JakobsHavn','jakobshavn'}
                    %Load polygon that outlines the JakobsHavn glacier
                    load(fullfile(PathDATA,'Geography','jakobshavn_basin.mat'));
                    IDX_mask = inpolygon(CS1a.GEO.LON(:),CS1a.GEO.LAT(:),polyLon,polyLat);
                otherwise
                    error('DOMain of interest not reckognized')
            end
        else
            IDX_mask = ingeoquad(CS1a.GEO.LAT(:),CS1a.GEO.LON(:),DOM(1,:),DOM(2,:));
        end
    end
end

%% Data editing
if is_cs
    IDXfd = (CS1a.GEO.MCD_FLAG.Block_Degraded | CS1a.GEO.MCD_FLAG.Blank_Block | ...
        CS1a.GEO.MCD_FLAG.Datation_Degraded     | CS1a.GEO.MCD_FLAG.Orbit_Propag_Err | ...
        CS1a.GEO.MCD_FLAG.Orbit_File_Change     | CS1a.GEO.MCD_FLAG.Orbit_Discontinuity | ...
        CS1a.GEO.MCD_FLAG.Echo_Saturation       | CS1a.GEO.MCD_FLAG.Other_Echo_Err | ...
        CS1a.GEO.MCD_FLAG.Rx1_Err_SARin         | CS1a.GEO.MCD_FLAG.Rx2_Err_SARin | ...
        CS1a.GEO.MCD_FLAG.Wind_Delay_Incon      | CS1a.GEO.MCD_FLAG.AGC_Incon | ...
        CS1a.GEO.MCD_FLAG.CAL1_Corr_Miss        | CS1a.GEO.MCD_FLAG.CAL1_Corr_IPF | ...
        CS1a.GEO.MCD_FLAG.DORIS_USO_Corr        | CS1a.GEO.MCD_FLAG.Complex_CAL1_Corr_IPF | ...
        CS1a.GEO.MCD_FLAG.TRK_ECHO_Err          | CS1a.GEO.MCD_FLAG.RX1_ECHO_Err | ...
        CS1a.GEO.MCD_FLAG.RX2_ECHO_Err          | CS1a.GEO.MCD_FLAG.NPM_Incon | ...
        CS1a.GEO.MCD_FLAG.Attitude_Corr_Missing | CS1a.GEO.MCD_FLAG.CAL1_Corr_Type);

    CS1a.GEO.IDXfd = ~IDX_mask | IDXfd(:);
elseif is_s6
    CS1a.GEO.IDXfd = ~IDX_mask | ~logical(CS1a.mcd_flags); %consider any set flag (value~=0) as bad
else
    CS1a.GEO.IDXfd = ~IDX_mask;
end
clear('IDX_mask','IDX','IDXfd')

if nnz(~CS1a.GEO.IDXfd) <= 50, run_check=0; CS1b=1; return; end
run_check=1;

%% Calibration
if is_cs
    CS1a = Calibrate_CryoSat_FBRdata(CS1a,DDAcf);
elseif is_s3
    CS1a = Calibrate_S3_L1Adata(CS1a,DDAcf);
elseif is_s6
    CS1a = Calibrate_S6_L1Adata(CS1a,DDAcf);
end

%% Preliminaries
% Some settings for looks
IDX1=true(length(CS1a.GEO.LAT),1);
if is_cs
    IDX1=CS1a.FBR.N_pulses ~= 0;
end
lat=CS1a.GEO.LAT(IDX1);
lon=CS1a.GEO.LON(IDX1);

% T=0.886*c*alt_nom/(2*f_c*v_nom*res_a); % integration time [s]
if is_cs
    T=1.9;%0.886*DDAcf.c*alt_nom/(2*f_c*median(CS1a.GEO.V.V(:))*res_a); % integration time [s]
elseif is_s3
    T=2.1;
elseif is_s6
    T=0.886*DDAcf.c*DDAcf.alt_nom/(2*DDAcf.fc*DDAcf.v_nom*DDAcf.res_a); %3.6307;
end
L=T*median(CS1a.GEO.V.V(IDX1))/2;
H=median(CS1a.GEO.H(IDX1));
res_a=(((L^2+H^2)^0.5+DDAcf.lambda0/4)^2-H^2)^0.5-L;
Re_local=rsphere('euler',lat(1),lon(1),lat(end),lon(end),[DDAcf.RefEll.SemimajorAxis DDAcf.RefEll.Eccentricity]);
Nl=median(CS1a.GEO.V.V(:))/(res_a*((Re_local+median(CS1a.GEO.H(:)))/Re_local)); % number of looks per seconds [1/s]
clear lat lon L H res_a Re_local

t_res=1/DDAcf.B; % time resolution of the waveform [s]
r_res=t_res*DDAcf.c/2; % range resolution of the waveform [m]
r_res_zp=r_res/DDAcf.os_ZP;

if DDAcf.ApplyWindow
    if strcmpi(DDAcf.Window,'hamming')
        W = hamming(DDAcf.Np,'periodic')';
    elseif strcmpi(DDAcf.Window,'hanning')
        W = hanning(DDAcf.Np)';
    else
        error('Window not reckognized!');
    end
else
    W = ones(1,DDAcf.Np);
end

% define fast-time vector
if ~is_s6
    Ts=DDAcf.tau_u/DDAcf.Np;
else
    Ts=DDAcf.B/DDAcf.Bt*DDAcf.tau_u/DDAcf.Np;
end
t=0:Ts:(DDAcf.Np-1)*Ts; % fast time samples [s]
t=t-max(t)*(0.5);
if is_s6
    pri = CS1a.MEA.PRI';
else
    pri = DDAcf.PRI;
end
tp_vec=(-(DDAcf.Nb-1)/2:1:(DDAcf.Nb-1)/2)'*pri; % vector of pulse timing [s]    

%% Interpolate fields to ground locations
if is_cs
    t_burst=CS1a.GEO.Elapsed_Time(~CS1a.GEO.IDXfd); % original time vector
else
    t_burst=CS1a.TIME(~CS1a.GEO.IDXfd);
end
t_burst=t_burst-t_burst(1);
    
if is_s6
    BRI=CS;
else
    BRI=DDAcf.BRI;
end
t_burst_corr=(0:length(t_burst)-1)'*BRI; % corrected time vector

%Construct array of valid entries of fields LAT, LON, H, H_rate, Vx, Vy,
%Vz, V, TAI, win_delay, BaseLine, Beam
ALLb = [CS1a.GEO.LAT(:),CS1a.GEO.LON(:),CS1a.GEO.H(:),CS1a.GEO.H_rate(:),CS1a.GEO.V.Vx(:),CS1a.GEO.V.Vy(:),CS1a.GEO.V.Vz(:),CS1a.GEO.V.V(:),CS1a.GEO.TAI.days(:),CS1a.GEO.TAI.secs(:),CS1a.GEO.TAI.secs(:).*1E6,CS1a.MEA.win_delay(:),CS1a.GEO.BaseLine.X(:),CS1a.GEO.BaseLine.Y(:),CS1a.GEO.BaseLine.Z(:),CS1a.GEO.Beam.X(:),CS1a.GEO.Beam.Y(:),CS1a.GEO.Beam.Z(:)];
ALLb = ALLb(~CS1a.GEO.IDXfd,:);

%Compute the piecewise polynomial forms of the cubic spline interpolants to
%the data values ALLc at latitudes LATc
dt=1/Nl; % time between looks

% TODO check here why there is +T/2 -T/2
if is_s6
    t_l=t_burst(1):dt:t_burst(end); % time vector for analysis location interpolation    
else
    t_l=t_burst(1)+T/2:dt:t_burst(end)-T/2; % time vector for analysis location interpolation    
end
pp   = spline(t_burst,ALLb');
ALLl = ppval(pp,t_l);

% IN=inpolygon(ALLl(2,:),ALLl(1,:),poly(:,1),poly(:,2));
% ALLl =ALLl(:,IN);

% lon_l=ALLl(2,:)';
% lat_l=ALLl(1,:)';
% alt_l=ALLl(3,:)';

% read relevant data for FF-SAR
lon_b=CS1a.GEO.LON(~CS1a.GEO.IDXfd); 
lat_b=CS1a.GEO.LAT(~CS1a.GEO.IDXfd); 
alt_b=CS1a.GEO.H(~CS1a.GEO.IDXfd); 
VX_b=CS1a.GEO.V.Vx(~CS1a.GEO.IDXfd); 
VY_b=CS1a.GEO.V.Vy(~CS1a.GEO.IDXfd); 
VZ_b=CS1a.GEO.V.Vz(~CS1a.GEO.IDXfd); 
% roll_b=CS1a.GEO.BaseLine.X(~CS1a.GEO.IDXfd);
% pitch_b=-CS1a.GEO.Beam.Y(~CS1a.GEO.IDXfd);
LAI=CS1a.MEA.LAI(~CS1a.GEO.IDXfd); 
if is_cs
    wd_b=CS1a.MEA.win_delay(~CS1a.GEO.IDXfd);
    Rtrk_b=wd_b*DDAcf.c/2; % one-way window delay in meters
else
    Rtrk_b=CS1a.MEA.ref_range;
end
wav=double(CS1a.FBR.data(:,:,~CS1a.GEO.IDXfd)); % waveforms

% make models for latitude/longitude/range/altitude
% % dt=1/Nl; % time between looks
% % t_l=t_burst(1)+T/2:dt:t_burst(end)-T/2; % time vector for analysis location interpolation
% % lon_l  = polyval(polyfit(t_burst,lon_b,2),t_l);
% % lat_l  = polyval(polyfit(t_burst,lat_b,2),t_l);
% % alt_l  = polyval(polyfit(t_burst,alt_b,2),t_l);
% % Rtrk_l = polyval(polyfit(t_burst,Rtrk_b,2),t_l);
% % h_l=alt_l-Rtrk_l; % estimates of elevation
A=ones(length(t_burst),5); % matrix A is the design matrix for a n-order polynomial through bursts
A(:,2)=t_burst;
A(:,3)=t_burst.^2;
A(:,4)=t_burst.^3;
A(:,5)=t_burst.^4;
lon_mod=(A'*A)^-1*A'*lon_b;
lat_mod=(A'*A)^-1*A'*lat_b;
alt_mod=(A'*A)^-1*A'*alt_b;
Rtrk_mod=(A'*A)^-1*A'*Rtrk_b;

% estimate ground locations
B=ones(length(t_l),5); % matrix B evaluates the models at the 'look' location
B(:,2)=t_l;
B(:,3)=t_l.^2;
B(:,4)=t_l.^3;
B(:,5)=t_l.^4;
lon_l=B*lon_mod;
lat_l=B*lat_mod;
alt_l=B*alt_mod;                        
Rtrk_l=B*Rtrk_mod; % this will only work for relatively flat surfaces (peak-to-peak max 30-40 m)

% some alternative pchip interpolation for Rtrk_l (through the center times
% of the tracker range settings to mitigate the problems with relatively bumpy surfaces / steep slopes)
trk_change_flag = (diff(Rtrk_b) ~= 0);
if sum(trk_change_flag) > 4
    IDXtrk = find(trk_change_flag);
    t_b_trk = t_burst(IDXtrk);
    if ~is_s6
        t_b_center = 0.5 * (t_b_trk(1:end-1) + t_b_trk(2:end) + DDAcf.BRI); % BRI/2 needs to be added for centering because of the 'righthanded' diff operation
    else
        t_b_center = 0.5 * (t_b_trk(1:end-1) + t_b_trk(2:end) + BRI); % BRI/2 needs to be added for centering because of the 'righthanded' diff operation
    end
    Rtrk_center = Rtrk_b(IDXtrk(2:end));
    trk_interpolator   = pchip(t_b_center,Rtrk_center);
    Rtrk_l = ppval(trk_interpolator,t_l)';
    
%     h_l = alt_l - Rtrk_l;
%     figure; plot(lat_b,alt_b-Rtrk_b);hold on; 
%     %plot(lat_l,h_l);hold on; 
%     %plot(lat_l,h_l_new); hold on; 
%     %scatter(lat_b(qq),alt_b(qq)-Rtrk_b(qq),'filled'); hold on; 
%     plot(lat_l,h_l); hold on;
%     %legend('alt_b-Rtrk_b','Marcel: h_l=alt_l-Rtrk_l','Cornelis: spline');
%     legend('alt_b-Rtrk_b','pchip interpolation');
%     xlabel('latitude')
%     ylabel('height')
end


% check if ground locations are inside polygon
% IN=inpolygon(lon_l,lat_l,poly(:,1),poly(:,2));
IN=inpolygon(lon_l,lat_l,DOM(2,:),DOM(1,:)); % instead of polygon, use DOM
lon_l=lon_l(IN);
lat_l=lat_l(IN);
alt_l=alt_l(IN);
Rtrk_l=Rtrk_l(IN);
t_l=t_l(IN);
h_l=alt_l-Rtrk_l; % estimates of elevation
% h_l=h_l(IN);

%% Create empty struct L1b format
[~,CS1b] = Cryo_L1b_struct('SIR_SAR_L1.DBL',numel(h_l),DDAcf,1);
CS1b.SAR.data = squeeze(CS1b.SAR.data);
CS1b.SAR.data = zeros(DDAcf.Np*DDAcf.os_ZP, numel(h_l)); % required because S6 uses 256 bins per waveform, which is not considered in Cryo_L1b_struct func
CS1b.SAR.beam_param = squeeze(CS1b.SAR.beam_param);
CS1b.GEO.LAT       = ALLl(1,:);      CS1b.GEO.LON           = ALLl(2,:);
CS1b.GEO.H         = ALLl(3,:);      CS1b.MEA.win_delay     = Rtrk_l'/(DDAcf.c/2); %ATTENTION: We give values for Rtrk_l for each look, but below we stack looks together in bunches of 50, which means they are all referenced to the same, centered Rtrk_l. This causes potential errors.
CS1b.GEO.H_rate    = ALLl(4,:);      CS1b.GEO.V.Vx          = ALLl(5,:);
CS1b.GEO.V.Vy      = ALLl(6,:);      CS1b.GEO.V.Vz          = ALLl(7,:);
CS1b.GEO.V.V       = ALLl(8,:);      CS1b.GEO.TAI.days      = ALLl(9,:);
CS1b.GEO.TAI.secs  = ALLl(10,:);     CS1b.GEO.TAI.microsecs = ALLl(11,:);
CS1a.GEO.BaseLine.X= ALLl(13,:);     CS1a.GEO.BaseLine.Y    = ALLl(14,:);
CS1a.GEO.BaseLine.Z= ALLl(15,:);     CS1a.GEO.Beam.X        = ALLl(16,:);
CS1a.GEO.Beam.Y    = ALLl(17,:);     CS1a.GEO.Beam.Z        = ALLl(18,:);
CS1b.GEO.Start_Time    = CS1a.GEO.Start_Time;
CS1b.GEO.Baseline_ID   = 'OWN';
%Set flag CS1b.GEO.MCD_FLAG.Block_Degraded;
CS1b.GEO.MCD_FLAG.Block_Degraded = double(isnan(CS1b.GEO.LAT));
%Compute roll, pitch, and yaw angles of the antenna bench (p. 16 of the CryoSat product handbook)
CS1b.GEO.Antenna_Bench_Roll  = rad2deg(CS1b.GEO.BaseLine.X);
CS1b.GEO.Antenna_Bench_Yaw   = rad2deg(-CS1b.GEO.BaseLine.Y); 
CS1b.GEO.Antenna_Bench_Pitch = rad2deg(-CS1b.GEO.Beam.Y);

% %% Create an output structure
% [CS1b] = Cryo_FF_L1b_struct(length(lon_l));
% CS1b.MEA.trk=Rtrk_l;
% CS1b.GEO.alt=alt_l;
% CS1b.GEO.lon=lon_l;
% CS1b.GEO.lat=lat_l;

% %% Interpolate Geophysical corrections to ground data
% % interpolation from the 1hz latitude to the look location latitude
% lat_1hz=mean(CS1a.GEO.LAT);
% CS1b.COR.mod=interp1(lat_1hz,CS1a.COR.model_ion,lat_l,'linear');
% CS1b.COR.gim=interp1(lat_1hz,CS1a.COR.gim_ion,lat_l,'linear');
% CS1b.COR.wet=interp1(lat_1hz,CS1a.COR.wet_trop,lat_l,'linear');
% CS1b.COR.dry=interp1(lat_1hz,CS1a.COR.dry_trop,lat_l,'linear');
% CS1b.COR.inv=interp1(lat_1hz,CS1a.COR.inv_bar,lat_l,'linear');
% CS1b.COR.loa=interp1(lat_1hz,CS1a.COR.ocean_loading_tide,lat_l,'linear');
% CS1b.COR.lon=interp1(lat_1hz,CS1a.COR.ocean_longperiod_tide,lat_l,'linear');
% CS1b.COR.sol=interp1(lat_1hz,CS1a.COR.solidearth_tide,lat_l,'linear');
% CS1b.COR.equ=interp1(lat_1hz,CS1a.COR.ocean_equilibrium_tide,lat_l,'linear');
% CS1b.COR.pol=interp1(lat_1hz,CS1a.COR.geocentric_polar_tide,lat_l,'linear');
% CS1b.GEO.TIME=interp1(CS1a.GEO.LAT(CS1a.GEO.LAT ~= 0),CS1a.TIME(CS1a.GEO.LAT ~= 0),lat_l,'linear'); % 20 hz
% CS1b.GEO.pitch=interp1(lat_b,pitch_b,lat_l,'linear'); % 20 hz
% CS1b.GEO.roll=interp1(lat_b,roll_b,lat_l,'linear'); % 20 hz

% satellite location at the zenith of the look locations
% compute shortest distance between track and look location
[x_l,y_l,z_l]=geodetic2enu(lat_l,lon_l,alt_l,lat_l,lon_l,h_l,wgs84Ellipsoid,'degrees');
R0_i=sqrt(x_l.^2+y_l.^2+z_l.^2);
clear x_l y_l z_l

disp('Number of looks:')
disp(length(h_l))

%% Estimate ranges and velocities for the target locations
% find nearest sub-satellite point to the transponder
if ~is_cs
    d=distance(TR.crete.lat,TR.crete.lon,lat_l,lon_l,[CONST.a CONST.e]);
else
    d=distance(TR.svalbard.lat,TR.svalbard.lon,lat_l,lon_l,[CONST.a CONST.e]);
end
[~,Id2]=min(d);

% h_l(Id2) = TR.crete.h;

% go through the nadir locations
% q2=0;
tic
if Id2<25
    Id2 = 25;
end
for i=Id2
% for i=25:50:length(h_l)
    fprintf('processing nadir look location i=%d ...\n',i);
    
    %if mod(i,1000)== 0
    %    fprintf('%d percent\n',i/length(h_l)*100)
    %end
% profile on -detail builtin -history -timer real
%     % satellite locations at burst in local reference frame of location 'l(i)'
%     [x_b,y_b,z_b]=geodetic2enu(lat_b,lon_b,alt_b,lat_l(i),lon_l(i),h_l(i),wgs84Ellipsoid,'degrees');
    
    % selection of bursts (radius of T second around POCA)
    q=abs(t_burst-t_l(i))<T/2;
    %q=find((t_burst-t_l(i)>-T/2 & t_burst-t_l(i) < 0)==1);
    %q=find((t_burst-t_l(i)<T/2 & t_burst-t_l(i) > 0)==1);
%     if sum(q==q2) ~= length(q)
        t_select=t_burst(q);
        t_select_corr=t_burst_corr(q);
        LAI_select=LAI(q);
        VX_select=VX_b(q);
        VY_select=VY_b(q);
        VZ_select=VZ_b(q);
        Rtrk_select=Rtrk_b(q);
        wav_select=wav(:,:,q);
        % now if s6 is used, the pulse repetition can vary and so tp_vec
        % becomes a matrix which needs to be accounted for
        if ~isvector(tp_vec)
            tp_vec_select = tp_vec(:,q);
        else
            tp_vec_select = tp_vec;
        end
%     end
%     x_select=x_b(q);
%     y_select=y_b(q);
%     z_select=z_b(q);

    % satellite locations at burst in local reference frame of location 'l(i)'
    [x_select,y_select,z_select]=geodetic2enu(lat_b(q),lon_b(q),alt_b(q),lat_l(i),lon_l(i),h_l(i),wgs84Ellipsoid,'degrees');

    % convert velocity vector into local reference frame
    [vx_b,vy_b,vz_b]=ecef2enuv(VX_select,VY_select,VZ_select,lat_l(i),lon_l(i),'degrees');
    
    % reference the time to the POCA time
    t_select=t_select-t_l(i);
    t_select_corr=t_select_corr-t_l(i);
    
    % (corrected) pulse timing
%     t_pulse=repmat(t_select_corr',Np,1)+repmat(tp_vec',1,length(t_select));
    t_pulse=t_select_corr'+tp_vec_select;
    t_pulse=t_pulse(:);
%     t_pulse=reshape(t_pulse,1,DDAcf.Nb,[]);
   
    % rearrange waveforms (for all pulses/bursts consecutively)
    wav_select=transpose(reshape(wav_select, [DDAcf.Np, DDAcf.Nb*length(t_select)]));
    
    % Apply window
    if DDAcf.ApplyWindow
        wav_select=wav_select.*W;
    end
    
    % compute ranges at burst locations
    R0=sqrt(x_select.^2+y_select.^2+z_select.^2); % estimated ranges to center of strip (based on satellite altitude and window delay)
    
    % radial velocity variation in slow time at every burst (from target positive)
    vr_0=(vx_b.*x_select+vy_b.*y_select+vz_b.*z_select)./R0;
    
    % interpolate radii and radial velocity to pulse locations (this also
    % overcomes the datation issues)
    C=ones(length(t_select),5); % design matrix C to compute models at the burst locations
    C(:,2)=t_select;
    C(:,3)=t_select.^2;
    C(:,4)=t_select.^3;
    C(:,5)=t_select.^4;
    R_mod=(C'*C)^-1*C'*R0;
    vr_mod=(C'*C)^-1*C'*vr_0;
    D=ones(numel(t_pulse),5); % matrix D evaluates the models at the pulse locations
    D(:,2)=t_pulse(:);
    D(:,3)=t_pulse(:).^2;
    D(:,4)=t_pulse(:).^3;
    D(:,5)=t_pulse(:).^4;
    R0_p=D*R_mod; % ranges to center to strip at pulses
    vr0_p=D*vr_mod; % radial velocities at pulses
%     R0_p=reshape(D*R_mod,size(t_pulse)); % ranges to center to strip at pulses
%     vr0_p=reshape(D*vr_mod,size(t_pulse)); % radial velocities at pulses
%     % interpolate radii and radial velocity to pulse locations (this also
%     % overcomes the datation issues)
%     R_mod  = polyfit(t_select,R0,4);
%     R0_p   = polyval(R_mod,t_pulse);   % ranges to center to strip at pulses
%     vr_mod = polyfit(t_select,vr_0,4);
%     vr0_p  = polyval(vr_mod,t_pulse);  % radial velocities at pulses
    
    % the window delay (tracker range) is kept constant over 64 pulses
    Rtrk=repmat(Rtrk_select',DDAcf.Nb,1);
    Rtrk=Rtrk(:);
%     Rtrk=reshape(Rtrk_select,1,1,[]);

    %{
    figure;
    if is_s6
        imagesc(abs(fft(transpose(wav_select))));
    else
        imagesc(abs(fftshift(fft(transpose(wav_select)),1)));
    end
    xlabel('pulse number');
    ylabel('bin');
    %title('Be aware there is an additional shift of one bin, because of the fftshift function')
    %}
    
    %% Fully Focused SAR algowav_selectwav_selectrithm   
    % RCMC, takes care of relative motion between target and satellite
    % (Doppler and Range) (see Egido IIIA)
    switch mission
        case 'CS'
            nTs=2;
        case {'S3A','S3B'}
            nTs=0.3;
        case {'S6A'}
            nTs=-0.5;
    end
    f_D=2*DDAcf.fc*vr0_p/DDAcf.c;
    RCMC=exp(-2*pi*1i*(2*DDAcf.s*(R0_p)/DDAcf.c-f_D)*(t+Ts*nTs));
    CTRP=exp(2*pi*1i*2*DDAcf.s*(Rtrk)/DDAcf.c*(t+Ts*nTs));
    wav_select=wav_select.*RCMC.*CTRP;
    
    % zero-padding (set zp to 2 if you want to zero-pad)
    wav_select=padarray(wav_select,[0 (DDAcf.os_ZP-1)*DDAcf.Np/2],0,'both'); 
    
    % Fast-time fft, range compression (see Egido IIIB)
    wav_select=fft(wav_select,[],2);
    if ~is_s6
        wav_select=fftshift(wav_select,2);
    end
    
%     %{
    figure
    imagesc(transpose(abs(wav_select)));
    xlabel('pulse number')
    ylabel('bin')
%     %}

    % Remove 'known' phase jumps of 0.5 pi (after every burst and additionally after every LAI change 
    if ~is_s6
        diff_LAI=[0; Rtrk_select(2:end)-Rtrk_select(1:end-1)]; % absolute removed w.r.t. earlier versions
        phase_LAI=round(cumsum(diff_LAI)'/(1/80E6*DDAcf.c/2))*1.18*pi; %(1/80E6*c/2) is basically the tracker range jumps than can be made by the instrument
        
        intra_burst_corr=1;
        if is_cs
            intra_burst_corr=0.5;           
        end
        
        phase_jumps=intra_burst_corr*pi*(0:length(t_select)-1)+phase_LAI-mean(phase_LAI);
        
        phase_jumps=repmat(phase_jumps,DDAcf.Nb,1);
        phasor_jumps=exp(-1i*phase_jumps(:));
        wav_select=wav_select.*phasor_jumps;
%         wav_select=wav_select.*reshape(phasor_jumps,1,1,[]);
    end


    % RVP, takes care of a higher order effects related to the chirp (see Egido IIIC)
    % (This might give some issues if the range if much larger than the
    % reference range R0_p)
    Ns_zp=DDAcf.os_ZP*DDAcf.Np;
    R_i=ones(length(Rtrk),1)*((1:Ns_zp)-DDAcf.RefBin*DDAcf.os_ZP)*r_res_zp+R0_i(i)*ones(length(t_pulse),Ns_zp);
    R_p=((R0_p*ones(1,Ns_zp)).^2+(R_i.^2-R0_i(i)^2)).^0.5;

    tau_i=2*(R_p-Rtrk*ones(1,Ns_zp))/DDAcf.c;
    
    if ~is_s6
        RVP=exp(2*pi*1i*DDAcf.s/2*(tau_i).^2);
        wav_select=wav_select.*RVP;
    end
    
    % along-track summing and RRP (see Egido IIID)
    AF=exp(2*pi*1i*DDAcf.fc*(tau_i));
    wav_select=wav_select.*AF;
    %SW=sum(wav_select);
    %CS1b.DAT.sca(i)=max(abs(SW).^2);
    %CS1b.DAT.wav(i,:)=abs(SW).^2/CS1b.DAT.sca(i);
    
%     %{
    figure;
    subplot(1,2,1)
    %[~,I_wav]=max(abs(wav_mod(floor(length(wav_mod(:,1))/2),:)));
    %plot(t_pulse,unwrap(angle(wav_mod(:,I_wav))-mean(angle(wav_mod(:,I_wav)))),'.')
    %dlmwrite('/home/mkleinherenbrink/Algorithms/GMT/FFSAR1/Transponder_model_true.ascii',[t_pulse 11/pi/2*unwrap(angle(wav_mod(:,I_wav))-mean(angle(wav_mod(:,I_wav))))])
    %axis([min(t_select)*1.1 max(t_select)*1.1 -2*pi 2*pi])
    %axis([min(t_select)*1.1 max(t_select)*1.1 -2 2])
    %title('model')
    %xlabel('time [s]')
    %ylabel('unwrapped phase [rad]')
    %subplot(1,2,2)
    [~,I_wav]=max(abs(wav_select(floor(length(wav_select(:,1))/2),:)));
    uwphase = unwrap(angle(wav_select(:,I_wav))-mean(unwrap(angle(wav_select(:,I_wav)))));
    
    P = polyfit(t_pulse,uwphase,2);
    uwphase_fit = P(1)*t_pulse.^2 + P(2)*t_pulse + P(3);
    
    plot(t_pulse,uwphase,'.'); hold on
    plot(t_pulse,uwphase_fit);
    axis([min(t_select)*1.1 max(t_select)*1.1 -0.5*pi 0.5*pi])
    title(mission);
    xlabel('time [s]');
    ylabel('unwrapped phase [rad]');
    
    subplot(1,2,2)
    plot(t_pulse,uwphase - (P(2)*t_pulse + P(3)),'.'); hold on
    title('detrended');
    xlabel('time [s]');
    ylabel('unwrapped phase [rad]');
    axis([min(t_select)*1.1 max(t_select)*1.1 -0.5*pi 0.5*pi])
%     %}
    
    % process 50 consecutive waveforms
    f_b=1/T;
    j=(-24:25)';
    f_f=j*f_b;
    phasor=exp(-2*pi*1i*f_f*(t_pulse(:)'));
    SW=abs(phasor*wav_select(:,:)).^2;
    CS1b.SAR.data(:,i+j)=SW';
    % %     CS1b.DAT.sca(i+j)=max(SW,[],2);
% %     CS1b.DAT.wav(i+j,:)=SW./CS1b.DAT.sca(i+j);
%     CS1b.SAR.echo_scaling(i+j)=max(SW,[],2);
%     CS1b.SAR.data(:,i+j)=(SW./CS1b.SAR.echo_scaling(i+j)')';
% profile viewer,profile off
end
toc

%% Copy all ancillary and auxiliary data fields from CS1a to CS1b
if is_cs
    CS1b         = CopyAncAuxDataFields(CS1a,CS1b);
end

%% Wrap angle in degrees to [-180 180]
CS1a.GEO.LON = wrapTo180(CS1a.GEO.LON);
CS1b.GEO.LON = wrapTo180(CS1b.GEO.LON);

%% Normalize to range [0,65534]
FAC                             = repmat(range(CS1b.SAR.data),size(CS1b.SAR.data,1),1,1)./(65534*(1 - (repmat(min(CS1b.SAR.data),size(CS1b.SAR.data,1),1,1)./CS1b.SAR.data)));
[DUM,CS1b.SAR.echo_scale_power] = log2(squeeze(min(FAC)));
CS1b.SAR.echo_scaling           = DUM*1E9;
CS1b.SAR.data                   = CS1b.SAR.data.*(1./min(FAC));

%% Add field BeamAngle
CS1b.SAR.('BeamAngle') = cell(numel(CS1b.GEO.LAT),1);

end

