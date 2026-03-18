function [CS1b,CS1a] = CryoSat_SARIn_L1a_to_L1b(FName,DOM,DEM)

%CRYOSAT_SARIN_L1A_TO_L1B produces L1b data from CryoSat baseline C L1a SIN
%data.

% close all

%% Settings
%General settings
LoadCommonSettings

%Computation settings
defval('FName','2017/04/CS_OFFL_SIR_SIN_FR_20170404T003732_20170404T004142_C001') %*.DBL file that contains level 1a data
defval('DOM',[59 84;-74 -10])                                                     %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('UseDEM2compSL',false)                                                     %Use DEM in computation of surface locations
defval('DEMmodels',{'nsidc0715_MEASURES_gimp_dem_v1','DTU13MSL'})                 %DEM models that will be included in compiling the reference DEM
defval('DEM',[])                                                                  %Compiled DEM

%Set paths
PathL1a = fullfile(PathDATA,'RadAlt','CryoSat','SIR_SIN_FR'); %Path to L1a data

%Remaining settings
defval('MakePlots',false)                 %Make plots
defval('FntSz',14)                        %Set fontsize to be used in figures

%% Load DEM (used to compute surface locations)
if UseDEM2compSL && ~isequal(class(DEM),'griddedInterpolant'), DEM = DARTutils.CompileReferenceDEM(DEMmodels,DOM); end

%% Read & Crop CryoSat Level 1a SIN data
%Read data
[~,CS1a] = Cryo_L1b_read(fullfile(PathL1a,sprintf('%s.DBL',FName)));

%Load Parameter Settings (Reference ellipsoid, SIRAL Characteristics, etc.)
DDAcf = DDA_ConfigFile('CS','SIN',CS1a.GEO.Baseline_ID);

%If segment crosses -180 meridian, wrap angle in degrees to [0 360]
if max(abs(diff(CS1a.GEO.LON(CS1a.FBR.N_pulses(:) > 0)))) > 350
    CS1a.GEO.LON = wrapTo360(CS1a.GEO.LON);
end

%Correct window delay for DORIS USO drift (pp 20 CryoSat Product Handbook)
if isequal(CS1a.GEO.Baseline_ID,'C')
    CS1a.MEA.win_delay = CS1a.MEA.win_delay.*(CS1a.GEO.USO+1);
end

%Correct window delay for internal path delay
CS1a.MEA.win_delay = CS1a.MEA.win_delay + (2*CS1a.MEA.ins_range_corr_rx_tx/DDAcf.c);

%Interpolate DEM to burst satellite positions or compute h_DEM using window delay
if ~isempty(DEM)
    %CS1a.GEO.('h_DEM') = DEM(CS1a.GEO.LON,CS1a.GEO.LAT);
    CS1a.GEO.('h_DEM') = double(DARTutils.EvalDEMatLATLON(DEM,CS1a.GEO.LAT,CS1a.GEO.LON));
else
    CS1a.GEO.('h_DEM') = CS1a.GEO.H-((CS1a.MEA.win_delay*DDAcf.c)/2);
end

%Transform acquisition time to datenum format
CS1a.('TIME') = datenum('2000','yyyy') + CS1a.GEO.TAI.days(:) + CS1a.GEO.TAI.secs(:)./86400 + CS1a.GEO.TAI.microsecs(:)./1e6./86400;

%Crop data
if isempty(DOM)
    IDX = true(numel(CS1a.GEO.LAT),1);
else
    if ischar(DOM)
        switch DOM
            case {'JakobsHavn','jakobshavn'}
                %Load polygon that outlines the JakobsHavn glacier
                load(fullfile(PathDATA,'Geography','jakobshavn_basin.mat'))
                IDX = inpolygon(CS1a.GEO.LON(:),CS1a.GEO.LAT(:),polyLon,polyLat);
            otherwise
                error('DOMain of interest not reckognized')
        end
    else
        IDX = ingeoquad(CS1a.GEO.LAT(:),CS1a.GEO.LON(:),DOM(1,:),DOM(2,:));
    end
end

%% Data editing
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
CS1a.GEO.('IDXfd') = ~IDX | IDXfd(:);
clear('IDX','IDXfd')

%% Calibration
CS1a         = Calibrate_CryoSat_FBRdata(CS1a,DDAcf);

%% Create empty struct L1b format
[~,CS1b_TMP] = Cryo_L1b_struct('SIR_SIN_L1.DBL',size(CS1a.FBR.data_ch1,4)+20,DDAcf);

%% Step 1 Surface Location
CS1b_TMP           = Compute_Surface_Locations(CS1a,CS1b_TMP,DDAcf,DEM);
%Re-create empty struct L1b format (now with proper array dimensions)
[~,CS1b]           = Cryo_L1b_struct('SIR_SIN_L1.DBL',size(CS1b_TMP.GEO.LAT,2),DDAcf);
CS1b.GEO.LAT       = CS1b_TMP.GEO.LAT;      CS1b.GEO.LON           = CS1b_TMP.GEO.LON;
CS1b.GEO.H         = CS1b_TMP.GEO.H;        CS1b.MEA.win_delay     = CS1b_TMP.MEA.win_delay;
CS1b.GEO.H_rate    = CS1b_TMP.GEO.H_rate;   CS1b.GEO.V             = CS1b_TMP.GEO.V;
CS1b.GEO.TAI       = CS1b_TMP.GEO.TAI;      CS1b.GEO.BURST_CNT     = CS1b_TMP.GEO.BURST_CNT;
CS1b.GEO.('h_DEM') = CS1b_TMP.GEO.h_DEM;    CS1b.GEO.Start_Time    = CS1a.GEO.Start_Time;
CS1b.GEO.BaseLine  = CS1b_TMP.GEO.BaseLine; CS1b.GEO.Baseline_ID   = 'OWN';
CS1b.GEO.Beam      = CS1b_TMP.GEO.Beam;     CS1b.MEA.int_phase_corr= CS1b_TMP.MEA.int_phase_corr;
clear('CS1b_TMP')
%Set flag CS1b.GEO.MCD_FLAG.Block_Degraded;
CS1b.GEO.MCD_FLAG.Block_Degraded = double(isnan(CS1b.GEO.LAT));
%Compute roll, pitch, and yaw angles of the antenna bench (p. 16 of the CryoSat product handbook)
CS1b.GEO.Antenna_Bench_Roll  = rad2deg(CS1b.GEO.BaseLine.X);
CS1b.GEO.Antenna_Bench_Yaw   = rad2deg(-CS1b.GEO.BaseLine.Y); 
CS1b.GEO.Antenna_Bench_Pitch = rad2deg(-CS1b.GEO.Beam.Y);

%% Step 2 Beam Angle
CS1a         = Compute_Beam_Angle(CS1a,CS1b,DDAcf);

%% Step 3.a Azimuth FFT
CS1a         = Apply_Azimuth_Processing(CS1a,CS1b,DDAcf);

%% Step 3.b Generate Stack
CS1b         = Generate_Stack(CS1a,CS1b,DDAcf);

%% Step 4 Geometry corrections
[CS1a,CS1b]  = Apply_Geometric_Corrections(CS1a,CS1b,DDAcf);

%% Step 5 Range compression
CS1b         = Apply_Range_Compression(CS1a,CS1b,DDAcf);

%% Antenna Compensation
if DDAcf.ACP
    CS1b = antenna_weighting(CS1a,CS1b,DDAcf); %correction for pitch not included
end

%% Step 6 Multi-looking
CS1b         = Compute_MultiLooked_Waveforms(CS1b);

% %% Step 7 Sigma0 scaling factor
% CS1b         = Compute_Sigma0_ScalingFactor(CS1a,CS1b,DDAcf);

%% Step 8 Copy all ancillary and auxiliary data fields from CS1a to CS1b
CS1b         = CopyAncAuxDataFields(CS1a,CS1b);

% CS1b = compute_stack_characterization_params(CS1b,DDAcf);

%% Wrap angle in degrees to [-180 180]
CS1a.GEO.LON = wrapTo180(CS1a.GEO.LON);
CS1b.GEO.LON = wrapTo180(CS1b.GEO.LON);

%% Normalize to range [0,65534]
FAC                             = repmat(range(CS1b.SIN.data),size(CS1b.SIN.data,1),1,1)./(65534*(1 - (repmat(min(CS1b.SIN.data),size(CS1b.SIN.data,1),1,1)./CS1b.SIN.data)));
[DUM,CS1b.SIN.echo_scale_power] = log2(squeeze(min(FAC)));
CS1b.SIN.echo_scaling           = DUM*1E9;
CS1b.SIN.data                   = CS1b.SIN.data.*(1/min(FAC));

%% Copy beam angles and geometry mask to CS1b
BA = reshape(CS1a.GEO.BeamAngle,size(CS1a.GEO.LAT,1),size(CS1a.GEO.LAT,2),size(CS1a.GEO.BeamAngle,2));
BA = permute(BA,[3,1,2]);
[CS1b.SIN.('BeamAngle'),CS1b.SIN.('Hmask')] = deal(cell(numel(find(~isnan(CS1b.GEO.LAT))),1));
for i=find(~isnan(CS1b.GEO.LAT))'
    CS1b.SIN.BeamAngle{i} = BA(CS1b.GEO.IDXstack{i});
    CS1b.SIN.Hmask{i}     = CS1a.FBR.Hmask(:,CS1b.GEO.IDXstack{i});
end

end
