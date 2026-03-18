function DATA = READ_CryoSat_SARIn_L2_ESA(YY,MM,FNamesPreFix)

%READ_CryoSat_SARIn_L2 reads CryoSat L2(I) data provided by ESA.

%% Settings
%General settings
LoadCommonSettings

%Computation settings
defval('YY',[])                  %Year for which you want to process data
defval('MM',[])                  %Month for which you want to process data
defval('FNamesPreFix','CS')      %Get DBL filenames
defval('RunID','GreenlandCoast') %Label that identifies experiment
defval('DEMmodels',{'nsidc0715_MEASURES_gimp_dem_v1','DTU13MSL'}) %DEM models that will be included in compiling the reference DEM
defval('DOM',[59 84;-74 -10])    %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]

%Set paths
PathL2     = fullfile(PathDATA,'RadAlt','CryoSat','SIR_SIN_L2');                       %Path to data

%% Load reference DEM (used in validation)
DEM   = DARTutils.CompileReferenceDEM(DEMmodels,DOM);

%% Load land/sea mask
MSK   = DARTutils.LoadICEmsk;

%% Process level 2 data
%Get filenames level 2 data
if isempty(YY)
    %Get years for which data are available
    YY = dir(PathL2);
    YY = {YY([YY.isdir]).name}.';
    YY = YY(~strcmp(YY,'.') & ~strcmp(YY,'..'));
    
    FILES = cell(numel(YY),12);
    for i=1:numel(YY)
        %Get filenames level 1b data for each month
        for j=1:12
            MM = sprintf('%02i',j);
            if ~exist(fullfile(PathL2,YY{i},MM),'dir'), continue, end
            DUM = dir(fullfile(PathL2,YY{i},MM,sprintf('%s*.DBL',FNamesPreFix)));
            for k = 1:numel(DUM)
                DUM(k).name = strcat(YY{i},filesep,MM,filesep,DUM(k).name);
            end
            FILES{i,j} = DUM;
        end
    end
    FILES = reshape(FILES',numel(YY)*12,1);
    FILES = vertcat(FILES{:});
else
    FILES = dir(fullfile(PathL2,YY,MM,sprintf('%s*.DBL',FNamesPreFix)));
    for k = 1:numel(FILES)
        FILES(k).name = strcat(YY,filesep,MM,filesep,FILES(k).name);
    end
end

%Return if no FILES are available
if isempty(FILES), DATA = []; return, end

%% Read ESA data and copy fields to struct DATA
DUM   = cell(numel(FILES),1);
for i=1:numel(FILES)
    %Read DBL file
    [~,CS]        = Cryo_L2_read(fullfile(PathL2,FILES(i).name));

    %Select SARIn data
    IDX           = CS.GEO.MEAS_MODE(:) == 3;

    %Use quality flags to edit data
    IDXfd = CS.MEA.Qual_flag_20Hz.record_degraded(:) == 1 | CS.MEA.Qual_flag_20Hz.orbit_err(:) == 1 | ...
        CS.MEA.Qual_flag_20Hz.height_err_r1(:) == 1 | CS.MEA.Qual_flag_20Hz.height_err_r2(:) == 1 | ...
        CS.MEA.Qual_flag_20Hz.height_err_r3(:) == 1 | ... %CS.MEA.Qual_flag_20Hz.SARin_Xtrack_angle_err(:) == 1 | ...
        CS.MEA.Qual_flag_20Hz.SARin_ch1_err(:) == 1 | CS.MEA.Qual_flag_20Hz.SARin_ch2_err(:) == 1 | ...
        CS.MEA.Qual_flag_20Hz.misp_err(:) == 1 | CS.MEA.Qual_flag_20Hz.delta_time_err(:) == 1 | ...
        CS.MEA.Qual_flag_20Hz.SIN_Baseline_bad(:) == 1 | CS.MEA.Qual_flag_20Hz.SIN_Out_of_Range(:) == 1 | ...
        CS.MEA.Qual_flag_20Hz.SIN_Velocity_bad(:) == 1 | CS.MEA.Qual_flag_20Hz.cal_wrng(:) == 1;
    IDX(IDXfd) = false;
    
    %Copy clean data to struct DATA
    DATA          = struct;
    TIME          = repmat(datenum('2000','yyyy') + CS.GEO.TAI.days(:)' + CS.GEO.TAI.secs(:)'./86400 + ...
        CS.GEO.TAI.microsecs(:)'./1e6./86400,20,1)+(CS.MEA.delta_time_20Hz./86400);
    DATA.('TIME') = TIME(IDX);
    DATA.('LAT')  = CS.MEA.LAT_20Hz(IDX);
    DATA.('LON')  = CS.MEA.LON_20Hz(IDX);
    DATA.('HEI1') = CS.MEA.surf_height_r1_20Hz(IDX);
    DATA.('HEI2') = CS.MEA.surf_height_r2_20Hz(IDX);
    DATA.('HEI3') = CS.MEA.surf_height_r3_20Hz(IDX);
    DATA.('ST')   = CS.COR.SURF_TYPE(IDX);
    DATA.('SARin_Xtrack_angle_err') = CS.MEA.Qual_flag_20Hz.SARin_Xtrack_angle_err(IDX);

    %Apply ice mask
    if ~isempty(DATA.LAT), DATA   = ApplyIceMask(DATA,MSK); end
    
    %Collect in DUM
    DUM{i}        = DATA;
end

%Concatenate fields
DUM   = vertcat(DUM{:});
names = fieldnames(DUM);
for i=1:numel(names)
    DATA.(names{i}) = vertcat(DUM(:).(names{i}));
end

%Interpolate DEM to measurement locations
DATA.('DEM') = double(DARTutils.EvalDEMatLATLON(DEM,DATA.LAT,DATA.LON));

%Save ALL data
save(fullfile(PathDATA,'RadAlt','CryoSat','SIR_SIN_L2','AllDATA.mat'),'DATA','-v7.3');
