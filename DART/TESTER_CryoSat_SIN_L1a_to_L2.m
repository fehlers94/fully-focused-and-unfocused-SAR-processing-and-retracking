% close all
clear CS1b DATA Retracker

%% Settings
LoadCommonSettings

% defval('FName','2017/04/CS_OFFL_SIR_SIN_FR_20170404T003732_20170404T004142_C001') %*.DBL file that contains level 1a data
defval('FName','2017/01/CS_OFFL_SIR_SIN_FR_20170101T000522_20170101T000543_C001') %*.DBL file that contains level 1a data
defval('DOM',[59 84;-74 -10])                                                     %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','OCOG')                                                        %Retracker to be used
defval('DEMmodels',{'nsidc0715_MEASURES_gimp_dem_v1','DTU13MSL'})                 %DEM models that will be included in compiling the reference DEM

%Load DEM (used to compute surface locations and/or validation)
if ~exist('DEM','var')
    DEM = DARTutils.CompileReferenceDEM(DEMmodels,DOM);
end

%Load land/sea mask
if ~exist('MSK','var')
    MSK   = DARTutils.LoadICEmsk;
end

%Set paths
PathL1bESA = fullfile(PathDATA,'RadAlt','CryoSat','SIR_SIN_L1'); %Path to data

%% L1a --> L1b --> L2
[CS1b,~]     = CryoSat_SARIn_L1a_to_L1b(FName,DOM);
CS1b.SIN.data(:,isnan(CS1b.GEO.LAT)) = NaN;
CreateRadarEchogram(CS1b)
[DATA,~]     = CryoSat_SARIn_L1b_to_L2(CS1b,DOM,Retracker,{'MaxNrPeaksPerWD',1},DEM);

%% ESA L1b --> L2
[DATA2,ESAL1] = CryoSat_SARIn_L1b_to_L2(sprintf('%s',regexprep(FName,'SIR_SIN_FR_','SIR_SIN_1B_')),DOM,Retracker,{'MaxNrPeaksPerWD',1},DEM);
CreateRadarEchogram(ESAL1)

%% Compare to DEM
DATA.('DEM')  = double(DARTutils.EvalDEMatLATLON(DEM,DATA.LAT,DATA.LON));
DATA2.('DEM') = double(DARTutils.EvalDEMatLATLON(DEM,DATA2.LAT,DATA2.LON));

%% Apply ice mask
DATA         = ApplyIceMask(DATA,MSK);
DATA2        = ApplyIceMask(DATA2,MSK);

%% Plot height profiles
figure,plot(DATA.LAT,DATA.HEI,'.')
hold on
plot(DATA2.LAT,DATA2.HEI,'g.')
plot(DATA.LAT,DATA.DEM,'r.')

%% Plot histograms
Z = double(DATA.HEI)-DATA.DEM;
Histogram(Z(abs(Z) <= 100),1,100)
Z = double(DATA2.HEI)-DATA2.DEM;
Histogram(Z(abs(Z) <= 100),1,100)
