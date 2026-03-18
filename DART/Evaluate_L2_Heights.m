function DEMxy = Evaluate_L2_Heights(DATA,DEM)

%Evaluate_L2_Heights returns retracked heights.

%% Settings
%General settings
LoadCommonSettings

%Computation settings
defval('DEMmodels',{'nsidc0715_MEASURES_gimp_dem_v1','DTU13MSL'}) %DEM models that will be included in compiling the reference DEM
defval('DOM',[59 84;-74 -10])                                     %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]

%Remaining settings
WGS84      = referenceEllipsoid('wgs84','m');

%% Load DEM
if ~exist('DEM','var')
    DEM = DARTutils.CompileReferenceDEM(DEMmodels,DOM);
end

%% Interpolate DEM to data points
DEMxy = double(DARTutils.EvalDEMatLATLON(DEM,DATA.LAT,DATA.LON));

%% Evaluate results
try
    Histogram(double(DATA.HEI)-DEMxy,1,100,'Differences (meters)')
catch exception
    warning('%s\n',exception.message)
    Histogram(DATA.HEI1-DEMxy,1,100,'Differences (meters)')
end

end
