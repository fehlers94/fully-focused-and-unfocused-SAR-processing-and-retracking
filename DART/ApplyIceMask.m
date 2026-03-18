function DATA = ApplyIceMask(DATA,MSK)

%ApplyIceMask returns all data points on ice. 

%The Greenland Ice Mapping Project (GIMP) ice mask is a raster binary land
%classification mask with 1 for glacier ice and 0 for all other terrain or
%water. Ice cover was mapped from a combination of orthorectified Landsat 7
%panchromatic band imagery, distributed by the USGS, and RADARSAT-1
%Synthetic Amplitude Radar (SAR) amplitude images produced by I. Joughin at
%the Applied Physics Laboratory, Univ. of Washington. Landsat imagery was
%from the months of July through September in 1999, 2000 and 2001 (mostly
%2000). RADARSAT images were from fall of 2000.

%Data Citation: Howat I.M., A. Negrete and B.E. Smith, The Greenland Ice
%Mapping Project (GIMP) land classification and surface elevation datasets,
%submitted to The Cryosphere on December 30, 2013

%% Settings
%General settings
LoadCommonSettings

%Remaining settings
WGS84      = referenceEllipsoid('wgs84','m');

%% Load mask
if ~exist('MSK','var')
    MSK = DARTutils.LoadICEmsk;
end

%% Interpolate MSK to data points
[x,y] = polarstereo_fwd(DATA.LAT,DATA.LON,WGS84.SemimajorAxis,WGS84.Eccentricity,70,-45);
IDXx  = MSK.x >= min(x)-5E3 & MSK.x <= max(x)+5E3;
IDXy  = MSK.y >= min(y)-5E3 & MSK.y <= max(y)+5E3;
if any(IDXx) && any(IDXy)
    MSKxy = reshape(interp2(MSK.x(IDXx),MSK.y(IDXy),MSK.z(IDXy,IDXx),x(:),y(:),'nearest',0),size(x));
else
    MSKxy = zeros(numel(DATA.LAT),1);
end

%% Apply mask
DATA = structfun(@(x) (x(MSKxy==1,:)), DATA, 'UniformOutput', false);

end