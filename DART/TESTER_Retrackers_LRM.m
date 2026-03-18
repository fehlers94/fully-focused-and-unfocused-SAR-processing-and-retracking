% close all
clear CS DATA FName
clc

%% Settings
LoadCommonSettings

defval('DOM',[59 84;-74 -10])                                                     %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('DEMmodels',{'ArcticDEM','DTU13MSL'})                                      %DEM models that will be included in compiling the reference DEM
% defval('DEMmodels',{'nsidc0715_MEASURES_gimp_dem_v1','DTU13MSL'})                 %DEM models that will be included in compiling the reference DEM
% defval('DEMmodels',{'GLAS_ICESat_Greenland_v1','DTU13MSL'})
% defval('FName','2018/12/CS_OFFL_SIR_LRM_1B_20181201T194252_20181201T194605_C001.DBL') %*.DBL/*.nc file that contains level 1b data
defval('FName','2018/12/CS_OFFL_SIR_LRM_1B_20181221T075239_20181221T075629_C001.DBL') %*.DBL/*.nc file that contains level 1b data
defval('SetRetr',{'SCmthd','''Slobbe2'''})

%Load DEM used to apply slope corrections
if ~exist('DEM','var')
    DEM     = DARTutils.CompileReferenceDEM(DEMmodels,DOM);
%     BlckSz  = 9;
%     DEM.z   = BlockMean(double(DEM.z),BlckSz,BlckSz);
%     DEM.x   = DEM.x(ceil(BlckSz/2):BlckSz:end-floor(BlckSz/2));
%     DEM.y   = DEM.y(ceil(BlckSz/2):BlckSz:end-floor(BlckSz/2));
%     DEM.GSP = DEM.GSP*BlckSz;
end

% %Compute aspect and slope
% if ~all(isfield(DEM,{'aA','sA'}))
%     DEM = DARTutils.Slope_AND_Aspect_OnLLgrd(DOM,DEMmodels,DEM);
% end

%Load DEM to be used in validation
if ~exist('DEMval','var')
    DEMval = DARTutils.CompileReferenceDEM({'ArcticDEM','DTU13MSL'},DOM);
end

%Select points in track
if strcmp(FName,'2018/12/CS_OFFL_SIR_LRM_1B_20181201T194252_20181201T194605_C001.DBL')
    IDXpoi = 1:4109;
%     IDXpoi = 1145;
%     IDXpoi = 2000:2150;
else
    IDXpoi = [];
%     IDXpoi = 3657;
end

%% Apply retracking
% [DATA,CS] = LRM_L1b_to_L2('CS',FName,DOM,'BetaX',[],DEM,[],IDXpoi);
% [DATA,CS] = LRM_L1b_to_L2('CS',FName,DOM,'BetaX_expTE',[],DEM,[],IDXpoi);
% [DATA,CS] = LRM_L1b_to_L2('CS',FName,DOM,'Brown',[],DEM,[],IDXpoi);
% [DATA,CS] = LRM_L1b_to_L2('CS',FName,DOM,'BrownHayne',[],DEM,[],IDXpoi);
% [DATA,CS] = LRM_L1b_to_L2('CS',FName,DOM,'D2P',[],DEM,[],IDXpoi);
% [DATA,CS] = LRM_L1b_to_L2('CS',FName,DOM,'FunctionFit',[],DEM,[],IDXpoi);
[DATA,CS] = LRM_L1b_to_L2('CS',FName,DOM,'ICE1',SetRetr,DEM,[],IDXpoi);
% [DATA,CS] = LRM_L1b_to_L2('CS',FName,DOM,'ICE2',[],DEM,[],IDXpoi);
% [DATA,CS] = LRM_L1b_to_L2('CS',FName,DOM,'ICE3',[],DEM,[],IDXpoi);
% [DATA,CS] = LRM_L1b_to_L2('CS',FName,DOM,'LIRT',[],DEM,[],IDXpoi);
%     [DATA,CS] = LRM_L1b_to_L2('CS',FName,DOM,'ICE1',[],[],[],IDXpoi);
% [DATA,CS] = LRM_L1b_to_L2('CS',FName,DOM,'OCOG',[],DEM,[],IDXpoi);
% [DATA,CS] = LRM_L1b_to_L2('CS',FName,DOM,'Threshold',[],DEM,[],IDXpoi);

%Interpolate DEM to corrected LAT/LON
CS.('DEM')      = DARTutils.EvalDEMatLATLON(DEMval,CS.LATc,CS.LONc);
DATA.('DEMval') = DARTutils.EvalDEMatLATLON(DEMval,DATA.LAT,DATA.LON);

if numel(IDXpoi) == 1, return, end

%Set IDXpoi if it is still empty
if isempty(IDXpoi)
    IDXpoi = find(~isnan(CS.HEI));
end

%% Compare to ESA level-2 data
try
    %Read level-2 file
    [~,ESAL2]     = Cryo_L2_read(fullfile(PathDATA,'RadAlt/CryoSat/SIR_LRM_L2',sprintf('%s',regexprep(FName,'SIR_LRM_1B_','SIR_LRM_2__'))));
    IDXl2         = ESAL2.MEA.Qual_flag_20Hz.height_err_r1(:) == 0 & ESAL2.MEA.LAT_20Hz(:) ~= 0;
    ESAL2.('DEM') = DARTutils.EvalDEMatLATLON(DEMval,ESAL2.MEA.LAT_20Hz,ESAL2.MEA.LON_20Hz);
catch
    ESAL2         = [];
end

%Plot tracks
figure ('Position',get (0,'Screensize')); subplot(2,3,1), plot(CS.GEO.LON(IDXpoi),CS.GEO.LAT(IDXpoi),'r.'), hold on
plot(CS.LONc(IDXpoi),CS.LATc(IDXpoi),'g.')
if ~isempty(ESAL2), plot(ESAL2.MEA.LON_20Hz(IDXl2),ESAL2.MEA.LAT_20Hz(IDXl2),'b.'), end

%Plot profiles
subplot(2,3,2),yyaxis left; plot(CS.LATc(IDXpoi),CS.DEM(IDXpoi),'r.-'), hold on
plot(CS.LATc(IDXpoi),CS.HEI(IDXpoi),'g.-')
if ~isempty(ESAL2), plot(ESAL2.MEA.LAT_20Hz(IDXl2),ESAL2.MEA.surf_height_r3_20Hz(IDXl2),'b.-'), end
yyaxis right; plot(CS.LATc(IDXpoi),CS.HEI(IDXpoi)-CS.DEM(IDXpoi),'g-'), hold on
if ~isempty(ESAL2), plot(ESAL2.MEA.LAT_20Hz(IDXl2),ESAL2.MEA.surf_height_r3_20Hz(IDXl2)-ESAL2.DEM(IDXl2),'b-'), end

%Plot histograms
subplot(2,3,3),Histogram(CS.HEI(IDXpoi)-CS.DEM(IDXpoi),1,100,'HEI - DEM [m]',[],false);
if ~isempty(ESAL2)
    IDXpoi = find(~isnan(CS.HEI));
    subplot(2,3,4),Histogram(ESAL2.MEA.surf_height_r1_20Hz(IDXpoi)-ESAL2.DEM(IDXpoi),1,100,'Surf height r1 - DEM [m]',[],false);
    subplot(2,3,5),Histogram(ESAL2.MEA.surf_height_r2_20Hz(IDXpoi)-ESAL2.DEM(IDXpoi),1,100,'Surf height r2 - DEM [m]',[],false);
    subplot(2,3,6),Histogram(ESAL2.MEA.surf_height_r3_20Hz(IDXpoi)-ESAL2.DEM(IDXpoi),1,100,'Surf height r3 - DEM [m]',[],false);
end
% figure ('Position',get (0,'Screensize')); subplot(2,2,1),plot(ESAL2.MEA.LAT_20Hz(:),CS.h_c(:),'.-'),hold on
% plot(ESAL2.MEA.LAT_20Hz(IDXl2),ESAL2.MEA.surf_height_r1_20Hz(IDXl2),'g-')
% plot(ESAL2.MEA.LAT_20Hz(IDXl2),ESAL2.MEA.surf_height_r2_20Hz(IDXl2),'r-')
% plot(ESAL2.MEA.LAT_20Hz(IDXl2),ESAL2.MEA.surf_height_r3_20Hz(IDXl2),'k-')
% 
% subplot(2,2,2),Histogram(CS.h_c-ESAL2.MEA.surf_height_r1_20Hz(:),1,100,'Differences w.r.t. surf height r1 [m]',[],false);
% subplot(2,2,3),Histogram(CS.h_c-ESAL2.MEA.surf_height_r2_20Hz(:),1,100,'Differences w.r.t. surf height r2 [m]',[],false);
% subplot(2,2,4),Histogram(CS.h_c-ESAL2.MEA.surf_height_r3_20Hz(:),1,100,'Differences w.r.t. surf height r3 [m]',[],false);
% 



% %Plot histogram of differences
% IDXtmp = DATA.SurfT == 2;
% Histogram(double(DATA.HEI(IDXtmp))-DATA.h_MSL(IDXtmp),1,100)
% 
% IDXtmp = IDXl2 & ESAL2.COR.SURF_TYPE(:) == 0;
% Histogram(ESAL2.MEA.surf_height_r1_20Hz(IDXtmp)-ESAL2.h_MSL(IDXtmp),1,100)
% 
% %Compare backscatter values
% figure,plot(IDXpoi,CS.sigma0(IDXpoi)),hold on
% plot(IDXpoi,ESAL2.MEA.backsc_sig_r1_20Hz(IDXpoi),'g-')
