% close all
clear FName CS DATA

%% Settings
LoadCommonSettings

defval('DOM',[59 84;-74 -10])                                                     %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('DEMmodels',{'ArcticDEM','DTU13MSL'})                                      %DEM models that will be included in compiling the reference DEM
% defval('FName','2013/05/CS_LTA__SIR_SIN_1B_20130523T130013_20130523T130223_C001.DBL') %*.DBL/*.nc file that contains level 1b data
defval('FName','2013/05/CS_LTA__SIR_SIN_1B_20130523T130013_20130523T130223_D001.nc') %*.DBL/*.nc file that contains level 1b data

%Load DEM (used to compute surface locations and/or validation)
if ~exist('DEM','var')
    DEM = DARTutils.CompileReferenceDEM(DEMmodels,DOM);
end
% IDXpoi = 1000;
% IDXpoi = 1000:1500;
% IDXpoi = 2800:2880;
IDXpoi = [];
% IDXpoi = 208;

% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'BetaX',{'MaxNrPeaksPerWD',1},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'BetaX',{'MaxNrPeaksPerWD',2},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'BetaX',{'MaxNrPeaksPerWD',8},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'BetaX_expTE',{'MaxNrPeaksPerWD',1},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'BetaX_expTE',{'MaxNrPeaksPerWD',2},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'BetaX_expTE',{'MaxNrPeaksPerWD',8},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'Brown',{'MaxNrPeaksPerWD',1},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'BrownHayne',{'MaxNrPeaksPerWD',1},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'D2P',{'MaxNrPeaksPerWD',1},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'D2P',{'MaxNrPeaksPerWD',8},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'FunctionFit',{'MaxNrPeaksPerWD',1},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'FunctionFit',{'MaxNrPeaksPerWD',8,'SolverAlgo','LM'},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'ICE1',{'MaxNrPeaksPerWD',1},DEM,[],IDXpoi);
% % [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'ICE1',{'MaxNrPeaksPerWD',8},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'OCOG',{'MaxNrPeaksPerWD',1},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'OCOG',{'MaxNrPeaksPerWD',8},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'Threshold',{'MaxNrPeaksPerWD',1},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'Threshold',{'MaxNrPeaksPerWD',8},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'Peak',{'MaxNrPeaksPerWD',1},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'Peak',{'MaxNrPeaksPerWD',8},DEM,[],IDXpoi);
% [DATA,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'SAMOSA2',{'MaxNrPeaksPerWD',1},DEM,[],IDXpoi);
[DATA,CS] = CryoSat_SARIn_L1b_to_L2(FName,DOM,'Swath',[],DEM,[],IDXpoi);

[~,ESAL2]         = Cryo_L2_read_nc(fullfile('/data/Projects/CryoSat2_Greenland/Data/RadAlt/CryoSat/SIR_SIN_L2',regexprep(FName,'SIR_SIN_1B_','SIR_SIN_2__')));
ESAL2.MEA.('DEM') = double(DARTutils.EvalDEMatLATLON(DEM,ESAL2.MEA.LAT_20Hz,ESAL2.MEA.LON_20Hz));
ESAL2.MEA.surf_height_r1_20Hz(ESAL2.MEA.surf_height_r1_20Hz > 5000) = NaN;

% figure,plot(IDXpoi,CS.sigma0(IDXpoi)),hold on
% plot(IDXpoi,ESAL2.MEA.backsc_sig_r1_20Hz(IDXpoi),'g-')

figure ('Position',get (0,'Screensize')); subplot(2,3,1), plot(CS.LON(:),CS.LAT(:),'.'), hold on
plot(ESAL2.MEA.LON_20Hz(:),ESAL2.MEA.LAT_20Hz(:),'g.')

subplot(2,3,2),plot(CS.HEI(IDXpoi,:)),hold on
plot(CS.DEM(IDXpoi,:),'r-')

subplot(2,3,3),plot(IDXpoi,ESAL2.MEA.surf_height_r1_20Hz(IDXpoi),'g-'),hold on
plot(IDXpoi,ESAL2.MEA.DEM(IDXpoi),'r-')

DUM = CS.HEI(DATA.IDXpoi)-CS.DEM(DATA.IDXpoi);
DUM(abs(DUM) >= 15) = NaN;
subplot(2,3,4),Histogram(DUM,1,100,'Differences w.r.t. surf height r1 [m]',[],false);
DUM = ESAL2.MEA.surf_height_r1_20Hz(DATA.IDXpoi)-ESAL2.MEA.DEM(DATA.IDXpoi);
DUM(abs(DUM) >= 15) = NaN;
subplot(2,3,5),Histogram(DUM,1,100,'Differences w.r.t. surf height r1 [m]',[],false);

PlotMap_Grid_Scat([DATA.LAT,DATA.LON,double(DATA.HEI-DATA.DEM)],[],[],[],[],[],DOM)

figure,imagescnan(double(CS.HEI-CS.DEM),'NanColor',[.5 .5 .5]),caxis([-25 25])
