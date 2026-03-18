% close all
clear CS DATA FName
clc

%% Settings
LoadCommonSettings

defval('DOM',[52 56;3 8])                                                         %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('DEMmodels',{'EuroDEM','SRTM','ASTGTM2','GEBCO'})                          %DEM models that will be included in compiling the reference DEM
defval('FName','2016/03/CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C001.DBL') %*.DBL file that contains level 1b data
% defval('FName','2016/03/CS_LTA__SIR_SAR_1B_20160302T183427_20160302T183612_D001.nc') %*.nc file that contains level 1b data
% defval('FName','2018/01/CS_OFFL_SIR_SAR_1B_20180101T000445_20180101T001018_C001.DBL') %*.DBL file that contains level 1b data
% defval('FName','2018/03/CS_OFFL_SIR_SAR_1B_20180301T162110_20180301T162555_C001.DBL') %*.DBL file that contains level 1b data
% defval('FName','2018/07/CS_OFFL_SIR_SAR_1B_20180713T120440_20180713T120614_C001.DBL') %*.DBL file that contains level 1b data

%Load DEM (used to compute surface locations and/or validation)
if ~exist('DEM','var')
    DEM = DARTutils.CompileReferenceDEM(DEMmodels,DOM);
end

%To retrack all points in the test files, DOM needs to be enlarged.
DOM = [52 59;3 8];

%Select points in track
if strcmp(FName,'2016/03/CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C001.DBL') || strcmp(FName,'2016/03/CS_LTA__SIR_SAR_1B_20160302T183427_20160302T183612_D001.nc')
%     IDXpoi = 1:2353; 
    % IDXpoi = 800;
    IDXpoi = 800:1000;
elseif any(strcmp(FName,{'2018/01/CS_OFFL_SIR_SAR_1B_20180101T000445_20180101T001018_C001.DBL','2018/03/CS_OFFL_SIR_SAR_1B_20180301T162110_20180301T162555_C001.DBL'}))
    DEM    = [];
    DOM    = [62 88;-17 140];
%     IDXpoi = [17 365 443 928 1082 1083 2625 2947 3966 3972 4079 4090 4109 4692 4693 4831 5516 5809 5872 6389 6537 6539 6541 6862 7126 7325 7365];
    IDXpoi = 928;
elseif strcmp(FName,'2018/07/CS_OFFL_SIR_SAR_1B_20180713T120440_20180713T120614_C001.DBL')
%     IDXpoi = 1:2100;
    IDXpoi = 1:200;
else
    error('Please specify!')
end

%% Apply retracking
% [DATA,CS] = SAR_L1b_to_L2('CS',FName,DOM,'BetaX',[],DEM,[],IDXpoi);
% [DATA,CS] = SAR_L1b_to_L2('CS',FName,DOM,'BetaX_expTE',[],DEM,[],IDXpoi);
% [DATA,CS] = SAR_L1b_to_L2('CS',FName,DOM,'Brown',[],DEM,[],IDXpoi);
% [DATA,CS] = SAR_L1b_to_L2('CS',FName,DOM,'BrownHayne',[],DEM,[],IDXpoi);
% [DATA,CS] = SAR_L1b_to_L2('CS',FName,DOM,'D2P',[],DEM,[],IDXpoi);
% [DATA,CS] = SAR_L1b_to_L2('CS',FName,DOM,'FunctionFit',[],DEM,[],IDXpoi);
% [DATA,CS] = SAR_L1b_to_L2('CS',FName,DOM,'OCOG',[],DEM,[],IDXpoi);
% [DATA,CS] = SAR_L1b_to_L2('CS',FName,DOM,'Threshold',[],DEM,[],IDXpoi);
% [DATA,CS] = SAR_L1b_to_L2('CS',FName,DOM,'SAMOSA2FF',[],DEM,[],IDXpoi);
[DATA,CS] = SAR_L1b_to_L2('CS',FName,DOM,'SAMOSA2',[],DEM,[],IDXpoi);
% load TEST_CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C001_20181221_Est3param_1step_InitCondAllParam_TRR1Emin4_LUTbesseli_t0ns_new_values_ChirpS_BW.mat
% load TEST_CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C001_20181214_Est4param_2steps_InitCondAllParam_MovAvg_TRR1Emin4.mat
% load TEST_CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C001_20181214_Est3param_1step_InitCondAllParam_TRR1Emin4.mat
% % load TEST_CS_DATA_20181122_InitCondAllParam_MovAvg_2step_TRR1E0_1Emin5.mat
% % load TEST_CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C001_20181210_Est3param_InitCondAllParam_MovAvg_1step_TRR1Emin4.mat
% % load TEST_CS_OFFL_SIR_SAR_1B_20180713T120440_20180713T120614_C001_20181210_Est3param_InitCondAllParam_MovAvg_1step_TRR1Emin4.mat
% [DATA,CS] = SAR_L1b_to_L2('CS',FName,DOM,'SAMOSA3FF',[],DEM,[],IDXpoi);
% [DATA,CS] = SAR_L1b_to_L2('CS',FName,DOM,'SAMOSA3',[],DEM,[],IDXpoi);

%% Compare estimated paramters to the ones provided by Dinardo Salvatore
if strcmp(FName,'2016/03/CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C001.DBL') || strcmp(FName,'2016/03/CS_LTA__SIR_SAR_1B_20160302T183427_20160302T183612_D001.nc')
    load RES_CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C002.txt
    figure
    subplot(4,1,1)
    plot(IDXpoi,RES_CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C002(IDXpoi,1)),hold on
    dt     = 1/(2*320E6);
    plot(IDXpoi,(CS.n_ret(IDXpoi)-(129))*dt*1E9,'r.-')
    title('Epoch [nanosec]')
    axis([min(IDXpoi) max(IDXpoi) floor(min((CS.n_ret(IDXpoi)-(129))*dt*1E9)) ceil(max((CS.n_ret(IDXpoi)-(129))*dt*1E9))])
    subplot(4,1,2)
    plot(IDXpoi,RES_CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C002(IDXpoi,3)),hold on
    plot(IDXpoi,CS.Pu(IDXpoi),'r.-')
    title('Pu')
    axis([min(IDXpoi) max(IDXpoi) 0 4])
    subplot(4,1,3)
    plot(IDXpoi,RES_CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C002(IDXpoi,2)),hold on
    plot(IDXpoi,CS.SWH(IDXpoi),'r.-')
    title('SWH [m]')
    axis([min(IDXpoi) max(IDXpoi) -1 3])
    subplot(4,1,4)
    plot(IDXpoi,RES_CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C002(IDXpoi,4)),hold on
    plot(IDXpoi,CS.sigma0(IDXpoi),'r.-')
    title('sigma0 [dB]')
    axis([min(IDXpoi) max(IDXpoi) min(CS.sigma0(IDXpoi)) max(CS.sigma0(IDXpoi))])
%     export_fig('figure_Comparison_1step_sub.png','-png','-r300'),close all

    figure
    subplot(2,2,1)
    Histogram(((CS.n_ret(IDXpoi)-(129))*dt*1E9)-RES_CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C002(IDXpoi,1),1,100,'Differences [ns]','Epoch',false,false)
    subplot(2,2,2)
    Histogram(CS.Pu(IDXpoi)-RES_CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C002(IDXpoi,3),1,100,'Differences','Pu',false,false)
    subplot(2,2,3)
    Histogram(CS.SWH(IDXpoi)-RES_CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C002(IDXpoi,2),1,100,'Differences [m]','SWH',false,false)
    subplot(2,2,4)
    Histogram(CS.sigma0(IDXpoi)-RES_CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C002(IDXpoi,4),1,100,'Differences [dB]','sigma0',false,false)
%     export_fig('histogram_Comparison_1step_sub.png','-png','-r300'),close all
    
    c      = 299792458; %Light velocity [m/s]
    NrBins = 256;       %Nr of bins/samples in any waveform
    R1     = (0.5*c*(CS.n_ret(IDXpoi)*dt - ((NrBins/2)+1)*dt));
    R2     = (0.5*c*(RES_CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C002(IDXpoi,1)*1E-9));
    Histogram(R1-R2,1,100)
%     export_fig('histogram_Comparison_1step_sub_1WayRange.png','-png','-r300')
%     TMP    = round((CS.n_ret(IDXpoi)*dt*1E9 - ((NrBins/2)+1)*dt*1E9),2);
%     R1b    = (0.5*c*(TMP*1E-9));
%     Histogram(R1b-R2,1,100)
end

% LON       = ncread('/data/Projects/CryoSat2_Greenland/Data/SWH/SWAN_DCSM_20160302_1800_2300.nc','x');
% LAT       = ncread('/data/Projects/CryoSat2_Greenland/Data/SWH/SWAN_DCSM_20160302_1800_2300.nc','y');
% SWH       = ncread('/data/Projects/CryoSat2_Greenland/Data/SWH/SWAN_DCSM_20160302_1800_2300.nc','wave_height_hm0');
% [LON,LAT] = meshgrid(LON,LAT);
% SWHint    = interp2(LON,LAT,squeeze(mean(SWH(:,:,1:2),3))',CS.GEO.LON(IDXpoi),CS.GEO.LAT(IDXpoi),'linear');
% plot(IDXpoi,SWHint(IDXpoi),'.-','Color',[.5 .5 .5])
% % SWHint    = interp2(LON,LAT,squeeze(SWH(:,:,2))',CS.GEO.LON(IDXpoi),CS.GEO.LAT(IDXpoi),'linear');
% % plot(IDXpoi,SWHint(IDXpoi),'--','Color',[.5 .5 .5])
% 
% LON    = ncread('/data/Projects/CryoSat2_Greenland/Data/SWH/SWAN_ZUNO_20160302_1800_2300.nc','x');
% LAT    = ncread('/data/Projects/CryoSat2_Greenland/Data/SWH/SWAN_ZUNO_20160302_1800_2300.nc','y');
% SWH    = ncread('/data/Projects/CryoSat2_Greenland/Data/SWH/SWAN_ZUNO_20160302_1800_2300.nc','wave_height_hm0');
% SWHt   = squeeze(mean(SWH(:,:,1:2),3));
% % SWHt   = squeeze(SWH(:,:,1));
% IDX    = ~isnan(SWHt);
% SWHint = griddata(LON(IDX(:)),LAT(IDX(:)),SWHt(IDX(:)),CS.GEO.LON(IDXpoi),CS.GEO.LAT(IDXpoi),'cubic');
% plot(IDXpoi,SWHint(IDXpoi),'g.-')
% % SWHt   = squeeze(SWH(:,:,1));
% % IDX    = ~isnan(SWHt);
% % SWHint = griddata(LON(IDX(:)),LAT(IDX(:)),SWHt(IDX(:)),CS.GEO.LON(IDXpoi),CS.GEO.LAT(IDXpoi),'cubic');
% % plot(IDXpoi,SWHint(IDXpoi),'g--')

%% Compare to ESA level-2 data
%Read level-2 file
[~,~,Fext] = fileparts(FName);
if isequal(Fext,'.DBL')
    [~,ESAL2] = Cryo_L2_read(fullfile(PathDATA,'RadAlt/CryoSat/SIR_SAR_L2',sprintf('%s',regexprep(FName,'SIR_SAR_1B_','SIR_SAR_2__'))));
else
    [~,ESAL2] = Cryo_L2_read_nc(fullfile(PathDATA,'RadAlt/CryoSat/SIR_SAR_L2',sprintf('%s',regexprep(FName,'SIR_SAR_1B_','SIR_SAR_2__'))));
end
IDXl2     = ESAL2.MEA.Qual_flag_20Hz.height_err_r1(:) == 0 & ESAL2.MEA.LAT_20Hz(:) ~= 0;

%The following corrections have been applied: dry_tropo + wet_tropo +
%inv_bar + iono_gim + ocean_tide + lpe_tide + ocean_loading + se_tide +
%geo_polar_tide + mode_window_offset_applied + SAR_retrack_applied +
%SAR_ocean_bias_applied + SAR_ice_bias_applied
% ListFields = fieldnames(ESAL2.MEA.Corr_Applic_flag_20Hz);
% for i=1:numel(ListFields)
%     DUM = ESAL2.MEA.Corr_Applic_flag_20Hz.(ListFields{i});
%     if any(DUM(:) == 1)
%         fprintf('%s\n',ListFields{i});
%     end
% end

%Restore applied inverted barometer correction and apply DAC. The latter
%should be obtained from the CS file
CorrST01                      = repmat(-ESAL2.COR.inv_baro'+CS.COR.dac,20,1);
ESAL2.MEA.surf_height_r1_20Hz = ESAL2.MEA.surf_height_r1_20Hz-CorrST01;

%So far, only the contributions due to tides and wind/mean sea level
%pressure variations have been removed. The main contributions that are
%left are the baroclinic and steric water level variations. Here, we
%compare the resulting sea surface heights to an altimeter-derived mean sea
%surface. This way, we account at least for the time-averaged baroclinic
%and steric contributions
DATA.('h_MSL')  = DeriveData_from_ModelsDTU('MSL','DTU18',[DATA.LAT,DATA.LON,zeros(size(DATA.LAT))],1,2,3,'TEST',[],[],datenum([2003 01 01]),false,false,false,false);
ESAL2.('h_MSL') = DeriveData_from_ModelsDTU('MSL','DTU18',[ESAL2.MEA.LAT_20Hz(:),ESAL2.MEA.LON_20Hz(:),zeros(size(ESAL2.MEA.LAT_20Hz(:)))],1,2,3,'TEST_ESAL2',[],[],datenum([2003 01 01]),false,false,false,false);

%Plot profiles
figure,plot(DATA.LAT,double(DATA.HEI)),hold on
plot(ESAL2.MEA.LAT_20Hz(IDXl2),ESAL2.MEA.surf_height_r1_20Hz(IDXl2),'g-')
plot(DATA.LAT,DATA.h_MSL,'-','Color',[.5 .5 .5])

%Plot histogram of differences
IDXtmp = DATA.SurfT == 0;
Histogram(double(DATA.HEI(IDXtmp))-DATA.h_MSL(IDXtmp),1,100)

IDXtmp = IDXl2 & ESAL2.COR.SURF_TYPE(:) == 0;
Histogram(ESAL2.MEA.surf_height_r1_20Hz(IDXtmp)-ESAL2.h_MSL(IDXtmp),1,100)

%Compare backscatter values
figure,plot(IDXpoi,CS.sigma0(IDXpoi)),hold on
plot(IDXpoi,ESAL2.MEA.backsc_sig_r1_20Hz(IDXpoi),'g-')

%% Compare to Geophysical Ocean Product (GOP)
if strcmp(FName,'2018/07/CS_OFFL_SIR_SAR_1B_20180713T120440_20180713T120614_C001.DBL')
    %ssha_20_ku:comment = "altitude of satellite [alt_20_ku] - ku band
    %corrected ocean altimeter range [range_ocean_20_ku] - gim ionospheric
    %correction [iono_cor_gim_01] - model dry tropospheric correction
    %[mod_dry_tropo_cor_01] - wet tropospheric correction
    %([mod_wet_tropo_cor_01] for NOP/IOP products, [gpd_wet_tropo_cor_01]
    %for GOP off line products) - sea state bias correction in Ku band
    %[sea_state_bias_01_ku] - solid earth tide height [solid_earth_tide_01]
    %- geocentric ocean tide height solution 1 [ocean_tide_sol1_01] -
    %geocentric pole tide height [pole_tide_01] - inverted barometer height
    %correction [inv_bar_cor_01] - high frequency fluctuations of the sea
    %surface topography ([hf_fluct_cor_01] for IOP/GOP off line products
    %only) - mean sea surface solution 1 [mean_sea_surf_sol1_01]";
    LAT_GOP  = ncread(fullfile(PathDATA,'RadAlt/CryoSat/SIR_GOP_L2',sprintf('%s.nc',regexprep(FName,'SIR_SAR_1B_','SIR_GOPR_2_'))),'lat_20_ku');
    LON_GOP  = ncread(fullfile(PathDATA,'RadAlt/CryoSat/SIR_GOP_L2',sprintf('%s.nc',regexprep(FName,'SIR_SAR_1B_','SIR_GOPR_2_'))),'lon_20_ku');
    hSat_GOP = ncread(fullfile(PathDATA,'RadAlt/CryoSat/SIR_GOP_L2',sprintf('%s.nc',regexprep(FName,'SIR_SAR_1B_','SIR_GOPR_2_'))),'alt_20_ku');
    R_GOP    = ncread(fullfile(PathDATA,'RadAlt/CryoSat/SIR_GOP_L2',sprintf('%s.nc',regexprep(FName,'SIR_SAR_1B_','SIR_GOPR_2_'))),'range_ocean_20_ku');
    SIG0_GOP = ncread(fullfile(PathDATA,'RadAlt/CryoSat/SIR_GOP_L2',sprintf('%s.nc',regexprep(FName,'SIR_SAR_1B_','SIR_GOPR_2_'))),'sig0_ocean_20_ku');
    SWH_GOP  = ncread(fullfile(PathDATA,'RadAlt/CryoSat/SIR_GOP_L2',sprintf('%s.nc',regexprep(FName,'SIR_SAR_1B_','SIR_GOPR_2_'))),'swh_ocean_20_ku');
    LAT1_GOP = ncread(fullfile(PathDATA,'RadAlt/CryoSat/SIR_GOP_L2',sprintf('%s.nc',regexprep(FName,'SIR_SAR_1B_','SIR_GOPR_2_'))),'lat_01');
    LON1_GOP = ncread(fullfile(PathDATA,'RadAlt/CryoSat/SIR_GOP_L2',sprintf('%s.nc',regexprep(FName,'SIR_SAR_1B_','SIR_GOPR_2_'))),'lon_01');
    
    %Corrections to be applied in case (i) surface == open oceans or
    %semi‐enclosed seas OR (ii) surface == enclosed seas or lakes (SSB corr. is
    %still lacking). Results in SSH as defined in Eq. 1 (Dinardo et al. 2018)
    CorrST01  = CS.COR.ocean_equilibrium_tide+CS.COR.ocean_longperiod_tide+CS.COR.ocean_loading_tide+...
        CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide+CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.gim_ion+...
        CS.COR.dac;
    %Corrections to be applied in case the instantaneous SSHs (SSHi, see Eq. 4
    %Dinardo et al. 2018) is needed. Note that the SSB correction still has to
    %be applied.
    CorrST01i = CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.ocean_loading_tide+CS.COR.gim_ion+...
        CS.COR.solidearth_tide+0.468*CS.COR.geocentric_polar_tide;
    %Corrections to be applied in case (i) surface == continental ice OR (ii)
    %surface == land
    CorrST23  = CS.COR.ocean_loading_tide+CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide+...
        CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.gim_ion;
    
    %Apply corrections by nearest neighbour interpolation
    ST                     = reshape(repmat(CS.COR.surf_type,20,1),numel(CS.COR.surf_type)*20,1);
    CorrST01               = reshape(repmat(CorrST01,20,1),numel(CorrST01)*20,1);
    CorrST23               = reshape(repmat(CorrST23,20,1),numel(CorrST23)*20,1);
    ST                     = ST(CS.GEO.H ~= 0); CorrST01 = CorrST01(CS.GEO.H ~= 0); CorrST23 = CorrST23(CS.GEO.H ~= 0);
    R_GOP(ST <= 1)         = R_GOP(ST <= 1) + CorrST01(ST <= 1);
    R_GOP(ST >= 2)         = R_GOP(ST >= 2) + CorrST23(ST >= 2);
    
    %Plot profiles
    figure,plot(DATA.LAT,double(DATA.HEI)),hold on
    plot(LAT_GOP,hSat_GOP-R_GOP,'g-')
    plot(DATA.LAT,DATA.h_MSL,'-','Color',[.5 .5 .5])
    
    DUM_HEI = CS.HEI(CS.GEO.H ~= 0);
    IDXtmp  = ST == 0;
    DUM     = DUM_HEI(IDXtmp)-(hSat_GOP(IDXtmp)-R_GOP(IDXtmp));
    Histogram(DUM(abs(DUM) <= 0.15),1,100)

    figure,plot(DATA.LAT(IDXtmp),DUM),axis([53 58 -.1 .1]), hold on
    plot(DATA.LAT(IDXtmp),movmean(DUM,20),'g','LineWidth',2)

    figure,plot(CS.GEO.LAT(IDXpoi),CS.sigma0(IDXpoi)),hold on
    plot(LAT_GOP,SIG0_GOP,'g-')
    
    figure,plot(CS.GEO.LAT(IDXpoi),CS.SWH(IDXpoi)),hold on
    plot(LAT_GOP,SWH_GOP,'g-')
end