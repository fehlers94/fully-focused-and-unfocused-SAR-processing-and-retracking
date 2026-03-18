% close all
clear CS1b DATA FName Retracker
clc

%% Settings
LoadCommonSettings

defval('DOM',[52 56;3 8])                                                         %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('DEMmodels',{'EuroDEM','SRTM','ASTGTM2','GEBCO'})                          %DEM models that will be included in compiling the reference DEM
defval('FName','2016/03/CS_OFFL_SIR1SAR_FR_20160302T183427_20160302T183612_C001') %*.DBL file that contains level 1b data
defval('Retracker','SAMOSA2')                                                     %Retracker to be used

%Load DEM (used to compute surface locations and/or validation)
if ~exist('DEM','var')
    DEM = DARTutils.CompileReferenceDEM(DEMmodels,DOM);
end

%Select points in track
if strcmp(FName,'2016/03/CS_OFFL_SIR1SAR_FR_20160302T183427_20160302T183612_C001')
    IDXpoi = 800:900;
else
    error('Please specify!')
end

%Set paths
PathL1bESA = fullfile(PathDATA,'RadAlt','CryoSat','SIR_SAR_L1'); %Path to data

%% L1a --> L1b --> L2
[CS1b,~]     = CryoSat_SAR_L1a_to_L1b(FName,DOM);
CS1b.SAR.data(:,isnan(CS1b.GEO.LAT)) = NaN;
% figure,imagesc(1:numel(CS1b.GEO.LAT),1:size(CS1b.SAR.data,1),reshape(CS1b.SAR.data,size(CS1b.SAR.data,1),numel(CS1b.GEO.LAT)))
% figure,plot((1:size(CS1b.SAR.data,1))',reshape(CS1b.SAR.data(:,700),size(CS1b.SAR.data,1),1))
% figure,plot((1:size(CS1b.SAR.data,1))',reshape(CS1b.SAR.data(:,700:1120),size(CS1b.SAR.data,1),421))
[DATA,CS1b]  = SAR_L1b_to_L2('CS',CS1b,DOM,Retracker,[],DEM,IDXpoi);

if strcmp(FName,'2016/03/CS_OFFL_SIR1SAR_FR_20160302T183427_20160302T183612_C001')
    CS = CS1b;
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

%% Load ESA L1b and --> L2
[~,ESAL1]    = Cryo_L1b_read(fullfile(PathL1bESA,sprintf('%s.DBL',regexprep(FName,'SIR1SAR_FR_','SIR_SAR_1B_'))));
IDX          = ingeoquad(ESAL1.GEO.LAT(:),ESAL1.GEO.LON(:),DOM(1,:),DOM(2,:));
% figure,imagesc(1:nnz(IDX),1:size(ESAL1.SAR.data,1),ESAL1.SAR.data(:,IDX))
% ESASAR = ESAL1.SAR.data(:,IDX); figure,plot((1:size(ESASAR,1))',reshape(ESASAR(:,700),size(ESASAR,1),1))
% ESASAR = ESAL1.SAR.data(:,IDX); figure,plot((1:size(ESASAR,1))',reshape(ESASAR(:,700:end),size(ESASAR,1),409))
ESASARwatt   = ESAL1.SAR.data .* repmat(reshape(1E-9 * ESAL1.SAR.echo_scaling .* 2.^(ESAL1.SAR.echo_scale_power),1,size(ESAL1.SAR.data,2),size(ESAL1.SAR.data,3)),256,1,1);
ESASARwatt   = ESASARwatt(:,IDX); figure,plot((1:size(ESASARwatt,1))',reshape(ESASARwatt(:,700:end),size(ESASARwatt,1),409)), clear ESASARwatt
[DATA2,~]    = SAR_L1b_to_L2('CS',ESAL1,DOM,Retracker,[],DEM,IDXpoi);

%% Compare to DEM
DATA.('DEM')  = DEM(DATA.LON,DATA.LAT);
DATA2.('DEM') = DEM(DATA2.LON,DATA2.LAT);

%Plot profiles
figure,plot(DATA.LAT,DATA.HEI,'.')
hold on
plot(DATA2.LAT,DATA2.HEI,'g.')
plot(DATA.LAT,DATA.DEM,'r.-')

% Histogram(double(DATA.HEI)-DATA.DEM,1,100)
% Histogram(double(DATA2.HEI)-DATA2.DEM,1,100)
% 
% IDX = DATA.SurfT == 0;
% Histogram(double(DATA.HEI(IDX))-DATA.DEM(IDX),1,100)
% IDX = DATA2.SurfT == 0;
% Histogram(double(DATA2.HEI(IDX))-DATA2.DEM(IDX),1,100)

IDX = DATA.SurfT == 0 & DATA.LAT >= 53.5;
Histogram(double(DATA.HEI(IDX))-DATA.DEM(IDX),1,100)
IDX = DATA2.SurfT == 0 & DATA.LAT >= 53.5;
Histogram(double(DATA2.HEI(IDX))-DATA2.DEM(IDX),1,100)

% ESAL2 = READ_CryoSat_SAR_L2_ESA('2016','03');
% ESAL2.('DEM')  = DEM(ESAL2.LON,ESAL2.LAT);
% Histogram(double(ESAL2.HEI1)-ESAL2.DEM,1,100)
% 
% figure,
% plot(ESAL2.LAT,ESAL2.HEI1,'.'),hold on
% plot(DATA.LAT,DATA.HEI,'g.')

%% Altitude of satellite (H) is the same!
% IDX = ESAL1.GEO.LAT(:) ~= 0;
% Z   = interp1(ESAL1.GEO.LAT(IDX),ESAL1.GEO.H(IDX),CS1b.GEO.LAT(:),'spline');
% plot(CS1b.GEO.LAT(:),CS1b.GEO.H(:)-Z,'.-')

%% Compare window delays
% IDX   = ESAL1.GEO.LAT(:) ~= 0;
% Z     = interp1(ESAL1.GEO.LAT(IDX),ESAL1.MEA.win_delay(IDX),CS1b.GEO.LAT(:),'spline');
% DDAcf = DDA_ConfigFile('CS','SAR',CS1b.GEO.Baseline_ID);
% figure,plot(CS1b.GEO.LAT(:),(CS1b.MEA.win_delay(:)-Z)*DDAcf.c,'.-')
