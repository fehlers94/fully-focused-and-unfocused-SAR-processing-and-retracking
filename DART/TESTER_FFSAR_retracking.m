clc;clear all;
close all;

mission = 'S3A';

% test FFSAR_Processor
LoadCommonSettings
DDAcf = DDA_ConfigFile(mission,'SAR');
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
work_dir = fullfile([getenv('HOME') '/TUDTUM/']);
data_dir = [work_dir '/ffsar_data/'];

% preconditions
if strcmp(mission,'S3A')
    FName_l1a = [data_dir 's3a/l1a-l1b-l2_sets/l1a/S3A_SR_1_SRA_A__20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement_l1a.nc'];
    FName_l1b = [data_dir 's3a/l1a-l1b-l2_sets/l1b/S3A_SR_1_SRA____20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement.nc'];
    FName_l2 = [data_dir 's3a/l1a-l1b-l2_sets/l2/S3A_SR_2_WAT____20190325T193107_20190325T201537_20191218T195020_2670_043_013______MR1_R_NT_004.SEN3/enhanced_measurement.nc'];
    %DOM = [33.8 37.8; 22.4 25.4]; %crete transponder, lat 35.3379, lon 23.7795
%     DOM = [34.35 34.6; 22.4 25.4]; %crete transponder, lat 35.3379, lon 23.7795
    DOM = [34.2 34.7; 22.4 25.4]; %crete transponder, lat 35.3379, lon 23.7795
elseif strcmp(mission,'S3B')
    FName_l1a = [data_dir 's3b/l1a-l1b-l2_sets/l1a/S3B_SR_1_SRA_A__20190620T083403_20190620T092432_20191107T073220_3029_026_335______MR1_R_NT_004.SEN3/measurement_l1a.nc'];
    FName_l1b = [data_dir 's3b/l1a-l1b-l2_sets/l1b/S3B_SR_1_SRA____20190620T083403_20190620T092432_20191107T073220_3029_026_335______MR1_R_NT_004.SEN3/measurement.nc'];
    FName_l2 = [data_dir 's3b/l1a-l1b-l2_sets/l2/S3B_SR_2_WAT____20190620T083403_20190620T092212_20191225T135656_2889_026_335______MR1_R_NT_004.SEN3/enhanced_measurement.nc'];
    %DOM = [33.8 37.8; 22.4 25.4]; %crete transponder, lat 35.3379, lon 23.7795
    DOM = [34.3 34.6; 22.4 25.4]; %crete transponder, lat 35.3379, lon 23.7795
elseif strcmp(mission,'CS')
    FName_l1a = [data_dir 'cs/l1a-l1b-l2_sets/l1a/CS_OFFL_SIR1SAR_FR_20160617T121101_20160617T121140_C001/CS_OFFL_SIR1SAR_FR_20160617T121101_20160617T121140_C001.DBL'];    
    FName_l1b = [data_dir 'cs/l1a-l1b-l2_sets/l1b/CS_OFFL_SIR_SAR_1B_20160617T121101_20160617T121140_C001.DBL'];    
    FName_l2 = [data_dir 'cs/l1a-l1b-l2_sets/l2/CS_OFFL_SIR_SAR_2__20160617T121101_20160617T121140_C001.DBL'];    
    DOM = [36.05 36.25; 22.4 25.4]; %crete transponder, lat 35.3379, lon 23.7795
end

%%

% L1b processing
hash_run = DataHash(struct('sets', FFSAR_processing_settings, 'fname', FName_l1a, 'DOM', DOM, 'DDAcf', DDAcf));
tmp_dir = fullfile(getenv('HOME'), '/temp/');
if ~exist(tmp_dir,'dir')
    mkdir(tmp_dir)
end

l1b_temp_filepath = [tmp_dir hash_run '.mat'];

if ~isfile(l1b_temp_filepath)
    tic
    
    ffproc = FFSAR_Processor(FName_l1a,DOM,FFSAR_processing_settings);
    ffproc.setup_proc();
    ffproc.proc();
    CS1b = ffproc.CS1b;
    save(l1b_temp_filepath, 'CS1b');
    
    toc
else
    fprintf('Cached L1b data loaded.')
    CS1b = load(l1b_temp_filepath).CS1b;
end


%% retracking ESA SAR L1b waveforms
[DATA3_1b,CS3_1b] = SAR_L1b_to_L2(mission,FName_l1b,DOM,'SAMOSA2',{'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseOwnPath','true'});
%% some illustrations of reference ranges
figure;
subplot(2,2,1)
plot(CS1b.GEO.LAT, CS1b.GEO.H - CS1b.MEA.win_delay*CONST.c/2,'o');hold on
%plot(CS1b.GEO.LAT, CS1b.MEA.win_delay*CONST.c/2);hold on
%plot(CS1b.GEO.LAT, detrend(CS1b.GEO.H));hold on
%plot(CS3_1b.GEO.LAT, CS3_1b.GEO.H - CS3_1b.MEA.ref_range,'.r', 'MarkerSize',25);
xlim([34.47,34.495])
ylim([8.3,9])
xlabel('Latitude')
ylabel('Sat. Alt. - tracker distance [m]')
legend('FF-SAR waveforms','UF-SAR waveforms')
subplot(2,2,2)
plot(CS1b.GEO.LAT, CS1b.MEA.win_delay*CONST.c/2,'o');hold on
xlabel('Latitude')
ylabel('tracker distance [m]')
xlim([34.47,34.495])

% activate if you want to align the tracker distances before FF-SAR
% waveform averaging
CS1b = waveform_shift_dummy(CS1b,600,mission);

subplot(2,2,3)
plot(CS1b.GEO.LAT, CS1b.GEO.H - CS1b.MEA.win_delay*CONST.c/2,'o');hold on
%plot(CS1b.GEO.LAT, CS1b.MEA.win_delay*CONST.c/2);hold on
%plot(CS1b.GEO.LAT, detrend(CS1b.GEO.H));hold on
%plot(CS3_1b.GEO.LAT, CS3_1b.GEO.H - CS3_1b.MEA.ref_range,'.r', 'MarkerSize',25);
xlim([34.47,34.495])
ylim([8.3,9])
xlabel('Latitude')
ylabel('Sat. Alt. - tracker distance [m]')
legend('FF-SAR waveforms','UF-SAR waveforms')
subplot(2,2,4)
plot(CS1b.GEO.LAT, CS1b.MEA.win_delay*CONST.c/2,'o');hold on
xlabel('Latitude')
ylabel('tracker distance [m]')
xlim([34.47,34.495])
%% averaging waveforms 
%  ATTENTION! This averaging averages waveforms with different reference
%  distances (win_delay piecewise constant) and so mixes 600/50 = 12
%  different reference ranges, that lie 4cm apart (in this S3 test case), which either makes 
%  or indirectly corrects for 12*4cm = 48 cm potential misaligment within the averaging!
CS1ba = FF_SAR_Average_Waveforms(CS1b,600,mission_from_fname(FName_l1a));
CS1ba.SAR.BeamAngle = zeros(size(CS1ba.GEO.LAT));
CS1ba.MEA.ref_range = CS1ba.MEA.win_delay*CONST.c/2;
CS1ba.SAR.scale_factor_ku = ones(size(CS1ba.GEO.LAT));

% TO-DO: 
% - CS1ba.SAR.BeamAngle is not averaged within FF_SAR_Average... as it is a cell array and can hence be set 0
% - resolve reference distance misalignment
% - CS1a.surf_type apparently does not get copied into CS1b
% - CS1a.MEA.ref_range must be provided
% - CS1a.SAR.scale_factor_ku does not get copied into CS1b
% - CS1a.MEA.ref_range and n_ret have transposed dimensions!? -> needed to be fixed by insertion (:), but is this universal?

%   The lines
%     CS.COR.cog_cor_01                = ncread(fullfile(FNameL2.folder,FNameL2.name),'cog_cor_01');
%     CS.COR.mod_instr_cor_range_01_ku = ncread(fullfile(FNameL2.folder,FNameL2.name),'mod_instr_cor_range_01_ku');
%   give no results/throw an exception when I try to read the L2 file along
%   with the L1b one, did ESA format changed?

%% plotting averaged waveforms
figure;
imagesc(1:256,CS1ba.GEO.LAT,CS1ba.SAR.data(:,1:end)')
ylabel('LAT')
xlabel('range bin')
title('averaged FF-SAR waveforms')

%% retracking of FF-SAR
[DATA3_FF,CS3_FF]                           = SAR_L1b_to_L2(mission,CS1ba,DOM,'SAMOSA2',{'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','true'});

%% read SWH and range (/ height) from L2 file
if strcmp(mission,'S3A')|strcmp(mission,'S3B')
    L2.SWH                    = ncread(FName_l2,'swh_ocean_20_ku');
    L2.LON                    = ncread(FName_l2,'lon_20_ku');
    L2.LAT                    = ncread(FName_l2,'lat_20_ku');
    L2.SWHcor                 = ncread(FName_l2,'net_instr_cor_swh_20_ku');
    L2.HEI                    = ncread(FName_l2,'alt_20_ku') - (ncread(FName_l2,'range_ocean_20_ku') - 0.55590); %reverse, that cog_cor has been added to range already to reproduce approximately 'raw' retracking result
elseif strcmp(mission,'CS')
    [~,DUM] = Cryo_L2_read(fullfile(FName_l2));
    L2.LON = DUM.MEA.LON_20Hz(:);
    L2.LAT = DUM.MEA.LAT_20Hz(:);
    L2.HEI = DUM.MEA.surf_height_r1_20Hz(:); % here, geophysical corrections have been applied already...
end
% undo SWH corrections for comparison
%L2.SWH = L2.SWH - L2.SWHcor;

%mask to DOM
L2mask = (L2.LAT > DOM(1))&(L2.LAT < DOM(3));
L1mask = (CS3_1b.GEO.LAT > DOM(1))&(CS3_1b.GEO.LAT < DOM(3));

%% plotting retracked height for both SAMOSA
figure;
subplot(1,2,1)
plot(DATA3_FF.LAT,DATA3_FF.HEI,'b.');hold on
plot(DATA3_1b.LAT,DATA3_1b.HEI,'g-');hold on
try
    plot(L2.LAT(L2mask),L2.HEI(L2mask),'r-');hold on
end
ylabel('Height [m]')
xlabel('LAT')
title('retracked distance')
legend('FF-SAR L1a->L2','SAR L1b->L2','ESA SAR L2')
grid on

subplot(1,2,2)
lat = L2.LAT(L2mask);
hei = L2.HEI(L2mask);
h_mask2 = (lat > min(DATA3_FF.LAT))&(lat < max(DATA3_FF.LAT)); 
h_mask = (DATA3_1b.LAT > min(DATA3_FF.LAT))&(DATA3_1b.LAT < max(DATA3_FF.LAT)); 

plot(DATA3_FF.LAT,polyval(polyfit(DATA3_FF.LAT,DATA3_FF.HEI,2),DATA3_FF.LAT)  -  polyval(polyfit(lat(h_mask2),hei(h_mask2),2),DATA3_FF.LAT),'bo');hold on

plot(DATA3_FF.LAT,polyval(polyfit(DATA3_1b.LAT(h_mask),DATA3_1b.HEI(h_mask),2),DATA3_FF.LAT)  -  polyval(polyfit(lat(h_mask2),hei(h_mask2),2),DATA3_FF.LAT),'go');hold on

ylim([-0.025,0.025])
ylabel('\Delta Height [m]')
xlabel('LAT')
title('retracked distance 2nd order fits')
legend('FF-SAR - EUM L2','SAR L1b->L2 - EUM L2')
grid on
%% plotting significant wave heights
%plot(DATA2.LAT,DATA2.SWH,'b.'); hold on
figure;
% subplot(1,2,1)
% plot(DATA3_FF.LAT,DATA3_FF.SWH,'bo'); hold on
% plot(DATA3_1b.LAT,DATA3_1b.SWH,'g-'); hold on
% try
%     plot(L2.LAT(L2mask),L2.SWH(L2mask),'r-')
% end
% ylabel('SWH [m]')
% xlabel('LAT')
% legend('FF-SAR L1a->L2','SAR L1b->L2','ESA SAR L2')
% 
% subplot(1,2,2)

n_mean = 30;
errorbar( movmean(DATA3_FF.LAT,n_mean),...
        movmean(DATA3_FF.SWH,n_mean),...
        movstd(DATA3_FF.SWH,n_mean), 'b.'); hold on
    
n_mean = 30;
errorbar( movmean(DATA3_1b.LAT,n_mean),...
        movmean(DATA3_1b.SWH,n_mean),...
        movstd(DATA3_1b.SWH,n_mean), 'g.'); hold on
    
n_mean = 30;
try
errorbar( movmean(L2.LAT(L2mask),n_mean),...
        movmean(L2.SWH(L2mask),n_mean),...
        movstd(L2.SWH(L2mask),n_mean), 'r.'); hold on
end
ylabel('SWH [m]')
xlabel('LAT')
legend('FF-SAR L1a->L2','SAR L1b->L2','ESA SAR L2')


%% H_rate causes bias?

% H_rate=median(CS3_1b.GEO.H_rate(L2mask));
% V = median(CS3_1b.GEO.V.V(L2mask));
% 300/V*H_rate;

%% compare an arbitrary retracked FF-SAR waveform with model
% n_wf = 1;
% figure;
% plot(CS3_FF.SAR.data(:,n_wf)/max(CS3_FF.SAR.data(:,n_wf))); hold on
% plot(CS3_FF.WDrecon(:,n_wf)); hold on


%% get the corrections from the L2 file (from S3_L1b_read)
if contains(mission,'S3')

    CS.COR.TIME1Hz                   = datenum('2000','yyyy') + double(ncread(FName_l2,'UTC_day_01')) + ncread(FName_l2,'UTC_sec_01')/86400;
    CS.COR.TIME20Hz                  = datenum('2000','yyyy') + double(ncread(FName_l2,'UTC_day_20_ku')) + ncread(FName_l2,'UTC_sec_20_ku')/86400;
    CS.COR.LAT1Hz                    = ncread(FName_l2,'lat_01');
    CS.COR.LON1Hz                    = ncread(FName_l2,'lon_01');
    CS.COR.LAT20Hz                   = ncread(FName_l2,'lat_20_ku');
    CS.COR.LON20Hz                   = ncread(FName_l2,'lon_20_ku');
    CS.COR.iono_cor_alt_20_ku        = ncread(FName_l2,'iono_cor_alt_20_ku');
    CS.COR.iono_cor_gim_01_ku        = ncread(FName_l2,'iono_cor_gim_01_ku');
    CS.COR.dry_trop                  = ncread(FName_l2,'mod_dry_tropo_cor_zero_altitude_01');
    CS.COR.wet_trop                  = ncread(FName_l2,'rad_wet_tropo_cor_01_ku');
    CS.COR.SSB                       = ncread(FName_l2,'sea_state_bias_01_ku');
    CS.COR.solidearth_tide           = ncread(FName_l2,'solid_earth_tide_01');
    CS.COR.ocean_equilibrium_tide    = ncread(FName_l2,'ocean_tide_sol1_01');
    CS.COR.ocean_longperiod_tide     = zeros(size(CS.COR.ocean_equilibrium_tide)); %This value has already been added to the two geocentric ocean tide height values recorded in the product (ocean_tide_sol1_01 and ocean_tide_sol2_01).
    CS.COR.ocean_loading_tide        = zeros(size(CS.COR.ocean_equilibrium_tide)); %This value has already been added to the corresponding ocean tide height value recorded in the product (ocean_tide_sol2_01)."
    CS.COR.geocentric_polar_tide     = ncread(FName_l2,'pole_tide_01');
    CS.COR.inv_bar                   = ncread(FName_l2,'inv_bar_cor_01');
    CS.COR.hf_fluct_cor              = ncread(FName_l2,'hf_fluct_cor_01');
    % CS.COR.cog_cor_01                = ncread(FName_l2,'cog_cor_01');
    % CS.COR.mod_instr_cor_range_01_ku = ncread(FName_l2,'mod_instr_cor_range_01_ku');

    %Compute sum of instrumental, propagation and geophysical corrections
    %Corrections to be applied in case (i) surface == open oceans or
    %semi‐enclosed seas OR (ii) surface == enclosed seas or lakes (SSB corr. is
    %still lacking). Results in SSH as defined in Eq. 1 (Dinardo et al. 2018)
    CorrST01           = CS.COR.ocean_equilibrium_tide+CS.COR.ocean_longperiod_tide+CS.COR.ocean_loading_tide+...
        CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide+CS.COR.dry_trop+CS.COR.wet_trop+...
        CS.COR.inv_bar+CS.COR.hf_fluct_cor+CS.COR.iono_cor_gim_01_ku;%+...
        %CS.COR.cog_cor_01+CS.COR.mod_instr_cor_range_01_ku;
    %Corrections to be applied in case the instantaneous SSHs (SSHi, see Eq. 4
    %Dinardo et al. 2018) is needed. Note that the SSB correction still has to
    %be applied.
    CorrST01i          = CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.ocean_loading_tide+...
                    CS.COR.solidearth_tide+0.468*CS.COR.geocentric_polar_tide+CS.COR.iono_cor_gim_01_ku;%+...
                    %CS.COR.cog_cor_01+CS.COR.mod_instr_cor_range_01_ku;

    % DATA3_FF.HEI       = DATA3_FF.HEI + interp1(CS.COR.LAT1Hz(~isnan(CorrST01)),CorrST01(~isnan(CorrST01)),DATA3_FF.LAT,'pchip',NaN);
    % DATA3_1b.HEI       = DATA3_1b.HEI + interp1(CS.COR.LAT1Hz(~isnan(CorrST01)),CorrST01(~isnan(CorrST01)),DATA3_1b.LAT,'pchip',NaN);
end
