function DATA = ALL_CryoSat_SARIn_L1b_to_L2(YY,MM,FNamesPreFix,RunID,DOM,Retracker,MaxNrPeaksPerWD)

%CRYOSAT_SARIN_L1B_TO_L2 computes elevations from CryoSat baseline B/C
%level 1b SARIn data.

% YYYY = 2010:2016;
% sprintf('qsub-eval -l %i%02i.log "mlab -eval \\"ALL_CryoSat_SARIn_L1b_to_L2(''%i'',''%02i'',[],''MaxNrPks8_MinPkProm0_1_ModTrRetr'',[],''Threshold'')\\""\nsleep 1s\n\n',repmat([reshape(repmat(YYYY,12,1),numel(YYYY)*12,1),reshape(repmat(1:12,numel(YYYY),1)',numel(YYYY)*12,1)],1,2)')

%% Settings
%General settings
LoadCommonSettings

%Computation settings
defval('YY',[])                  %Year for which you want to process data
defval('MM',[])                  %Month for which you want to process data
defval('FNamesPreFix','CS')      %Get DBL filenames
defval('RunID','GreenlandCoast') %Label that identifies experiment
defval('DOM',[59 84;-74 -10])    %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','Threshold')  %Retracker to be used
defval('MaxNrPeaksPerWD',8)      %Max number of peaks to be extracted per waveform
defval('DEMmodels',{'nsidc0715_MEASURES_gimp_dem_v1','DTU13MSL'}) %DEM models that will be included in compiling the reference DEM
% defval('DEMmodels',{'GLAS_ICESat_Greenland_v1','DTU13MSL'}) %DEM models that will be included in compiling the reference DEM

%Copy retracking settings to cell array
SetRetr = {'MaxNrPeaksPerWD',MaxNrPeaksPerWD};

%Set paths
PathL1b = fullfile(PathDATA,'RadAlt','CryoSat','SIR_SIN_L1');                                        %Path to data
% PathL1b = fullfile(PathDATA,'RadAlt','CryoSat','BaselineB');                                         %Path to data
PathOUT = fullfile(PathRSLT,sprintf('%s_%s_%s',RunID,DARTutils.GenerateStrDOM(DOM),Retracker),'L2'); %Path to output
if ~isdir(PathOUT); mkdir(PathOUT); end

%% Load reference DEM (used to geolocate scatterers)
DEM   = DARTutils.CompileReferenceDEM(DEMmodels,DOM);

%% Load land/sea mask
MSK   = DARTutils.LoadICEmsk;

%% Process level 1b data
%Get filenames level 1b data
if isempty(YY)
    %Get years for which data are available
    YY = dir(PathL1b);
    YY = {YY([YY.isdir]).name}.';
    YY = YY(~strcmp(YY,'.') & ~strcmp(YY,'..'));
    
    FILES = cell(numel(YY),12);
    for i=1:numel(YY)
        %Get filenames level 1b data for each month
        for j=1:12
            MM = sprintf('%02i',j);
            if ~exist(fullfile(PathL1b,YY{i},MM),'dir'), continue, end
            DUM = dir(fullfile(PathL1b,YY{i},MM,sprintf('%s*.DBL',FNamesPreFix)));
            for k = 1:numel(DUM)
                DUM(k).name = strcat(YY{i},filesep,MM,filesep,DUM(k).name);
            end
            FILES{i,j} = DUM;
        end
    end
    FILES = reshape(FILES',numel(YY)*12,1);
    FILES = vertcat(FILES{:});
else
    FILES = dir(fullfile(PathL1b,YY,MM,sprintf('%s*.DBL',FNamesPreFix)));
    for k = 1:numel(FILES)
        FILES(k).name = strcat(YY,filesep,MM,filesep,FILES(k).name);
    end
end

%Return if no FILES are available
if isempty(FILES), DATA = []; return, end

%Load data
DUM = cell(numel(FILES),1);
for i=1:numel(FILES)
    [~,FName,~] = fileparts(FILES(i).name);
    try
        load(fullfile(PathOUT,sprintf('%s.mat',FName)),'DATA');
    catch exception
        fprintf('%s\n',exception.message)
        DATA = CryoSat_SARIn_L1b_to_L2(FILES(i).name(1:end-4),DOM,Retracker,SetRetr,DEM);
        save(fullfile(PathOUT,sprintf('%s.mat',FName)),'DATA');
    end

    %Just for backward compatibility, can be removed in future releases
    if ~isfield(DATA,'kx2pi'), DATA.('kx2pi') = []; end
    if ~isfield(DATA,'Coh'), DATA.('Coh') = []; end
    
    %Apply ice mask
    if ~isempty(DATA.LAT), DATA   = ApplyIceMask(DATA,MSK); end
    
    %Collect in DUM
    DUM{i} = DATA;
end

%Concatenate fields
DUM   = vertcat(DUM{:});
names = fieldnames(DUM);
for i=1:numel(names)
    DATA.(names{i}) = vertcat(DUM(:).(names{i}));
end
