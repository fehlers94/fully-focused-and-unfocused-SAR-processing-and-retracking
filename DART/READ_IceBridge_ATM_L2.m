function DATA = READ_IceBridge_ATM_L2(YMD)

%READ_ICEBRIDGE_ATM_L2 reads Icebridge ATM L2 data as provided by ESA.

%% Settings
%General settings
LoadCommonSettings

%Computation settings
defval('YMD',[])                                                  %'yyyy.mm.dd' of which you want to read data
defval('DEMmodels',{'nsidc0715_MEASURES_gimp_dem_v1','DTU13MSL'}) %DEM models that will be included in compiling the reference DEM
defval('DOM',[59 84;-74 -10])                                     %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]

%Set paths
PathL2     = fullfile(PathDATA,'IceBridge','ATM','ILATM2'); %Path to data

%% Load/read IceBridge data
try
    load(fullfile(PathL2,'AllDATA.mat'),'DATA');
catch exception
    fprintf('%s\n',exception.message)
    
    %% Load reference DEM (used in validation)
    DEM   = DARTutils.CompileReferenceDEM(DEMmodels,DOM);
    
    %% Load land/sea mask
    MSK   = DARTutils.LoadICEmsk;
    
    %% Process level 2 data
    %Get filenames level 2 data
    if isempty(YMD)
        %Get years for which data are available
        YMD = dir(PathL2);
        YMD = {YMD([YMD.isdir]).name}.';
        YMD = YMD(~strcmp(YMD,'.') & ~strcmp(YMD,'..'));
        
        FILES = cell(numel(YMD),1);
        for i=1:numel(YMD)
            %Get filenames data for each day
            DUM = dir(fullfile(PathL2,YMD{i},'*.csv'));
            for k = 1:numel(DUM)
                DUM(k).name = strcat(YMD{i},filesep,DUM(k).name);
            end
            FILES{i} = DUM;
        end
        FILES = vertcat(FILES{:});
    else
        FILES = dir(fullfile(PathL2,YMD,'*.csv'));
        for k = 1:numel(FILES)
            FILES(k).name = strcat(YMD,filesep,FILES(k).name);
        end
    end
    
    %Return if no FILES are available
    if isempty(FILES), DATA = []; return, end
    
    %% Read ESA data and copy fields to struct DATA
    DUM   = cell(numel(FILES),1);
    for i=1:numel(FILES)
        %Read csv file
        CS            = csvread(fullfile(PathL2,FILES(i).name),10,0);
        DATA          = struct;
        DATA.('TIME') = datenum(FILES(i).name(1:10),'yyyy.mm.dd') + CS(:,1)./86400;
        DATA.('LAT')  = CS(:,2);
        DATA.('LON')  = CS(:,3);
        DATA.('HEI')  = CS(:,4);
        DATA.('RMS')  = CS(:,7);
        DATA.('NrO')  = CS(:,8);
        DATA.('TrID') = CS(:,11);
        
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
    save(fullfile(PathL2,'AllDATA.mat'),'DATA','-v7.3');
end
