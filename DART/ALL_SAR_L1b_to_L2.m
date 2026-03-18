function DATA = ALL_SAR_L1b_to_L2(SAT,YY,MM,FNamesPreFix,RunID,DOM,Retracker,SetRetr)

%ALL_SAR_L1b_to_L2 produces L2 data from SAR L1b data.

% YYYY = 2010:2019;
% sprintf('qsub-eval -l CS_%i%02i.log "mlab -eval \\"ALL_SAR_L1b_to_L2(''CS'',''%i'',''%02i'',[],''ROFI'',[51 55;2 10],''SAMOSA2'',{''UseDEM2CompInitRP'',''true'';''UseLSmskInWfC'',''true''})\\""\nsleep 1s\n\n',repmat([reshape(repmat(YYYY,12,1),numel(YYYY)*12,1),reshape(repmat(1:12,numel(YYYY),1)',numel(YYYY)*12,1)],1,2)')
% YYYY = 2016:2019;
% sprintf('qsub-eval -l S3A_%i%02i.log "mlab -eval \\"ALL_SAR_L1b_to_L2(''S3A'',''%i'',''%02i'',''S3A_SR_1_SRA'',''ROFI'',[51 55;2 10],''SAMOSA2'',{''UseDEM2CompInitRP'',''true'';''UseLSmskInWfC'',''true''})\\""\nsleep 1s\n\n',repmat([reshape(repmat(YYYY,12,1),numel(YYYY)*12,1),reshape(repmat(1:12,numel(YYYY),1)',numel(YYYY)*12,1)],1,2)')
% YYYY = 2018:2019;
% sprintf('qsub-eval -l S3B_%i%02i.log "mlab -eval \\"ALL_SAR_L1b_to_L2(''S3B'',''%i'',''%02i'',''S3B_SR_1_SRA'',''ROFI'',[51 55;2 10],''SAMOSA2'',{''UseDEM2CompInitRP'',''true'';''UseLSmskInWfC'',''true''})\\""\nsleep 1s\n\n',repmat([reshape(repmat(YYYY,12,1),numel(YYYY)*12,1),reshape(repmat(1:12,numel(YYYY),1)',numel(YYYY)*12,1)],1,2)')

%% Settings
%General settings
LoadCommonSettings

%Computation settings
defval('SAT','CS')                                                                %Satellite mission from which data are processed ('CS'/'S3A'/'S3B')
defval('YY',[])                                                                   %Year for which you want to process data
defval('MM',[])                                                                   %Month for which you want to process data
defval('FNamesPreFix','CS')                                                       %Set prefix of filenames that contain L1b data
defval('RunID','NorthSea')                                                        %Label that identifies experiment
defval('DOM',[51 55;2 10])                                                        %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','SAMOSA2')                                                     %Retracker to be used
defval('SetRetr',{})
for i=1:size(SetRetr,1), eval(sprintf('%s = %s;',SetRetr{i,1},SetRetr{i,2})); end
defval('UseDEM2CompInitRP',false)                                                 %Use DEM to compute initial retracking point [true/false].
defval('DEMmodels',{'EuroDEM','SRTM','ASTGTM2','GEBCO'})
defval('UseLSmskInWfC',false)                                                     %Use high-resolution land/sea mask in waveform type classification [true/false].

%Set paths
switch SAT
    case {'CS','CryoSat'}
        PathL1b = fullfile(PathDATA,'RadAlt','CryoSat','SIR_SAR_L1');             %Path to data
        FNExt   = 'DBL';
    case {'S3A','Sentinel-3A'}
        PathL1b = fullfile(PathDATA,'RadAlt','Sentinel3A','SR_1_SRA');            %Path to data
        FNExt   = 'SEN3';
    case {'S3B','Sentinel-3B'}
        PathL1b = fullfile(PathDATA,'RadAlt','Sentinel3B','SR_1_SRA');            %Path to data
        FNExt   = 'SEN3';
    otherwise
        error('SAT: %s not recognized',SAT)
end
PathOUT = fullfile(PathRSLT,sprintf('%s_%s_%s_%s',RunID,DARTutils.GenerateStrDOM(DOM),SAT,Retracker),'L1b'); %Path to output
if ~isdir(PathOUT); mkdir(PathOUT); end

%% Load reference DEM (used to identify peak that has to be retracked)
if UseDEM2CompInitRP
    DEM = DARTutils.CompileReferenceDEM(DEMmodels,DOM);
else
    DEM = [];
end

%% Load land/sea mask (1=LAND, 0=OCEAN)
if UseLSmskInWfC
    MSK = DARTutils.LoadLSmsk(DOM);
else
    MSK = [];
end

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
            DUM = dir(fullfile(PathL1b,YY{i},MM,sprintf('%s*.%s',FNamesPreFix,FNExt)));
            for k = 1:numel(DUM)
                DUM(k).name = strcat(YY{i},filesep,MM,filesep,DUM(k).name);
            end
            FILES{i,j} = DUM;
        end
    end
    FILES = reshape(FILES',numel(YY)*12,1);
    FILES = vertcat(FILES{:});
else
    FILES = dir(fullfile(PathL1b,YY,MM,sprintf('%s*.%s',FNamesPreFix,FNExt)));
    for k = 1:numel(FILES)
        FILES(k).name = strcat(YY,filesep,MM,filesep,FILES(k).name);
    end
end

%Return if no FILES are available
if isempty(FILES), DATA = []; return, end

%Apply L1b-L2 processing
DUM = cell(numel(FILES),1);
for i=1:numel(FILES)
    [~,FName,~] = fileparts(FILES(i).name);
    try
        load(fullfile(PathOUT,sprintf('%s.mat',FName)),'DATA');
    catch exception
        fprintf('%s\n',exception.message)
        DATA = SAR_L1b_to_L2(SAT,FILES(i).name,DOM,Retracker,SetRetr,DEM,MSK);
        save(fullfile(PathOUT,sprintf('%s.mat',FName)),'DATA');
    end
    
    if isempty(DATA.LAT), continue, end
    
    %Add file ID
    DATA.('FileID') = int16(ones(size(DATA.LAT))*i);
    if any(strcmp(SAT,{'S3A','Sentinel-3A','S3B','Sentinel-3B'}))
        DATA.('CCC') = int16(ones(size(DATA.LAT))*str2double(FILES(i).name(end-29:end-27)));
        DATA.('LLL') = int16(ones(size(DATA.LAT))*str2double(FILES(i).name(end-25:end-23)));
    end
    
    %Collect in DUM
    DUM{i} = DATA;
end

%Concatenate fields
DUM   = vertcat(DUM{:});
names = fieldnames(DUM);
for i=1:numel(names)
    DATA.(names{i}) = vertcat(DUM(:).(names{i}));
end

end