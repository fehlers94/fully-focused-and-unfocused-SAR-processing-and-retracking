function DATA = ALL_CryoSat_SAR_L1a_to_L1b(YY,MM,FNamesPreFix,RunID,DOM,UseDEM2compSL)

%CRYOSAT_SAR_L1A_TO_L1B produces L1b data from CryoSat baseline C L1a SAR
%data.

% YYYY = 2010:2016;
% sprintf('qsub-eval -l %i%02i.log "mlab -eval \\"ALL_CryoSat_SAR_L1a_to_L1b(''%i'',''%02i'',[],''TEST'',[])\\""\nsleep 1s\n\n',repmat([reshape(repmat(YYYY,12,1),numel(YYYY)*12,1),reshape(repmat(1:12,numel(YYYY),1)',numel(YYYY)*12,1)],1,2)')

%% Settings
%General settings
LoadCommonSettings

%Computation settings
defval('YY',[])                  %Year for which you want to process data
defval('MM',[])                  %Month for which you want to process data
defval('FNamesPreFix','CS')      %Get DBL filenames
defval('RunID','NorthSea')       %Label that identifies experiment
defval('DOM',[52 55;3 8])        %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('UseDEM2compSL',true)     %Use DEM in computation of surface locations
defval('DEMmodels',{'EuroDEM','SRTM','ASTGTM2','GEBCO'}) %DEM models that will be included in compiling the reference DEM

%Set paths
PathL1a = fullfile(PathDATA,'RadAlt','CryoSat','SIR_SAR_FR');                            %Path to data
PathOUT = fullfile(PathRSLT,sprintf('%s_%s',RunID,DARTutils.GenerateStrDOM(DOM)),'L2');  %Path to output
if ~isdir(PathOUT); mkdir(PathOUT); end

%% Load DEM (used to compute surface locations)
if UseDEM2compSL, DEM = DARTutils.CompileReferenceDEM(DEMmodels,DOM); else DEM = []; end

%% Process level 1a data
%Get filenames level 1a data
if isempty(YY)
    %Get years for which data are available
    YY = dir(PathL1a);
    YY = {YY([YY.isdir]).name}.';
    YY = YY(~strcmp(YY,'.') & ~strcmp(YY,'..'));
    
    FILES = cell(numel(YY),12);
    for i=1:numel(YY)
        %Get filenames level 1a data for each month
        for j=1:12
            MM = sprintf('%02i',j);
            if ~exist(fullfile(PathL1a,YY{i},MM),'dir'), continue, end
            DUM = dir(fullfile(PathL1a,YY{i},MM,sprintf('%s*.DBL',FNamesPreFix)));
            for k = 1:numel(DUM)
                DUM(k).name = strcat(YY{i},filesep,MM,filesep,DUM(k).name);
            end
            FILES{i,j} = DUM;
        end
    end
    FILES = reshape(FILES',numel(YY)*12,1);
    FILES = vertcat(FILES{:});
else
    FILES = dir(fullfile(PathL1a,YY,MM,sprintf('%s*.DBL',FNamesPreFix)));
    for k = 1:numel(FILES)
        FILES(k).name = strcat(YY,filesep,MM,filesep,FILES(k).name);
    end
end

%Return if no FILES are available
if isempty(FILES), DATA = []; return, end

%Apply L1a2L1b processing
for i=1:numel(FILES)
    [~,FName,~] = fileparts(FILES(i).name);
    try
        load(fullfile(PathOUT,sprintf('%s.mat',FName)),'DATA');
    catch exception
        fprintf('%s\n',exception.message)
        CS = CryoSat_SAR_L1a_to_L1b(FILES(i).name(1:end-4),DOM,DEM);
        save(fullfile(PathOUT,sprintf('%s.mat',FName)),'CS');
    end
end

end
