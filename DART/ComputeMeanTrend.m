function ComputeMeanTrend(RunID,DOM,Retracker,TileSize,ppn)

%COMPUTEMEANTREND computes the average trend per tile

%qsub-eval -nn 1 -np 12 -l ComputeMeanTrend.log "mlab -eval2 \"ComputeMeanTrend\""

clc

%% Settings
%General settings
LoadCommonSettings

%Computation settings
defval('RunID','MaxNrPks8_MinPkProm0_1') %Label that identifies experiment
defval('DOM',[59 84;-74 -10])            %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','Threshold')          %Retracker to be used
defval('TileSize',1)                     %Tile size [km]
defval('ppn',12)                         %#ppn's available

%Set paths
% PathOUT = fullfile(PathRSLT,sprintf('%s_%s_%s',RunID,DARTutils.GenerateStrDOM(DOM),Retracker)); %Path to output
PathOUT = fullfile(PathDATA,'RadAlt','CryoSat','SIR_SIN_L2'); %Path to output
if ~isdir(PathOUT); mkdir(PathOUT); end

%% Create parallel pool on cluster 
if numel(gcp('nocreate')) == 0,parpool(ppn),end

try
    load(fullfile(PathOUT,'StatsTILES.mat'));
catch exception
    warning('%s\n',exception.message)
    
    %% Load estimated trends
    DATA = Estimate_dhdt(RunID,DOM,Retracker,[],[],[],1);
    
    %% Divide data set into tiles and sort DATA by TileID
    EDGEx         = floor((min(DATA.x)-TileSize)/TileSize)*TileSize:TileSize:ceil((max(DATA.x)+TileSize)/TileSize)*TileSize;
    EDGEy         = floor((min(DATA.y)-TileSize)/TileSize)*TileSize:TileSize:ceil((max(DATA.y)+TileSize)/TileSize)*TileSize;
    [~,IDXx]      = histc(DATA.x,EDGEx);
    [~,IDXy]      = histc(DATA.y,EDGEy);
    DATA.('IDXp') = sub2ind([numel(EDGEy)-1 numel(EDGEx)-1],IDXy,IDXx);
    [~,IDX]       = sort(DATA.IDXp);
    DATA          = structfun(@(x) (x(IDX,:)), DATA, 'UniformOutput', false);
    
    %% For each tile, select indices of first and last data point in DATA and identify TileIDs of surrounding tiles
    [TILE.('C'),TILE.('R')]       = meshgrid(1:numel(EDGEx)-1,1:numel(EDGEy)-1);
    [TILE.('x'),TILE.('y')]       = meshgrid(EDGEx(1:end-1)+(.5*TileSize),EDGEy(1:end-1)+(.5*TileSize));
    TILE.('IDXp')                 = sub2ind([numel(EDGEy)-1 numel(EDGEx)-1],TILE.R,TILE.C);
    [TILE.('IDXb'),TILE.('IDXe')] = deal(zeros(size(TILE.R)));
    BE                            = [[1;find(diff(DATA.IDXp))+1],[find(diff(DATA.IDXp));size(DATA.IDXp,1)]];
    TILE.IDXb(DATA.IDXp(BE(:,1))) = BE(:,1);
    TILE.IDXe(DATA.IDXp(BE(:,2))) = BE(:,2);
    VALIDp                        = find(TILE.IDXb ~= 0);
    clear('EDGEx','EDGEy','IDXx','IDXy','IDX','BE')
    
    %% Compute statistics of estimated trends per tile
    [AVGdhdt,STDdhdt,NRdhdt,MINdhdt,MAXdhdt] = deal(nan(numel(TILE.IDXp),1));
    ALLdhdt                                  = DATA.dhdt;
    parfor i=1:numel(TILE.IDXp)
        if TILE.IDXb(i) == 0, continue, end
        IDX        = TILE.IDXb(i):TILE.IDXe(i);
        AVGdhdt(i) = nanmean(ALLdhdt(IDX));
        STDdhdt(i) = nanstd(ALLdhdt(IDX));
        NRdhdt(i)  = numel(IDX);
        MINdhdt(i) = nanmin(ALLdhdt(IDX));
        MAXdhdt(i) = nanmax(ALLdhdt(IDX));
    end
    TILE.('AVGdhdt') = AVGdhdt;
    TILE.('STDdhdt') = STDdhdt;
    TILE.('NRdhdt')  = NRdhdt;
    TILE.('MINdhdt') = MINdhdt;
    TILE.('MAXdhdt') = MAXdhdt;
    
    save(fullfile(PathOUT,'StatsTILES.mat'),'TILE')
end

WGS84           = referenceEllipsoid('wgs84','m');
[LONgrd,LATgrd] = meshgrid(linspace(DOM(2,1),DOM(2,2),2000),linspace(DOM(1,1),DOM(1,2),3000));
[Xgrd,Ygrd]     = polarstereo_fwd(LATgrd(:),LONgrd(:),WGS84.SemimajorAxis,WGS84.Eccentricity,70,-45);
AVGdhdt_grd     = interp2(TILE.x.*1000,TILE.y.*1000,reshape(TILE.AVGdhdt,size(TILE.x)),Xgrd(:),Ygrd(:),'linear');
PlotMap_Grid_Scat({LATgrd,LONgrd,reshape(AVGdhdt_grd,size(LATgrd))},1,2,3,-3,3,[59 84;-74 -10],'GMT_polar',[],[],[],[],[],[],[],'m/y');
figdisp(fullfile(PathOUT,sprintf('figure_AVGdhdt.%s',Set.PrintFormat)),[],Set.PrintOpt,Set.PrintFigures,Set.PrintFormat)

end
