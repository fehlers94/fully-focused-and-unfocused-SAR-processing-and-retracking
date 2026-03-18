function [DATA,TILE] = Evaluate_L2_Heights_IceBridge(Pth2AllDATA,StartTILE,EndTILE,DATA)

%EVALUATE_L2_HEIGHTS_IceBridge evaluates L2 ellipsoidal heights by 
%IceBridge ATM L2 Icessn Elevation, Slope, and Roughness, Version 2

% defval('NumNodes',36); %Nr of nodes available on CEASAR
% BE   = round(linspace(1,13542,NumNodes+1)');
% BE   = round(linspace(1,13537,NumNodes+1)');
% BE   = round(linspace(1,13915,NumNodes+1)'); %ESA L2 data
% BE   = [[1;BE(2:end-1)+1],BE(2:end)];
% sprintf('qsub-eval -nn 1 -np 4 -l Evaluate_L2_Heights_IceBridge%i.log "mlab -eval2 \\"Evaluate_L2_Heights_IceBridge([],%i,%i)\\""\nsleep 30s\n',[[1:size(BE,1)]',BE]')
% sprintf('qsub-eval -nn 1 -np 4 -l Evaluate_L2_Heights_IceBridge%i.log "mlab -eval2 \\"Evaluate_L2_Heights_IceBridge(''/home/cslobbe/Projects/CryoSat2_Greenland/Data/RadAlt/CryoSat/SIR_SIN_L2'',%i,%i)\\""\nsleep 30s\n',[[1:size(BE,1)]',BE]')

clc

%% Settings
%General settings
LoadCommonSettings

defval('RunID','MaxNrPks8_MinPkProm0_1') %Label that identifies experiment
defval('DOM',[59 84;-74 -10])            %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','Threshold')          %Retracker to be used
defval('Pth2AllDATA',fullfile(PathRSLT,sprintf('%s_%s_%s',RunID,DARTutils.GenerateStrDOM(DOM),Retracker))); %Path to AllDATA.mat
defval('TileSize',10)                    %Tile size [km]
defval('MaxDist_CS_ATM',0.1)             %Maximum allowed distance between CryoSat data point and ATM data point [km]

%Set paths
PathOUT = fullfile(Pth2AllDATA,'Val_IceBridge'); %Path to output
if ~isdir(PathOUT); mkdir(PathOUT); end

%Remaining settings
WGS84 = referenceEllipsoid('wgs84','m');

%% Read/pre-process RA data set
if ~exist('DATA','var'), load(fullfile(Pth2AllDATA,'AllDATA.mat')); end
if ~isfield(DATA,'HEI') && isfield(DATA,'HEI1'), DATA.('HEI') = DATA.HEI1; end
[DATA.('x'),DATA.('y')]       = polarstereo_fwd(DATA.LAT,DATA.LON,WGS84.SemimajorAxis,WGS84.Eccentricity,70,-45);
DATA.x                        = DATA.x./1E3; %m --> km
DATA.y                        = DATA.y./1E3; %m --> km

%Label each point in DATA_IB with a unique ID
DATA.('IDorg')                = int32((1:numel(DATA.LAT))');

%Divide data set into tiles
EDGEx                         = floor((min(DATA.x)-TileSize)/TileSize)*TileSize:TileSize:ceil((max(DATA.x)+TileSize)/TileSize)*TileSize;
EDGEy                         = floor((min(DATA.y)-TileSize)/TileSize)*TileSize:TileSize:ceil((max(DATA.y)+TileSize)/TileSize)*TileSize;
[DATA,TILE,VALIDp]            = DARTutils.DivideDATAinTILES(DATA,TileSize,EDGEx,EDGEy);

%% Read/pre-process IceBridge data
DATA_IB                       = READ_IceBridge_ATM_L2;
[DATA_IB.('x'),DATA_IB.('y')] = polarstereo_fwd(DATA_IB.LAT,DATA_IB.LON,WGS84.SemimajorAxis,WGS84.Eccentricity,70,-45);
DATA_IB.x                     = DATA_IB.x./1E3; %m --> km
DATA_IB.y                     = DATA_IB.y./1E3; %m --> km

%Remove duplicates
[~,IDX]                       = unique([DATA_IB.x,DATA_IB.y],'rows','stable');
DATA_IB                       = structfun(@(x) (x(IDX,:)), DATA_IB, 'UniformOutput', false);

%Label each point in DATA_IB with a unique ID
DATA_IB.('IDorg')             = int32((1:numel(DATA_IB.LAT))');

%Divide data set into tiles
[DATA_IB,TILE_IB,~]           = DARTutils.DivideDATAinTILES(DATA_IB,TileSize,EDGEx,EDGEy);

%% Interpolate IceBridge data to RA data locations
%Set indices of first and last tile to be processed by node
defval('StartTILE',1)
defval('EndTILE',numel(VALIDp))
[DATA.('HEI_IB'),DATA.('DIcp_IB')] = deal(nan(size(DATA.HEI)));
for i=VALIDp(StartTILE:EndTILE)'
    try
        load(fullfile(PathOUT,sprintf('Tile%i.mat',i)));
    catch exception
        warning('%s\n',exception.message)
    
        %Create vector with indices of all points within tile i or in one
        %of the eight surrounding tiles
        TwD = ~isnan(TILE_IB.SURp(:,i));
        BE  = [TILE_IB.IDXb(TILE_IB.SURp(TwD,i)),TILE_IB.IDXe(TILE_IB.SURp(TwD,i))];
        BE  = BE(BE(:,1) ~=0,:);
        DUM = cell(1,size(BE,1));
        for j=1:size(BE,1), DUM{j} = BE(j,1):BE(j,2); end
        IDX = horzcat(DUM{:});
        if isempty(IDX), continue,end
        
        %Crop DATA_IB for current tile
        DATA_IBp  = structfun(@(x) (x(IDX,:)), DATA_IB, 'UniformOutput', false);
        
        %Crop DATA for current tile
        DATAp     = structfun(@(x) (x(TILE.IDXb(i):TILE.IDXe(i),:)), DATA, 'UniformOutput', false);

        %Apply nearest neighbor interpolation
        dt        = delaunayTriangulation([DATA_IBp.x,DATA_IBp.y]);
        [xi,MinD] = nearestNeighbor(dt,[DATAp.x,DATAp.y]);
        
        %Save interpolated heights
        DUM       = [double(DATAp.IDsrt),DATA_IBp.HEI(xi),MinD];
        save(fullfile(PathOUT,sprintf('Tile%i.mat',i)),'DUM')

    end
    DATA.HEI_IB(DUM(:,1))  = DUM(:,2);
    DATA.DIcp_IB(DUM(:,1)) = DUM(:,3);
end

%% Compute statistics of differences per tile
if StartTILE == 1 && EndTILE == numel(VALIDp)
    [AVGdiff,STDdiff,NRdiff,MINdiff,MAXdiff] = deal(nan(numel(TILE.IDXp),1));
    ALLdiff                                  = double(DATA.HEI)-DATA.HEI_IB;
    ALLdiff(DATA.DIcp_IB > MaxDist_CS_ATM)   = NaN;
    % ALLdiff(DATA.SARin_Xtrack_angle_err == 1) = NaN;
    for i=1:numel(TILE.IDXp)
        if TILE.IDXb(i) == 0, continue, end
        IDX        = TILE.IDXb(i):TILE.IDXe(i);
        AVGdiff(i) = nanmean(ALLdiff(IDX));
        STDdiff(i) = nanstd(ALLdiff(IDX));
        NRdiff(i)  = numel(IDX);
        MINdiff(i) = nanmin(ALLdiff(IDX));
        MAXdiff(i) = nanmax(ALLdiff(IDX));
    end
    TILE.('AVGdiff') = AVGdiff;
    TILE.('STDdiff') = STDdiff;
    TILE.('NRdiff')  = NRdiff;
    TILE.('MINdiff') = MINdiff;
    TILE.('MAXdiff') = MAXdiff;
    save(fullfile(Pth2AllDATA,sprintf('StatsDiffIceBridgePerTile%ikm_MaxDist%im.mat',TileSize,MaxDist_CS_ATM*1000)),'EDGEx','EDGEy','TILE')
%     save(fullfile(Pth2AllDATA,sprintf('StatsDiffIceBridgePerTile%ikm_MaxDist%im_COHge0_7.mat',TileSize,MaxDist_CS_ATM*1000)),'EDGEx','EDGEy','TILE')
%     save(fullfile(Pth2AllDATA,sprintf('StatsDiffIceBridgePerTile%ikm_MaxDist%im_Appl_Xtrack_angle_err.mat',TileSize,MaxDist_CS_ATM*1000)),'EDGEx','EDGEy','TILE')
end

end
