function [DATA,TILE,DATA_IS] = Evaluate_L2_Heights_ICESat2(Pth2AllDATA,StartTILE,EndTILE,DATA)

%EVALUATE_L2_HEIGHTS_ICESat2 evaluates L2 ellipsoidal heights by 
%ICESat2 data

% defval('NumNodes',36); %Nr of nodes available on CEASAR
% BE   = round(linspace(1,13542,NumNodes+1)');
% BE   = round(linspace(1,13537,NumNodes+1)');
% BE   = round(linspace(1,13915,NumNodes+1)'); %ESA L2 data
% BE   = [[1;BE(2:end-1)+1],BE(2:end)];
% sprintf('qsub-eval -nn 1 -np 4 -l Evaluate_L2_Heights_ICESat2%i.log "mlab -eval2 \\"Evaluate_L2_Heights_ICESat2([],%i,%i)\\""\nsleep 30s\n',[[1:size(BE,1)]',BE]')
% sprintf('qsub-eval -nn 1 -np 4 -l Evaluate_L2_Heights_ICESat2%i.log "mlab -eval2 \\"Evaluate_L2_Heights_ICESat2(''/home/cslobbe/Projects/CryoSat2_Greenland/Data/RadAlt/CryoSat/SIR_SIN_L2'',%i,%i)\\""\nsleep 30s\n',[[1:size(BE,1)]',BE]')

clc

%% Settings
%General settings
LoadCommonSettings

defval('RunID','MaxNrPks8_MinPkProm0_1') %Label that identifies experiment
defval('DOM',[59 84;-74 -10])            %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','Threshold')          %Retracker to be used
defval('Pth2AllDATA',fullfile(PathRSLT,sprintf('%s_%s_%s',RunID,DARTutils.GenerateStrDOM(DOM),Retracker))); %Path to AllDATA.mat
defval('TileSize',50)                    %Tile size [km]
defval('MaxDist_CS_ATM',0.1)             %Maximum allowed distance between CryoSat data point and ATM data point [km]

%Set paths
PathOUT = fullfile(Pth2AllDATA,'Val_ICESat2'); %Path to output
if ~isdir(PathOUT); mkdir(PathOUT); end

%Remaining settings
WGS84 = referenceEllipsoid('wgs84','m');

%% Read/pre-process RA data set
if ~exist('DATA','var'), load(fullfile(Pth2AllDATA,'AllDATA.mat')); end
if ~isfield(DATA,'HEI') && isfield(DATA,'HEI1'), DATA.('HEI') = DATA.HEI1; end
[DATA.('x'),DATA.('y')]       = polarstereo_fwd(DATA.LAT,DATA.LON,WGS84.SemimajorAxis,WGS84.Eccentricity,70,-45);
DATA.x                        = DATA.x./1E3; %m --> km
DATA.y                        = DATA.y./1E3; %m --> km

%Label each point in DATA_IS with a unique ID
DATA.('IDorg')                = int32((1:numel(DATA.LAT))');

%Divide data set into tiles
EDGEx                         = floor((min(DATA.x)-2*TileSize)/TileSize)*TileSize:TileSize:ceil((max(DATA.x)+2*TileSize)/TileSize)*TileSize;
EDGEy                         = floor((min(DATA.y)-2*TileSize)/TileSize)*TileSize:TileSize:ceil((max(DATA.y)+2*TileSize)/TileSize)*TileSize;
[DATA,TILE,VALIDp]            = DARTutils.DivideDATAinTILES(DATA,TileSize,EDGEx,EDGEy);

%% Read/pre-process ICESat-2 data
DUM1 = load(fullfile(PathDATA,'LaserAlt','ICESat2','GRE_2018_12_ATLAS_1l.mat'));
DUM2 = load(fullfile(PathDATA,'LaserAlt','ICESat2','GRE_2018_12_ATLAS_1r.mat'));
DUM3 = load(fullfile(PathDATA,'LaserAlt','ICESat2','GRE_2018_12_ATLAS_2l.mat'));
DUM4 = load(fullfile(PathDATA,'LaserAlt','ICESat2','GRE_2018_12_ATLAS_2r.mat'));
DUM5 = load(fullfile(PathDATA,'LaserAlt','ICESat2','GRE_2018_12_ATLAS_3l.mat'));
DUM6 = load(fullfile(PathDATA,'LaserAlt','ICESat2','GRE_2018_12_ATLAS_3r.mat'));
DATA_IS.('TIME') = vertcat(DUM1.ICE.t,DUM2.ICE.t,DUM3.ICE.t,DUM4.ICE.t,DUM5.ICE.t,DUM6.ICE.t);
DATA_IS.('LAT')  = vertcat(DUM1.ICE.lat,DUM2.ICE.lat,DUM3.ICE.lat,DUM4.ICE.lat,DUM5.ICE.lat,DUM6.ICE.lat);
DATA_IS.('LON')  = vertcat(DUM1.ICE.lon,DUM2.ICE.lon,DUM3.ICE.lon,DUM4.ICE.lon,DUM5.ICE.lon,DUM6.ICE.lon);
DATA_IS.('HEI')  = vertcat(DUM1.ICE.h,DUM2.ICE.h,DUM3.ICE.h,DUM4.ICE.h,DUM5.ICE.h,DUM6.ICE.h);
clear DUM1 DUM2 DUM3 DUM4 DUM5 DUM6
[DATA_IS.('x'),DATA_IS.('y')] = polarstereo_fwd(DATA_IS.LAT,DATA_IS.LON,WGS84.SemimajorAxis,WGS84.Eccentricity,70,-45);
DATA_IS.x                     = DATA_IS.x./1E3; %m --> km
DATA_IS.y                     = DATA_IS.y./1E3; %m --> km

%Crop data outside Greenland
IDX                           = (DATA_IS.x > EDGEx(2) & DATA_IS.x < EDGEx(end-1)) & (DATA_IS.y > EDGEy(2) & DATA_IS.y < EDGEy(end-1));
DATA_IS                       = structfun(@(x) (x(IDX,:)), DATA_IS, 'UniformOutput', false);

%Remove duplicates
[~,IDX]                       = unique([DATA_IS.x,DATA_IS.y],'rows','stable');
DATA_IS                       = structfun(@(x) (x(IDX,:)), DATA_IS, 'UniformOutput', false);

%Label each point in DATA_IS with a unique ID
DATA_IS.('IDorg')             = int32((1:numel(DATA_IS.LAT))');

%Divide data set into tiles
[DATA_IS,TILE_IS,~]           = DARTutils.DivideDATAinTILES(DATA_IS,TileSize,EDGEx,EDGEy);

%% Interpolate IceBridge data to RA data locations
%Set indices of first and last tile to be processed by node
defval('StartTILE',1)
defval('EndTILE',numel(VALIDp))
[DATA.('HEI_IS'),DATA.('DIcp_IS')] = deal(nan(size(DATA.HEI)));
DATA.('IntMthd')                   = int8(nan(size(DATA.HEI)));
for i=VALIDp(StartTILE:EndTILE)'
    try
        load(fullfile(PathOUT,sprintf('Tile%i.mat',i)));
    catch exception
        warning('%s\n',exception.message)
    
        %Create vector with indices of all points within tile i or in one
        %of the eight surrounding tiles
        TwD = ~isnan(TILE_IS.SURp(:,i));
        BE  = [TILE_IS.IDXb(TILE_IS.SURp(TwD,i)),TILE_IS.IDXe(TILE_IS.SURp(TwD,i))];
        BE  = BE(BE(:,1) ~=0,:);
        DUM = cell(1,size(BE,1));
        for j=1:size(BE,1), DUM{j} = BE(j,1):BE(j,2); end
        IDX = horzcat(DUM{:});
        if isempty(IDX), continue,end
        
        %Crop DATA_IS for current tile
        DATA_ISp  = structfun(@(x) (x(IDX,:)), DATA_IS, 'UniformOutput', false);
        
        %Crop DATA for current tile
        DATAp     = structfun(@(x) (x(TILE.IDXb(i):TILE.IDXe(i),:)), DATA, 'UniformOutput', false);

        %Apply nearest neighbor interpolation
        dt        = delaunayTriangulation([DATA_ISp.x,DATA_ISp.y]);
        if size(dt,1) == 0, continue,end
        [xi,MinD] = nearestNeighbor(dt,[DATAp.x,DATAp.y]);

        %In case the CryoSat point is surrounded by ICESat2 points in all
        %four quadrants (within a few kilometers) a natural neighbour
        %interpolation is applied. In all other cases we use nearest
        %neighbour interpolation
        [HEIint,NrQD] = deal(nan(numel(DATAp.x),1));
        for j=1:numel(DATAp.x)
            %Identify ICESat2 points within neighborhood of CryoSat point
            IDX_l   = abs(DATA_ISp.x-DATAp.x(j)) <= 4 & abs(DATA_ISp.y-DATAp.y(j)) <= 4;
            %Evaluate in which quadrant each ICESat2 point is located and
            %compute total number of unique quadrants
            NrQD(j) = numel(unique(floor(wrapTo2Pi(atan2(DATA_ISp.y(IDX_l)-DATAp.y(j), DATA_ISp.x(IDX_l)-DATAp.x(j)))/(pi/2))));
            if NrQD(j) == 4
                %Apply natural neighbor interpolation
                F         = scatteredInterpolant(DATA_ISp.x(IDX_l),DATA_ISp.y(IDX_l),double(DATA_ISp.HEI(IDX_l)),'natural');
                HEIint(j) = F(DATAp.x(j),DATAp.y(j));
            else
                %Apply nearest neighbor interpolation
                HEIint(j) = double(DATA_ISp.HEI(xi(j)));
            end
        end
        
        %Save interpolated heights
%         DUM       = [double(DATAp.IDsrt),DATA_ISp.HEI(xi),MinD];
        DUM       = [double(DATAp.IDsrt),HEIint,MinD,NrQD];
        save(fullfile(PathOUT,sprintf('Tile%i.mat',i)),'DUM')

    end
    DATA.HEI_IS(DUM(:,1))  = DUM(:,2);
    DATA.DIcp_IS(DUM(:,1)) = DUM(:,3);
    DATA.IntMthd(DUM(:,1)) = int8(DUM(:,4));
end

%% Compute statistics of differences per tile
if StartTILE == 1 && EndTILE == numel(VALIDp)
    [TILE.('AVG'),TILE.('MAD'),TILE.('MED'),TILE.('MIN'),...
        TILE.('MAX'),TILE.('NR'),TILE.('RMS'),TILE.('SD')]   = deal(nan(size(TILE.C)));
    ALLdiff                                                  = double(DATA.HEI)-DATA.HEI_IS;
    ALLdiff(DATA.DIcp_IS > MaxDist_CS_ATM)                   = NaN;
    
    [TILE.AVG(VALIDp),TILE.MAD(VALIDp),TILE.MED(VALIDp),...
        TILE.MIN(VALIDp),TILE.MAX(VALIDp),TILE.NR(VALIDp),...
        TILE.RMS(VALIDp),TILE.SD(VALIDp)]                    = grpstats(ALLdiff,DATA.IDXp,{'mean',@(x)mad(x,1),'median','min','max','numel','rms','std'});
  
    TILE.('X') = mean([EDGEx(1:end-1);EDGEx(2:end)]);
    TILE.('Y') = mean([EDGEy(1:end-1);EDGEy(2:end)]);
    save(fullfile(Pth2AllDATA,sprintf('StatsDiffIceBridgePerTile%ikm_MaxDist%im.mat',TileSize,MaxDist_CS_ATM*1000)),'TILE')
end

end
