function DATA = Estimate_dhdt(RunID,DOM,Retracker,MaxDist2CP,StartTILE,EndTILE,ppn)

%ESTIMATE_DHDT estimates trends in time series of ellipsoidal heights

% defval('NumNodes',46); %Nr of nodes available on CEASAR
% BE   = round(linspace(1,13537,NumNodes+1)');
% BE   = round(linspace(1,13152,NumNodes+1)'); %ESA L2 data WITH data editing
% BE   = [[1;BE(2:end-1)+1],BE(2:end)];
% sprintf('qsub-eval -nn 1 -np 4 -l ComputeElevChange%i.log "mlab -eval2 \\"Estimate_dhdt([],[],[],[],%i,%i)\\""\nsleep 30s\n',[[1:size(BE,1)]',BE]')

clc

%% Settings
%General settings
LoadCommonSettings

%Computation settings
defval('RunID','MaxNrPks8_MinPkProm0_1') %Label that identifies experiment
defval('DOM',[59 84;-74 -10])            %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','Threshold')          %Retracker to be used
defval('TileSize',10)                    %Tile size [km]
defval('MaxDist2CP',1)                   %Maximum distance to computational point [km]
defval('ppn',4)                          %#ppn's available

%Set paths
% PathOUT = fullfile(PathRSLT,sprintf('%s_%s_%s',RunID,DARTutils.GenerateStrDOM(DOM),Retracker),'dhdt'); %Path to output
PathOUT = fullfile(PathDATA,'RadAlt','CryoSat','SIR_SIN_L2','dhdt'); %Path to output
if ~isdir(PathOUT); mkdir(PathOUT); end

%% Create parallel pool on cluster 
if numel(gcp('nocreate')) == 0,parpool(ppn),end
spmd, warning off, end

%% Load DATA
load(fullfile(PathDATA,'RadAlt','CryoSat','SIR_SIN_L2','AllDATA.mat'));
outliers         = DATA.SARin_Xtrack_angle_err == 1;
DATA             = structfun(@(x) (x(~outliers,:)), DATA, 'UniformOutput', false);
if ~isfield(DATA,'HEI') && isfield(DATA,'HEI1'), DATA.('HEI') = DATA.HEI1; end
% load(fullfile(PathRSLT,sprintf('%s_%s_%s',RunID,DARTutils.GenerateStrDOM(DOM),Retracker),'AllDATA.mat'));
% load(fullfile(PathRSLT,sprintf('%s_%s_%s',RunID,DARTutils.GenerateStrDOM(DOM),Retracker),'AllDATA_Jakobshvn.mat'));
% IDX = ingeoquad(DATA.LAT,DATA.LON,[67 72],[-52 -38]);
% DATA = structfun(@(x) (x(IDX,:)), DATA, 'UniformOutput', false);
% save('DATAJH.mat','DATA')
% load DATAJH.mat

%Label each point in DATA with a unique ID
DATA.('IDorg') = int32((1:numel(DATA.LAT))');

%Transform times from units of days to years relative to REFdate
DATA.TIME      = (DATA.TIME-Set.REFdate)./Set.Period;

%% Transform lat/lon data to map coordinates for a polar stereographic system
WGS84                   = referenceEllipsoid('wgs84','m');
[DATA.('x'),DATA.('y')] = polarstereo_fwd(DATA.LAT,DATA.LON,WGS84.SemimajorAxis,WGS84.Eccentricity,70,-45);
DATA.x                  = DATA.x./1E3; %m --> km
DATA.y                  = DATA.y./1E3; %m --> km

%% Divide data set into tiles and sort DATA by TileID
[DATA,TILE,VALIDp] = DARTutils.DivideDATAinTILES(DATA,TileSize);

%Set indices of first and last tile to be processed by node
defval('StartTILE',1)
defval('EndTILE',numel(VALIDp))

%% Compute trends
DATA.('dhdt') = nan(numel(DATA.LAT),1);
for i=VALIDp(StartTILE:EndTILE)'
    try
        load(fullfile(PathOUT,sprintf('Tile%i.mat',i)));
        DATA.dhdt(DUM(:,1)) = DUM(:,2);
    catch exception
        warning('%s\n',exception.message)
    
        %Create vector with indices of all points within tile i or in one
        %of the eight surrounding tiles
        BE  = [TILE.IDXb(TILE.SURp(:,i)),TILE.IDXe(TILE.SURp(:,i))];
        BE  = BE(BE(:,1) ~=0,:);
        DUM = cell(1,size(BE,1));
        for j=1:size(BE,1), DUM{j} = BE(j,1):BE(j,2); end
        IDX = horzcat(DUM{:});
        
        %Crop DATA for current tile
        DATAp = structfun(@(x) (x(IDX,:)), DATA, 'UniformOutput', false);
        
        %For each point in DATAp, compute dh/dt
        [dhdt,MSEfit] = deal(nan(numel(DATAp.LAT),1));
        parfor j=1:numel(DATAp.LAT)
            %Select points in DATAp within 1 km from point j
            IDXp = sqrt((DATAp.x-DATAp.x(j)).^2 + (DATAp.y-DATAp.y(j)).^2) < MaxDist2CP;
            if nnz(IDXp) < 14 || max(DATAp.TIME(IDXp)) - min(DATAp.TIME(IDXp)) < 1 || numel(unique(round(DATAp.TIME))) < 4, continue, end
            
            %Model Wouters et al., 2015
            A          = [DATAp.x(IDXp)-DATAp.x(j) DATAp.y(IDXp)-DATAp.y(j) DATAp.TIME(IDXp)];
            FunSurfMod = @(b,A) (b(1) + b(2)*A(:,1) + b(3)*A(:,2) + b(4)*A(:,1).^2 + b(5)*A(:,2).^2 + b(6)*A(:,1).*A(:,2)...
                + b(7)*A(:,1).^2.*A(:,2) + b(8)*A(:,1).*A(:,2).^2 + b(9)*A(:,1).^2.*A(:,2).^2 + b(10).*A(:,3));
            [xhat,~,~,~,MSEtmp] = nlinfit(A,DATAp.HEI(IDXp)-DATAp.DEM(IDXp),FunSurfMod,[DATAp.HEI(j)-DATAp.DEM(j) 0 0 0 0 0 0 0 0 0]);
            dhdt(j)    = xhat(end);
            MSEfit(j)  = MSEtmp;

            %Linear model, only accounting for tilts in plane that
            %is used to parameterize terrain
            % A               = [ones(nnz(IDXp),1),DATAp.x(IDXp)-DATAp.x(j) DATAp.y(IDXp)-DATAp.y(j) DATAp.TIME(IDXp)];
            % [xhat,~,MSEtmp] = lscov(A,DATAp.HEI(IDXp)-DATAp.DEM(IDXp));
            % dhdt(j)         = xhat(end);
            % MSEfit(j)       = MSEtmp;

            %Linear model similar to previous one, except that you rely on
            %DEM for representation terrain. What you estimate here is
            %change in tilt over time
            % A               = [ones(nnz(IDXp),1),(DATAp.x(IDXp)-DATAp.x(j))./DATAp.TIME(IDXp) (DATAp.y(IDXp)-DATAp.y(j))./DATAp.TIME(IDXp) DATAp.TIME(IDXp)];
            % [xhat,~,MSEtmp] = lscov(A,DATAp.HEI(IDXp)-DATAp.DEM(IDXp));
            % dhdt(j,:)       = xhat(:);
            % MSEfit(j)       = MSEtmp;

            % PlotResults(DATAp,IDXp,j,FunSurfMod,A,xhat,DOM)
            % close all
        end
        
        %Save computed trends
        DUM = [double(DATAp.IDsrt),dhdt,MSEfit];
        save(fullfile(PathOUT,sprintf('Tile%i.mat',i)),'DUM')
    end
end

end

function PlotResults(DATAp,IDXp,j,FunSurfMod,A,xhat,DOM)

figure

subplot(2,3,1)
load coast
plot(long,lat,'k')
hold on
plot(DATAp.LON(IDXp),DATAp.LAT(IDXp),'.')
axis([DOM(2,:) DOM(1,:)]);

if ~isempty(FunSurfMod)
    subplot(2,3,3)
    plot3(DATAp.x(IDXp)-DATAp.x(j),DATAp.y(IDXp)-DATAp.y(j),DATAp.HEI(IDXp),'.')
    hold on
    [dx,dy]      = meshgrid(linspace(min(DATAp.x(IDXp)-DATAp.x(j)),max(DATAp.x(IDXp)-DATAp.x(j)),100),linspace(min(DATAp.y(IDXp)-DATAp.y(j)),max(DATAp.y(IDXp)-DATAp.y(j)),100));
    TMPxhat      = xhat;
    TMPxhat(end) = 0;
    Atmp         = [dx(:) dy(:) zeros(numel(dx),1)];
    SurfMod      = FunSurfMod(TMPxhat,Atmp);
    surf(dx,dy,reshape(SurfMod,size(dx)),'edgecolor','none')
    colorbar('horizontal')
end

subplot(2,1,2)
plot(DATAp.TIME(IDXp),DATAp.HEI(IDXp),'.-')
hold on
plot(DATAp.TIME(IDXp),FunSurfMod(xhat,A),'rx')
TMPtimes = (min(DATAp.TIME(IDXp)):1/10:max(DATAp.TIME(IDXp)))';
plot(TMPtimes,[ones(nnz(TMPtimes),1),TMPtimes]*xhat([1,end])','g-')

end

