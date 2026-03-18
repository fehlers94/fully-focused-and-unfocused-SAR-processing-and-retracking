classdef DARTutils
    methods (Static)
        
        function DEM = CompileReferenceDEM(DEMmodels,DOM)
            %CompileDEM compiles a DEM covering both land and sea from
            %existing DEM models
            
            %% Settings
            %General settings
            LoadCommonSettings
            
            %Computation settings
            defval('DEMmodels',{'nsidc0715_MEASURES_gimp_dem_v1','DTU13MSL'}) %DEM models that will be included in compiling the reference DEM

            %Remaining settings
            WGS84      = referenceEllipsoid('wgs84','m');
            
            try
                load(fullfile(PathRSLT,'DEM',strcat('DEM',sprintf('_%s',DEMmodels{:}),sprintf('_%s.mat',DARTutils.GenerateStrDOM(DOM)))));
            catch exception
                warning('%s\n',exception.message)
            
                %Read/copy ArcticDEM to DEM
                if any(strcmp(DEMmodels,'ArcticDEM'))
                    %Read DEM
                    PathDEM          = fullfile(PathDATA,'DEM','ArcticDEM');
                    GTinfo           = geotiffinfo(fullfile(PathDEM,'arcticdem_mosaic_100m_v3.0.tif'));
                    X                = linspace(GTinfo.BoundingBox(1,1),GTinfo.BoundingBox(2,1),GTinfo.Width);
                    Y                = linspace(GTinfo.BoundingBox(2,2),GTinfo.BoundingBox(1,2),GTinfo.Height);
                    IDX_Xmin         = find(X < -639955,1,'last');
                    IDX_Xmax         = find(X > 855845,1,'first');
                    IDX_Ymin         = find(Y > -3355595,1,'last');
                    IDX_Ymax         = find(Y < -655595,1,'first');
                    DEM.('x')        = X(IDX_Xmin:IDX_Xmax);
                    DEM.('y')        = Y(IDX_Ymax:IDX_Ymin);
                    DEM.('z')        = single(imread(fullfile(PathDEM,'arcticdem_mosaic_100m_v3.0.tif'),'PixelRegion',{[IDX_Ymax IDX_Ymin],[IDX_Xmin IDX_Xmax]}));
                    DEM.('GSP')      = 100;
                    DEM.('CTP')      = 'CT_PolarStereographic';
                    DEM.('phi_c')    = 70;
                    DEM.('lambda_0') = -45;
                    DEM.('GCS')      = 'WGS84';

                    %Fill some gaps. Note, this DEM is not to be used
                    %outside Greenland. Be careful in using this DEM in
                    %Greenland's coastal waters
                    [DEMx,DEMy]      = meshgrid(DEM.x,DEM.y);
                    [LAT,LON]        = polarstereo_inv(DEMx,DEMy,WGS84.SemimajorAxis,WGS84.Eccentricity,DEM.phi_c,DEM.lambda_0);
                    D2C              = DARTutils.GetDist2Coast(DOM,LAT,LON);
                    IDXv             = D2C == 0 & DEM.z ~= -9999;
                    IDXiv            = D2C == 0 & DEM.z == -9999;
                    F                = scatteredInterpolant(DEMx(IDXv),DEMy(IDXv),double(DEM.z(IDXv)),'linear','none');
                    DEM.z(IDXiv)     = single(F(DEMx(IDXiv),DEMy(IDXiv)));
                    
                    %Identify sea-points
                    IDX              = D2C > 0;
                    clear('PathDEM','GTinfo','X','Y','IDX_Xmin','IDX_Xmax','IDX_Ymin','IDX_Ymax')
                end
                
                %Read/copy nsidc0715_MEASURES_gimp_dem_v1 to DEM
                if any(strcmp(DEMmodels,'nsidc0715_MEASURES_gimp_dem_v1'))
                    %Read DEM
                    PathDEM          = fullfile(PathDATA,'DEM','nsidc0715_MEASURES_gimp_dem_v1','90');
                    GTinfo           = geotiffinfo(fullfile(PathDEM,'gimpdem_90m_v1.1.tif'));
                    DEM.('x')        = linspace(GTinfo.BoundingBox(1,1),GTinfo.BoundingBox(2,1),GTinfo.Width);
                    DEM.('y')        = linspace(GTinfo.BoundingBox(2,2),GTinfo.BoundingBox(1,2),GTinfo.Height);
                    DEM.('z')        = single(geotiffread(fullfile(PathDEM,'gimpdem_90m_v1.1.tif')));
                    DEM.('GSP')      = 90;
                    DEM.('CTP')      = 'CT_PolarStereographic';
                    DEM.('phi_c')    = 70;
                    DEM.('lambda_0') = -45;
                    DEM.('GCS')      = 'WGS84';
                    %Identify sea-points
                    MSK.('z')        = geotiffread(fullfile(PathDATA,'Geography','GIMP','90','GimpOceanMask_90m.tif'));
                    IDX              = MSK.z == 1;
                    clear('PathDEM','GTinfo','MSK')
                end
                
                %Read/copy GLAS_ICESat_Greenland_v1 to DEM
                if any(strcmp(DEMmodels,'GLAS_ICESat_Greenland_v1'))
                    DEM              = enviread(fullfile(PathDATA,'DEM','GLAS_ICESat_Greenland_v1','NSIDC_Grn1km_wgs84_elev_cm.dat'));
                    DEM.z            = single(DEM.z)./100;
                    DEM.('GSP')      = 1000;
                    DEM.('CTP')      = 'CT_PolarStereographic';
                    DEM.('phi_c')    = 70;
                    DEM.('lambda_0') = -45;
                    DEM.('GCS')      = 'WGS84';
                    
                    %Fill some gaps. Note, this DEM is not to be used
                    %outside Greenland. Be careful in using this DEM in
                    %Greenland's coastal waters
                    [DEMx,DEMy]      = meshgrid(DEM.x,DEM.y);
                    [LAT,LON]        = polarstereo_inv(DEMx,DEMy,WGS84.SemimajorAxis,WGS84.Eccentricity,DEM.phi_c,DEM.lambda_0);
                    D2C              = DARTutils.GetDist2Coast(DOM,LAT,LON);
                    IDXv             = D2C == 0 & DEM.z ~= 0;
                    IDXiv            = D2C == 0 & DEM.z == 0;
                    F                = scatteredInterpolant(DEMx(IDXv),DEMy(IDXv),double(DEM.z(IDXv)),'linear','none');
                    DEM.z(IDXiv)     = single(F(DEMx(IDXiv),DEMy(IDXiv)));
                    
                    %Identify sea-points
                    IDX              = D2C > 0;
                end
                
                %Read/copy DTU13MSL to DEM
                if any(strcmp(DEMmodels,'DTU13MSL'))
                    [LOg,LAg]   = meshgrid(DOM(2,1)-1:1/60:DOM(2,2)+1,DOM(1,1)-1:1/60:DOM(1,2)+1);
                    MSL         = DeriveData_from_ModelsDTU('MSL','DTU13',[LAg(:),LOg(:),zeros(size(LAg(:)))],1,2,3,DARTutils.GenerateStrDOM(DOM));
                    if strcmp(DEM.CTP,'CT_PolarStereographic')
                        [x,y]     = meshgrid(DEM.x,DEM.y);
                        [lat,lon] = polarstereo_inv(x(IDX),y(IDX),WGS84.SemimajorAxis,WGS84.Eccentricity,DEM.phi_c,DEM.lambda_0);
                    else
                        error('Not implemented')
                    end
                    DEM.z(IDX)  = interp2(LOg,LAg,reshape(MSL,size(LOg)),lon,lat,'cubic');
                    clear('IDX','x','y','lat','lon','LOg','LAg','MSL')
                end
                
                if all(ismember(DEMmodels,{'EuroDEM','SRTM','ASTGTM2','GEBCO'}))
                    PathDEM          = fullfile(PathDATA,'DEM','MAT_EuroDEM_SRTM_ASTGTM2_GEBCO2014'); %Path to DEM
                    IDX              = 0;
                    [LATc,LONc,DEMc] = deal(cell(numel(DOM(1,1):DOM(1,2))+2,numel(DOM(2,1):DOM(2,2))+2));
                    for ew = DOM(2,1)-1:DOM(2,2)+1
                        if ew >= 0, ew_str = 'e'; else ew_str = 'w'; end
                        for ns = DOM(1,1)-1:DOM(1,2)+1
                            if ns >= 0, ns_str = 'n'; else ns_str = 's'; end
                            load(fullfile(PathDEM,sprintf('%s%02i%s%02i.mat',ew_str,ew,ns_str,ns)),'LAT','LON','DEMgrd','HAgrd','SIDgrd');
                            IDX                 = IDX+1;
                            DEMgrd(SIDgrd == 0) = 0;
                            LATc{IDX}           = LAT'; LONc{IDX} = LON; DEMc{IDX} = DEMgrd+HAgrd;
                        end
                    end
                    LATc = cell2mat(LATc); LONc = cell2mat(LONc); DEMc = cell2mat(DEMc);
                    DEM.('x') = LONc(1,:); DEM.('y') = LATc(:,1); DEM.('z') = DEMc;
                    clear('IDX','LATc','LONc','DEMc','LATgrd','LONgrd','DEMgrd','SIDgrd')
                    
                    %Create 2-D interpolant
                    [X,Y] = ndgrid(DEM.x,DEM.y);
                    DEM   = griddedInterpolant(X,Y,DEM.z','linear');
                end
                
                %Save compiled DEM
                save(fullfile(PathRSLT,'DEM',strcat('DEM',sprintf('_%s',DEMmodels{:}),sprintf('_%s.mat',DARTutils.GenerateStrDOM(DOM)))),'DEM');
            end
        end
        
        function STATS = ComputeStats(DATA,IDX)
            
            if ~any(IDX), STATS = []; return, end
            
            [S1,S2,S3,S4,S5,S6,S7,S8] = grpstats(double(DATA.HEI(IDX)-DATA.DEMval(IDX)),DATA.FileID(IDX),{'mean','numel','gname','std','min','max','median',@(x)mad(x,1)});
            STATS                     = [S1,S2,str2double(S3),S4,S5,S6,S7,S8*1.4826];
            for i=1:size(STATS,1)
                STATS(i,9) = DATA.TIME(find(DATA.FileID == STATS(i,3),1,'first'));
            end
            STATS = sortrows(STATS,9);
            
        end
        
        function [DATA,TILE,VALIDp] = DivideDATAinTILES(DATA,TileSize,EDGEx,EDGEy)
            defval('EDGEx',floor((min(DATA.x)-TileSize)/TileSize)*TileSize:TileSize:ceil((max(DATA.x)+TileSize)/TileSize)*TileSize)
            defval('EDGEy',floor((min(DATA.y)-TileSize)/TileSize)*TileSize:TileSize:ceil((max(DATA.y)+TileSize)/TileSize)*TileSize)
            
            %% Divide data set into tiles and sort DATA by TileID
            [~,IDXx]      = histc(DATA.x,EDGEx);
            [~,IDXy]      = histc(DATA.y,EDGEy);
            DATA.('IDXp') = sub2ind([numel(EDGEy)-1 numel(EDGEx)-1],IDXy,IDXx);
            [~,IDX]       = sort(DATA.IDXp);
            DATA          = structfun(@(x) (x(IDX,:)), DATA, 'UniformOutput', false);
            
            %Add new label to each point in DATA
            DATA.('IDsrt') = int32((1:numel(DATA.LAT))');
            
            %% For each tile, select indices of first and last data point in DATA and identify TileIDs of surrounding tiles
            [TILE.('C'),TILE.('R')]       = meshgrid(1:numel(EDGEx)-1,1:numel(EDGEy)-1);
            TILE.('IDXp')                 = sub2ind([numel(EDGEy)-1 numel(EDGEx)-1],TILE.R,TILE.C);
            [TILE.('IDXb'),TILE.('IDXe')] = deal(zeros(size(TILE.R)));
            BE                            = [[1;find(diff(DATA.IDXp))+1],[find(diff(DATA.IDXp));size(DATA.IDXp,1)]];
            TILE.IDXb(DATA.IDXp(BE(:,1))) = BE(:,1);
            TILE.IDXe(DATA.IDXp(BE(:,2))) = BE(:,2);
            TILE.('SURp')                 = nan(9,numel(TILE.IDXp));
            VALIDp                        = find(TILE.IDXb ~= 0);
            for i=VALIDp'
                TMPr           = TILE.R(i)+[-1 -1 -1 0 0 0 1 1 1];
                TMPc           = TILE.C(i)+[-1 0 1 -1 0 1 -1 0 1];
                TILE.SURp(:,i) = sub2ind([numel(EDGEy)-1 numel(EDGEx)-1],TMPr,TMPc);
            end
        end

        function h_DEM = EvalDEMatLATLON(DEM,LAT,LON)
            if strcmp(DEM.CTP,'CT_PolarStereographic')
                RefEll = referenceEllipsoid(DEM.GCS,'m');
                [x,y]  = polarstereo_fwd(LAT,LON,RefEll.SemimajorAxis,RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
                InDEMr = (x > min(DEM.x) & x < max(DEM.x)) & (y > min(DEM.y) & y < max(DEM.y));
                h_DEM  = nan(size(x));
                if any(InDEMr(:))
                    IDXx          = DEM.x >= min(x(InDEMr))-5*DEM.GSP & DEM.x <= max(x(InDEMr))+5*DEM.GSP;
                    IDXy          = DEM.y >= min(y(InDEMr))-5*DEM.GSP & DEM.y <= max(y(InDEMr))+5*DEM.GSP;
                    h_DEM(InDEMr) = interp2(DEM.x(IDXx),DEM.y(IDXy),DEM.z(IDXy,IDXx),x(InDEMr),y(InDEMr),'linear');
                end
            else
                error('Not implemented')
            end
        end
        
        function DEM = Slope_AND_Aspect_OnLLgrd(DOM,DEMmodels,DEM)
            %
            %SLOPE_AND_ASPECT_ONLLGRD computes aspect (:= direction of
            %steepest descent expressed as an azimuth measured clockwise
            %from north) and slope (:= magnitude of the gradient vector
            %converted to an angle) of data grid Z using a finite
            %difference method. Z has to be given on a regular lat/lon grid
            %defined with respect to the WGS84/GRS80 reference ellipsoid.
            %The output is stored in DEM.

            %% Settings
            %General settings
            LoadCommonSettings
            
            %Computation settings
            defval('DEMmodels',{'GLAS_ICESat_Greenland_v1','DTU13MSL'}) %DEM models that will be included in compiling the reference DEM

            %Remaining settings
            WGS84 = referenceEllipsoid('wgs84','m');
            
            %Compute slope and aspect
            try
                load(fullfile(PathRSLT,'DEM',strcat('DEM_SLOPE',sprintf('_%s',DEMmodels{:}),sprintf('_%s.mat',DARTutils.GenerateStrDOM(DOM)))));
            catch exception
                warning('%s\n',exception.message)
                
                if isequal(DEMmodels,{'GLAS_ICESat_Greenland_v1','DTU13MSL'})
                    %Check whether DEM is consistent with stored DEM
                    DUM = load(fullfile(PathRSLT,'DEM',strcat('DEM',sprintf('_%s',DEMmodels{:}),sprintf('_%s.mat',DARTutils.GenerateStrDOM(DOM)))));
                    if ~isequaln(DUM.DEM,DEM)
                        error('DEM not equal to stored DEM!')
                    end
                    
                    %Transform x,y coordinates to lat,lon coordinates
                    [DEMx,DEMy]               = meshgrid(DEM.x,DEM.y);
                    [la,lo]                   = polarstereo_inv(DEMx,DEMy,WGS84.SemimajorAxis,WGS84.Eccentricity,DEM.phi_c,DEM.lambda_0);
                    
                    %Grid DEM on regular lat,lon grid
                    [DEM.('lon'),DEM.('lat')] = meshgrid(DOM(2,1):km2deg(1)*cosd(DEM.phi_c):DOM(2,2),DOM(1,1):km2deg(1):DOM(1,2));
                    DEM.('z_LL')              = griddata(lo(:),la(:),double(DEM.z(:)),DEM.lon,DEM.lat,'linear');
                
                    %Compute aspect and slope
                    [DEM.('aA'),DEM.('sA')]   = gradientm(DEM.lat,DEM.lon,DEM.z_LL,WGS84);

                    %Save DEM including aspect and slope
                    save(fullfile(PathRSLT,'DEM',strcat('DEM_SLOPE',sprintf('_%s',DEMmodels{:}),sprintf('_%s.mat',DARTutils.GenerateStrDOM(DOM)))),'DEM','-v7.3');
                else
                    error('DEMmodels not recognized!')
                end
            end
            
        end

        function [aspectAngle,slopeAngle,dZdx,dZdy] = Slope_AND_Aspect(Z,GSP)
            %
            %SLOPE_AND_ASPECT computes aspect (:= direction of steepest
            %descent expressed as an azimuth measured clockwise from north)
            %and slope (:= magnitude of the gradient vector converted to an
            %angle) of data grid Z using a finite difference method. NOTE
            %THAT IT IS ASSUMED THAT Z IS GIVEN on a equi-distant
            %polar-stereographic grid with sampling interval GSP. The
            %ordering of the points is as follows:
            %
            % o---o---o     N  7---8---9 (row 1)
            % |   |   |     ^  |   |   |
            % o---o---o  =  |  4---5---6 (row 2)
            % |   |   |        |   |   |
            % o---o---o        1---2---3 (row 3)
            %           
            %                    --> E
            %
            %TEST:
            %[aspectAngle,slopeAngle,dZdx,dZdy] = DARTutils.Slope_AND_Aspect(rot90([29 29 29;30 30 30;32 32 32],0),5)
            %rot90([29 29 29;30 30 30;32 32 32],0)
            %[aspectAngle,slopeAngle,dZdx,dZdy] = DARTutils.Slope_AND_Aspect(rot90([29 29 29;30 30 30;32 32 32],2),5)
            %rot90([29 29 29;30 30 30;32 32 32],2)
            %[aspectAngle,slopeAngle,dZdx,dZdy] = DARTutils.Slope_AND_Aspect(rot90([29 29 29;30 30 30;32 32 32],1),5)
            %rot90([29 29 29;30 30 30;32 32 32],1)
            %[aspectAngle,slopeAngle,dZdx,dZdy] = DARTutils.Slope_AND_Aspect(rot90([29 29 29;30 30 30;32 32 32],3),5)
            %rot90([29 29 29;30 30 30;32 32 32],3)

            %For each grid cell, store values of 8 surrounding grid cells.
            %Note that in Matlab, row 1 is north from row 2
            NrPnts      = numel(Z(2:end-1,2:end-1));
            z1_9        = nan(NrPnts,9);
            z1_9(:,1)   = reshape(Z(3:end,1:end-2),NrPnts,1);
            z1_9(:,2)   = reshape(Z(3:end,2:end-1),NrPnts,1);
            z1_9(:,3)   = reshape(Z(3:end,3:end),NrPnts,1);
            z1_9(:,4)   = reshape(Z(2:end-1,1:end-2),NrPnts,1);
            z1_9(:,5)   = reshape(Z(2:end-1,2:end-1),NrPnts,1);
            z1_9(:,6)   = reshape(Z(2:end-1,3:end),NrPnts,1);
            z1_9(:,7)   = reshape(Z(1:end-2,1:end-2),NrPnts,1);
            z1_9(:,8)   = reshape(Z(1:end-2,2:end-1),NrPnts,1);
            z1_9(:,9)   = reshape(Z(1:end-2,3:end),NrPnts,1);

            %Compute gradient, slope and aspect
            dZdx        = reshape((z1_9(:,3) + 2*z1_9(:,6) + z1_9(:,9) - z1_9(:,1) - 2*z1_9(:,4) - z1_9(:,7))/(8*GSP),size(Z(2:end-1,2:end-1),1),size(Z(2:end-1,2:end-1),2));
            dZdy        = reshape((z1_9(:,7) + 2*z1_9(:,8) + z1_9(:,9) - z1_9(:,1) - 2*z1_9(:,2) - z1_9(:,3))/(8*GSP),size(Z(2:end-1,2:end-1),1),size(Z(2:end-1,2:end-1),2));
            slopeAngle  = atand(hypot(dZdy,dZdx));
            aspectAngle = wrapTo360(atan2d(-dZdx,-dZdy));
            
            %The aspect angle is indeterminate in regions of uniform Z, so
            %set it to NaN when both components of the gradient vector
            %vanish
            aspectAngle(dZdy == 0 & dZdx == 0) = NaN;
        end

        function DOMstr = GenerateStrDOM(DOM)
            %GenerateStrDOM generates string from geoquad defined by DOM
            if DOM(1,1) >= 0, TMP(1) = 'N'; else TMP(1) = 'S'; end
            if DOM(2,1) >= 0, TMP(2) = 'E'; else TMP(2) = 'W'; end
            if DOM(1,2) >= 0, TMP(3) = 'N'; else TMP(3) = 'S'; end
            if DOM(2,2) >= 0, TMP(4) = 'E'; else TMP(4) = 'W'; end
            DOMstr = sprintf('%.0f%s%.0f%s%.0f%s%.0f%s',abs(DOM(1,1)),TMP(1),abs(DOM(2,1)),TMP(2),abs(DOM(1,2)),TMP(3),abs(DOM(2,2)),TMP(4));
        end
        
        function D2C = GetDist2Coast(DOM,LAT,LON)
            %General settings
            LoadCommonSettings

            %Load cropped distance to coast grid
            D2C_file     = fullfile(PathDATA,'Geography','dist2coast','dist2coast_1deg_ocean.nc');
            LATd         = ncread(D2C_file,'lat');
            LONd         = ncread(D2C_file,'lon');
            STRTlon      = find(LONd <= DOM(2,1),1,'last'); STRTlat = find(LATd >= DOM(1,2),1,'last');
            ENDlon       = find(LONd >= DOM(2,2),1,'first'); ENDlat = find(LATd <= DOM(1,1),1,'first');
            if isempty(ENDlon), ENDlon = numel(LONd); end; if isempty(ENDlat), ENDlat = numel(LATd); end
            CNTlon       = numel(STRTlon:ENDlon);
            CNTlat       = numel(STRTlat:ENDlat);
            LATd         = ncread(D2C_file,'lat',STRTlat,CNTlat);
            LONd         = ncread(D2C_file,'lon',STRTlon,CNTlon);
            DIST         = ncread(D2C_file,'dist',[STRTlon STRTlat],[CNTlon CNTlat])';
            
            %Set distance to 0 in case point is on land
            IDXnan       = isnan(DIST);
            DIST(IDXnan) = 0;
            
            %Interpolate distance to data locations
            [LONd,LATd]  = meshgrid(double(LONd),double(LATd));
            D2C          = interp2(LONd,LATd,DIST,LON,LAT,'linear');
        end
        
        function MSK = LoadICEmsk
            %% Settings
            %General settings
            LoadCommonSettings

            %Load mask
            PathMSK    = fullfile(PathDATA,'Geography','GIMP','90');
            GTinfo     = geotiffinfo(fullfile(PathMSK,'GimpIceMask_90m.tif'));
            MSK.('x')  = linspace(GTinfo.BoundingBox(1,1),GTinfo.BoundingBox(2,1),GTinfo.Width);
            MSK.('y')  = linspace(GTinfo.BoundingBox(2,2),GTinfo.BoundingBox(1,2),GTinfo.Height);
            MSK.('z')  = geotiffread(fullfile(PathMSK,'GimpIceMask_90m.tif'));
            clear('GTinfo')
        end
        
        function MSK = LoadLSmsk(DOM)
            %LoadLSmsk compiles a land/sea mask
            
            %% Settings
            %General settings
            LoadCommonSettings
            
            try
                load(fullfile(PathRSLT,'LSmsk',strcat('MSK',sprintf('_%s.mat',DARTutils.GenerateStrDOM(DOM)))));
            catch exception
                warning('%s\n',exception.message)
                
                PathMSK          = fullfile(PathDATA,'DEM','MAT_EuroDEM_SRTM_ASTGTM2_GEBCO2014'); %Path to LSmask
                IDX              = 0;
                [LATc,LONc,MSKc] = deal(cell(numel(DOM(1,1):DOM(1,2))+2,numel(DOM(2,1):DOM(2,2))+2));
                for ew = DOM(2,1)-1:DOM(2,2)+1
                    if ew >= 0, ew_str = 'e'; else ew_str = 'w'; end
                    for ns = DOM(1,1)-1:DOM(1,2)+1
                        if ns >= 0, ns_str = 'n'; else ns_str = 's'; end
                        load(fullfile(PathMSK,sprintf('%s%02i%s%02i.mat',ew_str,ew,ns_str,ns)),'LAT','LON','LSmask');
                        IDX                 = IDX+1;
                        LATc{IDX}           = LAT'; LONc{IDX} = LON; MSKc{IDX} = LSmask;
                    end
                end
                LATc = cell2mat(LATc); LONc = cell2mat(LONc); MSKc = cell2mat(MSKc);
                MSK.('x') = LONc(1,:); MSK.('y') = LATc(:,1); MSK.('z') = MSKc;
                clear('IDX','LATc','LONc','MSKc','LSmask')
                
                %Create 2-D interpolant
                [X,Y] = ndgrid(MSK.x,MSK.y);
                MSK   = griddedInterpolant(X,Y,single(MSK.z)','nearest');
                
                %Save compiled DEM
                save(fullfile(PathRSLT,'LSmsk',strcat('MSK',sprintf('_%s.mat',DARTutils.GenerateStrDOM(DOM)))),'MSK');
            end
        end
        
        function [fit_param,Resnorm,exitflag] = SolveNLLS(fun,x0,XDATA,YDATA,lb,ub,options)
            [fit_param,Resnorm,~,exitflag] = lsqcurvefit(fun,x0,XDATA,YDATA,lb,ub,options);
        end
        
        function [fit_param,Resnorm,exitflag] = SolveTRRms(fun,x0,XDATA,YDATA,lb,ub,ms)
            problem                        = createOptimProblem('lsqcurvefit','x0',x0,'objective',fun,'lb',lb,'ub',ub,'xdata',XDATA,'ydata',YDATA);
            [fit_param,Resnorm,~,exitflag] = run(ms,problem,50);
        end
        
    end
end
