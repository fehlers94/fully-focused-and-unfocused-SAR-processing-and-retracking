function CS1b = Compute_Surface_Locations(CS1a,CS1b,DDAcf,DEM)

%COMPUTE_SURFACE_LOCATIONS computes the surface locations (and their
%corresponding datation and orbit parameters) defined by the intersection
%of the Doppler beams and the estimated surface positions along the
%satellite track

%Input:
%CS1a: Structure containing data fields read from the SAR FBR .DBL file
%CS1b: Empty structure containing data fields read from a SAR L1B .DBL file

%Output:
%CS1b.GEO.LAT, CS1b.GEO.LON: computed surface locations
%The following fields are interpolated from the CS1a data fields at the
%computed surface locations: CS1b.GEO.H, CS1b.MEA.win_delay,
%CS1b.GEO.H_rate, CS1b.GEO.V.Vx, CS1b.GEO.V.Vy, CS1b.GEO.V.Vz,
%CS1b.GEO.V.V, CS1b.GEO.TAI.days, CS1b.GEO.TAI.secs,
%CS1b.GEO.TAI.microsecs, CS1b.GEO.BURST_CNT, CS1b.GEO.h_DEM (if DEM is set
%it is interpolated from the DEM)

%% Preliminaries
%Select valid burst satellite positions
LATc = CS1a.GEO.LAT(~CS1a.GEO.IDXfd); LONc = CS1a.GEO.LON(~CS1a.GEO.IDXfd);

%Construct array of valid entries of fields H, h_DEM, H_rate, Vx, Vy, Vz,
%V, and TAI
ALLc = [CS1a.GEO.H(:),CS1a.GEO.h_DEM(:),CS1a.GEO.H_rate(:),CS1a.GEO.V.Vx(:),CS1a.GEO.V.Vy(:),CS1a.GEO.V.Vz(:),CS1a.GEO.V.V(:),CS1a.GEO.TAI.days(:),CS1a.GEO.TAI.secs(:),CS1a.GEO.TAI.microsecs(:),CS1a.MEA.win_delay(:),CS1a.GEO.BaseLine.X(:),CS1a.GEO.BaseLine.Y(:),CS1a.GEO.BaseLine.Z(:),CS1a.GEO.Beam.X(:),CS1a.GEO.Beam.Y(:),CS1a.GEO.Beam.Z(:)];
ALLc = ALLc(~CS1a.GEO.IDXfd,:);

%Compute the piecewise polynomial forms of the cubic spline interpolants to
%the data values ALLc at latitudes LATc
pp   = spline(LATc,ALLc');

%% Compute along-track sampling interval & fill gaps
%Compute along-track sampling interval
DIST    = distance(LATc(1:end-1),LONc(1:end-1),LATc(2:end),LONc(2:end),DDAcf.RefEll)/1000; %[km]
SmplInt = median(DIST);

%Fill gaps
if any(DIST > 1.9*SmplInt)
    %Densify latitude-longitude sampling in case of gaps
    [LATc,LONc] = interpm(LATc,LONc,1.9*km2deg(SmplInt),'gc');
    
    %Interpolate H, H_rate, win_delay, Vx, Vy, Vz, V, and TAI
    ALLc = ppval(pp,LATc)';
    
    %In case a DEM is exploited interpolate DEM to (LATc,LONc)
    if ~isempty(DEM)
        %ALLc(:,2) = DEM(LONc,LATc);
        ALLc(:,2) = double(DARTutils.EvalDEMatLATLON(DEM,LATc,LONc));
    end
end
clear('DIST','SmplInt')

%% Compute surface locations
%Transform geodetic coordinates of the burst surface positions (see Figure
%4.2-4 of SCOOP ATBD V1.1) to Cartesian coordinates
[Xnc,Ync,Znc] = geodetic2ecef(DDAcf.RefEll,LATc,LONc,ALLc(:,2));     %burst surface positions
Enc           = ones(size(Xnc));

%Initialization
IDX1       = 1; ii=0; CS1b.GEO.LAT(:) = NaN;
%Transform geodetic coordinates of FIRST burst satellite position and FIRST
%burst surface position == FIRST surface location to Cartesian coordinates
[Xp,Yp,Zp] = geodetic2ecef(DDAcf.RefEll,LATc(1),LONc(1),ALLc(1,1));  %platform = satellite
[Xs,Ys,Zs] = geodetic2ecef(DDAcf.RefEll,LATc(1),LONc(1),ALLc(1,2));  %first surface location
Vtmp       = ALLc(1,7);

%Start coarse and fine intersection loops
while IDX1 < numel(LATc)
    %Compute angular azimuth beam resolution
    AABR     = asind(DDAcf.lambda0./(2*DDAcf.Nb*DDAcf.PRI.*Vtmp));
    
    %Compute the angle of sight to each burst surface position
    %[LATc,LONc,h_DEM] = angle between nadir direction and vector from
    %satellite position to each burst surface position
    Vsp      = [Xs;Ys;Zs]-[Xp;Yp;Zp];
    Vncp     = [Xnc,Ync,Znc]-[Enc*Xp,Enc*Yp,Enc*Zp];
    % Alpha1nc = real(acosd( (Vncp*Vsp) ./ (norm(Vsp)*sqrt(sum(Vncp.^2,2))) ));
    Alpha1nc = atan2d(sqrt(sum((cross(Vncp,[Enc*Vsp(1),Enc*Vsp(2),Enc*Vsp(3)])).^2,2)),(Vncp*Vsp));
    
    %Find index of point for which applies Alpha1nc(IDX1) < AABR < Alpha1nc(IDX1+1)
    IDX1     = find(Alpha1nc < AABR,1,'last');

    %Step out if you are at end of track
    if IDX1+1 > numel(LATc), break, end
    
    %Densify latitude-longitude sampling between points (LATc(IDX1:IDX1+1),LONc(IDX1:IDX1+1))
    [LATf,LONf] = interpm(LATc(IDX1:IDX1+1),LONc(IDX1:IDX1+1),km2deg(1E-3),'gc');

    %Evaluate H, h_DEM, H_rate, Vx, Vy, Vz, and V at locations (LATf,LONf)
    ALLf     = ppval(pp,LATf)';
    
    %In case a DEM is exploited interpolate DEM to (LATf,LONf)
    if ~isempty(DEM)
        %ALLf(:,2) = DEM(LONf,LATf);
        ALLf(:,2) = double(DARTutils.EvalDEMatLATLON(DEM,LATf,LONf));
    end
    
    %Compute the angle of sight to the densified burst surface positions [LATf,LONf,h_DEM]
    [Xnf,Ynf,Znf] = geodetic2ecef(DDAcf.RefEll,LATf,LONf,ALLf(:,2)); %nadir locations
    Enf           = ones(size(Xnf));
    Vnfp          = [Xnf,Ynf,Znf]-[Enf*Xp,Enf*Yp,Enf*Zp];
    % Alpha1nf      = real(acosd( (Vnfp*Vsp) ./ (norm(Vsp)*sqrt(sum(Vnfp.^2,2))) ));
    Alpha1nf      = atan2d(sqrt(sum((cross(Vnfp,[Enf*Vsp(1),Enf*Vsp(2),Enf*Vsp(3)])).^2,2)),(Vnfp*Vsp));

    %Compute latitude and longitude of surface location + evaluate H,
    %h_DEM, H_rate, Vx, Vy, Vz, and V at this location
    LATs = interp1(Alpha1nf,LATf,AABR,'linear');
    LONs = interp1(Alpha1nf,LONf,AABR,'linear');
    ALLs = ppval(pp,LATs);
    
    %In case a DEM is exploited interpolate DEM to (LATs,LONs)
    if ~isempty(DEM)
        %ALLs(:,2) = DEM(LONs,LATs);
        ALLs(:,2) = double(DARTutils.EvalDEMatLATLON(DEM,LATs,LONs));
    end
    
    %Transform geodetic coordinates of the computed surface location and
    %its projection at satellite altitude to Cartesian coordinates
    [Xp,Yp,Zp] = geodetic2ecef(DDAcf.RefEll,LATs,LONs,ALLs(1)); %surface location at platform altitude
    [Xs,Ys,Zs] = geodetic2ecef(DDAcf.RefEll,LATs,LONs,ALLs(2)); %surface location
    Vtmp       = ALLs(7);

    %Copy info to CS1b
    ii                         = ii+1;
    CS1b.GEO.LAT(ii)           = LATs;
    CS1b.GEO.LON(ii)           = LONs;
    CS1b.GEO.H(ii)             = ALLs(1);
    CS1b.GEO.H_rate(ii)        = ALLs(3);
    CS1b.GEO.V.Vx(ii)          = ALLs(4);
    CS1b.GEO.V.Vy(ii)          = ALLs(5);
    CS1b.GEO.V.Vz(ii)          = ALLs(6);
    CS1b.GEO.V.V(ii)           = ALLs(7);
    CS1b.GEO.TAI.days(ii)      = ALLs(8);
    CS1b.GEO.TAI.secs(ii)      = ALLs(9);
    CS1b.GEO.TAI.microsecs(ii) = ALLs(10);
    CS1b.MEA.win_delay(ii)     = ALLs(11);
    CS1b.GEO.BaseLine.X(ii)    = ALLs(12);
    CS1b.GEO.BaseLine.Y(ii)    = ALLs(13);
    CS1b.GEO.BaseLine.Z(ii)    = ALLs(14);
    CS1b.GEO.Beam.X(ii)        = ALLs(15);
    CS1b.GEO.Beam.Y(ii)        = ALLs(16);
    CS1b.GEO.Beam.Z(ii)        = ALLs(17);
    CS1b.GEO.BURST_CNT(ii)     = ii;
end

%Interpolate DEM to computed surface locations or compute h_DEM using window delay
if ~isempty(DEM)
    %CS1b.GEO.('h_DEM') = DEM(CS1b.GEO.LON,CS1b.GEO.LAT);
    CS1b.GEO.('h_DEM') = double(DARTutils.EvalDEMatLATLON(DEM,CS1b.GEO.LAT,CS1b.GEO.LON));
else
    CS1b.GEO.('h_DEM') = CS1b.GEO.H-((CS1b.MEA.win_delay*DDAcf.c)/2);
end

%% Edit arrays
IDX                               = find(isnan(CS1b.GEO.LAT(1,:)),1,'first');
CS1b.GEO.LAT(:,IDX:end)           = [];
CS1b.GEO.LON(:,IDX:end)           = [];
CS1b.GEO.H(:,IDX:end)             = [];
CS1b.GEO.H_rate(:,IDX:end)        = [];
CS1b.GEO.V.Vx(:,IDX:end)          = [];
CS1b.GEO.V.Vy(:,IDX:end)          = [];
CS1b.GEO.V.Vz(:,IDX:end)          = [];
CS1b.GEO.V.V(:,IDX:end)           = [];
CS1b.GEO.TAI.days(:,IDX:end)      = [];
CS1b.GEO.TAI.secs(:,IDX:end)      = [];
CS1b.GEO.TAI.microsecs(:,IDX:end) = [];
CS1b.GEO.BURST_CNT(:,IDX:end)     = [];
CS1b.GEO.h_DEM(:,IDX:end)         = [];
CS1b.MEA.win_delay(:,IDX:end)     = [];
CS1b.GEO.BaseLine.X(:,IDX:end)    = [];
CS1b.GEO.BaseLine.Y(:,IDX:end)    = [];
CS1b.GEO.BaseLine.Z(:,IDX:end)    = [];
CS1b.GEO.Beam.X(:,IDX:end)        = [];
CS1b.GEO.Beam.Y(:,IDX:end)        = [];
CS1b.GEO.Beam.Z(:,IDX:end)        = [];

%Compute the piecewise polynomial forms of the cubic spline interpolants to
%the Internal Phase Correction Values at latitudes for which Cal4_Flag == 1
if isequal(CS1a.GEO.OPERATION_MODE,'SIR_FBR_SARIN')
    CS1b.MEA.int_phase_corr       = spline(CS1a.GEO.LAT(CS1a.GEO.MODE_ID.Cal4_Flag == 1),CS1a.MEA.int_phase_corr,CS1b.GEO.LAT);
end

end

