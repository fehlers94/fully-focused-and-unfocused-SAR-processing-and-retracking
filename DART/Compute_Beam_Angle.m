function CS1a = Compute_Beam_Angle(CS1a,CS1b,DDAcf)

%COMPUTE_BEAM_ANGLE computes, for every burst, the angles between the
%satellite velocity vector and the directions defined by the satellite
%location and the computed surface locations under the satellite's
%boresight. Moreover, it computes the angle at which the surfaces samples
%are seen w.r.t. the nadir direction of the satellite (look angle).

%Input:
%CS1a: Structure containing data fields read from the SAR FBR .DBL file
%CS1b: Structure containing data fields read from the SAR L1B .DBL file

%Output:
%CS1a.GEO.SurfLocIDs: indices of surface locations "observed" at each burst
%                     position. Note that columns in which the surface
%                     location indices/beam angles are stored correspond to
%                     the particular beams that should be steered to
%                     identified surface location.
%CS1a.GEO.BeamAngle:  beam angles to "observed" surface locations
%CS1a.GEO.LookAngle:  look angles

%% Preliminaries
%Allocate arrays to store IDs of surface locations "observed" by the
%at each satellite burst position, the associated beam angle, and the ID
%and beam angle to the nearest surface location
CS1a.GEO.('SurfLocIDs') = nan(numel(CS1a.GEO.LAT),DDAcf.Nb);
CS1a.GEO.('BeamAngle')  = nan(numel(CS1a.GEO.LAT),DDAcf.Nb);
CS1a.GEO.('LookAngle')  = nan(numel(CS1a.GEO.LAT),DDAcf.Nb);

%Transform geodetic coordinates of the computed surface locations to
%Cartesian coordinates
SurfLoc = lla2ecef([CS1b.GEO.LAT(:),CS1b.GEO.LON(:),CS1b.GEO.h_DEM(:)]);

%Transform geodetic coordinates of the burst satellite and burst surface
%locations to Cartesian coordinates
BSatPos = lla2ecef([CS1a.GEO.LAT(:),CS1a.GEO.LON(:),CS1a.GEO.H(:)]);
BSurPos = lla2ecef([CS1a.GEO.LAT(:),CS1a.GEO.LON(:),CS1a.GEO.h_DEM(:)]);

%Compute angular limits
%(http://www.radartutorial.eu/01.basics/Doppler%20Dilemma.en.html):
%v_{sat}*cos(theta) < c*f_{PRF} / 4*f_{TX}
q_min = acos(((DDAcf.c*DDAcf.fp) / (4*DDAcf.fc)) ./ CS1a.GEO.V.V);
q_max = pi - q_min;

%Determine whether track is ascending/descending
if CS1a.GEO.LAT(find(~CS1a.GEO.IDXfd,1,'first')) < CS1a.GEO.LAT(find(~CS1a.GEO.IDXfd,1,'last'))
    track_sign = -1; %ascending track
else
    track_sign = 1;  %descending track
end

%% For each satellite burst position find Np nearest surface locations + compute beam angles
for i=find(~CS1a.GEO.IDXfd)'
%     %Find the Np nearest surface locations 
%     [~,IDXs]   = sort(abs(CS1b.GEO.LAT(:)-CS1a.GEO.LAT(i)));
%     IDXnearest = IDXs(1);
%     IDXs       = sort(IDXs(1:DDAcf.Nb));
    
    %Find the Np nearest surface locations (nearest defined by Beam Angle).
    %Is NOT correct as terrain might be sloped.
    Surf_BSat  = bsxfun(@minus, SurfLoc, BSatPos(i,:));
    BeamAngle  = acos( (Surf_BSat*[CS1a.GEO.V.Vx(i);CS1a.GEO.V.Vy(i);CS1a.GEO.V.Vz(i)]) ./ (CS1a.GEO.V.V(i)*sqrt(sum(Surf_BSat.^2,2))));
    [~,IDXs]   = sort(abs(BeamAngle-(.5*pi)));
    IDXs       = sort(IDXs(1:DDAcf.Nb));

    %Create vectors from current burst satellite position to identified
    %surface locations
    Surf_BSat  = bsxfun(@minus, SurfLoc(IDXs,:), BSatPos(i,:));
    
    %Compute beam angles
    BeamAngle  = acos( (Surf_BSat*[CS1a.GEO.V.Vx(i);CS1a.GEO.V.Vy(i);CS1a.GEO.V.Vz(i)]) ./ (CS1a.GEO.V.V(i)*sqrt(sum(Surf_BSat.^2,2))));

    %Create vectors from current burst satellite position to burst surface
    %locations
    Bsur_BSat  = bsxfun(@minus, BSurPos(i,:), BSatPos(i,:))';

    %Compute look angles (angle = atan2(norm(cross(a,b)), dot(a,b)))
    Enf               = ones(size(Surf_BSat(:,1)));
    LookAngle         = atan2(sqrt(sum((cross(Surf_BSat,[Enf*Bsur_BSat(1),Enf*Bsur_BSat(2),Enf*Bsur_BSat(3)])).^2,2)),(Surf_BSat*Bsur_BSat));
    IDXbsf            = track_sign*CS1b.GEO.LAT(IDXs) > track_sign*CS1a.GEO.LAT(i);
    LookAngle(IDXbsf) = -LookAngle(IDXbsf);
    
    %Evaluate whether beam angles are within angular limits
    IDXb       = BeamAngle >= q_min(i) & BeamAngle <= q_max(i);
    
    %Store IDX of surface locations and beam/look angles. ColumnID in which
    %IDX and beam/look angles are stored, identifies that the beam
    %associated to ColumnID will be steered to the surface location IDX
    CS1a.GEO.SurfLocIDs(i,DDAcf.Nb-flipud(find(IDXb))+1) = IDXs(IDXb);
    CS1a.GEO.BeamAngle(i,DDAcf.Nb-flipud(find(IDXb))+1)  = BeamAngle(IDXb);
    CS1a.GEO.LookAngle(i,DDAcf.Nb-flipud(find(IDXb))+1)  = LookAngle(IDXb);
end

end
