function [CS1a,CS1b] = Apply_Geometric_Corrections(CS1a,CS1b,DDAcf)

%APPLY_GEOMETRIC_CORRECTIONS computes and applies all the corrections
%associated with the geometry of the scenario. These are the Doppler, slant
%range and window delay misalignments corrections. As the stack has already
%been generated, these compensations are performed for each stack.

%The Doppler correction is needed to remove the echoes' frequency shifts
%due to the sensor-­target velocity (compensate for the Doppler offset
%induced by movement of the platform while transmitting/receiving the
%pulse, see Prats--Iraola et al. (2014)). The correction is applied to the
%echoes in time domain, before the FFT step of the range compression.

%In addition, the geometry mask is computed

%Input:
%CS1a: Burst satellite positions, velocities, beam angles to surface
%      locations, and window delays
%CS1b: Computed surface locations, burst/beam indices contributing to
%      surface location, interpolated window delays, and stack indices

%Output:
%CS1a.FBR.data:      geocorrected beams (SAR mode)
%CS1a.FBR.data_ch1/2 geocorrected beams (SARIn mode)
%CS1a.FBR.Hmask:     geometry mask
%CS1b.MEA.win_delay: updated window delays

%% Preliminaries
%Transform geodetic coordinates of the computed surface locations to
%Cartesian coordinates
SurfLoc = lla2ecef([CS1b.GEO.LAT(:),CS1b.GEO.LON(:),CS1b.GEO.h_DEM(:)]);

%Transform geodetic coordinates of the burst satellite locations to
%Cartesian coordinates
BSatPos = lla2ecef([CS1a.GEO.LAT(:),CS1a.GEO.LON(:),CS1a.GEO.H(:)]);

%Set Hmask
switch CS1a.GEO.OPERATION_MODE
    case 'SIR_FBR_SAR'
        SizeDATA   = size(CS1a.FBR.data);
    case 'SIR_FBR_SARIN'
        SizeDATA   = size(CS1a.FBR.data_ch1);
    otherwise
        error('OPERATION_MODE: %s not recognized')
end
CS1a.FBR.('Hmask') = zeros(SizeDATA(1)*DDAcf.os_ZP,SizeDATA(2),SizeDATA(3),SizeDATA(4));
Lngth              = DDAcf.Np*DDAcf.os_ZP;

%% Apply all geometric corrections
DCfact = (-2*DDAcf.tau_u)/(DDAcf.lambda0);                       %Doppler Correction in range bins (i.e., sample units)
% DCfact = (-2*DDAcf.tau_u)/(DDAcf.B*DDAcf.lambda0) * (DDAcf.c/2); %Doppler Correction in meters
% DCfact = (-2*DDAcf.tau_u)/(DDAcf.B*DDAcf.lambda0);               %Doppler Correction in seconds
for i=find(~isnan(CS1b.GEO.LAT))'
    %Identified bursts that contribute to particular surface location
    IDXb      = CS1b.GEO.IDXbursts{i};

    %Doppler Correction in range bin units (i.e., sample units)
    DC        = DCfact*CS1a.GEO.V.V(IDXb).*cos(CS1a.GEO.BeamAngle(sub2ind(size(CS1a.GEO.BeamAngle),IDXb,CS1b.GEO.IDXbeams{i})));

    %Slant Range Correction in range bin units (i.e., sample units)
    Surf_BSat = sqrt(sum(bsxfun(@minus, SurfLoc(i,:), BSatPos(IDXb,:)).^2,2));
%     SRC       = (CS1b.MEA.win_delay(i) - (Surf_BSat * 2 / DDAcf.c)) * DDAcf.B; %DDOP
    SRC       = (((CS1b.GEO.H(i)-CS1b.GEO.h_DEM(i)) * 2 / DDAcf.c) - (Surf_BSat * 2 / DDAcf.c)) * DDAcf.B; %DDOP

    %Window delay misalignments in range bin units (i.e., sample units)
    DWD       = -(CS1b.MEA.win_delay(i)-CS1a.MEA.win_delay(IDXb))*DDAcf.B; %DDOP

    %Correct for saw-tooth behavior in echogram
    CS1b.MEA.win_delay(i) = CS1b.MEA.win_delay(i) + (nanmean(DWD)/DDAcf.B);
    DWD                   = DWD - nanmean(DWD);

    %Apply total correction
    CORRsum   = exp((1:DDAcf.Np)'*( sqrt(-1)*2*pi/DDAcf.Np.*(DC+SRC+DWD)'));
    switch CS1a.GEO.OPERATION_MODE
        case 'SIR_FBR_SAR'
            CS1a.FBR.data(:,CS1b.GEO.IDXstack{i})     = CS1a.FBR.data(:,CS1b.GEO.IDXstack{i}).*CORRsum;
        case 'SIR_FBR_SARIN'
            CS1a.FBR.data_ch1(:,CS1b.GEO.IDXstack{i}) = CS1a.FBR.data_ch1(:,CS1b.GEO.IDXstack{i}).*CORRsum;
            CS1a.FBR.data_ch2(:,CS1b.GEO.IDXstack{i}) = CS1a.FBR.data_ch2(:,CS1b.GEO.IDXstack{i}).*CORRsum;
        otherwise
            error('OPERATION_MODE: %s not recognized')
    end
    
    %Compute geometry mask
    %As the geometry corrections have been applied through an exponential
    %(that is equivalent to a circular shift in the other domain as
    %finite-length signals are considered), some samples may have suffered
    %a wrapping. This has to be solved through a mask and force these
    %samples to zero. This mask is computed using the sum of all the
    %three geometry corrections together
    sample_shift                           = round(ones(Lngth,1)*(DC+SRC+DWD)')*DDAcf.os_ZP;
    geom_mask                              = (1:Lngth)'*ones(1,numel(IDXb)) - sample_shift;
    geom_mask(geom_mask < 0)               = NaN;
    geom_mask(geom_mask > Lngth)           = NaN;
    geom_mask(~isnan(geom_mask))           = 1;
    geom_mask(isnan(geom_mask))            = 0;
    geom_mask(:,all(isnan(geom_mask)))     = NaN;
    CS1a.FBR.Hmask(:,CS1b.GEO.IDXstack{i}) = geom_mask;
end

end