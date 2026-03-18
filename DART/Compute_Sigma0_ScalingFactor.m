function CS1b = Compute_Sigma0_ScalingFactor(CS1a,CS1b,DDAcf)

%COMPUTE_SIGMA0_SCALINGFACTOR computes the scaling factor that allows to
%convert the power of the multi­looked waveform into sigma0 values
%(normalized radar cross section values).

%Based on the classical radar equation of the received power, the scaling
%factor per each beam pointing a specific surface can be obtained using Eq.
%4.2-17 [Scoop document], where the first fractional term is the so called
%instrumental scaling factor and the second one represents the external
%effects.

%Input:
%CS1a:  L1a waveforms
%CS1b:  Multilooked L1b waveforms
%DDAcf: structure containing universal constants, sensor parameters, and
%       control data values

%Output:
%CS1b.GEO.sigma0_scaling_factor: Scaling factor to compute the backscatter
%coefficient

error('This script is not correct!')

%% Preliminaries
%Set field sigma0_scaling_factor
switch CS1a.GEO.OPERATION_MODE
    case 'SIR_FBR_SAR'
        CS1b.GEO.('sigma0_scaling_factor') = zeros(size(CS1b.SAR.N_averaged_echoes));
    case 'SIR_FBR_SARIN'
        CS1b.GEO.('sigma0_scaling_factor') = zeros(size(CS1b.SIN.N_averaged_echoes));
    otherwise
        error('OPERATION_MODE: %s not recognized')
end

%Transform geodetic coordinates of the computed surface locations to
%Cartesian coordinates
SurfLoc = lla2ecef([CS1b.GEO.LAT(:),CS1b.GEO.LON(:),CS1b.GEO.h_DEM(:)]);

%Transform geodetic coordinates of the burst satellite locations to
%Cartesian coordinates
BSatPos = lla2ecef([CS1a.GEO.LAT(:),CS1a.GEO.LON(:),CS1a.GEO.H(:)]);

%Constants taken from https://www.researchgate.net/publication/262374503
%and Guidelines for reverting Waveform Power to Sigma Nought for CryoSat-2
%in SAR mode
rv          = 0.886;
wf          = 1.486*rv;

%% Compute sigma0 scaling factor
for i=find(~isnan(CS1b.GEO.LAT))'
    %Identify bursts that contribute to particular surface location
    IDXb      = CS1b.GEO.IDXbursts{i};
    Surf_BSat = sqrt(sum(bsxfun(@minus, SurfLoc(i,:), BSatPos(IDXb,:)).^2,2));

    %Compute alpha to account for Earth curvature
    N         = DDAcf.RefEll.SemimajorAxis./sqrt(1-DDAcf.RefEll.Eccentricity^2*sind(CS1b.GEO.LAT(i)).^2);
    alpha     = (Surf_BSat+N)/N;
   
    %Compute length and width of rectangular region illuminated by a beam
    Lx        = (DDAcf.lambda0*N)./(2*(DDAcf.Nb/DDAcf.fp).*CS1a.GEO.V.V(IDXb));
    Ly        = sqrt(2*N*DDAcf.PTR_width./alpha);

    %Compute the surface area
    SurfAREA  = 2*Lx.*Ly*wf;
    
    %Compute sigma0 scaling factor per beam
    Sigma0SFB = 10*log10(64) + 30*log10(pi)+ ...
                40*log10(Surf_BSat) - 10*log10(ParamSet.power_tx_ant_ku) -...
                2*ParamSet.antenna_gain_ku -20*log10(DDAcf.lambda0)-...
                10 * log10(SurfAREA) ;
    
    %Implementation of Scoop/dedop
    % azimuth_distance = (((1 + Surf_BSat) *DDAcf.lambda0 .* Surf_BSat) / ParamSet.R )./ (DDAcf.PRI *2 * CS1a.GEO.V.V(IDXb) * DDAcf.Nb);
    % range_distance   = 2 * sqrt((DDAcf.c * Surf_BSat)./(ParamSet.pulse_length * ParamSet.chirp_slope_ku./alpha));
    % SurfAREA         = azimuth_distance .* range_distance;
    % Sigma0SFB        = 10*log10(64) + 30*log10(pi)+ ...
    %                    40*log10(Surf_BSat) - 10*log10(ParamSet.power_tx_ant_ku) -...
    %                    2*ParamSet.antenna_gain_ku -20*log10(DDAcf.lambda0)-...
    %                    10 * log10(SurfAREA)+10 * log10(DDAcf.Np * DDAcf.os_ZP) -...
    %                    10 * log10(ParamSet.pulse_length * ParamSet.pulse_length * DDAcf.s); %DEDOP
    
    %Compute average
    CS1b.GEO.sigma0_scaling_factor(i) = mean(Sigma0SFB);
end
end


