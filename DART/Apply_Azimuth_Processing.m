function CS1a = Apply_Azimuth_Processing(CS1a,CS1b,DDAcf)

%AZIMUTH_PROCESSING steer the beams to the different surface locations. Two
%different approaches are considered: the exact method (Figure 4.2­8,
%SCOOP), and the approximate one (Figure 4.2­9, SCOOP), which is a
%simplification of the former one.

%The exact method uses all the beam angles computed in Section 4.2.2 to
%steer the beams to the surfaces. This implies that there will be an FFT
%process for each one of the surface locations. On the other hand, the
%approximate method simply uses the beam angle that is closest to the nadir
%to spread the other beams and steer them to the other surface locations.
%This means that the approximate method only goes through one FFT process.

%To determine which method is used, the variability of the surface (in
%terms of altitude) is computed by use of the standard deviation of the
%altitude of the surface locations that are seen by the current burst.
%Depending on the standard deviation value, being above or below a pre-set
%threshold, the exact or approximate method will be used, respectively. In
%this sense, the exact approach results particularly interesting for areas
%with high dynamic topography variation as it can be the case of some
%coastal regions.

%Input:
%CS1a: Indices of identified surface locations "observed" from each burst
%position, the computed beam angles, and the calibrated waveforms
%CS1b.GEO.H and CS1b.GEO.h_DEM

%Output:
%CS1a.FBR.data: contains steered beams

%% Preliminaries
ElevWD = CS1b.GEO.h_DEM(:);

%Define weighting or window. Weighting (applied before beam forming
%procedure takes place) can be used to minimize the impact of side-lobe
%effects in the Doppler/azimuth PTR and, so, the related Doppler
%ambiguities at the edges of the spectrum (edge beams). This is especially
%important when operating close to the coast as high reflectivity land
%scattering can contaminate the signal of interest. We must take into
%account that this weighting leads to a degradation on the along-track
%resolution or footprint.
if DDAcf.ApplyWindow
    if strcmpi(DDAcf.Window,'hamming')
        W = hamming(DDAcf.Nb)';
    elseif strcmpi(DDAcf.Window,'hanning')
        W = hanning(DDAcf.Nb)';
    else
        error('Window not reckognized!');
    end
else
    W = ones(1,DDAcf.Nb);
end

%% Apply Azimuth Processing
switch CS1a.GEO.OPERATION_MODE
    case 'SIR_FBR_SAR'
        for i=find(~CS1a.GEO.IDXfd)'
            IDXs = CS1a.GEO.SurfLocIDs(i,~isnan(CS1a.GEO.SurfLocIDs(i,:)));
            
            %Get beam angles to all surface locations "observed" by the
            %satellite at current satellite burst position
            BeamAngles = CS1a.GEO.BeamAngle(i,:);
            
            %Use exact method if surface variability exceeds pre-defined threshold,
            %else approximate method is used
            if std(ElevWD(IDXs)) > DDAcf.ThSurfVar
                burst  = CS1a.FBR.data(:,:,i);
                for j = find(~isnan(BeamAngles))
                    %Use beam angle associated to pulse j to spread the other beams
                    %and steer them to the other surface locations
                    burst_ps             = Apply_PhaseShift(burst,BeamAngles(j),DDAcf,CS1a.GEO.V.V(i));
                    %Compute the azimuth FFTs
                    burst_ps_fft         = Apply_AlongTrackFFT(burst_ps,DDAcf,W);
                    %Copy center beam to CS1a.FBR.data
                    CS1a.FBR.data(:,j,i) = burst_ps_fft(:,(DDAcf.Nb/2)+1);
                end
            else
                %Use nadir beam angle to spread the other beams and steer them to
                %the other surface locations
                CS1a.FBR.data(:,:,i) = Apply_PhaseShift(CS1a.FBR.data(:,:,i),BeamAngles(33),DDAcf,CS1a.GEO.V.V(i));
                %Compute the azimuth FFTs
                CS1a.FBR.data(:,:,i) = Apply_AlongTrackFFT(CS1a.FBR.data(:,:,i),DDAcf,W);
            end
        end
    case 'SIR_FBR_SARIN'
        %Allocate field 'PhaseShift'
        CS1a.FBR.('PhaseShift') = nan(DDAcf.Nb,DDAcf.Nb,numel(CS1a.GEO.LAT));
        
        for i=find(~CS1a.GEO.IDXfd)'
            IDXs = CS1a.GEO.SurfLocIDs(i,~isnan(CS1a.GEO.SurfLocIDs(i,:)));
            
            %Get beam angles to all surface locations "observed" by the
            %satellite at current satellite burst position
            BeamAngles = CS1a.GEO.BeamAngle(i,:);
            
            %Use exact method if surface variability exceeds pre-defined threshold,
            %else approximate method is used
            if std(ElevWD(IDXs)) > DDAcf.ThSurfVar
                burst_ch1  = CS1a.FBR.data_ch1(:,:,i);
                burst_ch2  = CS1a.FBR.data_ch2(:,:,i);
                for j = find(~isnan(BeamAngles))
                    %Compute phase shift
                    CS1a.FBR.PhaseShift(j,:,i) = CS1a.GEO.V.V(i) * cos(BeamAngles(j)) *DDAcf.PRI .* (0:DDAcf.Nb-1);
                    %Use beam angle associated to pulse j to spread the other beams
                    %and steer them to the other surface locations
                    burst_ps_ch1               = Apply_PhaseShift(burst_ch1,BeamAngles(j),DDAcf,CS1a.GEO.V.V(i));
                    burst_ps_ch2               = Apply_PhaseShift(burst_ch2,BeamAngles(j),DDAcf,CS1a.GEO.V.V(i));
                    %Compute the azimuth FFTs
                    burst_ps_fft_ch1           = Apply_AlongTrackFFT(burst_ps_ch1,DDAcf,W);
                    burst_ps_fft_ch2           = Apply_AlongTrackFFT(burst_ps_ch2,DDAcf,W);
                    %Copy center beam to CS1a.FBR.data
                    CS1a.FBR.data_ch1(:,j,i)   = burst_ps_fft_ch1(:,(DDAcf.Nb/2)+1);
                    CS1a.FBR.data_ch2(:,j,i)   = burst_ps_fft_ch2(:,(DDAcf.Nb/2)+1);
                end
            else
                %Use nadir beam angle to spread the other beams and steer them to
                %the other surface locations
                CS1a.FBR.data_ch1(:,:,i) = Apply_PhaseShift(CS1a.FBR.data_ch1(:,:,i),BeamAngles(33),DDAcf,CS1a.GEO.V.V(i));
                CS1a.FBR.data_ch2(:,:,i) = Apply_PhaseShift(CS1a.FBR.data_ch2(:,:,i),BeamAngles(33),DDAcf,CS1a.GEO.V.V(i));
                %Compute the azimuth FFTs
                CS1a.FBR.data_ch1(:,:,i) = Apply_AlongTrackFFT(CS1a.FBR.data_ch1(:,:,i),DDAcf,W);
                CS1a.FBR.data_ch2(:,:,i) = Apply_AlongTrackFFT(CS1a.FBR.data_ch2(:,:,i),DDAcf,W);
            end
        end
    otherwise
        error('OPERATION_MODE: %s not recognized')
end

end

function burst = Apply_PhaseShift(burst,beam_angle,DDAcf,vel)

%APPLY_PHASESHIFT applies a phase shift to each pulse of a burst. The phase
%shift depends on the beam angle (angle made by the velocity vector at the
%burst location to the vector joining burst to surface location)

%Input:
%burst:      the 64 waveforms of a burst
%beam angle: angle made by the velocity vector at the burst location to the
%            vector joining burst to surface location
%DDAcf:      structure containing universal constants, sensor parameters,
%            and control data values
%vel:        magnitude of the satellite's velocity at the respective burst
%            location 

%Outputs
%burst:      the 64 waveforms after application of the phase shift

bm_ang_phs = exp(-4 * sqrt(-1) * pi / DDAcf.lambda0 * vel * cos(beam_angle) *DDAcf.PRI .* (0:DDAcf.Nb-1));
burst      = bsxfun(@times,burst,bm_ang_phs);

end


function burst = Apply_AlongTrackFFT(burst,DDAcf,W)

%APPLY_ALONGTRACKFFT applies the along-track FFTs

%Input:
%burst: phase shifted waveforms of particular burst
%DDAcf: structure containing universal constants, sensor parameters, and
%       control data values

%Output:
%burst: azimuth fft'd waveforms

%Apply window
burst = bsxfun(@times,burst,W);

%Apply FFT and fftshift
burst = fftshift(fft(burst,DDAcf.Nb,2),2);

end

