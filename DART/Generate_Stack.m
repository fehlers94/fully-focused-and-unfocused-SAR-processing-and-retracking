function CS1b = Generate_Stack(CS1a,CS1b,DDAcf)

%GENERATE_STACK identifies all beams within all bursts contributing to a
%particular surface location.

%Input:
%CS1a: Indices of identified surface locations "observed" from each burst
%position

%Output:
%CS1b.GEO.IDXbursts = indices of bursts contributing to a surface location
%CS1b.GEO.IDXbeams  = indices of beams within bursts contributing to a surface location
%CS1b.GEO.IDXstack  = indices of waveforms that belong to stack

%% Settings

%% Preliminaries
%Burst and beam indices will be stored in cell array
CS1b.GEO.('IDXbursts') = cell(size(CS1b.GEO.LAT));
CS1b.GEO.('IDXbeams')  = cell(size(CS1b.GEO.LAT));
CS1b.GEO.('IDXstack')  = cell(size(CS1b.GEO.LAT));

switch CS1a.GEO.OPERATION_MODE
    case 'SIR_FBR_SAR'
        SizeDATA = [size(CS1a.FBR.data,2),size(CS1a.FBR.data,3)*size(CS1a.FBR.data,4)];
    case 'SIR_FBR_SARIN'
        SizeDATA = [size(CS1a.FBR.data_ch1,2),size(CS1a.FBR.data_ch1,3)*size(CS1a.FBR.data_ch1,4)];
    otherwise
        error('OPERATION_MODE: %s not recognized')
end

%% Identify all beams within all bursts contributing to each particular surface location
for i=find(~isnan(CS1b.GEO.LAT))'
    %Identify beam and burst indices
    [CS1b.GEO.IDXbeams{i},CS1b.GEO.IDXbursts{i}] = find(ismember(CS1a.GEO.SurfLocIDs,i)');
    
%     %Remove IDX of beams for which beam angle is outside acceptable range
%     if DDAcf.beam_selection
%         BeamAngle = CS1a.GEO.BeamAngle(sub2ind(size(CS1a.GEO.BeamAngle),CS1b.GEO.IDXbursts{i},CS1b.GEO.IDXbeams{i}));
%         IDXb      = BeamAngle >= DDAcf.Lwlt & BeamAngle <= DDAcf.Uplt;
%         BeamID    = CS1b.GEO.IDXbeams{i};  CS1b.GEO.IDXbeams{i}  = BeamID(IDXb);
%         BurstID   = CS1b.GEO.IDXbursts{i}; CS1b.GEO.IDXbursts{i} = BurstID(IDXb);
%     end
   
    %Store linear index of beams
    CS1b.GEO.IDXstack{i} = sub2ind(SizeDATA,CS1b.GEO.IDXbeams{i},CS1b.GEO.IDXbursts{i});
end

end





