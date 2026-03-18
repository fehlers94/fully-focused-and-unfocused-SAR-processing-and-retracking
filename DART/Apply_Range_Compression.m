function CS1b = Apply_Range_Compression(CS1a,CS1b,DDAcf)

%APPLY_RANGE_COMPRESSION conducts the range compression of the input bursts
%and then generates the power waveforms

%Input:
%CS1a: geocorrected beams
%CS1b: stack indices

%Output:
%CS1b.S**.Stack: range compressed stacks

switch CS1a.GEO.OPERATION_MODE
    case 'SIR_FBR_SAR'
        %Preliminaries
        CS1b.SAR.('Stack') = nan(DDAcf.os_ZP*DDAcf.Np,size(CS1a.FBR.data,2),size(CS1a.FBR.data,3),size(CS1a.FBR.data,4));
        
        %Apply range compression for each stack
        for i=find(~isnan(CS1b.GEO.LAT))'
            %Apply range FFT (apply zero padding)
            TMP = fftshift(fft(CS1a.FBR.data(:,CS1b.GEO.IDXstack{i}),DDAcf.os_ZP*DDAcf.Np,1),1);
            
            %Generate power waveforms (here you apply the geometry mask)
            CS1b.SAR.Stack(:,CS1b.GEO.IDXstack{i}) = abs(TMP.*CS1a.FBR.Hmask(:,CS1b.GEO.IDXstack{i})).^2;
        end
    case 'SIR_FBR_SARIN'
        %Preliminaries
        CS1b.SIN.('Stack_ch1')          = nan(DDAcf.os_ZP*DDAcf.Np,size(CS1a.FBR.data_ch1,2),size(CS1a.FBR.data_ch1,3),size(CS1a.FBR.data_ch1,4));
        CS1b.SIN.('Stack_ch2')          = nan(DDAcf.os_ZP*DDAcf.Np,size(CS1a.FBR.data_ch2,2),size(CS1a.FBR.data_ch2,3),size(CS1a.FBR.data_ch2,4));
        CS1b.SIN.('phase_difference')   = nan(DDAcf.os_ZP*DDAcf.Np,size(CS1b.GEO.IDXstack,1),size(CS1b.GEO.IDXstack,2));
        CS1b.SIN.('coherence')          = nan(DDAcf.os_ZP*DDAcf.Np,size(CS1b.GEO.IDXstack,1),size(CS1b.GEO.IDXstack,2));
        
        %Apply range compression for each stack
        for i=find(~isnan(CS1b.GEO.LAT))'
            %Apply range FFT (apply zero padding)
            TMP_ch1 = fftshift(fft(CS1a.FBR.data_ch1(:,CS1b.GEO.IDXstack{i}),DDAcf.os_ZP*DDAcf.Np,1),1);
            TMP_ch2 = fftshift(fft(CS1a.FBR.data_ch2(:,CS1b.GEO.IDXstack{i}),DDAcf.os_ZP*DDAcf.Np,1),1);
            
            %Compute Cross-product (See Wingham et al.: CryoSat: A mission to determine the fluctuations in Earth's land and marine ice fields )
            CrossProduct_ml                       = nanmean(TMP_ch1.*conj(TMP_ch2),2);
            
            %Compute Phase_difference and correct for Cal4 phase difference
            %correction as calculated in "Apply Calibration"
            CS1b.SIN.phase_difference(:,i)        = wrapToPi(angle(CrossProduct_ml)-CS1b.MEA.int_phase_corr(i));
            
            %Compute Multi-looked coherence
            CS1b.SIN.coherence(:,i)               = (abs(CrossProduct_ml))./(nanmean((TMP_ch2.*conj(TMP_ch2)),2));
            
            %Generate power waveforms (here you apply the geometry mask)
            CS1b.SIN.Stack_ch1(:,CS1b.GEO.IDXstack{i}) = abs(TMP_ch1.*CS1a.FBR.Hmask(:,CS1b.GEO.IDXstack{i})).^2;
            CS1b.SIN.Stack_ch2(:,CS1b.GEO.IDXstack{i}) = abs(TMP_ch2.*CS1a.FBR.Hmask(:,CS1b.GEO.IDXstack{i})).^2;
            
            %Different method of calculated multilooked Power
            %CS1b.SIN.data(:,i)                   = abs(CrossProduct_ml)
        end
    otherwise
        error('OPERATION_MODE: %s not recognized')
end

end
