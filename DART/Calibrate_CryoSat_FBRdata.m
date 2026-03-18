function CS1a = Calibrate_CryoSat_FBRdata(CS1a,DDAcf)

%CALIBRATE_CRYOSAT_FBRDATA applies power calibration, the interpulse
%amplitude and phase corrections, and the Low Pass Filter correction to
%CryoSat FBR data. See Dinardo (2016) and Scagliola & Fornari (2016). Both
%documents are available at
%https://wiki.services.eoportal.org/tiki-index.php?page=CryoSat+Technical+Notes
%
%Dinardo (2016), Guidelines for reverting Waveform Power to Sigma Nought
%for CryoSat-2 in SAR mode, Issue 2, Revision 2.
%NOTE THAT Eq. 15 IN THIS DOCUMENT IS NOT CORRECT (confirmed by Salvatore
%Dinardo, 20-6-2017, 12:18)
%
%Scagliola and Fornari (2016), CryoSat Characterization for FBR users,
%Issue 2.0.
%
%In Dinardo (2016), it is assumed that that a not unitary Fast Fourier
%Transform FFT implementation is used. Its definition, given in Annex C
%ibid. is in line with Matlab's FFT implementation.

%Input:
%CS1a: Structure containing data fields read from the SAR FBR .DBL file

%Output:
%CS1a.FBR.data: contains calibrated FBR data

%% FBR INDIVIDUAL ECHOES CALIBRATION IN POWER
%Evaluation of the Static Part for the Power Gain (Dinardo 2016, Section 2.1)
GAIN_ADC               = 10*log10(DDAcf.ADC_MULT^2);
GAIN_PROC_RANGE        = 10*log10(DDAcf.Np^2);                    %DDAcf.Np = Number of SAR/SARIn FBR Echo Samples = 128/512
GAIN_PROC_DOPPLER      = 10*log10(DDAcf.Nb^2);
GAIN_DIGITAL           = GAIN_ADC+GAIN_PROC_RANGE+GAIN_PROC_DOPPLER;
GAIN_STATIC_RX1        = GAIN_DIGITAL+CS1a.MEA.Gain_Rx1;
GAIN_STATIC_RX2        = GAIN_DIGITAL+CS1a.MEA.Gain_Rx2;

%Power Gain: Dynamic Part for Baseline C FBR products (Dinardo 2016, Section 2.3)
AGC_SETTING_CORRECTION = CS1a.MEA.AGC_1 + CS1a.MEA.AGC_2;                  %Automatic Gain Control Setting Correction
GAIN_DYNAMIC_RX1       = AGC_SETTING_CORRECTION + CS1a.MEA.ins_gain_corr_rx_tx;
GAIN_DYNAMIC_RX2       = AGC_SETTING_CORRECTION + CS1a.MEA.ins_gain_corr_rx;

%Total Power Gain (Dinardo 2016, Section 2.4)
GAIN_RX1               = GAIN_DYNAMIC_RX1 + GAIN_STATIC_RX1;
GAIN_RX2               = GAIN_DYNAMIC_RX2 + GAIN_STATIC_RX2;

%Since we do not divide the output of the FFT by sqrt(# samples), no
%compensation for zero-padding scale factor is required as described in
%Dinardo (2016), section 3.1.1.

%Compensation for Weighting Application (Note that Eq. 15 provided in Dinardo (2016) is not correct!)
if DDAcf.ApplyWindow
    if strcmpi(DDAcf.Window,'hamming')
        SCALE_WIN = 10*log10(mean(hamming(DDAcf.Nb).^2));
    elseif strcmpi(DDAcf.Window,'hanning')
        SCALE_WIN = 10*log10(mean(hanning(DDAcf.Nb).^2));
    else
        error('Window not reckognized!');
    end
else
    SCALE_WIN = 0;
end

%% Apply calibration corrections/compensation
switch CS1a.GEO.OPERATION_MODE
    case 'SIR_FBR_SAR'
        %Apply power calibration.
        TMP = bsxfun(@times,double(CS1a.FBR.data),reshape(sqrt(10.^(-(GAIN_RX1/10))),1,1,size(GAIN_RX1,1),size(GAIN_RX1,2)));
        
        %Apply ground processing gain correction (applied to echoes, hence sqrt(*))
        TMP = TMP*sqrt(10.^(-SCALE_WIN/10));
        
        %Apply pulse-to-pulse amplitude and pulse-to-pulse phase corrections (CryoSat Characterization for FBR users)
        TMP = bsxfun(@times,TMP,DDAcf.p2p_amp'.*exp(sqrt(-1)*DDAcf.p2p_phs'));
        
        %Apply Low Pass Filter amplitude gain correction (CryoSat Characterization for FBR users)
        CS1a.FBR.data = ifft(ifftshift(bsxfun(@times,fftshift(fft(TMP,[],1),1),DDAcf.lpf_sar),1),[],1);
        
    case 'SIR_FBR_SARIN'
        %Internal Phase Calibration (Cal4 mode)
        Cal4FBR_ch1             = double(CS1a.FBR.data_ch1(:,:,CS1a.GEO.MODE_ID.Cal4_Flag == 1));
        Cal4FBR_ch2             = double(CS1a.FBR.data_ch2(:,:,CS1a.GEO.MODE_ID.Cal4_Flag == 1));
        %Compute cross-product
        Cal4_CP                 = (Cal4FBR_ch1.*conj(Cal4FBR_ch2));
        %Compute Phase_difference
        CS1a.MEA.int_phase_corr = squeeze(nanmean(nanmean(angle(Cal4_CP))));
        
        %Apply power calibration. Note that size of CH1 is not 4D
        TMPRX1 = bsxfun(@times,double(reshape(CS1a.FBR.data_ch1,size(CS1a.FBR.data_ch2))),reshape(sqrt(10.^(-(GAIN_RX1/10))),1,1,size(GAIN_RX1,1),size(GAIN_RX1,2)));
        TMPRX2 = bsxfun(@times,double(CS1a.FBR.data_ch2),reshape(sqrt(10.^(-(GAIN_RX2/10))),1,1,size(GAIN_RX2,1),size(GAIN_RX2,2)));

        %Apply ground processing gain correction (applied to echoes, hence sqrt(*))
        TMPRX1 = TMPRX1*sqrt(10.^(-SCALE_WIN/10));
        TMPRX2 = TMPRX2*sqrt(10.^(-SCALE_WIN/10));
        
        %Apply pulse-to-pulse amplitude and pulse-to-pulse phase corrections (CryoSat Characterization for FBR users)
        TMPRX1 = bsxfun(@times,TMPRX1,DDAcf.p2p_amp_rx1'.*exp(sqrt(-1)*DDAcf.p2p_phs_rx1'));
        TMPRX2 = bsxfun(@times,TMPRX2,DDAcf.p2p_amp_rx2'.*exp(sqrt(-1)*DDAcf.p2p_phs_rx2'));

        %Apply Low Pass Filter amplitude gain correction (CryoSat Characterization for FBR users)
        CS1a.FBR.data_ch1 = ifft(ifftshift(bsxfun(@times,fftshift(fft(TMPRX1,[],1),1),DDAcf.lpf_sarin_rx1),1),[],1);
        CS1a.FBR.data_ch2 = ifft(ifftshift(bsxfun(@times,fftshift(fft(TMPRX2,[],1),1),DDAcf.lpf_sarin_rx2),1),[],1);
        
        %SET CAL4 mode measurements to zero
        CS1a.FBR.data_ch1(:,:,CS1a.GEO.MODE_ID.Cal4_Flag == 1) = 0;
        CS1a.FBR.data_ch2(:,:,CS1a.GEO.MODE_ID.Cal4_Flag == 1) = 0;
        
    otherwise
        error('OPERATION_MODE: %s not recognized')
end

end

