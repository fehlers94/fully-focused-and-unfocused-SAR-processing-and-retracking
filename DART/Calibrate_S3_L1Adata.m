function CS1a = Calibrate_S3_L1Adata(CS1a,DDAcf)
% pendant to Calibrate_CryoSat_FBRdata
% to be completed

% undo burst phase correction
%CS1a.FBR.data = bsxfun(@times,CS1a.FBR.data,exp(-1i*CS1a.MEA.burst_phase_cor(:,1))');
% old: testing with TESTER_FFSAR_Processor_TransponderTest.m indicates in the
% phase plots, that this correction must be applied this way (-i),
% different from the ATBD, which states (+i)
CS1a.FBR.data = bsxfun(@times,CS1a.FBR.data,transpose(exp(+1i*CS1a.MEA.burst_phase_cor(:,1))));
% update: the transpose in matlab is the complex transpose that inverts the
% phase, that caused the weird sign behaviour. This suggests, that the CAL1
% still needed to be applied!

% apply burst power correction to all pulses as suggested in ATBD
%CS1a.FBR.data = bsxfun(@times,CS1a.FBR.data,CS1a.MEA.burst_power_cor(:,1)');
% visualisation shows, that this makes the power alignment over a
% transponder worse again!

% apply agc_ku (in dB, mentioned directly in the netcdf, that it is not yet applied)
% make vector 3 dimensional
agc(1,1,:) = CS1a.SAR.agc_ku';
% and apply it
CS1a.FBR.data = CS1a.FBR.data./(sqrt(10.^(-agc/10)));


end

%Notes:
% >> Level-1A products: geo-located bursts of echoes  with all calibrations
%   applied. << source: https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-altimetry/product-types
%   
% >> The L1A is an intermediate output of the SAR processor. 
%   L1A complex waveforms should be fully calibrated (including both 
%   instrumental gains and calibration corrections) and aligned in range 
%   within each burst. The time tag is given at the surface (that is when 
%   the middle of the burst reaches the surface). L1A is the starting point 
%   for the SAR processing which provides high resolution products.<< source:
%   https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-altimetry/processing-levels
% 
% >>  There are three types of calibrations:
%
%    - calibration due to internal calibration of the instrument (CAL1)
%    - calibration due to the gain profile range window (CAL2)
%    - calibration due to the calibration of the Automatic Gain Control (AGC) correction tables.
%    
%     The main processing steps of the SAR-Ku Level-1 measurement chain are:
%
%    1. Determine surface type.
%    2. Correct tracker ranges for USO frequency drift.
%    3. Compute internal path correction.
%    4. Correct the AGC for instrumental errors.
%    5. Correct waveforms by CAL1 and CAL2.
%    6. Compute surface locations.
%    7. Apply Doppler correction.
%    8. Apply slant range correction.
%    9. Align the waveforms.
%    10. Perform the Doppler beams stack multi-looking.
%    11. Compute Sigma0 scaling factor. <<
%
%       source: https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-altimetry/processing-levels/level-1
%   
%   >>The reference ellipsoid model in the SENTINEL-3 mission is WGS84. The
%   geodetic coordinates (longitude, latitude and altitude) of each point
%   measurement are referenced to WGS84 reference ellipsoid.<<
%   Time and space variable conventions: https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-altimetry/definitions/conventions
%
%
%   Excerpt of the ATBD:
%   2.20 - To apply phase calibration (CAL1) corrections (SAR mode)
%   2.21 - To apply power calibration (CAL1) corrections (SAR mode)
%   2.22 - To corect the waveforms for GPRW (SAR mode)   [GPRW is basically CAL2 correction]
%   2.23 - To correct waveforms from AGC (SAR mode)
%
%
%   >> Level-1A ECHO_SAR_Ku Product: Waveforms;
%   Echo waveform corrected for the Gain Profile Range Window (GPRW)
%   effect<< Source: https://sentinel.esa.int/web/sentinel/level-1a-echo_sar_ku-waveforms
%   
%   This collection suggests, that only the Automatic Gain Control
%   correction needs to be applied before FF-SAR processing.