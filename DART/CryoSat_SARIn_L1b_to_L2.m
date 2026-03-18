function [DATA,CS] = CryoSat_SARIn_L1b_to_L2(FName,DOM,Retracker,SetRetr,DEM,MSK,IDXpoi)

%CRYOSAT_SARIN_L1B_TO_L2 processes level 2 data from CryoSat level 1b SARIn
%data

%% Settings
%General settings
LoadCommonSettings

%Computation settings
defval('FName','2013/05/CS_LTA__SIR_SIN_1B_20130523T130013_20130523T130223_D001.nc') %*.DBL/*.nc file that contains level 1b data
defval('DOM',[59 84;-74 -10])                                                     %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','Threshold')                                                   %Retracker to be used
defval('SetRetr',{})
for i=1:size(SetRetr,1), eval(sprintf('%s = %f;',SetRetr{i,1},SetRetr{i,2})); end
defval('SAMOSAimpl','S+')                                                         %SAMOSA implementation to be used ('S+' = 'SAMOSA+' as described by Salvatore Dinardo et al. 2018; 'DPM_v2_5_2' = Gommenginger et al. 2017)
defval('SolverAlgo','TRR')                                                        %Solver Algorithm to be applied in case Retracker is an analytical retracker ('LM' = Levenberg-Marquardt; 'TRR' = trust-region-reflective)
defval('MaxNrPeaksPerWD',8)                                                       %Max # of peaks to be extracted per waveform
defval('ThresMinPP',0.1);                                                         %MinPeakProminence; threshold for selecting significant peaks (see findpeaks function)
defval('max_waveform_energy',150)                                                 %Max allowed energy of NORMALIZED waveform (i.e., area below the curve)
defval('ApplyWfClass',false)                                                      %Apply surface type classification based on waveform parameters [true/false].
defval('UseLSmskInWfC',false)                                                     %Use high-resolution land/sea mask in waveform type classification [true/false].
defval('MSK',[])                                                                  %Land/sea mask used in waveform classification
defval('minpeakdist',10);                                                         %Minimum distance between peaks. Default is 20.
defval('MINbin',50);                                                              %Minimum bin index of interval in which peaks are detected for NON-zero padded waveform
defval('MAXbin',487);                                                             %Maximum bin index of interval in which peaks are detected for NON-zero padded waveform
defval('ThreshDispersion',0.97);                                                  %Min value R (= some measure of dispersion which = 0 in case dispersion is large and 1 in case dispersion is low) should have (used in case swath processing is applied)
defval('MaxDoN',30E3);                                                            %Max allowed distance-of-nadir
defval('MAfromSTR',true)                                                          %Obtain mispointing angle from star tracker data
defval('Wf_Norm_Aw',1)                                                            %Width (in gates) of the sliding window for waveform normalization
defval('Wf_TN_First',3);                                                          %First gate of the NON-zero padded waveform to estimate the amplitude of the thermal noise floor
defval('Wf_TN_Last',8);                                                           %Last gate of the NON-zero padded waveform to estimate the amplitude of the thermal noise floor
defval('DEMmodels',{'ArcticDEM','DTU13MSL'})                                      %DEM models that will be included in compiling the reference DEM

%Create default optimization options for SolverAlgo with the named
%parameters altered with the specified values
if ~any(strcmp(Retracker,{'ICE1','OCOG','Threshold','Peak','Swath'}))
    switch SolverAlgo
        case 'LM'
            %Levenberg-marquardt
            options = optimoptions('lsqcurvefit','Display','off','Algorithm','levenberg-marquardt','StepTolerance',1E-2,'FunctionTolerance',1E-2,'OptimalityTolerance',1E-2);
        case 'TRR'
            %Trust-region-reflective
            options = optimoptions('lsqcurvefit','Display','off','Algorithm','trust-region-reflective','StepTolerance',1E-2,'FunctionTolerance',1E-2,'OptimalityTolerance',1E-2);
        otherwise
            error('Solver not reckognized!')
    end
end

%Remaining settings
defval('MakePlots',exist('IDXpoi','var')) %Make plots
if MakePlots; if ~isscalar(IDXpoi); MakePlots = false; end; end
defval('FntSz',14)                        %Set fontsize to be used in figures

%% Load reference DEM (used in geolocation of scatterers)
if ~exist('DEM','var')
    DEM = DARTutils.CompileReferenceDEM(DEMmodels,DOM);
end

%% Read & Crop CryoSat Level 1b SARIn data
%Read data
if isstruct(FName)
    CS = FName;
else
    [~,~,Fext] = fileparts(FName);
    if isequal(Fext,'.DBL')
        [~,CS] = Cryo_L1b_read(fullfile(PathDATA,'RadAlt','CryoSat','SIR_SIN_L1',FName));

        %Apply scaling
        CS.SIN.echo_scaling     = CS.SIN.echo_scaling.*1E-9;
        CS.SIN.phase_difference = CS.SIN.phase_difference*1E-6;
        CS.SIN.coherence        = CS.SIN.coherence*1E-3;
    else
        [~,CS] = Cryo_L1b_read_nc(fullfile(PathDATA,'RadAlt','CryoSat','SIR_SIN_L1',FName));
    end
end

%Physical constants & instrument characteristics
DDAcf = DDA_ConfigFile('CS','SIN',CS.GEO.Baseline_ID);

%Wrap angle in degrees to [-180 180]
CS.GEO.LON = wrapTo180(CS.GEO.LON);
%If segment crosses -180 meridian, wrap angle in degrees to [0 360]
if max(abs(diff(CS.GEO.LON(CS.SIN.N_averaged_echoes(:) > 0)))) > 350
    CS.GEO.LON = wrapTo360(CS.GEO.LON);
end

%Remove field AVG to save RAM
CS = rmfield(CS,'AVG');

%If it does not exists, set field with beam angles
if ~isfield(CS.SIN,'BeamAngle')
    CS.SIN.('BeamAngle') = cell(numel(CS.GEO.LAT),1);
    for i=1:numel(CS.GEO.LAT)
        CS.SIN.BeamAngle{i} = linspace(CS.SIN.beam_param(8,i),CS.SIN.beam_param(9,i),CS.SIN.beam_param(12,i))';
    end
end

%Correct window delay for DORIS USO drift (pp 20 CryoSat Product Handbook)
if isequal(CS.GEO.Baseline_ID,'C')
    CS.MEA.win_delay = CS.MEA.win_delay.*(CS.GEO.USO+1);
end

%Compute reference range
CS.MEA.('ref_range') = 0.5*DDAcf.c*CS.MEA.win_delay(:);

%Transform acquisition time to datenum format
CS.('TIME') = datenum('2000','yyyy') + CS.GEO.TAI.days(:) + CS.GEO.TAI.secs(:)./86400 + CS.GEO.TAI.microsecs(:)./1e6./86400 - DDAcf.Biases.Datation/86400;

%Get Roll, Yaw, and Pitch angles from star tracker data
if (isequal(CS.GEO.Baseline_ID,'C') || isequal(CS.GEO.Baseline_ID,'OWN')) && MAfromSTR
    try
        tSTR          = datetime(CS.GEO.Start_Time./(24.*60.*60) + datenum('01-Jan-2000 00:00:00') - 1,'ConvertFrom','datenum');
        mpSTR         = READ_STR_Mispointing_Angles(tSTR);
        DUM           = interp1(mpSTR.Time,[mpSTR.Pitch mpSTR.Roll mpSTR.Yaw],CS.TIME+DDAcf.Biases.Datation/86400,'spline');
        IDXnan        = isnan(interp1(mpSTR.Time,mpSTR.Roll,CS.TIME+DDAcf.Biases.Datation/86400));
        DUM(IDXnan,:) = NaN;
        CS.GEO.Antenna_Bench_Roll(~isnan(DUM(:,2)))  = DUM(~isnan(DUM(:,2)),2);
        CS.GEO.Antenna_Bench_Roll(isnan(DUM(:,2)))   = CS.GEO.Antenna_Bench_Roll(isnan(DUM(:,2))) + DDAcf.Biases.Roll;
        %Note that yaw angle is defined in different sign convention!
        CS.GEO.Antenna_Bench_Yaw(~isnan(DUM(:,3)))   = -DUM(~isnan(DUM(:,3)),3);
        %Note that pitch angle is defined in different sign convention!
        CS.GEO.Antenna_Bench_Pitch(~isnan(DUM(:,1))) = -DUM(~isnan(DUM(:,1)),1);
    catch exception
        fprintf('%s\n',exception.message)
    end
else
    CS.GEO.Antenna_Bench_Roll = CS.GEO.Antenna_Bench_Roll + DDAcf.Biases.Roll;
end

%Convert Roll, Yaw & Pitch angles to radians
CS.GEO.Antenna_Bench_Roll  = deg2rad(CS.GEO.Antenna_Bench_Roll);
CS.GEO.Antenna_Bench_Yaw   = deg2rad(CS.GEO.Antenna_Bench_Yaw);
CS.GEO.Antenna_Bench_Pitch = deg2rad(CS.GEO.Antenna_Bench_Pitch);

%The power echo	sample values are all scaled to	fit	between	0 and 65535.
%The scaling factors can change	for	each waveform. To convert these back to
%values in Watts the following equation should be used (CryoSat Product
%Handbook - Baseline D 1.1, section 2.4.4):
%Power in Watts	= scaled value * (scale	factor) * 2^scale power
CS.SIN.data = bsxfun(@times,CS.SIN.data,reshape(CS.SIN.echo_scaling .* 2.^(CS.SIN.echo_scale_power),1,size(CS.SIN.data,2),size(CS.SIN.data,3)));

%Identify flagged data
IDXfd = CS.GEO.MCD_FLAG.Block_Degraded(:) == 1 | CS.GEO.MCD_FLAG.Blank_Block(:) == 1 | ...
    CS.GEO.MCD_FLAG.Datation_Degraded(:) == 1 | CS.GEO.MCD_FLAG.Orbit_Propag_Err(:) == 1 | ...
    CS.GEO.MCD_FLAG.Echo_Saturation(:) == 1 | CS.GEO.MCD_FLAG.Other_Echo_Err(:) == 1 | ...
    CS.GEO.MCD_FLAG.Rx1_Err_SARin(:) == 1 | CS.GEO.MCD_FLAG.Rx2_Err_SARin(:) == 1 | ...
    CS.GEO.MCD_FLAG.Wind_Delay_Incon(:) == 1 | CS.GEO.MCD_FLAG.AGC_Incon(:) == 1 | ...
    CS.GEO.MCD_FLAG.TRK_ECHO_Err(:) == 1 | CS.GEO.MCD_FLAG.RX1_ECHO_Err(:) == 1 | ...
    CS.GEO.MCD_FLAG.RX2_ECHO_Err(:) == 1 | CS.GEO.MCD_FLAG.NPM_Incon(:) == 1 | ...
    CS.GEO.MCD_FLAG.Power_Scaling_Err(:) == 1;

%Compute sum of instrumental, propagation and geophysical corrections
%Corrections to be applied in case (i) surface == open oceans or
%semi‐enclosed seas OR (ii) surface == enclosed seas or lakes (SSB corr. is
%still lacking). Results in SSH as defined in Eq. 1 (Dinardo et al. 2018)
CorrST01  = CS.COR.ocean_equilibrium_tide+CS.COR.ocean_longperiod_tide+CS.COR.ocean_loading_tide+...
    CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide+CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.gim_ion+...
    CS.COR.dac;
%Corrections to be applied in case the instantaneous SSHs (SSHi, see Eq. 4
%Dinardo et al. 2018) is needed. Note that the SSB correction still has to
%be applied.
CorrST01i = CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.ocean_loading_tide+CS.COR.gim_ion+...
    CS.COR.solidearth_tide+0.468*CS.COR.geocentric_polar_tide;
%Corrections to be applied in case (i) surface == continental ice OR (ii)
%surface == land
CorrST23  = CS.COR.ocean_loading_tide+CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide+...
    CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.gim_ion;
warning('In case of sea ice, different set of corrections should be applied!')

%Interpolate corrections to 20 Hz
if isequal(CS.GEO.Baseline_ID,'C') || isequal(CS.GEO.Baseline_ID,'OWN')
    CS.COR.time_cor = CS.GEO.TAI.days*86400 + CS.GEO.TAI.secs + CS.GEO.TAI.microsecs/1e6 - DDAcf.Biases.Datation;
    CS.COR.time_cor = CS.COR.time_cor(1,:);
end
TIME_COR         = datenum('2000','yyyy') + CS.COR.time_cor'./86400;
CS.('surf_type') = interp1(TIME_COR,CS.COR.surf_type',CS.TIME,'nearest','extrap');
CorrST01         = interp1(TIME_COR,CorrST01',CS.TIME,'linear','extrap');  CorrST01(IDXfd)  = NaN;
CorrST01i        = interp1(TIME_COR,CorrST01i',CS.TIME,'linear','extrap'); CorrST01i(IDXfd) = NaN;
CorrST23         = interp1(TIME_COR,CorrST23',CS.TIME,'linear','extrap');  CorrST23(IDXfd)  = NaN;

%Determine Zero-Padding Oversampling Factor and Waveform sampling interval.
%If os_ZP ~= 1, adjust reference bin.
DDAcf.Ns = size(CS.SIN.data,1);        %Nr of bins/samples in any waveform
os_ZP    = DDAcf.Ns/DDAcf.Np;            %Zero-Padding Oversampling Factor
RefBin   = (DDAcf.RefBin-1)*os_ZP + 1;
dt       = 1/(os_ZP*DDAcf.B);          %Waveform sampling interval [s]

%Compute the local radii of curvature of the Earth's surface (Maulik Jain,
%Improved sea level determination in the Arctic regions through development
%of tolerant altimetry retracking, Eq. 5.2, pp. 47)
CS.GEO.('Re')  = sqrt(DDAcf.RefEll.SemimajorAxis^2*cosd(geocentricLatitude(CS.GEO.LAT,DDAcf.RefEll.Flattening)).^2 + DDAcf.RefEll.SemiminorAxis^2*sind(geocentricLatitude(CS.GEO.LAT,DDAcf.RefEll.Flattening)).^2);

%Compute slope of orbit
if CS.GEO.LAT(find(~IDXfd,1,'first')) < CS.GEO.LAT(find(~IDXfd,1,'last'))
    track_sign = -1; %ascending track
else
    track_sign = 1;  %descending track
end
CS.GEO.('orbit_slope') = track_sign*((DDAcf.RefEll.SemimajorAxis^2 - DDAcf.RefEll.SemiminorAxis^2)./(2*CS.GEO.Re.^2)).*sin(2*deg2rad(CS.GEO.LAT)) - (-CS.GEO.H_rate./CS.GEO.V.V);

%Normalize each waveform by max value of waveform
NORMfactor = 1./max(movmean(CS.SIN.data,Wf_Norm_Aw,1),[],1);
NORMdata   = NORMfactor .* CS.SIN.data;

%Estimate the normalized thermal noise
CS.MEA.('TN') = mean(NORMdata(Wf_TN_First*os_ZP:Wf_TN_Last*os_ZP,:,:));

%Crop data to area of interest
if ischar(DOM)
    switch DOM
        case {'JakobsHavn','jakobshavn'}
            %Load polygon that outlines the JakobsHavn glacier
            load(fullfile(PathDATA,'Geography','jakobshavn_basin.mat'))
            IDX = inpolygon(CS.GEO.LON(:),CS.GEO.LAT(:),polyLon,polyLat);
        otherwise
            error('DOMain of interest not reckognized')
    end
else
    IDX = ingeoquad(CS.GEO.LAT(:),CS.GEO.LON(:),DOM(1,:),DOM(2,:));
end

%% Data editing
%Exclude flagged data
IDX(IDXfd) = false;

%Identify waveforms for which power == 0 for all entries
IDX(squeeze(all(CS.SIN.data == 0,1))) = false;

%Energy NORMALIZED waveform < max_waveform_energy && power at first bins of
%NORMALIZED waveform should be at noise level
IDX(squeeze(trapz(1:DDAcf.Ns,NORMdata,1) >= max_waveform_energy) | squeeze(any(NORMdata(1:MINbin*os_ZP,:,:) > .1))) = false;

%Exclude first/last waveforms acquired after/before mode switch
disp('To be implemented: Exclude first/last waveforms acquired after/before mode switch!')

%Set entries for which (data == 0 & phase_difference == 0), resulting in
%coherence > 1, to NaN
IDXivEntr = CS.SIN.data == 0 & CS.SIN.phase_difference == 0;

%Return if no data remain
if ~any(IDX)
    DATA = struct('TIME',[],'LAT',[],'LON',[],'HEI',[],'DEM',[],'IDXpoi',[],'kx2pi',[],'Coh',[]);
    return
end

%Select points of interest (by default entire track is retracked)
defval('IDXpoi',find(IDX)')

%% Classify waveforms
if ApplyWfClass
    %Preliminaries
    [CS.('n_ret'),CS.('Pu'),CS.('WFc')] = deal(nan(numel(CS.GEO.LAT),1));
    
    %Apply threshold retracker to obtain retracking point and Pu
    [CS.n_ret(IDXpoi),CS.Pu(IDXpoi)]    = RETRACKER.Threshold_mat(NORMdata(:,IDXpoi));
    
    %Compute range and backscatter coefficient (sigma0 / sigma-naught)
    CS.('range')   = CS.MEA.ref_range + (0.5*DDAcf.c*(CS.n_ret*dt - RefBin*dt)) - DDAcf.Biases.Range;
    CS.Pu          = CS.Pu ./ NORMfactor(:);
    CS.('sigma0')  = Compute_Sigma0(CS,DDAcf);
    CS.Pu          = CS.Pu .* NORMfactor(:);
    
    %Apply land/sea mask to find points wrongly classified as land
    if UseLSmskInWfC
        classIN              = CS.surf_type(IDXpoi);
        LorS                 = MSK(CS.GEO.LON(IDXpoi),CS.GEO.LAT(IDXpoi));
        IDXs                 = LorS == 0 & ismember(classIN,[2,3]);
        classIN(IDXs)        = 0;
        CS.surf_type(IDXpoi) = classIN;
        clear('classIN','LorS','IDXs')
    end
    
    %Classify waveforms
    CS.WFc(IDXpoi) = Classify_Waveforms('CS',NORMdata(:,IDXpoi),CS.sigma0(IDXpoi)',CS.surf_type(IDXpoi)');
else
    CS.WFc(IDXpoi) = CS.surf_type(IDXpoi);
end

%% Retrack waveforms
if ~isequal(Retracker,'Swath')
    %Preliminaries
    [CS.('n_ret'),CS.('Pu'),CS.('SWH')] = deal(nan(numel(CS.GEO.LAT),MaxNrPeaksPerWD));
    [CS.('nu'),CS.('ExitF'),CS.('MF')]  = deal(nan(numel(CS.GEO.LAT),MaxNrPeaksPerWD));
    BinIDs                              = (1:DDAcf.Ns)';
    CS.('WDrecon')                      = nan(size(CS.SIN.data));
    
    %Define/load look-up-tables
    if strncmp(Retracker,'SAMOSA',6)
        if strcmp(SAMOSAimpl,'S+')
            %Generate look-up-tables for fast evaluation of modified Bessel
            %functions of the first kind
            LUT_x    = logspace(-16,4,10000)';
            LUT_B14  = griddedInterpolant(LUT_x,besseli(1/4,LUT_x,1),'spline');
            LUT_Bm14 = griddedInterpolant(LUT_x,besseli(-1/4,LUT_x,1),'spline');
            LUT_B34  = griddedInterpolant(LUT_x,besseli(3/4,LUT_x,1),'spline');
            LUT_Bm34 = griddedInterpolant(LUT_x,besseli(-3/4,LUT_x,1),'spline');
        elseif strcmp(SAMOSAimpl,'DPM_v2_5_2')
            %Load look-up-tables alpha_p parameter and the F0 and F1 terms
            fid      = fopen('LUT_Alpha_P.txt'); LUT_AP = cell2mat(textscan(fid,'%f %f\n','HeaderLines',1)); fclose(fid);
            LUT_AP   = griddedInterpolant(LUT_AP(:,1),LUT_AP(:,2),'linear','nearest');
            fid      = fopen('F0.txt'); LUT_F0 = cell2mat(textscan(fid,'%f %f\n','HeaderLines',1)); fclose(fid);
            LUT_F0   = griddedInterpolant(LUT_F0(:,1),LUT_F0(:,2),'spline','nearest');
            fid      = fopen('F1.txt'); LUT_F1 = cell2mat(textscan(fid,'%f %f\n','HeaderLines',1)); fclose(fid);
            LUT_F1   = griddedInterpolant(LUT_F1(:,1),LUT_F1(:,2),'spline','nearest');
        else
            error('Which SAMOSA implementation should be used?')
        end
    end
    
    %Apply retracking
    for i = IDXpoi
        if ~IDX(i); continue; end
        
        %Copy normalized waveform i to vector WD
        WD = NORMdata(:,i);
        
        %Find MaxNrPeaksPerWD largests peaks in NORMALIZED waveform and return associated bin index
        if MaxNrPeaksPerWD > 1
            [Ypk,Xpk,Wpk] = findpeaks(WD(MINbin*os_ZP:MAXbin*os_ZP),'MinPeakProminence',ThresMinPP,'MINPEAKDISTANCE',minpeakdist,'NPEAKS',MaxNrPeaksPerWD,'SortStr','descend');
            [~,IDXocc]    = sort(Xpk);
            Ypk = Ypk(IDXocc); Xpk = Xpk(IDXocc); Wpk = Wpk(IDXocc);
        else
            [Ypk,Xpk,Wpk] = findpeaks(WD(MINbin*os_ZP:MAXbin*os_ZP),'NPEAKS',MaxNrPeaksPerWD,'SortStr','descend');
        end
        Xpk             = Xpk+(MINbin*os_ZP)-1;
        if isempty(Xpk); continue; end
        NrP             = numel(Xpk);
        
        %Set initial values for SAMOSA retracker
        if strncmp(Retracker,'SAMOSA',6)
            if NrP ~= 1, continue, end
            t0_0 = (Xpk-1 - RefBin)*dt*1E9;
            Pu0  = 1;
            SWH0 = 2;
            nu0  = 10;
        end
        
        %Apply retracking
        switch Retracker
            case 'BetaX'
                %X-parameter Beta-retracker with a linear trailing edge (Martin
                %et al., 1983)
                if NrP == 1
                    [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.Beta5(WD,BinIDs,Ypk,Xpk,'linear');
                else
                    [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.BetaX(WD,BinIDs,Ypk,Xpk,Wpk,'linear');
                end
            case 'BetaX_expTE'
                %X-parameter Beta-retracker with an exponential trailing edge
                %(Deng & Featherstone, 2006)
                if NrP == 1
                    [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.Beta5(WD,BinIDs,Ypk,Xpk,'exponential');
                else
                    [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.BetaX(WD,BinIDs,Ypk,Xpk,Wpk,'exponential');
                end
            case 'Brown'
                %Brown Theoretical Ocean Model (Passaro et al., 2014)
                if NrP ~= 1, continue, end
                [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.Brown(WD,BinIDs,Ypk(1),Xpk(1),CS.GEO.H(i),dt*1E9);
            case 'BrownHayne'
                %Brown-Hayne Theoretical Ocean Model (Gommenginger et al., 2011)
                if NrP ~= 1, continue, end
                [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.BrownHayne(WD,BinIDs,Ypk(1),Xpk(1),CS.GEO.H(i),dt*1E9);
            case 'D2P'
                %D2P (Delay/Doppler Phase-monopulse Radar Altimeter) retracker
                %(Giles et al., 2007)
                [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.D2P(WD,BinIDs,Ypk,Xpk,Wpk);
            case 'FunctionFit'
                %"Function Fit" retracker (Surface Topography Mission (STM)
                %SRAL/MWR L2 Algorithms Definition, Accuracy and Specification
                %[SD-03] [SD-07])
                [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.FunctionFit(WD,BinIDs,Ypk,Xpk,Wpk);
            case {'ICE1','OCOG','Threshold'}
                %Select parts of waveform around identified peaks
                if MaxNrPeaksPerWD > 1
                    sub_n  = Xpk-floor(minpeakdist/2) + [0:minpeakdist];
                    sub_WD = WD(sub_n)';
                else
                    sub_n  = BinIDs';
                    sub_WD = WD;
                end
                if isvector(sub_WD), if isrow(sub_WD), sub_WD = sub_WD'; end; end
                if strcmp(Retracker,'ICE1')
                    %Offset Centre Of Gravity retracker (Bamber, 1994)
                    [CS.n_ret(i,1:NrP),CS.Pu(i,1:NrP)] = RETRACKER.ICE1(sub_WD,sub_n);
                elseif strcmp(Retracker,'OCOG')
                    %Offset Centre Of Gravity retracker (Wingham et al., 1986)
                    [CS.n_ret(i,1:NrP),CS.Pu(i,1:NrP)] = RETRACKER.OCOG(sub_WD,sub_n);
                else
                    %Threshold retracker (Davis, 1997).
                    [CS.n_ret(i,1:NrP),CS.Pu(i,1:NrP)] = RETRACKER.Threshold(WD,sub_n,sub_WD);
                end
            case 'Peak'
                CS.n_ret(i,1:NrP) = Xpk - .2*Wpk;
                CS.Pu(i,1:NrP)    = Ypk;
            case 'SAMOSA2'
                %SAMOSA2 retracker (Ray et al. 2015; Dinardo et al., 2018)
                [x0,lb,ub,XDATA,IDXnr,BeamIDX] = RETRACKER.SAMOSA(WD,BinIDs,DDAcf,CS.MEA.TN(i),Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),CS.GEO.H(i),CS.GEO.V.V(i),CS.GEO.orbit_slope(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),1,0,CS.SIN.BeamAngle{i}');
            case 'SAMOSA2FF'
                %SAMOSA2 retracker (Ray et al. 2015; Dinardo et al., 2018)
                [x0,lb,ub,XDATA,IDXnr,BeamIDX] = RETRACKER.SAMOSA(WD,BinIDs,DDAcf,CS.MEA.TN(i),Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),CS.GEO.H(i),CS.GEO.V.V(i),CS.GEO.orbit_slope(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),1,0,0);
            case 'SAMOSA3'
                %SAMOSA3 retracker (Ray et al. 2015; Dinardo et al., 2018)
                [x0,lb,ub,XDATA,IDXnr,BeamIDX] = RETRACKER.SAMOSA(WD,BinIDs,DDAcf,CS.MEA.TN(i),Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),CS.GEO.H(i),CS.GEO.V.V(i),CS.GEO.orbit_slope(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),0,0,CS.SIN.BeamAngle{i}');
            case 'SAMOSA3FF'
                %SAMOSA3 retracker (Ray et al. 2015; Dinardo et al., 2018)
                [x0,lb,ub,XDATA,IDXnr,BeamIDX] = RETRACKER.SAMOSA(WD,BinIDs,DDAcf,CS.MEA.TN(i),Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),CS.GEO.H(i),CS.GEO.V.V(i),CS.GEO.orbit_slope(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),0,0,0);
            otherwise
                error('Retracker not implemented!')
        end
        
        %Set function handle
        if strncmp(Retracker,'SAMOSA',6) && strcmp(SAMOSAimpl,'S+')
            fun = @(x,n) RETRACKER.SAMOSAfun(x,n,DDAcf,LUT_B14,LUT_Bm14,LUT_B34,LUT_Bm34,BeamIDX);
        elseif strncmp(Retracker,'SAMOSA',6) && strcmp(SAMOSAimpl,'DPM_v2_5_2')
            fun = @(x,n) RETRACKER.SAMOSA_DPM_v2_5_2_fun(x,n,DDAcf,LUT_AP,LUT_F0,LUT_F1,BeamIDX);
        end
        
        %In case of analytical retrackers, solve the non-linear obs. eq.
        if ~any(strcmp(Retracker,{'ICE1','OCOG','Threshold','Peak'}))
            if strncmp(Retracker,'SAMOSA',6)
                %In case of the SAMOSA retracker, the 4 unknown parameters are not
                %solved simultaneously. For ocean waveforms, nu (the inverse of the
                %mean-square slope of the sea surface) is set to 0 and not
                %estimated. In case of lead waveforms, the SWH is set to 0 and not
                %estimated. For land contaminated waveforms, Dinardo et al. (2018)
                %apply a dual step retracking where first Pu, t0, and SWH are
                %estimated. Thereafter, the SWH is set 0 and Pu, t0, and nu are
                %(re-)estimated.
                
                %Solve for parameters in case of ocean waveform
                if any(CS.WFc(i) == [0 1 5 99])
                    %Estimate Pu, t0, and SWH
                    XDATA(end)           = 3;
                    x0_tmp               = x0(1:3); lb_tmp = lb(1:3); ub_tmp = ub(1:3);
                    if strcmp(SolverAlgo,'LM'), lb_tmp = []; ub_tmp = []; end
                    [x,Resnorm,exitflag] = DARTutils.SolveNLLS(fun,x0_tmp,XDATA,WD,lb_tmp,ub_tmp,options);
                    x(4)                 = 0;
                    
                    %Verify reason why algorithm is terminated
                    CS.ExitF(i)  = exitflag;
                    if exitflag <= 0, continue, end
                    
                    %Assess whether or not ocean waveform is a land contaminated
                    %waveform
                    if 100*sqrt(Resnorm/DDAcf.Ns) > 4, CS.WFc(i) = 6; XDATA(end-8) = x(3); end
                end
                
                %Solve for parameters in case of lead or land contaminated waveform
                if any(CS.WFc(i) == [2 3 4 6])
                    if CS.WFc(i) == 4, XDATA(end-8) = 0; end
                    
                    %(Re-)estimate Pu, t0, and nu (inverse of the mean-square slope of
                    %the sea surface).
                    XDATA(end)           = 4;
                    x0_tmp               = x0([1 2 4]); lb_tmp = lb([1 2 4]); ub_tmp = ub([1 2 4]);
                    if strcmp(SolverAlgo,'LM'), lb_tmp = []; ub_tmp = []; end
                    [x,Resnorm,exitflag] = DARTutils.SolveNLLS(fun,x0_tmp,XDATA,WD,lb_tmp,ub_tmp,options);
                    x                    = [x(1:2) XDATA(end-8) x(3)];
                    
                    %Verify again the reason why algorithm is terminated
                    CS.ExitF(i) = exitflag;
                    if exitflag <= 0, continue, end
                end
            else
                if any(isnan(fun(x0,XDATA))), x0 = x0 + abs(randn(size(x0)))*.01; end
                
                if strcmp(SolverAlgo,'LM'), lb = []; ub = []; end
                [x,Resnorm,exitflag] = DARTutils.SolveNLLS(fun,x0,XDATA,WD,lb,ub,options);
                
                %Verify reason why algorithm is terminated
                CS.ExitF(i)  = exitflag;
                if exitflag <= 0 || isequal(x,x0), continue, end
            end
            
            %Copy retracking location [bins] and other estimated parameters to CS
            if any(strcmp(Retracker,{'Brown','BrownHayne'}))
                CS.n_ret(i,1:NrP) = x(IDXnr)./(dt*1E9);
                CS.Pu(i,1:NrP)    = x(IDXnr-1);
            elseif any(strcmp(Retracker,{'SAMOSA2','SAMOSA2FF','SAMOSA3','SAMOSA3FF'}))
                CS.n_ret(i,1:NrP) = (x(IDXnr)/(dt*1E9)) + RefBin;
                CS.Pu(i,1:NrP)    = x(IDXnr-1);
                CS.SWH(i,1:NrP)   = x(IDXnr+1);
                CS.nu(i,1:NrP)    = x(IDXnr+2);
                CS.MF(i,1:NrP)    = 100*sqrt(Resnorm/DDAcf.Ns);
            else
                CS.n_ret(i,1:NrP) = x(IDXnr);
                CS.Pu(i,1:NrP)    = x(IDXnr-1);
            end
            
            %Reconstruct waveform
            CS.WDrecon(:,i) = fun(x,XDATA);
        end
    end
    clear('WD','sub_n','sub_WD','Ypk','Xpk','IDXpks','i','j')
    
    %Data editing
    CS.n_ret(CS.n_ret < MINbin*os_ZP | CS.n_ret > MAXbin*os_ZP) = NaN;
else
    %Set window length and determine total number of waveforms stored in CS
    M1          = 64;
    M2          = M1/2;
    M3          = M2/2;
    NrWaveforms = numel(CS.GEO.LAT);

    %Divide phase difference vectors into overlapping segments with length M1 and an overlap of M1/2
    W       = 1-M3:M2:DDAcf.Ns-M2; W(2:M1,:) = 1; W = cumsum(W,1);
    W       = W + sum(W < 1); W = W - sum(W > DDAcf.Ns);
    PDseg   = CS.SIN.phase_difference(W(:),:);  % [M1*M1/2-1 x NrWaveforms]
    PDseg   = reshape(PDseg,M1,[],NrWaveforms); % [M1 x M2-1 x NrWaveforms]

    %Divide UNWRAPPED phase difference vectors into NON-overlapping segments with length M2
    DUM1    = unwrap(flipud(CS.SIN.phase_difference(1:DDAcf.Ns/2,:,:)),1.9*pi);
    DUM2    = unwrap(CS.SIN.phase_difference(DDAcf.Ns/2:end,:,:),1.9*pi);
    PDseg_u = reshape(cat(1,flipud(DUM1),DUM2(2:end,:,:)),M2,[],NrWaveforms);

    %The approach implemented here is based on Foresta, L., Gourmelen, N.,
    %Weissgerber, F., Nienow, P., Williams, J. J., Shepherd, A., ... &
    %Plummer, S. (2018). Heterogeneous and rapid ice loss over the
    %Patagonian Ice Fields revealed by CryoSat-2 swath radar altimetry.
    %Remote Sensing of Environment, 211, 441-455. They describe it as:
    %
    %Compute slope of phase differces per segment by applying a Fourier
    %transform on the normalized complex coherence. Compared to linear
    %regression, this approach is both more efficient computationally as
    %well as independent on phase wrapping. The Fourier transform enables
    %to test a large number of possible slopes and the one with the highest
    %correlation with the input data is selected. The signal is oversampled
    %to take into account that the slope of CS-2's phase difference can
    %represent non-integer frequencies. Thus, each overlapping section has
    %two possible slopes. A correlation is applied again to the data in
    %each overlapping section, this time using only its two estimated slope
    %values. Sections whose correlation is below a set threshold (for this
    %work, ThreshDispersion) are considered noisy and discarded. Finally,
    %the remaining segments are used to unwrap the phase difference. With
    %this procedure, no smoothing is applied to the phase difference and no
    %threshold is set on the power or coherence
    %
    %In our approach the windows are centered around the NON-overlapping
    %segments. So, we compute one slope per segement. We also skip the part
    %where the ``remaining segments are used...''

    %Apply FFT on normalized complex coherence
    Nfft        = DDAcf.Ns/2; %zero-padding applied
    pdDFT       = fft(exp(1i*PDseg),Nfft);

    %Compute one-sided power spectrum
    pdPS              = abs(pdDFT(1:Nfft/2+1,:,:)/Nfft);
    pdPS(2:end-1,:,:) = 2*pdPS(2:end-1,:,:);
    fPS               = 2*pi*(0:(Nfft/2))/Nfft;

    %Obtain index of max power spectrum & retrieve corresponding frequency
    [~,IDXmax]  = max(pdPS,[],1);
    PDslope     = fPS(IDXmax); % [1 x M2-1 x NrWaveforms]

    %Determine sign of slope for each of the NON-overlapping segments
    PDslope_sgn = sign(sum(diff(PDseg_u)));

    %Reconstruct PD for each NON-overlapping segment
    PDr     = PDslope_sgn.*(-M2/2:M2/2-1)'.*PDslope;

    %Compute R (see eq. 2, Weissgerber and Gourmelen, 2017)
    R       = reshape(repmat(abs(sum(exp(1i*((PDseg_u - mean(PDseg_u)) - (PDr - mean(PDr))))))/M2,M2,1,1),size(CS.SIN.phase_difference));

    %Select BinIDs that meet settings
    % IDbins            = find(CS.SIN.data(MINbin*os_ZP:MAXbin*os_ZP,:,:) >= ThreshPower & CS.SIN.coherence(MINbin*os_ZP:MAXbin*os_ZP,:,:) >= ThreshCoherence);
    % IDbins            = find(CS.SIN.data(MINbin*os_ZP:MAXbin*os_ZP,:,:) >= ThreshPower & CS.SIN.MSdPD(MINbin*os_ZP:MAXbin*os_ZP,:,:) <= ThreshCoherence);
    IDbins            = find(R(MINbin*os_ZP:MAXbin*os_ZP,:,:) >= ThreshDispersion & ~IDXivEntr(MINbin*os_ZP:MAXbin*os_ZP,:,:));

    %Determine the equivalent subscript values corresponding to IDbins into
    %array CS.SIN.data
    [I1,I2,I3]        = ind2sub(size(CS.SIN.data(MINbin*os_ZP:MAXbin*os_ZP,:,:)),IDbins);
    I1                = I1+(MINbin*os_ZP)-1;

    %Determine linear indeces of waveforms in CS.SIN.data
    IDwf              = sub2ind([size(CS.SIN.data,2),size(CS.SIN.data,3)],I2,I3);

    %Determine number of selected bins per waveform and set MaxNrPeaksPerWD
    [NrBpW,~,~]       = histcounts(IDwf,.5:(size(CS.SIN.data,2)*size(CS.SIN.data,3))+.5);
    MaxNrPeaksPerWD   = max(NrBpW);

    %Concatenate vectors with BinIDs per waveform. Note that these have
    %different lengths. Stolen from padcat.m
    SZ                = zeros(numel(CS.GEO.LAT),2);
    SZ(:,2)           = NrBpW;
    SZ(SZ(:,2)~=0,1)  = 1;
    M                 = zeros([MaxNrPeaksPerWD+1,numel(CS.GEO.LAT)]) ;
    M(sub2ind(size(M), SZ(:,2).'+1, 1:numel(CS.GEO.LAT))) = 1 ;
    M                 = cumsum(M(1:end-1,:),1) ; %remove last row
    M(M==1)           = NaN;                     %put the fillers in
    M(M==0)           = I1;                      %put the values in

    %Copy matrix M to CS.n_ret and set entries associated to waveforms not
    %belonging to IDXpoi to NaN
    CS.('n_ret')                                      = M';
    CS.n_ret(~ismember(1:numel(CS.GEO.LAT),IDXpoi),:) = NaN;
    CS.('Pu')                                         = ones(numel(CS.GEO.LAT),1); %Set amplitude of waveform power to 1
    clear('IDbins','I1','I2','I3','IDwf','NrBpW','SZ','M')
end

% %To reproduce ESA's results
% %Apply: [DATAma4,CS] = CryoSat_SARIn_L1b_to_L2([],DOM,'Threshold',{'MaxNrPeaksPerWD',1},DEM);
% [~,CS_ESA_I]               = Cryo_L2I_read(fullfile('/data/Projects/CryoSat2_Greenland/Data/RadAlt/CryoSat/SIR_SIN_L2I/2013/05','CS_LTA__SIR_SINI2__20130501T000133_20130501T000357_C001.DBL'));
% RBin1                      = CS_ESA_I.MEA.Retracked_Range_r1/(DDAcf.c/(4*DDAcf.B))+256;
% CS.n_ret(~isnan(CS.n_ret)) = RBin1(~isnan(CS.n_ret));
% CS.n_ret(find(CS_ESA_I.MEA.MEA_Quality_Flag.Height_err_r1 == 1)) = 512;
% RBin1(CS_ESA_I.MEA.MEA_Quality_Flag.Height_err_r1 == 1) = NaN;
% Histogram(double(CS.n_ret(1:3261)-RBin1),1,500,'Differences (bins)')

%% Compute corrected range
%Compute range Eqn 2.8‐1 CryoSat Product Handbook
CS.('range')      = repmat(CS.MEA.ref_range,1,MaxNrPeaksPerWD) + (0.5*DDAcf.c*(CS.n_ret*dt - RefBin*dt)) - DDAcf.Biases.Range;
CS.('range_SSHi') = CS.range;

%Apply propagation/geophysical corrections
ST                     = repmat(CS.surf_type,1,MaxNrPeaksPerWD);
CorrST01               = repmat(CorrST01,1,MaxNrPeaksPerWD);
CorrST01i              = repmat(CorrST01i,1,MaxNrPeaksPerWD);
CorrST23               = repmat(CorrST23,1,MaxNrPeaksPerWD);
CS.range(ST <= 1)      = CS.range(ST <= 1) + CorrST01(ST <= 1);
CS.range(ST >= 2)      = CS.range(ST >= 2) + CorrST23(ST >= 2);
CS.range_SSHi(ST <= 1) = CS.range_SSHi(ST <= 1) + CorrST01i(ST <= 1);
CS.range_SSHi(ST >= 2) = CS.range_SSHi(ST >= 2) + CorrST23(ST >= 2);
clear('CorrST01','CorrST01i','CorrST23','ST')

%% Interferometric processing
%Compute curvature radius N (required to compute surface height)
N        = repmat(DDAcf.RefEll.SemimajorAxis./sqrt(1-DDAcf.RefEll.Eccentricity^2*sind(CS.GEO.LAT(:)).^2),1,MaxNrPeaksPerWD);

%Compute along-track azimuth and obtain analytical representation of
%azimuth as function of latitude
LAT      = CS.GEO.LAT(CS.GEO.H ~= 0);
LON      = CS.GEO.LON(CS.GEO.H ~= 0);
AlngTrAz = azimuth(LAT(1:end-1),LON(1:end-1),LAT(2:end),LON(2:end),DDAcf.RefEll);
if numel(AlngTrAz) < 40
    [p_Az,~,mu] = polyfit(LAT(1:end-1),AlngTrAz,2);
    if mu(2) == 0, mu(2) = 1; end
elseif numel(AlngTrAz) < 60
    [p_Az,~,mu] = polyfit(LAT(1:end-1),AlngTrAz,3);
else
    [p_Az,~,mu] = polyfit(LAT(1:end-1),AlngTrAz,4);
end
AZIM     = polyval(p_Az,CS.GEO.LAT(:),[],mu);
clear('p_Az','mu','AlngTrAz','LAT','LON')

%Declare kTIMES2pi and arrays where locations/heights will be stored
kTIMES2pi         = (-8:8)*2*pi;
NumCand           = numel(kTIMES2pi);
[LAT,LON,HEI,DoN] = deal(nan([size(CS.n_ret,1)*MaxNrPeaksPerWD,NumCand]));

%Determine (linear) indices
IDXwf_a   = repmat((1:size(CS.n_ret,1))',1,MaxNrPeaksPerWD);
IDXvalid  = ~isnan(CS.n_ret);
IDXbs     = round(CS.n_ret(IDXvalid));
[IDXwf,~] = ind2sub(size(CS.n_ret),find(IDXvalid));
IDXbw     = sub2ind([DDAcf.Ns,size(CS.n_ret,1)],IDXbs,IDXwf);

%Smooth phase difference in complex domain by applying moving average filter
FWidth = 4;
if FWidth > 1
    SmoothPD = angle(filtfilt(ones(FWidth,1)/FWidth,1,exp(sqrt(-1)*CS.SIN.phase_difference)));
else
    SmoothPD = CS.SIN.phase_difference;
end
% %Gaussian
% % sigma         = 2;
% % error('The Gaussian filter is not properly implemented as CS.SIN.phase_difference contains NaNs!')
% % gaussFilter   = exp(-linspace(-10*sigma/2,10*sigma/2,10*sigma).^ 2 / (2*sigma^2));
% % gaussFilter   = gaussFilter / sum (gaussFilter);
% % SmoothPD      = angle(filtfilt(gaussFilter,1,exp(sqrt(-1)*CS.SIN.phase_difference)));

%Compute possible angles of arrival (phase difference might be
%wrapped). Here, the variable roll correction (pp. 17 CryoSat Product
%Handbook) need to be applied as well as the residual roll bias
%correction (Gourmelen et al. 2017, Garcia-Mondejar et al. 2017)
if ~all(CS.GEO.INS_CFG.SIRAL_ID(:) == 0), warning('SIRAL_ID ~= 0, what to do'), keyboard, end
theta    = asin(-((SmoothPD(IDXbw) + repmat(kTIMES2pi,numel(IDXbw),1)) .* (DDAcf.lambda0/(2*pi*DDAcf.ABL))) - CS.GEO.Antenna_Bench_Roll(IDXwf));

% %To reproduce ESA's results
% theta(1) = -1*(CS_ESA_I.MEA.Cross_Track_Angle(i));
% if CS_ESA_I.MEA.MEA_Quality_Flag.Height_err_r1(i) == 1, CS.GEO.BaseLine.X(i) = 0; end

%Compute associated surface heights relative to the reference ellipsoid
Rtmp              = repmat(CS.range(IDXvalid),1,NumCand);
Htmp              = repmat(CS.GEO.H(IDXwf),1,NumCand);
Ntmp              = repmat(N(IDXvalid),1,NumCand);
HN                = Htmp + Ntmp;
tp                = asin( sin(theta).*Rtmp ./ sqrt(Rtmp.^2 + HN.^2 - (2.*Rtmp.*HN.*cos(theta))) );
hP                = Htmp - Rtmp.*cos(theta) + Ntmp.*(1 - cos(tp))./cos(tp);
HEI(IDXvalid,:)   = hP;
% %Compute height using eq. 31 of Nanda's thesis
% HEI(IDXvalid,:) = (Htmp - Rtmp.*cos(theta) + Ntmp.*(1 - cos(Rtmp./Ntmp.*theta))).*cos(Rtmp./Ntmp.*theta);
% %Compute height using eq. 31 of Nanda's thesis
% HEI(IDXvalid,:) = (Htmp - Rtmp.*cos(theta) + 6378137.*(1 - cos(Rtmp./6378137.*theta))).*cos(Rtmp./6378137.*theta);

%Compute position. Note that a positive/negative theta corresponds to a
%signal coming from the right/left hand side of the observer sitting on
%the antenna bench with its dorsal spine aligned with the nadir
%direction and looking forward in the along‐track flight direction.
%Note: "... It is important to note that phase-wrapping can occur when
%the across track offset is great enough and this can have the effect
%of the echo appearing to come from the other side of the ground track"
[LAT(IDXvalid,:),LON(IDXvalid,:)] = reckon(repmat(CS.GEO.LAT(IDXwf),1,NumCand),repmat(CS.GEO.LON(IDXwf),1,NumCand),abs((Rtmp.*sin(theta)).*(Ntmp./(hP+Ntmp))),AZIM(IDXwf)+sign(theta).*90,DDAcf.RefEll);

%Store Distance of Nadir
DoN(IDXvalid,:) = Rtmp.*sin(theta);

clear('N','Rtmp','Htmp','HN','tp','kTIMES2pi')

%Exclude all computed heights for which the abs(Distance of Nadir) > width
%in the across‐track direction (is approximately equal to 15 km since it
%can be computed as 2*h*tan(theta/2) being h the altitude and theta = 1.2
%degrees the antenna beamwidth across‐track)
HEI(abs(DoN) > MaxDoN) = NaN;

%Interpolate DEM to all possible scatterer positions but exclude DEM values
%that are obviously outside the range window
DEMxy   = double(DARTutils.EvalDEMatLATLON(DEM,LAT,LON));
R0      = repmat(CS.MEA.ref_range(IDXwf_a(:)),1,NumCand) - (0.5*DDAcf.c*(RefBin*dt));
Rend    = repmat(CS.MEA.ref_range(IDXwf_a(:)),1,NumCand) + (0.5*DDAcf.c*(RefBin*dt));
XYZsat  = CoordinateTrans('ell2xyz',[CS.GEO.LON(IDXwf_a(:)),CS.GEO.LAT(IDXwf_a(:)),CS.GEO.H(IDXwf_a(:))],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');
XYZsat  = repmat(XYZsat,NumCand,1);
LoLaDEM = [LON(:),LAT(:),DEMxy(:)];
XYZdem  = CoordinateTrans('ell2xyz',LoLaDEM,1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');
DISTsd  = reshape(sqrt(sum((XYZsat-XYZdem).^2,2)),size(CS.n_ret,1)*MaxNrPeaksPerWD,NumCand);
IDXinRW = (DISTsd >= R0-100) & (DISTsd <= Rend+100);
DEMxy(~IDXinRW) = NaN;

%Select candidate height that best fits to reference DEM
[CS.('LAT'),CS.('LON'),CS.('HEI'),CS.('DEM'),CS.('DoN'),CS.('UnW'),CS.('theta')] = deal(nan(size(CS.n_ret)));
[~,IDXmd]        = min(abs(HEI-DEMxy),[],2);
CS.UnW(IDXvalid) = IDXmd(IDXvalid(:))-11;
IDXmd            = sub2ind(size(HEI),(1:size(HEI,1))',IDXmd);
CS.LAT(IDXvalid) = LAT(IDXmd(IDXvalid(:)));
CS.LON(IDXvalid) = LON(IDXmd(IDXvalid(:)));
CS.HEI(IDXvalid) = HEI(IDXmd(IDXvalid(:)));
CS.DEM(IDXvalid) = DEMxy(IDXmd(IDXvalid(:)));
CS.DoN(IDXvalid) = DoN(IDXmd(IDXvalid(:)));
CS.theta(IDXvalid) = rad2deg(theta(sub2ind(size(theta),(1:nnz(IDXvalid))',CS.UnW(IDXvalid)+11)));

%Copy coherence information to CS
CS.('Coh')       = nan(size(CS.n_ret));
CS.Coh(IDXvalid) = CS.SIN.coherence(IDXbw);

%Copy power information to CS
CS.('Pwr')       = nan(size(CS.n_ret));
% CS.Pwr(IDXvalid) = CS.SIN.data(IDXbw);
CS.Pwr(IDXvalid) = NORMdata(IDXbw);

%Copy R (see eq. 2, Weissgerber and Gourmelen, 2017)
if isequal(Retracker,'Swath')
    CS.('R')         = nan(size(CS.n_ret));
    CS.R(IDXvalid)   = R(IDXbw);
end

%% Compute backscatter coefficient (sigma0 / sigma-naught)
CS.Pu = CS.Pu ./ NORMfactor(:);
CS.('sigma0') = Compute_Sigma0(CS,DDAcf);
% CS.Pu = CS.Pu .* NORMfactor(:);

%% Analyze output
if MakePlots
    for i = IDXpoi
        if all(isnan(CS.n_ret(i,:))), continue, end
        
        figure('Position',get(0,'Screensize'));
        
        %Select valid entries
        IDXvalid = ~isnan(CS.n_ret(i,:));
        
        %Plot observed/reconstructed waveform
        WD        = CS.SIN.data(:,IDXpoi);
        WDnorm    = NORMdata(:,i);
        StrYlabel = 'Power (Watts)';
        if ~any(strcmp(Retracker,{'ICE1','OCOG','Threshold','Peak','Swath'})); WD = WDnorm; StrYlabel = 'Normalized power'; end
        subplot(2,3,1),plot(1:DDAcf.Ns,WD,'.-','LineWidth',2),hold on
        if isequal(Retracker,'Swath')
            plot(CS.n_ret(IDXpoi,IDXvalid),WD(CS.n_ret(IDXpoi,IDXvalid)),'kx','LineWidth',2)
        else
            if MaxNrPeaksPerWD > 1
                [~,Xpk,~,~] = findpeaks(NORMdata(MINbin*os_ZP:MAXbin*os_ZP,IDXpoi),'MinPeakProminence',ThresMinPP,'MINPEAKDISTANCE',minpeakdist,'NPEAKS',MaxNrPeaksPerWD,'SortStr','descend');
            else
                [~,Xpk,~,~] = findpeaks(NORMdata(MINbin*os_ZP:MAXbin*os_ZP,IDXpoi),'NPEAKS',MaxNrPeaksPerWD,'SortStr','descend');
            end
            plot(Xpk+(MINbin*os_ZP)-1,WD(Xpk+(MINbin*os_ZP)-1),'kx','LineWidth',2)
        end
        if isfield(CS,'WDrecon')
            if ~all(isnan(CS.WDrecon(:,i)))
                plot(1:DDAcf.Ns,CS.WDrecon(:,i),'g-','LineWidth',2)
                plot(CS.n_ret(i,1:nnz(~isnan(CS.n_ret(i,:)))),interp1(1:DDAcf.Ns,CS.WDrecon(:,i),CS.n_ret(i,1:nnz(~isnan(CS.n_ret(i,:))))),'ro','LineWidth',2)
            else
                plot(CS.n_ret(i,1:nnz(~isnan(CS.n_ret(i,:)))),interp1(1:DDAcf.Ns,WD,CS.n_ret(i,1:nnz(~isnan(CS.n_ret(i,:))))),'ro','LineWidth',2)
            end
        else
            plot(CS.n_ret(i,1:nnz(~isnan(CS.n_ret(i,:)))),interp1(1:DDAcf.Ns,WD,CS.n_ret(i,1:nnz(~isnan(CS.n_ret(i,:))))),'ro','LineWidth',2)
        end
        axis([0 DDAcf.Ns 0 max(WD)])
        title('Observed waveform','fontsize',FntSz)
        legend({'Observed waveform';'Peak locations';'Retracking points'},'box','off')
        xlabel('Bin','fontsize',FntSz)
        ylabel(StrYlabel,'fontsize',FntSz)
        set(gca,'fontsize',FntSz)
        
        %Plot coherence time series
        subplot(2,3,2),plot(1:DDAcf.Ns,CS.SIN.coherence(:,i),'.-','LineWidth',2),hold on
        plot(round(CS.n_ret(i,IDXvalid)),CS.SIN.coherence(round(CS.n_ret(i,IDXvalid)),i),'ko','LineWidth',2,'MarkerSize',8)
        b = axis; axis([0 DDAcf.Ns b(3) b(4)])
        title('Observed coherence','fontsize',FntSz)
        xlabel('Bin','fontsize',FntSz)
        ylabel('Coherence','fontsize',FntSz)
        set(gca,'fontsize',FntSz)
        
        %Plot phase difference time series
        subplot(2,3,3)
        if isequal(Retracker,'Swath'),yyaxis left; end
        plot(1:DDAcf.Ns,CS.SIN.phase_difference(:,i),'-','LineWidth',2),hold on
        %Smooth phase difference in complex domain
        %No smoothing
        % SmoothPD    = CS.SIN.phase_difference(:,i);
        %Moving average
        FWidth        = 4;
        SmoothPD      = movmean(exp(sqrt(-1)*CS.SIN.phase_difference(:,i)),FWidth,1);
        SmoothPD      = movmean(SmoothPD(end:-1:1,:,:),FWidth,1);
        SmoothPD      = angle(SmoothPD(end:-1:1,:,:));
        %Gaussian
        % sigma       = 2;
        % error('The Gaussian filter is not properly implemented as CS.SIN.phase_difference contains NaNs!')
        % gaussFilter = exp(-linspace(-10*sigma/2,10*sigma/2,10*sigma).^ 2 / (2*sigma^2));
        % gaussFilter = gaussFilter / sum (gaussFilter);
        % SmoothPD    = angle(filtfilt(gaussFilter,1,exp(sqrt(-1)*CS.SIN.phase_difference(:,i))));
        plot(1:DDAcf.Ns,SmoothPD,'g-','LineWidth',2)
        %Unwrap phase differences
        DUM1     = unwrap(flipud(SmoothPD(1:DDAcf.Ns/2)),1.9*pi);
        DUM2     = unwrap(SmoothPD(DDAcf.Ns/2:end),1.9*pi);
        SmoothPD = [flipud(DUM1);DUM2(2:end)];
        plot(1:DDAcf.Ns,SmoothPD,'r-','LineWidth',2),hold on
        b = axis; axis([0 DDAcf.Ns b(3) b(4)])
        ylabel('Phase difference','fontsize',FntSz)
        if isequal(Retracker,'Swath')
            plot(reshape(1:DDAcf.Ns,M2,M2),PDr(:,:,i)+mean(PDseg_u(:,:,i)),'--','Color',[0.5,0.5,0.5],'LineWidth',1);
            yyaxis right
            plot(1:DDAcf.Ns,R(:,i),'x');
            ylabel('R [-]'); ylim([0,1]);
        end
        title('Observed phase differences','fontsize',FntSz)
        xlabel('Bin','fontsize',FntSz)
        set(gca,'fontsize',FntSz)
        
        %Plot coastline
        load(fullfile(Pth_mat,'COASTS','GSHHS-coastlines-f-59-84-m74-m10.mat'))
        subplot(2,3,4),plot(LOLA(:,1),LOLA(:,2),'k-'),hold on
        plot(CS.GEO.LON(CS.GEO.LON~=0),CS.GEO.LAT(CS.GEO.LON~=0),'.')
        plot(CS.GEO.LON(i),CS.GEO.LAT(i),'ro','LineWidth',2)
        axis([minmax(LOLA(:,1)) minmax(LOLA(:,2))])
        xlabel('Longitude (deg)','fontsize',FntSz)
        ylabel('Latitude (deg)','fontsize',FntSz)
        set(gca,'fontsize',FntSz)
        
        %Plot DEM around point of interest including ground track
        [x,y]  = polarstereo_fwd(CS.GEO.LAT,CS.GEO.LON,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
        [~,Ix] = min(abs(DEM.x(:)-x(i))); [~,Iy] = min(abs(DEM.y(:)-y(i)));
        MDidx  = ceil(MaxDoN/DEM.GSP);
        subplot(2,3,5),imagesc(DEM.x(Ix+(-MDidx:MDidx))./1000,DEM.y(Iy+(-MDidx:MDidx))./1000,DEM.z((Iy+(-MDidx:MDidx)),(Ix+(-MDidx:MDidx)))),hold on
        axis square, colorbar, set(gca,'YDir','normal')
        set(get(findobj(gcf,'tag','Colorbar'),'ylabel'),'String','meters','FontSize',FntSz);
        %Plot direction of flight
        quiver(x(i)./1000,y(i)./1000,(x(i+1)-x(i))./1000,(y(i+1)-y(i))./1000,20,'Color','k','LineWidth',2,'MaxHeadSize',0.8);
        %Plot across-track points
        [lat_ac,lon_ac] = reckon(CS.GEO.LAT(i),CS.GEO.LON(i),-MaxDoN:DEM.GSP/2:MaxDoN,AZIM(i)+90,DDAcf.RefEll);
        [x_ac,y_ac]     = polarstereo_fwd(lat_ac,lon_ac,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
        plot(x_ac./1000,y_ac./1000,'k.')
        %Plot identified scatterers
        [x_rp,y_rp]     = polarstereo_fwd(CS.LAT(i,:),CS.LON(i,:),DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
        plot(x_rp./1000,y_rp./1000,'o','Color',[.5 .5 .5],'LineWidth',2,'MarkerSize',12)
        title('DEM','fontsize',FntSz)
        xlabel('X (km)','fontsize',FntSz)
        ylabel('Y (km)','fontsize',FntSz)
        set(gca,'fontsize',FntSz)
        
        %Plot DEM in across-track direction
        R0      = CS.MEA.ref_range(i) - (0.5*DDAcf.c*(RefBin*dt));
        Rend    = CS.MEA.ref_range(i) + (0.5*DDAcf.c*(RefBin*dt));
        XYZsat  = CoordinateTrans('ell2xyz',[CS.GEO.LON(i),CS.GEO.LAT(i),CS.GEO.H(i)],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');
        TmpDoN  = (-MaxDoN:DEM.GSP/2:MaxDoN)./1000;
        ShftAT  = linspace(0,300,ceil(300/DEM.GSP)+1)-150;
        subplot(2,3,6),hold on
        for j=1:numel(ShftAT)
            [lat_at,lon_at] = reckon(CS.GEO.LAT(i),CS.GEO.LON(i),ShftAT(j),AZIM(i),DDAcf.RefEll);
            [lat_ac,lon_ac] = reckon(lat_at,lon_at,-MaxDoN:DEM.GSP/2:MaxDoN,AZIM(i)+90,DDAcf.RefEll);
            DEMxy           = double(DARTutils.EvalDEMatLATLON(DEM,lat_ac,lon_ac));
            XYZdem          = CoordinateTrans('ell2xyz',[lon_ac(:),lat_ac(:),DEMxy(:)],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');
            IDXinRW         = (pdist2(XYZsat,XYZdem) >= R0-100) & (pdist2(XYZsat,XYZdem) <= Rend+100);
            plot(TmpDoN,DEMxy,'.-','LineWidth',2,'Color',0.5+abs(ShftAT(j)/max(ShftAT)).*[.25 .25 .25])
            DEMxy(~IDXinRW) = NaN;
            plot(TmpDoN,DEMxy,'k.-','LineWidth',2)
        end
        plot(CS.DoN(i,:)/1000,CS.HEI(i,:),'ro','LineWidth',2)
        legend({'Out of range window +/- 100 m','In range window +/- 100 m'},'location','northwest','box','off')
        grid on
        title('DEM in across-track direction','fontsize',FntSz)
        xlabel('Distance of nadir point (km)','fontsize',FntSz)
        ylabel('DEM (m)','fontsize',FntSz)
        set(gca,'fontsize',FntSz)
        
        %Set background color to white
        set(gcf, 'Color',[1 1 1])
        
        figure, yyaxis left, hold on
        plot(round(CS.n_ret(i,IDXvalid)),CS.HEI(i,IDXvalid),'.-','LineWidth',2,'MarkerSize',8), hold on
        plot(round(CS.n_ret(i,IDXvalid)),CS.DEM(i,IDXvalid),'g.-','LineWidth',2,'MarkerSize',8)
        b = axis; axis([MINbin*os_ZP MAXbin*os_ZP b(3) b(4)])
        ylabel('Elevation','fontsize',FntSz)
        yyaxis right
        plot(round(CS.n_ret(i,IDXvalid)),CS.HEI(i,IDXvalid)-CS.DEM(i,IDXvalid),'r.-','LineWidth',2,'MarkerSize',8), hold on
        legend({'Retracked heights','DEM heights','Difference'},'location','best','box','off')
        title('Retracked height versus DEM as a function of range bin','fontsize',FntSz)
        xlabel('Bin','fontsize',FntSz)
        ylabel('Difference','fontsize',FntSz)
        set(gca,'fontsize',FntSz)
    end
end

%% Save output
DATA            = struct;
IDX             = ~isnan(CS.HEI);
CS.TIME         = repmat(CS.TIME,1,MaxNrPeaksPerWD);
DATA.('TIME')   = CS.TIME(IDX);
DATA.('LAT')    = CS.LAT(IDX);
DATA.('LON')    = CS.LON(IDX);
DATA.('HEI')    = single(CS.HEI(IDX));
DATA.('DEM')    = single(CS.DEM(IDX));
IDXpoi          = repmat((1:numel(CS.GEO.LAT))',1,MaxNrPeaksPerWD);
DATA.('IDXpoi') = int32(IDXpoi(IDX));
DATA.('kx2pi')  = int32(CS.UnW(IDX));
DATA.('Coh')    = CS.Coh(IDX);
DATA.('DoN')    = CS.DoN(IDX);
DATA.('Pwr')    = CS.Pwr(IDX);
if isequal(Retracker,'Swath'), DATA.('R') = CS.R(IDX); end
DATA.('theta')  = CS.theta(IDX);
DATA.('sigma0') = single(CS.sigma0(IDX));

end
