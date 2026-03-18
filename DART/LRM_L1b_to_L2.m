function [DATA,CS] = LRM_L1b_to_L2(SAT,FName,DOM,Retracker,SetRetr,DEM,MSK,IDXpoi)

%LRM_L1B_TO_L2 processes level 2 data from CryoSat level 1b LRM data.

%% Settings
%General settings
LoadCommonSettings

%Computation settings
defval('SAT','CS')                                                                %Satellite mission from which data are processed ('CS')
defval('FName','2018/12/CS_OFFL_SIR_LRM_1B_20181201T194252_20181201T194605_C001.DBL') %*.DBL/*.nc file that contains level 1b data
defval('DOM',[59 84;-74 -10])                                                     %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','Threshold')                                                   %Retracker to be used
defval('SetRetr',{})
for i=1:size(SetRetr,1), eval(sprintf('%s = %s;',SetRetr{i,1},SetRetr{i,2})); end
defval('SolverAlgo','TRR')                                                        %Solver Algorithm to be applied in case Retracker is an analytical retracker ('LM' = Levenberg-Marquardt; 'TRR' = trust-region-reflective)
defval('MaxNrPeaksPerWD',1)                                                       %Max # of peaks to be extracted per waveform
defval('ThresMinPP',0.1);                                                         %MinPeakProminence; threshold for selecting significant peaks (see findpeaks function)
defval('minpeakdist',10);                                                         %Minimum distance between peaks.
defval('max_waveform_energy',150)                                                 %Maximum allowed energy of NORMALIZED waveform (i.e., area below the curve)
defval('ApplyWfClass',true)                                                       %Apply surface type classification based on waveform parameters [true/false].
defval('MAfromSTR',true)                                                          %Obtain mispointing angle from star tracker data
defval('MINbin',10);                                                              %Minimum bin index of interval in which retracking point is accepted for NON-zero padded waveform
defval('MAXbin',125);                                                             %Maximum bin index of interval in which retracking point is accepted for NON-zero padded waveform
defval('Wf_Norm_Aw',1)                                                            %Width (in gates) of the sliding window for waveform normalization
defval('Wf_TN_First',8);                                                          %First gate of the NON-zero padded waveform to estimate the amplitude of the thermal noise floor
defval('Wf_TN_Last',15);                                                          %Last gate of the NON-zero padded waveform to estimate the amplitude of the thermal noise floor
defval('ApplySubWFretr',false);                                                   %Apply sub-waveform retracking
defval('UseDEM2CompInitRP',false)                                                 %Use DEM to compute initial retracking point [true/false].
defval('DEM',[])                                                                  %DEM used to i) compute initial retracking point, ii) apply slope correction
defval('UseLSmskInWfC',false)                                                     %Use high-resolution land/sea mask in waveform type classification [true/false].
defval('MSK',[])                                                                  %Land/sea mask used in waveform classification
defval('ApplySlopeCorrection',true)                                               %Apply slope correction
defval('SCmthd','Roemer')                                                         %Slope correction method ['Bamber'/'Roemer'/'Slobbe'/'Slobbe2']
if (UseDEM2CompInitRP || ApplySlopeCorrection) && isempty(DEM), error('Provide DEM as input!'), end
if UseLSmskInWfC && isempty(MSK), error('Provide land/sea mask as input!'), end

%Create default optimization options for SolverAlgo with the named
%parameters altered with the specified values
if ~any(strcmp(Retracker,{'ICE1','ICE3','OCOG','Threshold'}))
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

%Physical constants & instrument characteristics
DDAcf   = DDA_ConfigFile(SAT,'LRM');

%Remaining settings
defval('MakePlots',exist('IDXpoi','var')) %Make plots
if MakePlots; if ~isscalar(IDXpoi); MakePlots = false; end; end
defval('FntSz',14)                        %Set fontsize to be used in figures

%% Read & Crop Level 1b LRM data
switch SAT
    case {'CS','CryoSat'}
        %Read data
        if isstruct(FName)
            CS     = FName;
        else
            [~,~,Fext] = fileparts(FName);
            if isequal(Fext,'.DBL')
                [~,CS] = Cryo_L1b_read(fullfile(PathDATA,'RadAlt','CryoSat','SIR_LRM_L1',FName));
            else
                [~,CS] = Cryo_L1b_read_nc(fullfile(PathDATA,'RadAlt','CryoSat','SIR_LRM_L1',FName));
            end
        end

        %Regenerate DDAcf - update Baseline_ID dependent biases
        DDAcf = DDA_ConfigFile(SAT,'SAR',CS.GEO.Baseline_ID);

        %Wrap angle in degrees to [-180 180]
        CS.GEO.LON = wrapTo180(CS.GEO.LON);
        %If segment crosses -180 meridian, wrap angle in degrees to [0 360]
        if max(abs(diff(CS.GEO.LON(CS.LRM.N_averaged_echoes(:) > 0)))) > 350
            CS.GEO.LON = wrapTo360(CS.GEO.LON);
        end
        
        %Remove field AVG to save RAM
        CS = rmfield(CS,'AVG');
        
        %Correct window delay for DORIS USO drift (pp 20 CryoSat Product Handbook)
        if isequal(CS.GEO.Baseline_ID,'C')
            CS.MEA.win_delay = CS.MEA.win_delay.*(CS.GEO.USO+1);
        end
        
        %Compute reference range
        CS.MEA.('ref_range') = 0.5*DDAcf.c*CS.MEA.win_delay(:);
        
        %Transform acquisition time to datenum format
        CS.('TIME') = datenum('2000','yyyy') + CS.GEO.TAI.days(:) + CS.GEO.TAI.secs(:)./86400 + CS.GEO.TAI.microsecs(:)./1e6./86400 - DDAcf.Biases.Datation/86400;
        
        %Get Roll, Yaw, and Pitch angles from star tracker data
        if isequal(CS.GEO.Baseline_ID,'C') && MAfromSTR
            try
                tSTR          = datetime(CS.GEO.Start_Time./(24.*60.*60) + datenum('01-Jan-2000 00:00:00') - 1,'ConvertFrom','datenum');
                mpSTR         = READ_STR_Mispointing_Angles(tSTR);
                DUM           = interp1(mpSTR.Time,[mpSTR.Pitch mpSTR.Roll mpSTR.Yaw],CS.TIME+DDAcf.Biases.Datation/86400,'spline');
                IDXnan        = isnan(interp1(mpSTR.Time,mpSTR.Roll,CS.TIME+DDAcf.Biases.Datation/86400));
                DUM(IDXnan,:) = NaN;
                CS.GEO.Antenna_Bench_Roll(~isnan(DUM(:,2)))  = DUM(~isnan(DUM(:,2)),2);
                %Note that yaw angle is defined in different sign convention!
                CS.GEO.Antenna_Bench_Yaw(~isnan(DUM(:,3)))   = -DUM(~isnan(DUM(:,3)),3);
                %Note that pitch angle is defined in different sign convention!
                CS.GEO.Antenna_Bench_Pitch(~isnan(DUM(:,1))) = -DUM(~isnan(DUM(:,1)),1);
            catch exception
                fprintf('%s\n',exception.message)
            end
        end
        
        % %Salvatore Dinardo applies the pitch and roll biases proposed by Remko
        % %Scharroo. In doing so, he undoes the biases applied in producing baseline
        % %C (Main evolutions and expected quality improvements in Baseline C Level
        % %1b products - ARESYS / ESA, V1.3, https://wiki.services.eoportal.org/tiki-
        % %index.php?page=CryoSat+Technical+Notes)
        % CS.GEO.Antenna_Bench_Roll  = CS.GEO.Antenna_Bench_Roll  + 0.1062 - 0.086;
        % CS.GEO.Antenna_Bench_Pitch = CS.GEO.Antenna_Bench_Pitch + 0.0550 - 0.096;
        
        %Convert Roll, Yaw & Pitch angles to radians
        CS.GEO.Antenna_Bench_Roll  = deg2rad(CS.GEO.Antenna_Bench_Roll);
        CS.GEO.Antenna_Bench_Yaw   = deg2rad(CS.GEO.Antenna_Bench_Yaw);
        CS.GEO.Antenna_Bench_Pitch = deg2rad(CS.GEO.Antenna_Bench_Pitch);
        
        %The power echo	sample values are all scaled to	fit	between	0 and 65535.
        %The scaling factors can change	for	each waveform. To convert these back to
        %values in Watts the following equation should be used (CryoSat Product
        %Handbook, Eqn 4.2‐1):
        %Power in Watts	= scaled value * (scale	factor * 10^‐9) * 2^scale power
        CS.LRM.data = bsxfun(@times,CS.LRM.data,reshape((CS.LRM.echo_scaling.*1E-9) .* 2.^(CS.LRM.echo_scale_power),1,size(CS.LRM.data,2),size(CS.LRM.data,3)));
        
        %Identify flagged data
        IDXfd = CS.GEO.MCD_FLAG.Block_Degraded(:) == 1 | CS.GEO.MCD_FLAG.Blank_Block(:) == 1 | ...
            CS.GEO.MCD_FLAG.Datation_Degraded(:) == 1 | CS.GEO.MCD_FLAG.Orbit_Propag_Err(:) == 1 | ...
            CS.GEO.MCD_FLAG.Echo_Saturation(:) == 1 | CS.GEO.MCD_FLAG.Other_Echo_Err(:) == 1 | ...
            CS.GEO.MCD_FLAG.Rx1_Err_SARin(:) == 1 | CS.GEO.MCD_FLAG.Rx2_Err_SARin(:) == 1 | ...
            CS.GEO.MCD_FLAG.Wind_Delay_Incon(:) == 1 | CS.GEO.MCD_FLAG.AGC_Incon(:) == 1 | ...
            CS.GEO.MCD_FLAG.TRK_ECHO_Err(:) == 1 | CS.GEO.MCD_FLAG.RX1_ECHO_Err(:) == 1 | ...
            CS.GEO.MCD_FLAG.RX2_ECHO_Err(:) == 1 | CS.GEO.MCD_FLAG.NPM_Incon(:) == 1 | ...
            CS.GEO.MCD_FLAG.Power_Scaling_Err(:) == 1 | ~strcmp(cellstr(CS.LRM.FLAG(:,:)'),'0000000000000000');
        
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
        if isequal(CS.GEO.Baseline_ID,'C')
            CS.COR.time_cor = CS.GEO.TAI.days*86400 + CS.GEO.TAI.secs + CS.GEO.TAI.microsecs/1e6 - DDAcf.Biases.Datation;
            CS.COR.time_cor = CS.COR.time_cor(1,:);
        end
        TIME_COR         = datenum('2000','yyyy') + CS.COR.time_cor'./86400;
        CS.('surf_type') = interp1(TIME_COR,CS.COR.surf_type',CS.TIME,'nearest','extrap');
        CorrST01         = interp1(TIME_COR,CorrST01',CS.TIME,'linear','extrap');  CorrST01(IDXfd)  = NaN;
        CorrST01i        = interp1(TIME_COR,CorrST01i',CS.TIME,'linear','extrap'); CorrST01i(IDXfd) = NaN;
        CorrST23         = interp1(TIME_COR,CorrST23',CS.TIME,'linear','extrap');  CorrST23(IDXfd)  = NaN;
        FlgC01           = isnan(CorrST01); FlgC01i = isnan(CorrST01i); FlgC23 = isnan(CorrST23);

    otherwise
        error('SAT: %s not recognized',SAT)
end

%Determine Zero-Padding Oversampling Factor and Waveform sampling interval.
%If os_ZP ~= 1, adjust reference bin.
DDAcf.Ns = size(CS.LRM.data,1);        %Nr of bins/samples in any waveform
os_ZP    = DDAcf.Ns/DDAcf.Np;          %Zero-Padding Oversampling Factor
RefBin   = (DDAcf.RefBin-1)*os_ZP + 1;
dt       = 1/(os_ZP*DDAcf.B);          %Waveform sampling interval [s]

%Compute mispointing angles
CS.GEO.('MPA') = acos(cos(CS.GEO.Antenna_Bench_Pitch)./sqrt(sum(cat(3,sin(CS.GEO.Antenna_Bench_Roll),sin(CS.GEO.Antenna_Bench_Yaw),cos(CS.GEO.Antenna_Bench_Pitch)).^2,3)));

%Normalize each waveform by max value of waveform
NORMfactor = 1./max(movmean(CS.LRM.data,Wf_Norm_Aw,1),[],1);
NORMdata   = NORMfactor .* CS.LRM.data;

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

%If a slope correction is applied, crop data to DOM covered by DEM (i.e.,
%DOM covered by DEM might be smaller)
if ApplySlopeCorrection
    [x_sat,y_sat]      = polarstereo_fwd(CS.GEO.LAT,CS.GEO.LON,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
    InDOM_DEM          = (x_sat > DEM.x(1)+10E3 & x_sat < DEM.x(end)-10E3) & (y_sat < DEM.y(1)-10E3 & y_sat > DEM.y(end)+10E3);
    IDX(~InDOM_DEM(:)) = false;
end

%% Data editing
%Exclude flagged data
IDX(IDXfd) = false;

%Identify waveforms for which power == 0 for all entries
IDX(squeeze(all(CS.LRM.data == 0,1))) = false;

%Energy NORMALIZED waveform < max_waveform_energy && power at first bins of
%NORMALIZED waveform should be at noise level
IDX(squeeze(trapz(1:DDAcf.Ns,NORMdata,1) >= max_waveform_energy) | squeeze(any(NORMdata(1:MINbin*os_ZP,:,:) > .2))) = false;

%Return if no data remain
if ~any(IDX)
    DATA = struct('TIME',[],'LAT',[],'LON',[],'HEI',[],'SurfT',[]);
    return
end

%Select points of interest (by default entire track is retracked)
defval('IDXpoi',find(IDX)')

%% Classify waveforms
if ApplyWfClass
%     %Preliminaries
%     [CS.('n_ret'),CS.('Pu'),CS.('WFc')] = deal(nan(numel(CS.GEO.LAT),1));
%     
%     %Apply threshold retracker to obtain retracking point and Pu
%     [CS.n_ret(IDXpoi),CS.Pu(IDXpoi)]    = RETRACKER.Threshold_mat(NORMdata(:,IDXpoi));
%     
%     %Compute range and backscatter coefficient (sigma0 / sigma-naught)
%     CS.('range')   = CS.MEA.ref_range + (0.5*DDAcf.c*(CS.n_ret*dt - RefBin*dt)) - DDAcf.Biases.Range;
%     switch SAT
%         case {'CS','CryoSat'}
%             CS.Pu          = CS.Pu ./ NORMfactor(:);
%             CS.('sigma0')  = Compute_Sigma0(CS,DDAcf);
%             CS.Pu          = CS.Pu .* NORMfactor(:);
%         case {'S3A','Sentinel-3A','S3B','Sentinel-3B'}
%             %Surface Topography Mission (STM) SRAL/MWR L2 Algorithms
%             %Definition, Accuracy and Specification, Section 2.15.3.3
%             CS.('sigma0')  = CS.LRM.scale_factor_ku_l1b_echo_sar_ku + 10*log10(CS.Pu) + DDAcf.Biases.sigma0;
%         otherwise
%             error('SAT: %s not recognized',SAT)
%     end
%     
%     %Apply land/sea mask to find points wrongly classified as land
%     if UseLSmskInWfC
%         classIN              = CS.surf_type(IDXpoi);
%         LorS                 = MSK(CS.GEO.LON(IDXpoi),CS.GEO.LAT(IDXpoi));
%         IDXs                 = LorS == 0 & ismember(classIN,[2,3]);
%         classIN(IDXs)        = 0;
%         CS.surf_type(IDXpoi) = classIN;
%         clear('classIN','LorS','IDXs')
%     end
%     
%     %Classify waveforms
%     CS.WFc(IDXpoi) = Classify_Waveforms(SAT,NORMdata(:,IDXpoi),CS.sigma0(IDXpoi)',CS.surf_type(IDXpoi)');
% else
    CS.('WFc')     = deal(nan(numel(CS.GEO.LAT),1));
    CS.WFc(IDXpoi) = CS.surf_type(IDXpoi);
end
    
%% Preliminaries
%Declare arrays
[CS.('n_ret'),CS.('Pu'),CS.('SWH')] = deal(nan(numel(CS.GEO.LAT),1));
[CS.('nu'),CS.('ExitF'),CS.('MF')]  = deal(nan(numel(CS.GEO.LAT),1));
CS.('PCorr')                        = nan(numel(CS.GEO.LAT),1);
BinIDs                              = (1:DDAcf.Ns)';
CS.('WDrecon')                      = nan(size(CS.LRM.data));

%Use DEM to compute initial retracking point
Xpk_ALL = nan(size(CS.GEO.H(:)));
if UseDEM2CompInitRP
    Xpk_ALL = ((CS.GEO.H(:) - DEM(CS.GEO.LON(:),CS.GEO.LAT(:))) - CorrST01(:) + DDAcf.Biases.Range - CS.MEA.ref_range(:) + (0.5*DDAcf.c*RefBin*dt))/(0.5*DDAcf.c*dt);
    Xpk_ALL(Xpk_ALL < 1 | Xpk_ALL > DDAcf.Ns) = NaN;
end

% profile on -detail builtin -history

%% Retrack waveforms
for i = IDXpoi
    %Copy normalized waveform i to vector WD
    WD        = NORMdata(:,i);
    
    if ~IDX(i) || any(isnan(WD)); continue; end

    %Find MaxNrPeaksPerWD largests peaks in NORMALIZED waveform and return associated bin index
    [Ypk,Xpk,Wpk] = findpeaks(movmean(WD(MINbin*os_ZP:MAXbin*os_ZP),5),'MinPeakProminence',ThresMinPP,'MINPEAKDISTANCE',minpeakdist,'NPEAKS',MaxNrPeaksPerWD,'SortStr','descend');
    [~,IDXocc]    = sort(Xpk);
    Ypk           = Ypk(IDXocc); Xpk = Xpk(IDXocc); Wpk = Wpk(IDXocc);
    Xpk           = Xpk+(MINbin*os_ZP)-1;
    if isempty(Xpk); continue; end
    if ~isnan(Xpk_ALL(i)), Xpk = Xpk_ALL(i); end
    
    %Apply retracking
    switch Retracker
        case 'BetaX'
            %X-parameter Beta-retracker with a linear trailing edge (Martin
            %et al., 1983)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.BetaX(WD,BinIDs,Ypk,Xpk,Wpk,'linear');
        case 'BetaX_expTE'
            %X-parameter Beta-retracker with an exponential trailing edge
            %(Deng & Featherstone, 2006)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.BetaX(WD,BinIDs,Ypk,Xpk,Wpk,'exponential');
        case 'Brown'
            %Brown Theoretical Ocean Model (Passaro et al., 2014)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.Brown(WD,BinIDs,Ypk,Xpk,CS.GEO.H(i),dt*1E9);
        case 'BrownHayne'
            %Brown-Hayne Theoretical Ocean Model (Gommenginger et al., 2011)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.BrownHayne(WD,BinIDs,Ypk,Xpk,CS.GEO.H(i),dt*1E9);
        case 'D2P'
            %D2P (Delay/Doppler Phase-monopulse Radar Altimeter) retracker
            %(Giles et al., 2007)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.D2P(WD,BinIDs,Ypk,Xpk,Wpk);
        case 'FunctionFit'
            %"Function Fit" retracker (Surface Topography Mission (STM)
            %SRAL/MWR L2 Algorithms Definition, Accuracy and Specification
            %[SD-03] [SD-07])
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.FunctionFit(WD,BinIDs,Ypk,Xpk,Wpk,CS.MEA.TN(i));
        case 'ICE1'
            %Offset Centre Of Gravity retracker (Bamber, 1994)
            [CS.n_ret(i),CS.Pu(i)] = RETRACKER.ICE1(WD(MINbin*os_ZP:MAXbin*os_ZP),BinIDs(MINbin*os_ZP:MAXbin*os_ZP)');
        case 'ICE2'
            %ICE-2 retracker (Surface Topography Mission (STM)
            %SRAL/MWR L2 Algorithms Definition, Accuracy and Specification
            %[SD-03] [SD-07])
            [fun,x0,lb,ub,XDATA,IDXnr,IDXew] = RETRACKER.ICE2(WD,BinIDs,Ypk,Xpk,CS.MEA.TN(i));
        case 'ICE3'
            %ICE-3 retracker (Coastal and Hydrology Altimetry product
            %(PISTACH) handbook)
            [CS.n_ret(i),CS.Pu(i)] = RETRACKER.ICE3(WD,BinIDs,Xpk,CS.MEA.TN(i));
        case 'LIRT'
            %LIRT retracker (D. J. Brockley (2019). CryoSat2 : L2 Design Summary Document)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.LIRT(WD,BinIDs,Ypk,Xpk,CS.MEA.TN(i),CS.GEO.MPA(i),CS.GEO.H(i),DDAcf.theta_y,DDAcf.theta_x,DDAcf.B,dt*1E9);
        case 'OCOG'
            %Offset Centre Of Gravity retracker (Wingham et al., 1986)
            [CS.n_ret(i),CS.Pu(i)] = RETRACKER.OCOG(WD(MINbin*os_ZP:MAXbin*os_ZP),BinIDs(MINbin*os_ZP:MAXbin*os_ZP)');
        case 'Threshold'
            %Threshold retracker (Davis, 1997).
            [CS.n_ret(i),CS.Pu(i)] = RETRACKER.Threshold(WD,BinIDs',WD);
        otherwise
            error('Retracker not implemented!')
    end

    %In case of analytical retrackers, solve the non-linear obs. eq.
    if ~any(strcmp(Retracker,{'ICE1','ICE3','OCOG','Threshold'}))

        %Apply sub-waveform retracking
        XDATAorg = XDATA;
        if ApplySubWFretr
            XDATA = [XDATA(IDXew);XDATA(DDAcf.Ns+1:end)];
            WD    = WD(IDXew);
        end
        
        %Solve for parameters
        x0_tmp               = x0; lb_tmp = lb; ub_tmp = ub;
        if strcmp(SolverAlgo,'LM'), lb_tmp = []; ub_tmp = []; end
        [x,Resnorm,exitflag] = DARTutils.SolveNLLS(fun,x0_tmp,XDATA,WD,lb_tmp,ub_tmp,options);
        
        %Verify again the reason why algorithm is terminated
        CS.ExitF(i) = exitflag;
        if exitflag <= 0, continue, end
        
        %Assess quality of fit
        CS.MF(i)    = 100*sqrt(Resnorm/DDAcf.Ns);
        CS.PCorr(i) = corr(WD,fun(x,XDATA),'type','Pearson');
               
        %Copy retracking location [bins] and other estimated parameters to CS
        if any(strcmp(Retracker,{'Brown','BrownHayne','LIRT'}))
            CS.n_ret(i) = x(IDXnr)/(dt*1E9);
            CS.Pu(i)    = x(IDXnr-1);
        else
            CS.n_ret(i) = x(IDXnr(1));
            CS.Pu(i)    = x(IDXnr(1)-1);
        end
        
        %Reconstruct waveform
        CS.WDrecon(:,i) = fun(x,XDATAorg);
    end
end
clear('WD','sub_n','sub_WD','Ypk','Xpk','IDXpks','i','j')

% profile viewer
% profile off

%% Data editing
CS.n_ret(CS.n_ret < MINbin*os_ZP | CS.n_ret > MAXbin*os_ZP) = NaN;

%% Compute corrected range (corrections include instrumental, propagation and geophysical corrections)
%Compute range Eqn 2.8‐1 CryoSat Product Handbook
CS.('range')      = CS.MEA.ref_range + (0.5*DDAcf.c*(CS.n_ret*dt - RefBin*dt)) - DDAcf.Biases.Range;
CS.('range_SSHi') = CS.range;

%Apply propagation/geophysical corrections
ST                     = CS.surf_type;
CS.('SumCorrST')       = nan(size(CorrST01)); CS.('FlgC')       = true(size(CorrST01));
CS.SumCorrST(ST <= 1)  = CorrST01(ST <= 1);   CS.FlgC(ST <= 1)  = FlgC01(ST <= 1);
CS.SumCorrST(ST >= 2)  = CorrST23(ST >= 2);   CS.FlgC(ST >= 2)  = FlgC23(ST >= 2);
CS.('SumCorrSTi')      = CS.SumCorrST;        CS.('FlgCi')      = CS.FlgC;
CS.SumCorrSTi(ST <= 1) = CorrST01i(ST <= 1);  CS.FlgCi(ST <= 1) = FlgC01i(ST <= 1);
CS.range               = CS.range + CS.SumCorrST;
CS.range_SSHi          = CS.range_SSHi + CS.SumCorrSTi;
clear('CorrST01','CorrST01i','CorrST23','ST','FlgC01','FlgC01i','FlgC23')

%Compute (instantaneous) (sea) surface heights relative to the reference ellipsoid
CS.('HEI')  = CS.GEO.H(:) - CS.range;
CS.('SSHi') = CS.GEO.H(:) - CS.range_SSHi;

%% Compute & Apply slope correction
if ApplySlopeCorrection
    %Set min/max range used to determine the part of the terrain visible by
    %the radar
    CS.('R0')   = CS.MEA.ref_range - (0.5*DDAcf.c*(RefBin*dt)) + CS.SumCorrST;
    CS.('Rend') = CS.MEA.ref_range + (0.5*DDAcf.c*(RefBin*dt)) + CS.SumCorrST;

    %Apply `slope' correction
    if strcmp(SCmthd,'Roemer')
        %Determine 3D position of approximate impact point
        [CS.('LATc'),CS.('LONc'),CS.('hr'),CS.('SI')] = deal(nan(numel(CS.GEO.LAT),1));
        CS          = Compute_Impact_Point_Roemer(CS,IDXpoi,DEM,DDAcf);
        
        %Compute terrain corrected height
        CS.HEI      = CS.HEI + (CS.SI(:) - (CS.GEO.H(:)-CS.hr(:)));
    elseif strcmp(SCmthd,'Slobbe2')
        %Set min/max range used to determine the part of the terrain
        %visible by the radar
        Fdata       = movmean(NORMdata(:,:),5,1);
        [~,IDXmax]  = max(Fdata,[],1);
        IDXmax      = IDXmax';
        IDXmin      = zeros(size(IDXmax));
        for i = IDXpoi
            if all(Fdata(1:IDXmax(i),i) > 0.05)
                IDXmin(i) = 1;
            else
                IDXmin(i) = IDXmax(i) - find(flipud(Fdata(1:IDXmax(i),i)) < 0.05,1,'first') + 1;
            end
        end
        CS.('R0')   = CS.MEA.ref_range + (0.5*DDAcf.c*(IDXmin*dt - RefBin*dt)) - DDAcf.Biases.Range + CS.SumCorrST;
        [n_ret,~]   = RETRACKER.Threshold_mat(NORMdata(:,:));
        CS.('Rend') = min([CS.range + 2.5,CS.MEA.ref_range + (0.5*DDAcf.c*(n_ret'*dt - RefBin*dt)) - DDAcf.Biases.Range + CS.SumCorrST],[],2);
        clear Fdata IDXmax IDXmin n_ret
        
        %Compute updated height and lat/lon coordinates. 
        [CS.('LATc'),CS.('LONc'),CS.('MinD2CDP')] = deal(nan(numel(CS.GEO.LAT),1));
        CS          = Compute_Impact_Point_Slobbe(CS,IDXpoi,DEM,DDAcf);
    else
        %Compute aspect and slope w.r.t the reference ellipsoid. As DEM
        %contains ellipsoidal heights, SLOPE from DARTutils.Slope_AND_Aspect is
        %correct. ASPECT from DARTutils.Slope_AND_Aspect is in the map
        %projection relative to which the DEM is given. Hence, we convert
        %ASPECT to ASPECT on the reference ellipsoid. Also check:
        %https://icesat4.gsfc.nasa.gov/cryo_data/icesat_elevation_slope_grids.php
        [aA,sA,~,~] = DARTutils.Slope_AND_Aspect(double(DEM.z),DEM.GSP);
        aA          = aA + 180; %Note that the aspect angle is the direction of the steepest descent
        [DEMx,DEMy] = meshgrid(DEM.x(2:end-1),DEM.y(2:end-1));
        DEMx2       = DEMx + sind(aA);
        DEMy2       = DEMy + cosd(aA);
        [la1,lo1]   = polarstereo_inv(DEMx,DEMy,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
        [la2,lo2]   = polarstereo_inv(DEMx2,DEMy2,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
        aA          = azimuth(la1,lo1,la2,lo2,DDAcf.RefEll);
        
        %Interpolate aspect and slope to sub-satellite points (SCmthd ==
        %'Bamber') or point of first reflection, i.e. point on DEM closest to
        %satellite (SCmthd == 'Slobbe')
        if strcmp(SCmthd,'Slobbe')
            CS.('LATr') = nan(size(CS.GEO.LAT)); CS.('LONr') = nan(size(CS.GEO.LON)); CS.('SI') = nan(size(CS.GEO.LON)); CS.('hr') = nan(size(CS.GEO.LON));
            CS          = ClosestPoint2DEM(CS,IDXpoi,DEM,DDAcf);
            [Xs,Ys]     = polarstereo_fwd(CS.LATr,CS.LONr,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
        else
            [Xs,Ys]     = polarstereo_fwd(CS.GEO.LAT,CS.GEO.LON,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
        end
        aAi         = cosd(aA)+1i*sind(aA);
        aA_s        = wrapTo360(rad2deg(angle(interp2(DEMx,DEMy,aAi,Xs,Ys,'linear'))));
        sA_s        = interp2(DEMx,DEMy,sA,Xs,Ys,'linear');
        clear aAi
        
        %Compute curvature radii and corrected height (see
        %/data/Projects/CryoSat2_Greenland/Literature/Slope_Correction/Slope_Correction_equations.odt)
        %for a deriviation of the equations [Bamber 1994, and Cooper 1989]
        R_N         = DDAcf.RefEll.SemimajorAxis./sqrt(1-DDAcf.RefEll.Eccentricity^2*sind(CS.GEO.LAT).^2);
        R_M         = (DDAcf.RefEll.SemimajorAxis*(1-DDAcf.RefEll.Eccentricity^2))./(1-DDAcf.RefEll.Eccentricity^2*sind(CS.GEO.LAT).^2).^1.5;
        R_alpha     = (R_N.*R_M)./(R_N.*cosd(aA_s).^2 + R_M.*sind(aA_s).^2);
        R_s         = R_alpha + CS.GEO.H;
        aGAMMA      = asind(CS.range.*sind(sA_s(:))./R_s(:));
        CS.HEI      = (R_s(:).*sind(sA_s(:)-aGAMMA)./sind(sA_s(:))) - R_alpha(:);
        
        %Compute associated lat/lon coordinates
        [CS.('LATc'),CS.('LONc')] = reckon(CS.GEO.LAT(:),CS.GEO.LON(:),R_alpha(:).*deg2rad(aGAMMA),aA_s(:),DDAcf.RefEll);
        
        %Corrected lat,lon according to Bamber 1994. So far, I have not been
        %able to derive these equations
        % XYZsat = CoordinateTrans('ell2xyz',[CS.GEO.LON(:),CS.GEO.LAT(:),CS.GEO.H(:)],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');
        % dX     = R_alpha(:).*deg2rad(aGAMMA).*cosd(aA_s(:));
        % dY     = R_alpha(:).*deg2rad(aGAMMA).*sind(aA_s(:));
        % LLHsat = CoordinateTrans('xyz2ell',[XYZsat(:,1) + dX,XYZsat(:,2) + dY,XYZsat(:,3)],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');
        % CS.('LONc') = LLHsat(:,1);
        % CS.('LATc') = LLHsat(:,2);
    end
else
    CS.('LONc') = CS.GEO.LON(:);
    CS.('LATc') = CS.GEO.LAT(:);
end

%% Compute backscatter coefficient (sigma0 / sigma-naught)
switch SAT
    case {'CS','CryoSat'}
        CS.Pu         = CS.Pu ./ NORMfactor(:);
        CS.('sigma0') = Compute_Sigma0(CS,DDAcf);
        CS.Pu         = CS.Pu .* NORMfactor(:);
    case {'S3A','Sentinel-3A','S3B','Sentinel-3B'}
        %Surface Topography Mission (STM) SRAL/MWR L2 Algorithms
        %Definition, Accuracy and Specification, Section 2.15.3.3
        CS.('sigma0') = CS.LRM.scale_factor_ku_l1b_echo_sar_ku + 10*log10(CS.Pu) + DDAcf.Biases.sigma0;
    otherwise
        error('SAT: %s not recognized',SAT)
end

%% Analyze output
if MakePlots
    for i = IDXpoi
        if all(isnan(CS.n_ret(i,:))), continue, end
    
        figure('Position',get(0,'Screensize'));

        %Select valid entries
        IDXvalid = ~isnan(CS.n_ret(i,:));

        %Plot observed/reconstructed waveform
        WD        = CS.LRM.data(:,IDXpoi);
        WDnorm    = NORMdata(:,i);
        StrYlabel = 'Power (Watts)';
        if ~any(strcmp(Retracker,{'ICE1','ICE3','OCOG','Threshold'})); WD = WDnorm; StrYlabel = 'Normalized power'; end
        subplot(2,3,1),plot(1:DDAcf.Ns,WD,'.-','LineWidth',2),hold on
        if ~all(isnan(CS.WDrecon(:,IDXpoi)))
            plot(1:DDAcf.Ns,CS.WDrecon(:,IDXpoi),'g-','LineWidth',2)
            plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:DDAcf.Ns,CS.WDrecon(:,IDXpoi),CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
            LegendEntries = {'Observed waveform','Reconstructed waveform','Retracking points'};
        else
            plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:DDAcf.Ns,WD,CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
            LegendEntries = {'Observed waveform','Retracking points'};
        end
        axis([0 DDAcf.Ns 0 max(WD)])
        title('Observed waveform','fontsize',FntSz)
        legend(LegendEntries,'box','off')
        xlabel('Bin','fontsize',FntSz)
        ylabel(StrYlabel,'fontsize',FntSz)
        set(gca,'fontsize',FntSz)

        WD        = CS.LRM.data(:,IDXpoi);
        subplot(2,3,2),plot(1:DDAcf.Ns,WD,'.-','LineWidth',2),hold on
        if ~all(isnan(CS.WDrecon(:,IDXpoi)))
            plot(1:DDAcf.Ns,fun(x,XDATAorg)/NORMfactor(IDXpoi),'g-','LineWidth',2)
            plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:DDAcf.Ns,fun(x,XDATAorg)/NORMfactor(IDXpoi),CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
            LegendEntries = {'Observed waveform','Reconstructed waveform','Retracking points'};
        else
            plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:DDAcf.Ns,WD,CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
            LegendEntries = {'Observed waveform','Retracking points'};
        end
        plot(1:DDAcf.Ns,ones(1,DDAcf.Ns)*CS.Pu(IDXpoi)/NORMfactor(IDXpoi),'k--');
        LegendEntries = horzcat(LegendEntries,{'Pu'});
        axis([0 DDAcf.Ns 0 max(WD)])
        title('Observed waveform','fontsize',FntSz)
        legend(LegendEntries,'box','off')
        xlabel('Bin','fontsize',FntSz)
        ylabel('Power (Watts)','fontsize',FntSz)
        set(gca,'fontsize',FntSz)

        %Plot coastline
        LOLA = Get_River_Border_data('GSHHS','coastlines','f',DOM(1,:),DOM(2,:));
        subplot(2,3,4),plot(LOLA(:,1),LOLA(:,2),'k-'),hold on
        plot(CS.GEO.LON(CS.GEO.LON~=0),CS.GEO.LAT(CS.GEO.LON~=0),'.')
        plot(CS.GEO.LON(i),CS.GEO.LAT(i),'ro','LineWidth',2)
        axis([minmax(LOLA(:,1)) minmax(LOLA(:,2))])
        xlabel('Longitude (deg)','fontsize',FntSz)
        ylabel('Latitude (deg)','fontsize',FntSz)
        set(gca,'fontsize',FntSz)

        %Plot DEM around point of interest including ground track
        if ~isempty(DEM)
        [x,y]  = polarstereo_fwd(CS.GEO.LAT,CS.GEO.LON,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
        [~,Ix] = min(abs(DEM.x(:)-x(i))); [~,Iy] = min(abs(DEM.y(:)-y(i)));
        MDidx  = ceil(20E3/DEM.GSP);
        subplot(2,3,5),imagesc(DEM.x(Ix+(-MDidx:MDidx))./1000,DEM.y(Iy+(-MDidx:MDidx))./1000,DEM.z((Iy+(-MDidx:MDidx)),(Ix+(-MDidx:MDidx)))),hold on
        axis square; cb = colorbar; set(gca,'YDir','normal')
        cb.Label.String = 'meters'; cb.Label.FontSize = FntSz;
        %Plot direction of flight
        quiver(x(i)./1000,y(i)./1000,(x(i+1)-x(i))./1000,(y(i+1)-y(i))./1000,20,'Color','k','LineWidth',2,'MaxHeadSize',0.8);
        %Plot DEM points within range window + mark closest point
        R0          = CS.MEA.ref_range(i) - (0.5*DDAcf.c*(RefBin*dt));
        Rend        = CS.MEA.ref_range(i) + (0.5*DDAcf.c*(RefBin*dt));
        XYZsat      = CoordinateTrans('ell2xyz',[CS.GEO.LON(i),CS.GEO.LAT(i),CS.GEO.H(i)],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');
        [DEMx,DEMy] = meshgrid(DEM.x(Ix+(-MDidx:MDidx)),DEM.y(Iy+(-MDidx:MDidx)));
        DEMz        = double(DEM.z((Iy+(-MDidx:MDidx)),(Ix+(-MDidx:MDidx))));
        [la,lo]     = polarstereo_inv(DEMx(:),DEMy(:),DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
        XYZdem      = CoordinateTrans('ell2xyz',[lo(:),la(:),DEMz(:)],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');
        RngSS       = pdist2(XYZsat,XYZdem);
        IDXinRW     = (RngSS >= R0-1) & (RngSS <= Rend+1);
        IDXorg      = 1:numel(RngSS); IDXorg = IDXorg(IDXinRW);
        [~,IDXsrt]  = sort(RngSS(IDXinRW));
        NrPnts      = floor(nnz(IDXinRW)./5);
        for j=1:5
            IDXtmp  = IDXsrt(((j-1)*NrPnts)+1:j*NrPnts);
            plot(DEMx(IDXorg(IDXtmp))./1000,DEMy(IDXorg(IDXtmp))./1000,'x','Color',[.2 .2 .2]*j)
        end
        plot(DEMx(IDXorg(IDXsrt(1)))./1000,DEMy(IDXorg(IDXsrt(1)))./1000,'rx','LineWidth',2)

        [xc,yc]     = polarstereo_fwd(CS.LATc,CS.LONc,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
        plot(xc/1000,yc/1000,'ro','LineWidth',2)
        
        mycmap = cell2mat(struct2cell(load(fullfile(Pth_mat,'ColormapStdev'))));
        mycmap = interp1((1:size(mycmap,1))',mycmap,linspace(1,size(mycmap,1),256)');
        ax2 = subplot(2,3,6); imagesc(DEM.x(Ix+(-MDidx:MDidx))./1000,DEM.y(Iy+(-MDidx:MDidx))./1000,abs(reshape(RngSS,numel(-MDidx:MDidx),numel(-MDidx:MDidx))-CS.range(IDXpoi))); hold on
        axis square; cb2 = colorbar; set(gca,'YDir','normal')
        cb2.Label.String = 'meters'; cb2.Label.FontSize = FntSz;
        colormap(ax2,mycmap);
        caxis([0 60])
        %Plot direction of flight
        quiver(x(i)./1000,y(i)./1000,(x(i+1)-x(i))./1000,(y(i+1)-y(i))./1000,20,'Color','k','LineWidth',2,'MaxHeadSize',0.8);
        [~,IDXmR] = min(abs(RngSS-CS.range(IDXpoi)));
        plot(DEMx(IDXmR)./1000,DEMy(IDXmR)./1000,'bx','LineWidth',2)
        plot(xc/1000,yc/1000,'bo','LineWidth',2)
        
        end
       
        %Set background color to white
        set(gcf, 'Color',[1 1 1])
    end
end

%% Save output
DATA              = struct;
IDX               = ~isnan(CS.HEI);
DATA.('TIME')     = CS.TIME(IDX);
DATA.('LAT')      = CS.LATc(IDX);
DATA.('LON')      = CS.LONc(IDX);
DATA.('HEI')      = single(CS.HEI(IDX));
if ApplySlopeCorrection && strcmp(SCmthd,'Slobbe2'), DATA.('MinD2CDP') = single(CS.MinD2CDP(IDX)); end
DATA.('sumCORR')  = single(CS.SumCorrST(IDX));
DATA.('FlgC')     = CS.FlgC(IDX);
DATA.('SSHi')     = single(CS.SSHi(IDX));
DATA.('sumCORRi') = single(CS.SumCorrSTi(IDX));
DATA.('FlgCi')    = CS.FlgCi(IDX);
DATA.('sigma0')   = single(CS.sigma0(IDX));
DATA.('SWH')      = single(CS.SWH(IDX));
DATA.('SurfT')    = int8(CS.surf_type(IDX));
DATA.('WFc')      = int8(CS.WFc(IDX));
DATA.('nu')       = single(CS.nu(IDX));
DATA.('ExitF')    = int8(CS.ExitF(IDX));
DATA.('MF')       = single(CS.MF(IDX));
DATA.('PCorr')    = single(CS.PCorr(IDX));

end


function CS = ClosestPoint2DEM(CS,IDXpoi,DEM,DDAcf)

%Set number of grid cells to be considered in all directions from grid
%point closest to satellite point
Didx          = ceil(6E3/DEM.GSP);

%Transform DEM to ECEF cartesian coordinates
[DEMx,DEMy]   = meshgrid(DEM.x,DEM.y);
[la,lo]       = polarstereo_inv(DEMx,DEMy,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
XYZdem        = CoordinateTrans('ell2xyz',[lo(:),la(:),double(DEM.z(:))],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');
clear la lo

%Transform DEM shifted half the grid spacing to ECEF cartesian coordinates.
%Each point of this shifted grid is in the center of four grid points. It
%is used to construct four triangles
[DEMxc,DEMyc] = meshgrid(DEM.x(1:end-1)+(.5*DEM.GSP),DEM.y(1:end-1)-(.5*DEM.GSP));
DEMzc         = interp2(DEMx,DEMy,double(DEM.z),DEMxc,DEMyc,'linear');
[lac,loc]     = polarstereo_inv(DEMxc,DEMyc,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
XYZdemc       = CoordinateTrans('ell2xyz',[loc(:),lac(:),DEMzc(:)],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');
clear lac loc

%Find the DEM grid point closest to the satellite positions
[x_sat,y_sat] = polarstereo_fwd(CS.GEO.LAT,CS.GEO.LON,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
[~,Ix]        = min(abs(DEMxc(1,:)-x_sat(:)),[],2); [~,Iy] = min(abs(DEMyc(:,1)'-y_sat(:)),[],2);
clear x_sat y_sat

%Transform the satellite positions to ECEF cartesian coordinates
XYZsat        = CoordinateTrans('ell2xyz',[CS.GEO.LON(:),CS.GEO.LAT(:),CS.GEO.H(:)],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');

%Compile matrix with GLOBAL VERTEX INDICES of the 4 triangles per grid cell:
% 1 - 2
% | 5 |
% 4 - 3
IDX           = nan(size(DEMx));  IDX(:)  = 1:numel(IDX);
IDXc          = nan(size(DEMxc)); IDXc(:) = 1:numel(IDXc);
IDXtr_g       = nan(numel(DEMxc),5);
IDXtr_g(:,1)  = reshape(IDX(1:end-1,1:end-1),numel(DEMxc),1);
IDXtr_g(:,2)  = reshape(IDX(1:end-1,2:end),numel(DEMxc),1);
IDXtr_g(:,3)  = reshape(IDX(2:end,2:end),numel(DEMxc),1);
IDXtr_g(:,4)  = reshape(IDX(2:end,1:end-1),numel(DEMxc),1);
IDXtr_g(:,5)  = reshape(IDXc,numel(DEMxc),1);

%Compile matrix with LOCAL VERTEX INDICES of the 4 triangles per grid cell
IDXc          = nan(2*Didx+1); IDXc(:) = 1:numel(IDXc);
IDX           = nan(2*Didx+2); IDX(:)  = (1:numel(IDX))+numel(IDXc);
IDXtr_l       = nan(numel(IDXc),5);
IDXtr_l(:,1)  = reshape(IDX(1:end-1,1:end-1),numel(IDXc),1);
IDXtr_l(:,2)  = reshape(IDX(1:end-1,2:end),numel(IDXc),1);
IDXtr_l(:,3)  = reshape(IDX(2:end,2:end),numel(IDXc),1);
IDXtr_l(:,4)  = reshape(IDX(2:end,1:end-1),numel(IDXc),1);
IDXtr_l(:,5)  = reshape(IDXc,numel(IDXc),1);
clear IDX IDXc

%Extract list of unique LOCAL vertex indices
[IDXlgrd_l,IA] = unique(IDXtr_l(:,1:4)','stable');
IDXlgrd_l      = IDXlgrd_l-size(IDXtr_l,1);

%Construct faces
ALLfaces       = reshape(IDXtr_l(:,[1 2 5;2 3 5;3 4 5;4 1 5]')',3,4*size(IDXtr_l,1))';

%For each point, construct a local triangulated surface and compute the
%shortest line connecting the satellite position and this triangulation in
%3D
for i = IDXpoi
    if isnan(CS.range(i))
        continue
    end
    
    %Get the GLOBAL vertex indices of the LOCAL grid. Note, we separately
    %get the indices of the original grid and those of the shifted grid
    %(the center points used to create the triangles)
    [IDXlgrd_g_c,IDXlgrd_g_r] = meshgrid((Ix(i)+(-Didx:Didx)),(Iy(i)+(-Didx:Didx)));
    IDXlgrd_gc                = sub2ind(size(DEMxc),IDXlgrd_g_r(:),IDXlgrd_g_c(:));
    IDXlgrd_g                 = IDXtr_g(IDXlgrd_gc,1:4)';
    
    %Copy vertex XYZ coordinates to XYZdem_lgrd and extract unique vertices
    DUM                       = XYZdem(IDXlgrd_g(:),:);
    DUM                       = DUM(IA,:);
    XYZdem_lgrd               = nan(size(DUM));
    XYZdem_lgrd(IDXlgrd_l,:)  = DUM;

    %Select faces/vertices in range window
    RngSS                     = pdist2(XYZsat(i,:),XYZdemc(IDXlgrd_gc,:));
    IDXinRW                   = (RngSS >= CS.R0(i)) & (RngSS <= CS.Rend(i));
    IDXFinRW                  = ismember(ALLfaces(:,3),find(IDXinRW));
    XYZdem_lgrd               = XYZdem_lgrd(unique(ALLfaces(IDXFinRW,1:2))-numel(IDXlgrd_gc),:);
    IDXlgrd_gc                = IDXlgrd_gc(IDXinRW);
    FinRW                     = ALLfaces(IDXFinRW,:);
    [VertIDXorg,IDXfaces]     = sort(FinRW(:));
    FinRW(IDXfaces)           = [1;1+cumsum(diff(VertIDXorg) ~= 0)];
    
    %Compile list of faces and vertices
    FV.faces                  = FinRW;
    FV.vertices               = [XYZdemc(IDXlgrd_gc,:);XYZdem_lgrd];
    
    %Compute shortest line between satellite and DEM
    [SI,SP]                   = point2trimesh(FV, 'QueryPoints', XYZsat(i,:));
    CS.SI(i)                  = abs(SI);
    
    %Transform ECEF cartesian coordinates back to ellipsoidal
    SP_llh                    = CoordinateTrans('xyz2ell',SP,1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');
    
    %Store lat/lon of point on DEM closest to satellite position
    CS.LATr(i)                = SP_llh(2); CS.LONr(i) = SP_llh(1); CS.hr(i) = SP_llh(3);
end

end


function CS = Compute_Impact_Point_Roemer(CS,IDXpoi,DEM,DDAcf)

%COMPUTE_IMPACT_POINT_ROEMER computes the impact point, i.e. the position
%to which the radar measurement actually refers, according to the method
%proposed by Roemer et al. 2007 (Refined analysis of radar altimetry data
%applied to the region of the subglacial Lake Vostok/Antarctica)

%Set number of grid cells to be considered in all directions from grid
%point closest to satellite point
Didx          = round(8E3/DEM.GSP)+2; %2 cells are extra to handle situation minimum range is on edge of grid
SzeLOCgrd     = numel(-Didx:Didx);

%Size PLF. The PLF is assumed to be a rectangular area. Here we assume the
%footprint has a size of 2 by 2 km
Sz_PLF        = 2*ceil(1E3/DEM.GSP) + 1;

%Find the indices of the DEM grid point closest to the satellite positions
[DEMx,DEMy]   = meshgrid(DEM.x,DEM.y);
[x_sat,y_sat] = polarstereo_fwd(CS.GEO.LAT,CS.GEO.LON,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
[~,Ix]        = min(abs(DEMx(1,:)-x_sat(:)),[],2); [~,Iy] = min(abs(DEMy(:,1)'-y_sat(:)),[],2);

%Get the GLOBAL vertex indices of the LOCAL grids
IDXlgrd_g     = repmat((-Didx:Didx)',1,SzeLOCgrd);
IDXlgrd_g     = IDXlgrd_g + repmat((-Didx:Didx)*size(DEMx,1),SzeLOCgrd,1);
IDXlgrd_g     = IDXlgrd_g(:)' + sub2ind(size(DEMx),Iy,Ix);
IDXlgrd_g     = IDXlgrd_g(IDXpoi,:);

%Transform DEM coordinates of local grids to ECEF cartesian coordinates
[la,lo]       = polarstereo_inv(DEMx(IDXlgrd_g),DEMy(IDXlgrd_g),DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
XYZdem        = CoordinateTrans('ell2xyz',[lo(:),la(:),double(DEM.z(IDXlgrd_g(:)))],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');

%Transform the satellite positions to ECEF cartesian coordinates
XYZsat        = CoordinateTrans('ell2xyz',[CS.GEO.LON(IDXpoi'),CS.GEO.LAT(IDXpoi'),CS.GEO.H(IDXpoi')],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');

%Compute ranges between satellite and DEM points
dX            = XYZsat(:,1) - reshape(XYZdem(:,1),size(IDXlgrd_g));
dY            = XYZsat(:,2) - reshape(XYZdem(:,2),size(IDXlgrd_g));
dZ            = XYZsat(:,3) - reshape(XYZdem(:,3),size(IDXlgrd_g));
RngSS         = reshape(vecnorm([dX(:),dY(:),dZ(:)],2,2),size(IDXlgrd_g));

%Calculate the averaging operators in both dimensions
[~,~,~,CTl,CTr] = blockmean(randn(SzeLOCgrd),[Sz_PLF Sz_PLF],[Sz_PLF-1 Sz_PLF-1]);

%Compute dX and dY based on which local grids are created
[dXloc,dYloc]   = meshgrid((-Didx+floor(Sz_PLF/2):Didx-floor(Sz_PLF/2)).*mode(diff(DEM.x)),(-Didx+floor(Sz_PLF/2):Didx-floor(Sz_PLF/2)).*mode(diff(DEM.y)));
% [dXfloc,dYfloc] = meshgrid(dXloc(1,1):10:dXloc(1,end),dYloc(1,1):-10:dYloc(end,1));
[dXfloc,dYfloc] = meshgrid(-90:10:90,90:-10:-90);
Xc              = DEMx(sub2ind(size(DEMx),Iy(IDXpoi),Ix(IDXpoi)));
Yc              = DEMy(sub2ind(size(DEMy),Iy(IDXpoi),Ix(IDXpoi)));

%For each point, construct local grid and computed location of grid point
%closest to satellite position
[Xmin,Ymin] = deal(nan(numel(IDXpoi),1));
for i = 1:numel(IDXpoi)
    if isnan(CS.range(IDXpoi(i)))
        continue
    end
    
    %Apply blockmean operator with overlap
    LPF_RngSS  = CTl'*reshape(RngSS(i,:),SzeLOCgrd,SzeLOCgrd)*CTr/prod([Sz_PLF Sz_PLF]);
    %LPF_RngSS  = blockmean(reshape(RngSS(i,:),SzeLOCgrd,SzeLOCgrd),[Sz_PLF Sz_PLF],[Sz_PLF-1 Sz_PLF-1]);
    
    %Densify grid
    LPF_RngSS_r   = LPF_RngSS(3:end-2,3:end-2);
    [~,IDXmin]    = min(LPF_RngSS_r(:));
    [IDr,IDc]     = ind2sub(size(LPF_RngSS_r),IDXmin);
    IDr           = IDr+2; IDc = IDc+2;
    LPF_RngSSf    = interp2(dXloc(IDr-2:IDr+2,IDc-2:IDc+2)+Xc(i),dYloc(IDr-2:IDr+2,IDc-2:IDc+2)+Yc(i),LPF_RngSS(IDr-2:IDr+2,IDc-2:IDc+2),dXfloc+dXloc(IDr,IDc)+Xc(i),dYfloc+dYloc(IDr,IDc)+Yc(i),'linear');
    
    % %Densify grid
    % LPF_RngSSf    = interp2(dXloc+Xc(i),dYloc+Yc(i),LPF_RngSS,dXfloc + Xc(i),dYfloc + Yc(i),'linear');
    
    %Identify location for which mean range is minimum
    [CS.SI(IDXpoi(i)),IDXminR] = min(LPF_RngSSf(:));
    Xmin(i) = dXfloc(IDXminR)+dXloc(IDr,IDc)+Xc(i);
    Ymin(i) = dYfloc(IDXminR)+dYloc(IDr,IDc)+Yc(i);
end

%Compute lat/lon coordinates of impact point
[CS.LATc(IDXpoi),CS.LONc(IDXpoi)] = polarstereo_inv(Xmin,Ymin,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);

%Evaluate DEM height @ impact point
CS.hr(IDXpoi) = DARTutils.EvalDEMatLATLON(DEM,CS.LATc(IDXpoi),CS.LONc(IDXpoi));

end


function CS = Compute_Impact_Point_Slobbe(CS,IDXpoi,DEM,DDAcf)

%COMPUTE_IMPACT_POINT_SLOBBE computes the impact point, i.e. the position
%to which the radar measurement actually refers, as the average of all
%DEM grid points in view of the satellite. These are the points for which
%the range between the satellite and DEM surface is in the domain where the
%leading edge is rising. Subsequently, Roemer's method is used to compute
%the slope correction 

%Set number of grid cells to be considered in all directions from grid
%point closest to satellite nadir point
Didx          = ceil(4E3/DEM.GSP);

%Find the indices of the DEM grid points closest to the satellite nadir
%positions
[DEMx,DEMy]   = meshgrid(DEM.x,DEM.y);
[x_sat,y_sat] = polarstereo_fwd(CS.GEO.LAT,CS.GEO.LON,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
[~,Ix]        = min(abs(DEMx(1,:)-x_sat(:)),[],2); [~,Iy] = min(abs(DEMy(:,1)'-y_sat(:)),[],2);

%Get the GLOBAL vertex indices of the LOCAL grids
IDXlgrd_g     = repmat((-Didx:Didx)',1,numel(-Didx:Didx));
IDXlgrd_g     = IDXlgrd_g + repmat((-Didx:Didx)*size(DEMx,1),numel(-Didx:Didx),1);
IDXlgrd_g     = IDXlgrd_g(:)' + sub2ind(size(DEMx),Iy,Ix);
IDXlgrd_g     = IDXlgrd_g(IDXpoi,:);

%Transform DEM coordinates of local grids to ECEF cartesian coordinates
X_DEM         = DEMx(IDXlgrd_g); Y_DEM = DEMy(IDXlgrd_g);
[la,lo]       = polarstereo_inv(X_DEM,Y_DEM,DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);
XYZdem        = CoordinateTrans('ell2xyz',[lo(:),la(:),double(DEM.z(IDXlgrd_g(:)))],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');

%Transform the satellite positions to ECEF cartesian coordinates
XYZsat        = CoordinateTrans('ell2xyz',[CS.GEO.LON(IDXpoi'),CS.GEO.LAT(IDXpoi'),CS.GEO.H(IDXpoi')],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');

%Compute ranges between satellite and DEM points
dX            = XYZsat(:,1) - reshape(XYZdem(:,1),size(IDXlgrd_g));
dY            = XYZsat(:,2) - reshape(XYZdem(:,2),size(IDXlgrd_g));
dZ            = XYZsat(:,3) - reshape(XYZdem(:,3),size(IDXlgrd_g));
RngSS         = reshape(vecnorm([dX(:),dY(:),dZ(:)],2,2),size(IDXlgrd_g));

%Select for each row whether range is between R0 and Rend (here we make use
%of implicit array expansion).
IDXinRW       = (RngSS >= CS.R0(IDXpoi)) & (RngSS <= CS.Rend(IDXpoi));

%Select entries for which for no DEM point the range is between R0 and
%Rend. Compute the differences between RngSS and the retracked range.
%Replace R0 and Rend by the minimum of these differences and the minimum of
%these differences + Rend-R0. Redo selection
IDXour            = ~any(IDXinRW,2);
R0                = CS.R0(IDXpoi);
Rend              = CS.Rend(IDXpoi);
RngSS_min_r       = RngSS-CS.range(IDXpoi);
RngSS_min_r       = RngSS_min_r(IDXour,:);
IDXinRW(IDXour,:) = (RngSS_min_r >= min(RngSS_min_r,[],2)) & (RngSS_min_r <= min(RngSS_min_r,[],2)+(Rend(IDXour)-R0(IDXour)));

%Use IDXinRW to put all values in RngSS, X_DEM, and Y_DEM  to NaN for which
%the range is outside the interval R0 to Rend
X_DEM(~IDXinRW) = NaN; Y_DEM(~IDXinRW) = NaN;

%Compute lat,lon coordinates of impact point as the mean of all DEM points
%in view of the satellite
[CS.LATc(IDXpoi),CS.LONc(IDXpoi)] = polarstereo_inv(nanmean(X_DEM,2),nanmean(Y_DEM,2),DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity,DEM.phi_c,DEM.lambda_0);

%Compute terrain corrected height (as in Roemer)
DEM_i          = DARTutils.EvalDEMatLATLON(DEM,CS.LATc(IDXpoi),CS.LONc(IDXpoi));
XYZdem_i       = CoordinateTrans('ell2xyz',[CS.LONc(IDXpoi),CS.LATc(IDXpoi),DEM_i],1,2,3,[DDAcf.RefEll.SemimajorAxis,DDAcf.RefEll.Eccentricity^2],'EarthSph');
Rng_i          = vecnorm(XYZsat-XYZdem_i,2,2);
CS.HEI(IDXpoi) = CS.HEI(IDXpoi) + nanmean((Rng_i - (CS.GEO.H(IDXpoi)'-DEM_i)),2);

%Compute minimum distance between impact point and closest DEM point
CS.MinD2CDP(IDXpoi) = min(sqrt((X_DEM-nanmean(X_DEM,2)).^2 + (Y_DEM-nanmean(Y_DEM,2)).^2),[],2);

% PntID = 10;
% figure,plot(x_sat(IDXpoi(PntID)),y_sat(IDXpoi(PntID)),'kx','LineWidth',2)
% hold on,plot(DEMx(IDXlgrd_g(PntID,:)),DEMy(IDXlgrd_g(PntID,:)),'b.'),axis square
% lgrd = IDXlgrd_g(PntID,:);
% hold on,plot(DEMx(lgrd(IDXinRW(PntID,:))),DEMy(lgrd(IDXinRW(PntID,:))),'g.')
% TMP_DEMx  = reshape(DEMx(IDXlgrd_g(PntID,:)),numel(-Didx:Didx),numel(-Didx:Didx));
% TMP_DEMy  = reshape(DEMy(IDXlgrd_g(PntID,:)),numel(-Didx:Didx),numel(-Didx:Didx));
% TMP_RngSS = reshape(RngSS(PntID,:),numel(-Didx:Didx),numel(-Didx:Didx));
% figure,imagesc(TMP_DEMx(1,:),TMP_DEMy(:,1),TMP_RngSS-CS.range(IDXpoi(PntID)));
% axis square; colorbar; set(gca,'YDir','normal')

end
