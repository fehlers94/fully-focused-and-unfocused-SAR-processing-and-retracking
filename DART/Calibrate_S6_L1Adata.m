function CS1a = Calibrate_S6_L1Adata(CS1a,DDAcf)

%cal2
% CS1a.FBR.data = bsxfun(@times,CS1a.FBR.data,exp(1i*CS1a.MEA.cal2_correction(:,1)));

%burst phase
% CS1a.FBR.data = bsxfun(@times,CS1a.FBR.data,exp(1i*CS1a.MEA.burst_phase_cor(1,:)));

%burst power
% CS1a.FBR.data = bsxfun(@times,CS1a.FBR.data,CS1a.MEA.burst_power_cor(1,:));

% power-to-antenna-scaling
power_to_antenna_scaling(1,1,:) = CS1a.SAR.scale_factor_ku';
% CS1a.FBR.data = CS1a.FBR.data ./ (sqrt(10.^(power_to_antenna_scaling/10)));

% agc
agc(1,1,:) = CS1a.SAR.agc_ku';
% CS1a.FBR.data = CS1a.FBR.data./(sqrt(10.^(-agc/10)));


end