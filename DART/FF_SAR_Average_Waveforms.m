function CS1b = FF_SAR_Average_Waveforms(CS1b,N_avg,mission,IDX)

defval('N_avg',600) % number of single-look waveforms to be averaged
defval('pl',5E-6)   % limit on power
defval('IDX',[])   % limit on power

% this function first applies a moving average with N_avg and then
% subsamples it on the indices IDX. This gives us the freedom to supply the
% IDX that correspond to the latitudes/longitudes closest to the EUMETSAT
% L2 files in the first place, while just doing the moving average that
% corresponds to only approximately the right spacing

if isempty(IDX)
    IDX = 1+N_avg/2:N_avg:numel(CS1b.GEO.LAT);
    IDX(IDX + N_avg/2-1 > numel(CS1b.GEO.LAT)) = [];
end 
    
%% GEO
FldN = fieldnames(CS1b.GEO);
for i=1:numel(FldN)
    if strcmp(FldN{i},'TAI') % do not average seconds and microseconds but take center values them instead to avoid decimals in days, seconds, etc.
        CS1b.GEO.(FldN{i}).days = CS1b.GEO.(FldN{i}).days(IDX);
        CS1b.GEO.(FldN{i}).secs = CS1b.GEO.(FldN{i}).secs(IDX);
        CS1b.GEO.(FldN{i}).microsecs = CS1b.GEO.(FldN{i}).microsecs(IDX);
    elseif isstruct(CS1b.GEO.(FldN{i}))
        CS1b.GEO.(FldN{i}) = ApplyAveragingStrct(CS1b.GEO.(FldN{i}),N_avg,IDX);
    elseif isnumeric(CS1b.GEO.(FldN{i})) && ~isscalar(CS1b.GEO.(FldN{i}))
        CS1b.GEO.(FldN{i}) = ApplyAveraging(CS1b.GEO.(FldN{i}),N_avg,IDX);
    end
end

%% MEA
CS1b.MEA = ApplyAveragingStrct(CS1b.MEA,N_avg,IDX);

%% COR
FldN = fieldnames(CS1b.COR);
for i=1:numel(FldN)
    if isstruct(CS1b.COR.(FldN{i}))
        CS1b.COR.(FldN{i}) = ApplyAveragingStrct(CS1b.COR.(FldN{i}),N_avg,IDX);
    elseif isnumeric(CS1b.COR.(FldN{i})) && ~isscalar(CS1b.COR.(FldN{i}))
        CS1b.COR.(FldN{i}) = ApplyAveraging(CS1b.COR.(FldN{i}),N_avg,IDX);
    end
end

%% SAR
%Undo normalization.
%The power echo	sample values are all scaled to	fit	between	0 and 65535.
%The scaling factors can change	for	each waveform. To convert these back to
%values in Watts the following equation should be used (CryoSat Product
%Handbook, Eqn 4.2‐1):
%Power in Watts	= scaled value * (scale	factor * 10^‐9) * 2^scale power
if strcmp(mission,'CS')
    CS1b.SAR.data = bsxfun(@times,CS1b.SAR.data,reshape((CS1b.SAR.echo_scaling.*1E-9) .* 2.^(CS1b.SAR.echo_scale_power),1,size(CS1b.SAR.data,2),size(CS1b.SAR.data,3)));
end

FldN = fieldnames(CS1b.SAR);
for i=1:numel(FldN)
    if strcmp(FldN{i},'data_pseudoDD')
        CS1b.SAR.(FldN{i}) = CS1b.SAR.data_pseudoDD(:,IDX); % take the center value but do not apply any average on this field
    elseif isstruct(CS1b.SAR.(FldN{i}))
        CS1b.SAR.(FldN{i}) = ApplyAveragingStrct(CS1b.SAR.(FldN{i}),N_avg,IDX);
    elseif isnumeric(CS1b.SAR.(FldN{i})) && ~isscalar(CS1b.SAR.(FldN{i}))
        CS1b.SAR.(FldN{i}) = ApplyAveraging(CS1b.SAR.(FldN{i}),N_avg,IDX);
    end
end

% SAR: fix averaged flags, if not equal 0.0 then set it to 1.0
FldN = fieldnames(CS1b.SAR.FLAG);
for i=1:numel(FldN)
    CS1b.SAR.FLAG.(FldN{i}) = double(logical(CS1b.SAR.FLAG.(FldN{i})));
end

if strcmp(mission,'CS')
    %Normalize to range [0,65534]
    FAC                             = repmat(range(CS1b.SAR.data),size(CS1b.SAR.data,1),1,1)./(65534*(1 - (repmat(min(CS1b.SAR.data),size(CS1b.SAR.data,1),1,1)./CS1b.SAR.data)));
    [DUM,CS1b.SAR.echo_scale_power] = log2(squeeze(min(FAC)));
    CS1b.SAR.echo_scaling           = DUM*1E9;
    CS1b.SAR.data                   = CS1b.SAR.data.*(1./min(FAC));
end
end

function Avg = ApplyAveragingStrct(FLD,N_avg,IDX)
    MAvg = structfun(@(x)movmean(x,N_avg,2), FLD, 'UniformOutput', false);
    Avg  = structfun(@(x) (x(:,IDX)), MAvg, 'UniformOutput', false);
end

function Avg = ApplyAveraging(FLD,N_avg,IDX)
    MAvg = movmean(FLD,N_avg,2);
    Avg  = MAvg(:,IDX);
end

