function [out_l2, out_l1b] = FFSAR_proc_and_retrack_l2target(dataset, ffsar_proc_sets, proc_retrack_sets)

out_l2 = [];
% %get l1a fileinfo
% [fpath_l1a, basename_l1a, fext] = fileparts(dataset.L1A);
% if contains(basename_l1a, 'measurement')
%     [~, basename_l1a] = fileparts(fpath_l1a);
% end
% 
% basename_l1a_parts = split(basename_l1a,'_');
% mission = basename_l1a_parts{1};

mission = mission_from_fname(dataset.L1A);

is_s6 = contains(mission, 'S6');
is_s3 = contains(mission, 'S3');
is_cs = contains(mission, 'CS');

% load constants
DDAcf = DDA_ConfigFile(mission,'SAR');
[CONST,~,~] = FFSAR_LoadCommonSettings;

%% L1b processing
hash_run = DataHash(struct('sets', ffsar_proc_sets, 'fname', dataset.L1A, 'DOM', dataset.DOM, 'DDAcf', DDAcf));
tmp_dir = fullfile(getenv('HOME'), '/temp/');
if ~exist(tmp_dir,'dir')
    mkdir(tmp_dir)
end

l1b_temp_filepath = [tmp_dir hash_run '.mat'];

if ~isfile(l1b_temp_filepath) || proc_retrack_sets.ignore_cashed_files
    tic
    ffproc = FFSAR_Processor(dataset.L1A,dataset.DOM,ffsar_proc_sets);
    ffproc.setup_proc();
    ffproc.proc();
    toc

    CS1b_ff = ffproc.CS1b;
    save(l1b_temp_filepath, 'CS1b_ff', '-v7.3');
else
    fprintf('Cached L1b data loaded.\n')
    CS1b_ff = load(l1b_temp_filepath).CS1b_ff;
end

%% apply IRF correction
if ffsar_proc_sets.do_correct_IRF && is_s6
    kernel = IRF_corr(CS1b_ff, ffsar_proc_sets);
    CS1b_ff.SAR.data = conv2(CS1b_ff.SAR.data,kernel,'same');

    fprintf('IRF correction applied.')
end

%% L2 processing
l2_exists = isfield(dataset,'L2');

% auxiliary data    
if l2_exists
    l2_eum.LAT = ncread(dataset.L2,[DDAcf.ncgrp_l1b_l2 DDAcf.var_mapping_l2.lat]);
    l2_eum.LON = ncread(dataset.L2,[DDAcf.ncgrp_l1b_l2 DDAcf.var_mapping_l2.lat]);
    l2_eum.DIST2COAST = ncread(dataset.L2,[DDAcf.ncgrp_l1b_l2 DDAcf.var_mapping_l2.dist2coast]);
    l2_eum.record_ind = ncread(dataset.L2,[DDAcf.ncgrp_l1b_l2 'l2_record_counter']);

    % read L2 and define startLoc,count indices
    eps_ = 0.001;
    IDX_mask = ingeoquad(l2_eum.LAT,l2_eum.LON,dataset.DOM(1,:) + [-eps_ eps_],dataset.DOM(2,:));
%     IDX_mask = ingeoquad(l2_eum.LAT,l2_eum.LON,dataset.DOM(1,:),dataset.DOM(2,:));
    IDX = find(IDX_mask);
    record_inds_l2 = l2_eum.record_ind(IDX) + 1;  %make it 1-based
    startLoc_l2 = min(record_inds_l2);
    count_l2 = length(record_inds_l2);
end

% do_cleanse = true;
if exist('do_cleanse')
    % take pointwise-product of normalised single-look radargram and get maximum
    data_norm = CS1b_ff.SAR.data ./ max(CS1b_ff.SAR.data);
    n_avg_sl = 60;
    data_movprod = movprod(data_norm, n_avg_sl, 2);
    [~, max_inds_movprod] = max(data_movprod, [], 1);

    % get bad single-look inds
    % check if maximum is within +-n_inds_range range around movprod maximum (pointwise product)
    n_inds_range = 30;
    max_inds_movprod_mean = floor(mean(max_inds_movprod));
    mask_inds_bad = max_inds_movprod < (max_inds_movprod_mean - n_inds_range) | (max_inds_movprod > (max_inds_movprod_mean + n_inds_range));
    
    % mask out range bins before movprod mean, whose power exceeds thr_pwr
    thr_pwr = 0.3;
    mask_inds_bad = mask_inds_bad | (any(data_norm(1:max_inds_movprod_mean - n_inds_range, :) > thr_pwr));

    %%
%     figure; plot(inds_20_intp, max_inds_movprod)
%     hold on;
%     plot(inds_20_intp(mask_inds_bad), max_inds_movprod(mask_inds_bad), '.r')
    
    %% cleanse waveform
    CS1b_ff.SAR.data(:, mask_inds_bad) = 0;
end

l1b_len = length(CS1b_ff.GEO.LAT);
num_skipped_before = 0;
num_skipped_after = 0;

if l2_exists && proc_retrack_sets.collocate_with_l2
    % try to get the right frequency (multilooking) to match the L2 latitudes
    lat_l2 = l2_eum.LAT(startLoc_l2:startLoc_l2+count_l2-1);
    lat_ff = CS1b_ff.GEO.LAT;

    [~,closestIndex] = min(abs(lat_l2 - lat_ff),[],2);

    % get the corrected number of multilooks necessary to agree approximately with the L2 spacing
    N_avg = round(mean(diff(closestIndex)));

    oversample = round(proc_retrack_sets.target_posting_rate / 20); % from 20Hz to 80Hz means 4 times 'oversampling'
    if oversample > 1
        if mod(oversample,2) ~= 0 % odd-numbered oversampling
            diff_inds_subsamples = round(N_avg/oversample*(-(oversample-1)/2:(oversample-1)/2));
        else % odd-numbered oversampling
            diff_inds_subsamples = round(N_avg/oversample*(-(oversample/2-1):oversample/2));
        end
    else
        diff_inds_subsamples = 0;
    end
    
    closestIndex = closestIndex + diff_inds_subsamples;
    closestIndex = reshape(closestIndex',[], 1);
    if any(closestIndex <= 0)
        num_skipped_before = sum(closestIndex <= 0);
%         startLoc_l2 = startLoc_l2 + num_skipped_before;
%         count_l2 = count_l2 - 1;
        
        closestIndex = closestIndex(closestIndex > 0);
    end
    if any(closestIndex > l1b_len)
        num_skipped_after = sum(closestIndex > l1b_len);
%         count_l2 = count_l2 - num_skipped_after;
        
        closestIndex = closestIndex(closestIndex <= l1b_len);
    end

    N_avg = round(N_avg / oversample);
    
    CS1b_ff_avg = FF_SAR_Average_Waveforms(CS1b_ff,N_avg,mission,closestIndex);
    
elseif proc_retrack_sets.target_posting_rate > 0
    current_posting_rate = 1/mean(diff(CS1b_ff.GEO.Elapsed_Time));

    N_avg = ceil(current_posting_rate / proc_retrack_sets.target_posting_rate);
    closestIndex = (ceil((l1b_len / N_avg) / 2):N_avg:l1b_len);

    CS1b_ff_avg = FF_SAR_Average_Waveforms(CS1b_ff,N_avg,mission,closestIndex);
else
    error('invalid proc_and_retrack configuration.')
end

% retracking
if proc_retrack_sets.enable_l2_proc
    [out_l2,~ ] = SAR_L1b_to_L2(mission,CS1b_ff_avg,dataset.DOM,'SAMOSA2',...
    {'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'; 'fit_zero_doppler', mat2str(proc_retrack_sets.fit_zero_doppler)});
end

out_l1b = CS1b_ff_avg;

%% add some data to l1b/l2
l2_len = numel(CS1b_ff_avg.GEO.LAT);
if l2_exists
    IDXpoi = startLoc_l2:startLoc_l2+count_l2-1;
    
    if oversample > 1
        IDXpoi = repelem(IDXpoi - 1, oversample);  %-1, make it 0-based
        IDXpoi = IDXpoi(1+num_skipped_before:end-num_skipped_after);
    end
    out_l1b.IDXpoi = IDXpoi;
    
    if numel(out_l1b.IDXpoi) ~= l2_len
        error('l2 mapping index error. ');
    end
    
    out_l1b.GEO.dist2coast = l2_eum.DIST2COAST;
        
    % add to L2 data
    if proc_retrack_sets.enable_l2_proc
        out_l2.IDXpoi = IDXpoi;
        out_l2.GEO.dist2coast = l2_eum.DIST2COAST;
    end
else
    out_l1b.IDXpoi = (1:length(CS1b_ff_avg.GEO.LAT));
end

% append non-averaged single-look data to out_l1b
out_l1b.non_average.data = CS1b_ff.SAR.data;
if isfield(CS1b_ff.SAR, 'data_opt')
    out_l1b.non_average.data_opt = CS1b_ff.SAR.data_opt;
end
if isfield(CS1b_ff.SAR, 'data_opt2')
    out_l1b.non_average.data_opt2 = CS1b_ff.SAR.data_opt2;
end

if ffsar_proc_sets.split_aperture
    out_l1b.non_average.data_T_0 = CS1b_ff.SAR.data_T_0;
    out_l1b.non_average.data_0_T = CS1b_ff.SAR.data_0_T;
end
if isfield(CS1b_ff.SAR.FLAG, 'cleanse_pdm')
    out_l1b.non_average.cleanse_pdm = CS1b_ff.SAR.FLAG.cleanse_pdm;
end

if isfield(CS1b_ff.SAR, 'data_pseudoDD')
    out_l1b.non_average.data_pseudoDD = CS1b_ff.SAR.data_pseudoDD;
end
out_l1b.non_average.closestIndex = closestIndex;
out_l1b.non_average.N_avg = N_avg;
out_l1b.non_average.tracker_range = CS1b_ff.MEA.tracker_range;
out_l1b.non_average.LAT = CS1b_ff.GEO.LAT;
out_l1b.non_average.LON = CS1b_ff.GEO.LON;

end
