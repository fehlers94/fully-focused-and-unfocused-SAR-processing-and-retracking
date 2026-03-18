%% script to batch-process FFSAR data for coastalFFSAR study
%settings
batch_T_postrate = [...
    2.1, 20;...    
    2.1, 80;...
    2.1, 60;...
    2.1, 100;...
    2.1, 140;...
    2.1, 160;...
    2.4, 20;...
    2.4, 80;...
    2.4, 100;...
    ]';
batch_T = unique(batch_T_postrate(1,:));

skip_if_exists = true;
reverse_filer_order = false;

%create illumination time T and posting rate combinations
map_ignore_cached_files = containers.Map(batch_T, true([numel(batch_T),1]));


for i=1:size(batch_T_postrate,2)
    T  = batch_T_postrate(1,i);
%     T  = batch_T(i);
    posting_rate = batch_T_postrate(2,i);
    
    ignore_cached_files = map_ignore_cached_files(T);
    map_ignore_cached_files(T) = false;
    
    fprintf('processing FF-SAR data with T=%.1f, posting_rate=%i (ignore_cached_files: %d) \n', T, posting_rate, ignore_cached_files);
    

    FFSAR_proc_and_retrack('coastal_ffsar', T, posting_rate, ignore_cached_files, skip_if_exists, reverse_filer_order)
end
