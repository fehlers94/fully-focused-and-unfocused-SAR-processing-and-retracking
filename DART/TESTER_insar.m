%% Settings
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;

FFSAR_processing_settings.integration_time = 2;
% FFSAR_processing_settings.integration_time = 3.413;
%     FFSAR_processing_settings.do_correct_IRF = true;
proc_retrack_sets.target_posting_rate = 0;

nc_base_dir = '/lfs/DGFI24/bigdata/s6_data/JASON_CS_S6A_L1A_ALT_HR_NTC_F/'; 

% DOM = [61.95 62.04; 0, 360];
% DOM = [62.01 62.04; 0, 360];
% DOM = [62.01 62.13; 0, 360];

% DOM = [62.04 62.07; 0, 360];
% DOM = [62.27 62.30; 0, 360];
% dataset_l1a_b_l2(1).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211104T104138_20211104T113755_20211128T071705_3377_036_122_061_EUM__OPE_NT_F04.SEN6/measurement.nc'];
% dataset_l1a_b_l2(1).DOM = DOM;
% dataset_l1a_b_l2(2).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211114T084010_20211114T093626_20211208T062830_3376_037_122_061_EUM__OPE_NT_F04.SEN6/measurement.nc'];
% dataset_l1a_b_l2(2).DOM = DOM;
% % 
% DOM = [64.00 64.10; 0, 360];
% dataset_l1a_b_l2(1).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211106T093222_20211106T102839_20211129T001017_3377_036_172_086_EUM__OPE_NT_F04.SEN6/measurement.nc'];
% dataset_l1a_b_l2(1).DOM = DOM;
% 
% dataset_l1a_b_l2(2).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211116T073054_20211116T082710_20211208T235954_3376_037_172_086_EUM__OPE_NT_F04.SEN6/measurement.nc'];
% dataset_l1a_b_l2(2).DOM = DOM;


DOM = [64.40 64.50; 0, 360];
dataset_l1a_b_l2(1).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211029T073556_20211029T083213_20211118T174959_3377_035_219_109_EUM__OPE_NT_F04.SEN6/measurement.nc'];
dataset_l1a_b_l2(1).DOM = DOM;

dataset_l1a_b_l2(2).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211029T073556_20211029T083213_20211118T174959_3377_035_219_109_EUM__OPE_NT_F04.SEN6/measurement.nc'];
dataset_l1a_b_l2(2).DOM = DOM;



%% L1b-L2-processing of each dataset
res = struct([]);
for i = 1:numel(dataset_l1a_b_l2)
    dataset = dataset_l1a_b_l2(i);
    
    is_l1a_set = isfield(dataset,'L1A') && ~isempty(dataset);
          
    %get l1a fileinfo
    [fpath_l1a, basename, ~] = fileparts(dataset.L1A);
    if contains(basename, 'measurement')
        [~, basename] = fileparts(fpath_l1a);
    end
   
    mission = mission_from_fname(basename);
    DDAcf = DDA_ConfigFile(mission,'SAR');

    fprintf(['Processing file ',num2str(i),' of ',num2str(numel(dataset_l1a_b_l2)),': ' basename '\n']);
    DOM = dataset.DOM;
    
    %% L1b processing
    hash_run = DataHash(struct('sets', FFSAR_processing_settings, 'fname', dataset.L1A, 'DOM', dataset.DOM, 'DDAcf', DDAcf));
    tmp_dir = fullfile(getenv('HOME'), '/temp/');
    if ~exist(tmp_dir,'dir')
        mkdir(tmp_dir)
    end

    l1b_temp_filepath = [tmp_dir hash_run '.mat'];

    if ~isfile(l1b_temp_filepath)
        if isempty(gcp('nocreate')) && (FFSAR_processing_settings.n_cores_parfor > 0)
            n_cores = FFSAR_processing_settings.n_cores_parfor;
            parpool(n_cores);
            fprintf(['parfor processing enabled, using ' num2str(n_cores) ' cores.'])
        end

        tic
        ffproc = FFSAR_Processor(dataset.L1A,dataset.DOM,FFSAR_processing_settings);
        ffproc.setup_proc();
        ffproc.proc();
        toc

        CS1b_ff = ffproc.CS1b;
        save(l1b_temp_filepath, 'CS1b_ff', '-v7.3');
    else
        fprintf('Cached L1b data loaded.\n')
        CS1b_ff = load(l1b_temp_filepath).CS1b_ff;
    end

    res(i).data = CS1b_ff;
end

%%
distance_km = deg2km(distance(res(1).data.GEO.LAT(1),res(1).data.GEO.LON(1),res(2).data.GEO.LAT(1),res(2).data.GEO.LON(1)));
disp(['distance_km between first latlon points: ' num2str(distance_km) 'km']);

phases1 = angle(res(1).data.SAR.dataIQ);
% phases1 = unwrap(angle(res(1).data.SAR.dataIQ),[],1);
% phases1 = unwrap(angle(res(1).data.SAR.dataIQ),[],2);
phases2 = angle(res(2).data.SAR.dataIQ);
% phases2 = unwrap(angle(res(2).data.SAR.dataIQ),[],1);
% phases2 = unwrap(angle(res(2).data.SAR.dataIQ),[],2);

N_avg = 10;
IDX = N_avg / 2:N_avg:size(phases1, 2);
phases1 = ApplyAveraging(phases1,N_avg,IDX);
phases2 = ApplyAveraging(phases2,N_avg,IDX);

% phase_diff = phases1 - phases2;
phase_diff = unwrap(phases1 - phases2,[],1);

row = 320;

figure;

subplot(1,4,1);
imagesc(log10(res(1).data.SAR.data))
title('log10(power1)')
xlabel('along-track [sl]'); ylabel('range gates [bin]');

subplot(1,4,2);
imagesc(log10(res(2).data.SAR.data))
title('log10(power2)')
xlabel('along-track [sl]'); ylabel('range gates [bin]');

subplot(1,4,3);
imagesc(phase_diff);

xlabel(['along-track [' num2str(N_avg) '-ml]']); ylabel('range gates [bin]');
colorbar;
% caxis([8,11.8]);
colormap(jet);
yline(row,'-','unwrap(phase\_diff)');

title('unwrap(phase\_diff)')

subplot(1,4,4);
plot(unwrap(phase_diff(row,:)))

xlabel(['along-track [' num2str(N_avg) '-ml]']); ylabel(['unwrap(phase\_diff(' num2str(row) '))'])

%%

function Avg = ApplyAveraging(FLD,N_avg,IDX)
    MAvg = movmean(FLD,N_avg,2);
    Avg  = MAvg(:,IDX);
end
