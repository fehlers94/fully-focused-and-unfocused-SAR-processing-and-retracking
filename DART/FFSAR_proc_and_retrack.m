function [d] = FFSAR_proc_and_retrack(scenario, T, posting_rate, ignore_cashed_files, skip_if_exists, reverse_file_order)

% clear all
% close all
% clc

LoadCommonSettings

% defaultl dirs
tudtum_dir = [getenv('HOME') '/TUDTUM/'];
tudtum_data_dir = [tudtum_dir '/ffsar_data/'];
l1b_dest_dir = [getenv('HOME') '/tmp/'];

%% Settings
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;

% proc-and-retrack_settings default params
proc_retrack_sets = struct;
proc_retrack_sets.collocate_with_l2 = true;
proc_retrack_sets.target_posting_rate = 20;
proc_retrack_sets.enable_l2_proc = true;
proc_retrack_sets.fit_zero_doppler = false;
defval('ignore_cashed_files',false)
proc_retrack_sets.ignore_cashed_files = ignore_cashed_files;

defval('skip_if_exists',false)
defval('reverse_file_order',false)

%% scenarios
if ~exist('scenario')
    scenario = 'default_S6A';
    % scenario = 'default_S3A';
%     scenario = 'coastal_ffsar';
    % scenario = 'inland_proc_campaign';
    
    show_plots = true;
else
    show_plots = false;
end

dataset_l1a_b_l2 = struct([]);

if contains(scenario, 'coastal_ffsar')
    if exist('T')
        FFSAR_processing_settings.integration_time = T;
    else
%         FFSAR_processing_settings.integration_time = 3.413;
        FFSAR_processing_settings.integration_time = 2.1;
    end
    FFSAR_processing_settings.do_correct_IRF = false;
    FFSAR_processing_settings.combine_n_look_locations = 50;
    FFSAR_processing_settings.output_pseudo_delay_doppler_processing = true;
    FFSAR_processing_settings.num_coherent_bursts = 1;
    FFSAR_processing_settings.simulate_range_walk = false;
    FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant';  % pchip_smooth
    FFSAR_processing_settings.split_aperture = true;
    
    proc_retrack_sets.collocate_with_l2 = true;
    
    if exist('posting_rate')
        proc_retrack_sets.target_posting_rate = posting_rate;
    else
        proc_retrack_sets.target_posting_rate = 160;
    end
    proc_retrack_sets.enable_l2_proc = false;
    proc_retrack_sets.fit_zero_doppler = ~FFSAR_processing_settings.do_correct_IRF;
    
    coastal_ffsar_dir = '/lfs/DGFI24/bigdata/coastal_ffsar';
%     coastal_ffsar_dir = '/nfs/DGFI145/C/work_flo/coastal_ffsar';
    orig_data_dir = [coastal_ffsar_dir '/orig_data/f06/'];
    dest_data_dir = [coastal_ffsar_dir '/processed/'];
    
    l1b_dest_dir = [dest_data_dir 'l1b/T_' strrep(num2str(FFSAR_processing_settings.integration_time),'.','-') '_posting_rate_' num2str(proc_retrack_sets.target_posting_rate) '/'];
    l2_dest_dir = [dest_data_dir 'l2/T_' strrep(num2str(FFSAR_processing_settings.integration_time),'.','-') '_posting_rate_' num2str(proc_retrack_sets.target_posting_rate) '/'];
      
    doms = {};
    
    cycle_minmax = [5, 42];
%     cycle_minmax = [40,40];
    pass_all = [213, 18, 196, 120, 44];
%     pass_all = [120];

    doms_all = containers.Map([213, 18, 196, 120, 44],{[53.73 53.95; 0, 360], [53.65 53.88; 0, 360], [53.13 53.375; 0, 360], [51.794 52.05; 0, 360], [50.99 51.406; 0, 360]});  % until d2c=30km
    
    ncs_l1a = {}; ncs_l1b = {}; ncs_l2 = {}; cycles = {}; passes = {};

    for i=1:numel(pass_all)
        sel_pass = pass_all(i);
        DOM = doms_all(sel_pass);
    
        [all_ncs_l1a, all_cycle_pass] = get_s6_files_from_dir([orig_data_dir 'P4_1A_HR_____/'], 0, sel_pass);
        all_ncs_l1b = get_s6_files_from_dir([orig_data_dir 'P4_1B_HR_____/'], 0, sel_pass);
        all_ncs_l2 = get_s6_files_from_dir([orig_data_dir 'P4_2__HR_____/'], 0, sel_pass);
        
        if length(all_ncs_l1a) ~= length(all_ncs_l2)
            error('number L1a and L2 is not the same!')
        end

        ncs_model = {};
        for i=1:numel(all_ncs_l1a)
            cycle = all_cycle_pass{i}(1);
            pass = all_cycle_pass{i}(2);
            if cycle >= cycle_minmax(1) && cycle <= cycle_minmax(2)
                ncs_l1a{end+1} = all_ncs_l1a{i};
                ncs_l1b{end+1} = all_ncs_l1b{i};
                ncs_l2{end+1} = all_ncs_l2{i};
                doms{end+1} = DOM;
                cycles{end+1} = cycle;
                passes{end+1} = pass;
            end
        end
    end

    if reverse_file_order
        dataset_l1a_b_l2 = struct('L1A', flip(ncs_l1a), 'L1B', flip(ncs_l1b), 'L2', flip(ncs_l2), 'DOM', flip(doms), 'cycle', flip(cycles), 'pass', flip(passes));
    else
        dataset_l1a_b_l2 = struct('L1A', ncs_l1a, 'L1B', ncs_l1b, 'L2', ncs_l2, 'DOM', doms, 'cycle', cycles, 'pass', passes);
    end


elseif contains(scenario, 'default_S6A')
%     FFSAR_processing_settings.integration_time = 2;
    FFSAR_processing_settings.integration_time = 2.1;
%     FFSAR_processing_settings.integration_time = 3.413;
    
    FFSAR_processing_settings.do_correct_IRF = false;
    FFSAR_processing_settings.combine_n_look_locations = 50;
    FFSAR_processing_settings.output_pseudo_delay_doppler_processing = true;
    FFSAR_processing_settings.num_coherent_bursts = 1;
    FFSAR_processing_settings.simulate_range_walk = false;
    FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant';  % pchip_smooth
    FFSAR_processing_settings.split_aperture = false;
    
    proc_retrack_sets.fit_zero_doppler = ~FFSAR_processing_settings.do_correct_IRF;
    proc_retrack_sets.collocate_with_l2 = true;
    proc_retrack_sets.target_posting_rate = 20;
    
    nc_base_dir = [tudtum_data_dir 's6a/l1a-l1b-l2/crete/']; 
%     dataset_l1a_b_l2(1).L1A = [nc_base_dir 'S6A_P4_1A_HR______20210901T212456_20210901T222027_20210902T134010_3331_030_018_009_EUM__OPE_ST_F03.nc'];
%     dataset_l1a_b_l2(1).L1B = [nc_base_dir 'S6A_P4_1B_HR______20210901T212456_20210901T222025_20210902T134009_3329_030_018_009_EUM__OPE_ST_F03.nc'];
%     dataset_l1a_b_l2(1).L2 = [nc_base_dir 'S6A_P4_2__HR_STD__ST_030_018_20210901T212456_20210901T222025_F03.nc'];
    
%     dataset_l1a_b_l2(1).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211120T051222_20211120T060838_20211121T100952_3376_038_018_009_EUM__OPE_ST_F04.nc'];
%     dataset_l1a_b_l2(1).L2 = [nc_base_dir 'S6A_P4_2__HR_STD__ST_038_018_20211120T051224_20211120T060836_F04.nc'];
%     dataset_l1a_b_l2(1).L1B = [nc_base_dir 'S6A_P4_1B_HR______20211120T051224_20211120T060836_20211121T100951_3372_038_018_009_EUM__OPE_ST_F04.nc'];
    
    dataset_l1a_b_l2(1).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211120T051222_20211120T060838_20220430T212624_3376_038_018_009_EUM__REP_NT_F06.nc'];
    dataset_l1a_b_l2(1).L2 = [nc_base_dir 'S6A_P4_2__HR_STD__NT_038_018_20211120T051224_20211120T060836_F06.nc'];
    dataset_l1a_b_l2(1).L1B = [nc_base_dir 'S6A_P4_1B_HR______20211120T051224_20211120T060836_20220430T212619_3372_038_018_009_EUM__REP_NT_F06.nc'];
    
%     dataset_l1a_b_l2(1).DOM = [33.59 33.60; 0, 360]; % open ocean
%     dataset_l1a_b_l2(1).DOM = [33.55 33.60; 0, 360]; % open ocean
%     dataset_l1a_b_l2(1).DOM = [-33.60 -33.55; 0, 360]; % open ocean
%     dataset_l1a_b_l2(1).DOM = [36.35 36.54; 0, 360]; % coast

%     dataset_l1a_b_l2(1).DOM = [33.5 33.6612; 0, 360]; % open ocean
    dataset_l1a_b_l2(1).DOM = [33.5 33.562; 0, 360]; % open ocean

elseif contains(scenario, 'default_S3A')
    nc_base_dir = [tudtum_data_dir 's3a/l1a-l1b-l2_sets/']; 
   
    dataset_l1a_b_l2(1).L1A = [nc_base_dir 'l1a/S3A_SR_1_SRA_A__20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement_l1a.nc'];
    dataset_l1a_b_l2(1).L1B = [nc_base_dir 'l1b/S3A_SR_1_SRA____20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement.nc'];
    dataset_l1a_b_l2(1).L2 = [nc_base_dir 'l2/S3A_SR_2_WAT____20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004_concat.nc'];

    dataset_l1a_b_l2(1).DOM = [33.45 33.60; 0, 360];
elseif contains(scenario, 'inland_proc_campaign')
    FFSAR_processing_settings.integration_time = 2;
    FFSAR_default_processing_settings.focal_point_height_interpolation = 'piecewise_constant';
    
    % proc-and-retrack_settings
    proc_retrack_sets = struct;
    proc_retrack_sets.collocate_with_l2 = false;
    proc_retrack_sets.target_posting_rate = 640;  % if collocate_with_l2 is set to false
    proc_retrack_sets.enable_l2_proc = false;
    
%     nc_base_dir = [data_dir 's6a/l1a/'];

%     dataset_l1a_b_l2(1).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211203T181851_20211203T191508_20211204T141044_3377_039_111_055_EUM__OPE_ST_F04.SEN6/measurement.nc'];
%     dataset_l1a_b_l2(1).L1B = [nc_base_dir 'l1b/S3A_SR_1_SRA____20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement.nc'];
%     dataset_l1a_b_l2(1).L2 = [nc_base_dir 'l2/S3A_SR_2_WAT____20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004_concat.nc'];
%     dataset_l1a_b_l2(1).DOM = [53.19 53.68; 12.28, 13.04];

%     dataset_l1a_b_l2(1).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211202T035403_20211202T045020_20211203T084909_3377_039_070_035_EUM__OPE_ST_F04.SEN6/measurement.nc'];
%     dataset_l1a_b_l2(1).DOM = [44.40 44.70; 0, 359];
    
    l1b_dest_dir = '/nfs/DGFI145/C/work_flo/ffsar_inland_proc/s6_processed_2/';
%
%     nc_base_dir = '/nfs/DGFI8/H/work_laura/FFSARdata/S6data/';
%     dataset_l1a_b_l2(1).mission = 'S6A';
%     dataset_l1a_b_l2(1).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211207T015706_20211207T025323_20211208T085426_3377_039_196_098_EUM__OPE_ST_F04.SEN6/measurement.nc'];
%     dataset_l1a_b_l2(1).DOM = [47.49 47.75; 0, 359];
%     
%     dataset_l1a_b_l2(2).mission = 'S6A';
%     dataset_l1a_b_l2(2).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211129T184457_20211129T194114_20211130T141042_3377_039_009_004_EUM__OPE_ST_F04.SEN6/measurement.nc'];
%     dataset_l1a_b_l2(2).DOM = [47.49 47.75; 0, 359];
%     
%     dataset_l1a_b_l2(3).mission = 'S6A';
%     dataset_l1a_b_l2(3).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211207T015706_20211207T025323_20211208T085426_3377_039_196_098_EUM__OPE_ST_F04.SEN6/measurement.nc'];
%     dataset_l1a_b_l2(3).DOM = [51.52 51.92; 0, 359];
%     
%     dataset_l1a_b_l2(4).mission = 'S6A';
%     dataset_l1a_b_l2(4).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211203T181851_20211203T191508_20211204T141044_3377_039_111_055_EUM__OPE_ST_F04.SEN6/measurement.nc'];
%     dataset_l1a_b_l2(4).DOM = [53.19 53.68; 0, 359];
%     
%     dataset_l1a_b_l2(5).mission = 'S6A';
%     dataset_l1a_b_l2(5).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211206T013531_20211206T023148_20211207T085346_3377_039_170_085_EUM__OPE_ST_F04.SEN6/measurement.nc'];
%     dataset_l1a_b_l2(5).DOM = [53.19 53.68; 0, 359];
% 
%     dataset_l1a_b_l2(6).mission = 'S6A';
%     dataset_l1a_b_l2(6).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211129T184457_20211129T194114_20211130T141042_3377_039_009_004_EUM__OPE_ST_F04.SEN6/measurement.nc'];
%     dataset_l1a_b_l2(6).DOM = [48.78 49.01; 0, 359];
% 
%     dataset_l1a_b_l2(7).mission = 'S6A';
%     dataset_l1a_b_l2(7).L1A = [nc_base_dir 'S6A_P4_1A_HR______20211130T031053_20211130T040710_20211201T085046_3377_039_018_009_EUM__OPE_ST_F04.SEN6/measurement.nc'];
%     dataset_l1a_b_l2(7).DOM = [52.33 52.70; 0, 359];
    
    % second processing
    nc_base_dir = '/DGFI8/A/original_data/sentinel6a/JASON_CS_S6A_L1A_ALT_HR_NTC_F/';
    
    dataset_l1a_b_l2(1).DOM = [48.78 49.01; 0, 359];
    
    dataset_l1a_b_l2 = repelem(dataset_l1a_b_l2,11);
    
%     dataset_l1a_b_l2(1).L1A = [nc_base_dir '027/S6A_P4_1A_HR______20210802T190240_20210802T195412_20210822T000518_3092_027_009_004_EUM__OPE_NT_F03.SEN6.measurement.nc'];
%     dataset_l1a_b_l2(2).L1A = [nc_base_dir '028/S6A_P4_1A_HR______20210812T170111_20210812T175727_20210903T011031_3376_028_009_004_EUM__OPE_NT_F03.SEN6.measurement.nc'];
%     dataset_l1a_b_l2(3).L1A = [nc_base_dir '029/S6A_P4_1A_HR______20210822T145942_20210822T154957_20210913T170503_3015_029_009_004_EUM__OPE_NT_F03.SEN6.measurement.nc']; % stripped inland meas
    dataset_l1a_b_l2(4).L1A = [nc_base_dir '030/S6A_P4_1A_HR______20210901T125814_20210901T134829_20210922T182709_3015_030_009_004_EUM__OPE_NT_F03.SEN6.measurement.nc']; 
    dataset_l1a_b_l2(5).L1A = [nc_base_dir '032/S6A_P4_1A_HR______20210921T085517_20210921T095134_20211010T102941_3377_032_009_004_EUM__OPE_NT_F03.SEN6.measurement.nc'];
    dataset_l1a_b_l2(6).L1A = [nc_base_dir '033/S6A_P4_1A_HR______20211001T065349_20211001T075006_20211029T012108_3377_033_009_004_EUM__OPE_NT_F03.SEN6.measurement.nc'];
    dataset_l1a_b_l2(7).L1A = [nc_base_dir '034/S6A_P4_1A_HR______20211011T045220_20211011T054837_20211108T011838_3377_034_009_004_EUM__OPE_NT_F03.SEN6.measurement.nc'];
    dataset_l1a_b_l2(8).L1A = [nc_base_dir '035/S6A_P4_1A_HR______20211021T025051_20211021T034708_20211115T175707_3377_035_009_004_EUM__OPE_NT_F04.SEN6.measurement.nc'];
    dataset_l1a_b_l2(9).L1A = [nc_base_dir '036/S6A_P4_1A_HR______20211031T004923_20211031T014540_20211119T081300_3377_036_009_004_EUM__OPE_NT_F04.SEN6.measurement.nc'];
    dataset_l1a_b_l2(10).L1A = [nc_base_dir '037/S6A_P4_1A_HR______20211109T224754_20211109T234411_20211130T064335_3377_037_009_004_EUM__OPE_NT_F04.SEN6.measurement.nc'];
    dataset_l1a_b_l2(11).L1A = [nc_base_dir '038/S6A_P4_1A_HR______20211119T204626_20211119T214242_20211210T055410_3376_038_009_004_EUM__OPE_NT_F04.SEN6.measurement.nc'];

end

%% L1b-L2-processing of each dataset
for i = 1:numel(dataset_l1a_b_l2)
    dataset = dataset_l1a_b_l2(i);
       
    is_l1a_set = isfield(dataset,'L1A') && ~isempty(dataset);
    is_l2_set = isfield(dataset,'L2') && ~isempty(dataset);
%     if ~isfield(dataset_l1a_b_l2, 'L1A') | isempty(dataset.L1A)
%         disp('no L1A file set, skipping dataset entry...')
%         continue;
%     end
   
    %get l1a fileinfo
    if is_l1a_set
        [fpath_l1a, basename_l1a, ~] = fileparts(dataset.L1A);
        if contains(basename_l1a, 'measurement')
            [~, basename_l1a] = fileparts(fpath_l1a);
        end
    end
    
    %get l2 fileinfo
    if is_l2_set
        [fpath_l2, basename_l2, ~] = fileparts(dataset.L2);
        if contains(basename_l2, 'measurement')
            [~, basename_l2] = fileparts(fpath_l2);
        end
    end
    
    if exist('basename_l1a')
        basename = basename_l1a;
    elseif exist('basename_l2')
        basename = basename_l2;
    end
    
	lats_min_max = sort(dataset.DOM(1,:));
    l1b_destfile = [fullfile(l1b_dest_dir, replace(replace(basename, '1A', '1B'), 'EUM', 'DAR')) '_' num2str(lats_min_max(1),'%.2f') '_' num2str(lats_min_max(2),'%.2f') '.nc'];
    
    if skip_if_exists && isfile(l1b_destfile)
        fprintf('%s already exists, skipping...\n', basename)
        continue;
    end
    
    mission = mission_from_fname(basename);
    DDAcf = DDA_ConfigFile(mission,'SAR');

    fprintf(['Processing file ',num2str(i),' of ',num2str(numel(dataset_l1a_b_l2)),': ' basename '\n']);
    DOM = dataset.DOM;
    
    % ff-sar-process EUM L1a->L1b
    if is_l1a_set
        [L2_ff, L1b_ff] = FFSAR_proc_and_retrack_l2target(dataset, FFSAR_processing_settings, proc_retrack_sets);
    end
        
    %% read baseline L2
    if isfield(dataset,'L2')
        l2_eum_lat_all = ncread(dataset.L2,[DDAcf.ncgrp_l1b_l2 DDAcf.var_mapping_l2.lat]);
        l2_eum_lon_all = ncread(dataset.L2,[DDAcf.ncgrp_l1b_l2 DDAcf.var_mapping_l2.lat]);
        eps_lat = 0.002;
        IDX_mask = ingeoquad(l2_eum_lat_all,l2_eum_lon_all,dataset.DOM(1,:)+[-eps_lat,eps_lat],dataset.DOM(2,:));
        IDX = find(IDX_mask);
        startLoc=min(IDX);
        count=max(IDX)-startLoc+1;

        L2_eum.LAT = ncread(dataset.L2,[DDAcf.ncgrp_l1b_l2 DDAcf.var_mapping_l2.lat],startLoc,count);
        L2_eum.LON = ncread(dataset.L2,[DDAcf.ncgrp_l1b_l2 DDAcf.var_mapping_l2.lon],startLoc,count);
        L2_eum.SWH = ncread(dataset.L2,[DDAcf.ncgrp_l1b_l2 DDAcf.var_mapping_l2.swh],startLoc,count);
        alt_l2_eum = ncread(dataset.L2,[DDAcf.ncgrp_l1b_l2 DDAcf.var_mapping_l2.alt],startLoc,count);
        range_l2_eum = ncread(dataset.L2,[DDAcf.ncgrp_l1b_l2 DDAcf.var_mapping_l2.range],startLoc,count);
        L2_eum.HEI = alt_l2_eum - range_l2_eum;
        L2_eum.TIME = ncread(dataset.L2,[DDAcf.ncgrp_l1b_l2 DDAcf.var_mapping_l2.time],startLoc,count);
        L2_eum.ALT = ncread(dataset.L2,[DDAcf.ncgrp_l1b_l2 DDAcf.var_mapping_l2.alt],startLoc,count)';
        L2_eum.sigma0 = ncread(dataset.L2,[DDAcf.ncgrp_l1b_l2 DDAcf.var_mapping_l2.alt],startLoc,count);

        % reverse cog_cor in case of S3
        if contains(mission, 'S3')
            L2_eum.HEI = L2_eum.HEI + 0.55590;
        end
    end
    
    %% read model
    if isfield(dataset,'MODEL')
        model_lat_all = ncread(dataset.MODEL, DDAcf.var_mapping_kuststrook.lat);
        model_lon_all = ncread(dataset.MODEL, DDAcf.var_mapping_kuststrook.lon);
        IDX_mask = ingeoquad(model_lat_all,model_lon_all,dataset.DOM(1,:),dataset.DOM(2,:));
        IDX = find(IDX_mask);
        startLoc = min(IDX);
        count = max(IDX)-startLoc+1;
        
        model.SWH = ncread(dataset.MODEL,DDAcf.var_mapping_kuststrook.swh, startLoc,count);
%         model.HEI = ncread(dataset.MODEL,DDAcf.var_mapping_kuststrook.range, startLoc,count);
        model.HEI = ncread(dataset.MODEL,'waterlevel', startLoc,count);
    end
    
    %% retrack baseline L1b product
    if proc_retrack_sets.enable_l2_proc && isfield(dataset,'L1B')
        [L2_retracked,CS1b_eum] = SAR_L1b_to_L2(mission,dataset.L1B,DOM,'SAMOSA2',...
        {'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false'});
    end

    %% plotting retracked results vs. L2 file
    if proc_retrack_sets.enable_l2_proc && show_plots
        figure;
        fontsize_legends = 8;
        fontsize_title = 8;

        if exist('L2_ff')
            inds_ff = linspace(L2_retracked.IDXpoi(1), L2_retracked.IDXpoi(end), length(L2_ff.TIME));
        end

        % SWH
        ax1 = subplot(2,1,1);
        plot(L2_retracked.IDXpoi, L2_retracked.SWH); hold on;
        plot(L2_retracked.IDXpoi, L2_eum.SWH); hold on;
        legend_entries = {['EUM L1b->L2, std=' num2str(std(L2_retracked.SWH)) ', median=' num2str(median(L2_retracked.SWH)) ', median\_bias_{EUM}=' num2str(median(L2_retracked.SWH - L2_eum.SWH)) ', RMSE_{EUM}=' num2str(sqrt(mean(L2_retracked.SWH - L2_eum.SWH)^2))],...
                          ['EUM L2, std=' num2str(std(L2_eum.SWH)) ', median=' num2str(median(L2_eum.SWH))]      
            };
        
        if exist('L2_ff')
            try
                legend_entries{end+1} = ['FF-SAR, std=' num2str(std(L2_ff.SWH)) ', median=' num2str(median(L2_ff.SWH)) ', median\_bias_{EUM}=' num2str(median(L2_ff.SWH - L2_eum.SWH)) ', RMSE_{EUM}=' num2str(sqrt(mean(L2_ff.SWH - L2_eum.SWH)^2))];
            catch
                legend_entries{end+1} = ['FF-SAR, std=' num2str(std(L2_ff.SWH)) ', median=' num2str(median(L2_ff.SWH))];
            end
            
            plot(inds_ff, L2_ff.SWH); hold on;
            
        end
        if isfield(dataset,'MODEL')
            plot(inds_ff, model.SWH); hold on;
            legend_entries{end+1} = ['Model, std=' num2str(std(model.SWH)) ', median=' num2str(median(model.SWH))];
        end      
     
        if FFSAR_processing_settings.output_pseudo_delay_doppler_processing
            L1b_uf = L1b_ff;
            L1b_uf.SAR.data = L1b_ff.SAR.data_pseudoDD;
            [L2_uf,~ ] = SAR_L1b_to_L2(mission,L1b_uf,dataset.DOM,'SAMOSA2',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false';});
        
            plot(L2_retracked.IDXpoi(1:length(L2_uf.SWH)), L2_uf.SWH); hold on;
            legend_entries{end+1} = ['UF-SAR, std=' num2str(std(L2_uf.SWH)) ', median=' num2str(median(L2_uf.SWH))];
        end
        
        legend(legend_entries, 'FontSize',fontsize_legends,'Location', 'southeast')
        
        grid on;
        xlabel('record index i');
        ylabel('SWH [m]');

        % range
        ax2 = subplot(2,1,2);
        plot(L2_retracked.IDXpoi, L2_retracked.HEI,'-');hold on;
        plot(L2_retracked.IDXpoi, L2_eum.HEI,'-');hold on;
        legend_entries = {['EUM L1b->L2, RMSE_{EUM}=' num2str(sqrt(mean(L2_retracked.HEI - L2_retracked.HEI)^2)) ', detrended std=',num2str(std(detrend(L2_retracked.HEI))),' m'],...
                          ['EUM L2, detrended std=',num2str(std(detrend(L2_eum.HEI))),' m']};
        
        if exist('L2_ff')
            plot(inds_ff, L2_ff.HEI,'-'); hold on
            legend_entries{end+1} = ['FF-SAR, detrended std=',num2str(std(detrend(L2_ff.HEI))),' m'];
        end
%         if isfield(dataset,'MODEL')
%             plot(inds_ff, model.HEI); hold on;
%             legend_entries{end+1} = ['Model, std=' num2str(std(model.HEI)) ', median=' num2str(median(model.HEI))];
%         end
        
        if FFSAR_processing_settings.output_pseudo_delay_doppler_processing
            plot(inds_ff, L2_uf.HEI,'-');hold on
            legend_entries{end+1} = ['UF-SAR, detrended std=',num2str(std(detrend(L2_uf.HEI))),' m'];
        end

        legend(legend_entries, 'FontSize',fontsize_legends,'Location', 'southeast')
        
        grid on
        xlabel('record index i')
        ylabel('uncorrected SSH [m]')

        linkaxes([ax1 ax2],'x')
        
        % plot(detrend(L2_eum.HEI),detrend(L2_ff.HEI),'b.')

        % add title
        if exist('L2_ff') && isfield(L1b_ff.GEO, 'dist2coast')
            d2c_km = median(L2_ff.GEO.dist2coast) / 1000;
            sgtitle({basename, ['dist2coast: ' num2str(d2c_km) 'km']}, 'FontSize', fontsize_title);
        else
            sgtitle(basename, 'FontSize', fontsize_title);
        end
    end
    
    %% plot echo radargram single
    if exist('L1b_ff') && show_plots && false
        max_range_gate = 150;
        figure;
        imagesc(L1b_ff.IDXpoi, (1:max_range_gate), log10(L1b_ff.SAR.data(1:max_range_gate,:)));

        xlabel('20-Hz record index'); ylabel('range gates [bin]');    
        colorbar; colormap(parula); caxis([8,11.8]);
        title([num2str(proc_retrack_sets.target_posting_rate) '-Hz multilooked (FFSAR)']);
    end
    
    %% plot echo radargram
    if exist('L1b_ff') && show_plots && false
        inds_20 = linspace(L1b_ff.IDXpoi(1)+0.5, L1b_ff.IDXpoi(end)+0.5, length(L1b_ff.IDXpoi));
        inds_20_intp = linspace(L1b_ff.IDXpoi(1)-0.5, L1b_ff.IDXpoi(end)-0.5, size(L1b_ff.non_average.data,2));
        
        max_range_gate = 150;
    %     max_range_gate = size(L1b_ff.SAR.data,1);
  
        figure;
        ax1 = subplot(2,2,1);
    %     imagesc(L1b_ff.GEO.LAT, (1:max_range_gate), log10(L1b_ff.non_average.data(1:max_range_gate,:)));
        imagesc(inds_20_intp, (1:max_range_gate), log10(L1b_ff.non_average.data_pseudoDD(1:max_range_gate,:)));
        xline(inds_20, '-k', 'LineWidth', 1)

        xlabel('20-Hz record index'); ylabel('range gates [bin]');
        colorbar; caxis([8,11.8]); colormap(parula);
        title('single-looks (UFSAR)');

        ax2 = subplot(2,2,2);
    %     imagesc(L1b_ff.GEO.LAT, (1:max_range_gate), log10(L1b_ff.SAR.data(1:max_range_gate,:)));
        imagesc(L1b_ff.IDXpoi, (1:max_range_gate), log10(L1b_ff.SAR.data_pseudoDD(1:max_range_gate,:)));
        xline(inds_20, '-k', 'LineWidth', 1);

        xlabel('20-Hz record index'); ylabel('range gates [bin]');    
        colorbar; colormap(parula); caxis([8,11.8]);
        title([num2str(proc_retrack_sets.target_posting_rate) '-Hz multilooked (UFSAR)']);

        ax3 = subplot(2,2,3);
        imagesc(inds_20_intp, (1:max_range_gate), log10(L1b_ff.non_average.data(1:max_range_gate,:)));
        xline(inds_20, '-k', 'LineWidth', 1)

        xlabel('20-Hz record index'); ylabel('range gates [bin]');
        colorbar; caxis([8,11.8]); colormap(parula);
        title('single-looks (FFSAR)');

        ax4 = subplot(2,2,4);
        imagesc(L1b_ff.IDXpoi, (1:max_range_gate), log10(L1b_ff.SAR.data(1:max_range_gate,:)));
        xline(inds_20, '-k', 'LineWidth', 1);

        xlabel('20-Hz record index'); ylabel('range gates [bin]');    
        colorbar; colormap(parula); caxis([8,11.8]);
        title([num2str(proc_retrack_sets.target_posting_rate) '-Hz multilooked (FFSAR)']);
        
        linkaxes([ax1,ax2,ax3,ax4],'xy');
    end
    
    %% single-look radargram (plot for paper)
    if exist('L1b_ff') && show_plots && false
        fig_h = figure;
        posting_rate_reduce_factor = int32(proc_retrack_sets.target_posting_rate / 20);
        last_n_inds = 150;
        total_len = size(L1b_ff.SAR.data,2);
        lat = L1b_ff.GEO.LAT;
        lon = L1b_ff.GEO.LAT;
        sel_lat = L1b_ff.GEO.LAT(end-last_n_inds:end);
        sel_lon = L1b_ff.GEO.LAT(end-last_n_inds:end);
        
        lat_interferer_centre = 51.7937;
        mean_distance = mean(diff(1e3*deg2km(distance(lat,lon,lat(1),lon(1)))*8));
        
        [d,pos_idx] = min(abs(lat-lat_interferer_centre));
%         x_vec = L1b_ff.GEO.LAT(end-last_n_inds:end);
%         x_vec = (total_len-last_n_inds:total_len);
        x_vec = 1e3*deg2km(distance(sel_lat,sel_lon,lat(pos_idx),lon(pos_idx)));
        x_vec(find(x_vec == 0):end) = -x_vec(find(x_vec == 0):end);
%         x_vec = x_vec(end-last_n_inds:end);
%           x_vec = (total_len-last_n_inds:total_len);
%         sel_radargram = L1b_ff.non_average.data_pseudoDD(1:max_range_gate,single_look_offset:end);
%         sel_radargram = L1b_ff.non_average.data(1:max_range_gate,single_look_offset:end);
%         sel_radargram = L1b_ff.non_average.data_opt(1:max_range_gate,single_look_offset:end);
%         sel_radargram = L1b_ff.non_average.data_opt(1:max_range_gate,single_look_offset:end);

%         ax1 = subplot(2,1,1);
        ha = tight_subplot(2,2,.0,[.1 .0],[.1 .05]);
        
        axes(ha(1))
        sel_radargram = L1b_ff.SAR.data(1:max_range_gate,end-last_n_inds:end);
%         imagesc((total_len-last_n_inds:total_len), (1:max_range_gate), 10*log10(sel_radargram));
        imagesc(x_vec, (1:max_range_gate), 10*log10(sel_radargram));
        ylabel('range gates [bin]');
        c = colorbar('Location', 'northoutside');
        caxis([80,118]);
        colormap(parula);
        annotation('textbox', [0.13 0.55 0.3 0.3],'String', '(a)','FitBoxToText','on', 'color', 'white', 'EdgeColor','none', 'FontSize', 8);
        set(gca,'xtick',[])

        for i=[0:3]
            xline(i*mean_distance,'color','r')
        end
        
%         ax2 = subplot(2,1,2);
        axes(ha(3))
        plot((total_len-last_n_inds:total_len), 10*log10(sum(sel_radargram,1)));
        plot(x_vec, 10*log10(sum(sel_radargram,1)));
        xlabel('distance-to-interferer-centre [m]');
        ylabel('across-track power [dBm]');
        grid on;
        annotation('textbox', [0.13 0.25 0.3 0.3],'String', '(c)','FitBoxToText','on', 'color', 'black', 'EdgeColor','none', 'FontSize', 8);
        xlim([min(x_vec) max(x_vec)]);

        for i=[0:3]
            xline(i*mean_distance,'color','r')
        end
        
        axes(ha(2))
        sel_radargram = L1b_ff.SAR.data_opt(1:max_range_gate,end-last_n_inds:end);
%         imagesc((total_len-last_n_inds:total_len), (1:max_range_gate), 10*log10(sel_radargram));
        imagesc(x_vec, (1:max_range_gate), 10*log10(sel_radargram));
        c = colorbar('Location', 'northoutside');
        caxis([80,118]);
        colormap(parula);
        annotation('textbox', [0.55 0.55 0.3 0.3],'String', '(b)','FitBoxToText','on', 'color', 'white', 'EdgeColor','none', 'FontSize', 8);
        set(gca,'xtick',[]); set(gca,'ytick',[]);

        for i=[0:3]
            xline(i*mean_distance,'color','r')
        end
        
%         ax2 = subplot(2,1,2);
        axes(ha(4))
%         plot((total_len-last_n_inds:total_len), 10*log10(sum(sel_radargram,1)));
        plot(x_vec, 10*log10(sum(sel_radargram,1)));
        grid on;
        annotation('textbox', [0.55 0.25 0.3 0.3],'String', '(d)','FitBoxToText','on', 'color', 'black', 'EdgeColor','none', 'FontSize', 8);
        xlabel('distance-to-interferer-centre [m]');
        xlim([min(x_vec) max(x_vec)]);
        set(gca,'ytick',[]);

        for i=[0:3]
            xline(i*mean_distance,'color','r')
        end
        
%         linkaxes([ax1,ax2],'x')
        linkaxes([ha(1),ha(2),ha(3),ha(4)],'x')   
        
        % export plot data to .mat file (for plottig in Python)
        if dataset.cycle == 40 && dataset.pass == 120
            data_export = struct;
            data_export.data = L1b_ff.SAR.data;
            data_export.data_pseudoDD = L1b_ff.SAR.data_pseudoDD;
            data_export.LAT = L1b_ff.GEO.LAT;
            data_export.LON = L1b_ff.GEO.LON;

            save([orig_data_dir '/sandbank_interferer_l1b_data.mat'], 'data_export');
        end
    end
   
    %%
%     [max_vals, max_inds] = max(L1b_ff.non_average.data, [], 1);
%     figure;
% %     plot(inds_20_intp, movmean(max_inds,L1b_ff.N_avg))
%     plot(inds_20_intp, movmean(max_inds,10))
% %     plot(inds_20_intp, max_inds)
    
%     %%
%     data_norm = L1b_ff.non_average.data ./ max(L1b_ff.non_average.data);
%     data_movprod = movprod(data_norm, 30, 2);
%     [~, max_inds_movprod] = max(data_movprod, [], 1);
%     figure;
%     plot(inds_20_intp, max_inds_movprod)
% 
%     % get bad sl_inds
%     n_inds_range = 30;
%     mask_inds_bad = max_inds_movprod < (mean(max_inds_movprod) - n_inds_range) | (max_inds_movprod > (mean(max_inds_movprod) + n_inds_range));
%     
%     hold on;
%     plot(inds_20_intp(mask_inds_bad), max_inds_movprod(mask_inds_bad), '.r')
%     
%     %cleanse waveform
%     data_cleansed = L1b_ff.non_average.data;
%     data_cleansed(:,mask_inds_bad) = 0;
%     
%     %average cleansed waveform
%     MAvg = movmean(data_cleansed,L1b_ff.non_average.N_avg,2);
%     data_cleansed_avg = MAvg(:,L1b_ff.non_average.closestIndex);
%     
%     % plot cleansed single-look radargram
%     figure;
%     imagesc(inds_20, (1:max_range_gate), log10(data_cleansed_avg(1:max_range_gate,:)));
%     xlabel('20-Hz record index')
%     xline(inds_20, '-k', 'LineWidth', 1)
% 
%     colorbar;
%     caxis([8,11.8]);
%     ylabel('range gates [bin]');
%     colormap(parula)
%     
    %% export to NC file
    % L1b
    if ~isdir(l1b_dest_dir); mkdir(l1b_dest_dir); end
    
    export_L1b_to_nc(L1b_ff, l1b_destfile, basename);
%     export_L1b_to_nc(L1b_uf, l1b_destfile, basename);
    
    % L2
    if proc_retrack_sets.enable_l2_proc
        l2_destfile = [fullfile(l1b_dest_dir, replace(replace(basename, '1A', '2_'), 'EUM', 'DAR')) '_' num2str(lats_min_max(1),'%.2f') '_' num2str(lats_min_max(2),'%.2f') '.nc'];
%         export_L2_to_nc(L1b_ff, L2_ff, l2_destfile, basename);
%         export_L2_to_nc(L1b_uf, L2_uf, l2_destfile, basename);
        export_L2_to_nc(L1b_uf, L2_eum, l2_destfile, basename);
    end
end

end