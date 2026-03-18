tudtum_data_dir = fullfile([getenv('HOME') '/TUDTUM/ffsar_data/']);
matlab_data_export_dir = fullfile([getenv('HOME') '/MEGA/thesis_matlab_export/']);

scenarios = {...
    {[tudtum_data_dir 's3a/l1a/crete/S3A_SR_1_SRA_A__20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement_l1a.nc'] 'tr_phase_s3' 'crete' [35.33 35.345; 23 24]},
    {[tudtum_data_dir 's3a/l1a/crete/S3A_SR_1_SRA_A__20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement_l1a.nc'] 'tr_radargram_s3' 'crete' [35.32 35.355; 23 24] },
    {[tudtum_data_dir 's6a/l1a-l1b-l2/crete/S6A_P4_1A_HR______20211120T051222_20211120T060838_20211121T100952_3376_038_018_009_EUM__OPE_ST_F04.nc'] 'tr_phase_s6' 'crete' [35.33 35.345; 23 24]},
    {[tudtum_data_dir 's6a/l1a-l1b-l2/crete/S6A_P4_1A_HR______20211120T051222_20211120T060838_20211121T100952_3376_038_018_009_EUM__OPE_ST_F04.nc'] 'tr_radargram_s6' 'crete' [35.33 35.345; 23 24] },
%     {'/nfs/DGFI145/C/work_flo/coastal_ffsar/orig_data/f06/P4_1A_HR_____/S6A_P4_1A_HR______20211127T195413_20211127T205030_20220509T170629_3377_038_213_106_EUM__REP_NT_F06.SEN6/measurement.nc' 'coastal_ffsar_radargram' '' [53.73 53.83; 0 360]},
    };

for i=1:length(scenarios)
    scen = scenarios{i};
    FName = scen{1};
    name_run = scen{2};
    name_run_ext = scen{3};
    DOM = scen{4};
    
    [CONST,TR,FFSAR_default_processing_settings] = FFSAR_LoadCommonSettings;
    FFSAR_default_processing_settings.output_pseudo_delay_doppler_processing = 1;

    is_s6 = contains(FName, 'S6A');

    % process only look location that is closest to the transponder
    proc_sets = FFSAR_default_processing_settings;
    proc_sets.integration_time = 2.1;
        
    % transponder scenario
    if contains(name_run,'tr')
        proc_sets.transponder_test_mode = contains(name_run,'tr_phase');
        tr_latlon = TR.(name_run_ext);
        proc_sets.transponder_latlonalt = [tr_latlon.lat, tr_latlon.lon, tr_latlon.h];
        proc_sets.split_aperture = false;
    end

    % init and set FFSAR_Processor object
    ffproc = FFSAR_Processor(FName,DOM,proc_sets);

    % run processing
    ffproc.setup_proc()
    ffproc.proc()


    %% plots
    range_bin_max = 256;
    
    if contains(name_run, 'tr_phase') 
        uw_phase_detrended = plot_phase(ffproc);
        assert(std(uw_phase_detrended)< 2.2);
        
        save(fullfile([matlab_data_export_dir name_run '.mat']), 'uw_phase_detrended')
    elseif contains(name_run, 'tr_radargram')
        CS1b = ffproc.CS1b;
        data = CS1b.SAR.data;
        lat = CS1b.GEO.LAT;
        lon = CS1b.GEO.LON;
        Lx = ffproc.res_a;
        data_pseudoDD = CS1b.SAR.data_pseudoDD;
        r_abs_pow_ffsar = abs(CS1b.SAR.data(1:range_bin_max,:));
        r_abs_pow_ufsar = abs(CS1b.SAR.data_pseudoDD(1:range_bin_max,:));

        figure;
        ax1 = subplot(1,2,1);
        imagesc(CS1b.GEO.LAT,[],log10(r_abs_pow_ffsar));
        title('complex radargram FF-SAR');
        
        ax2 = subplot(1,2,2);
        imagesc(CS1b.GEO.LAT,[],log10(r_abs_pow_ufsar));
        title('complex radargram UF-SAR');
        
        linkaxes([ax1 ax2], 'xy')
        
%         [~,max_idx] = max(abs(r_abs_pow(:)));
%         [row_max,col_max]=ind2sub(size(abs(r_abs_pow)),max_idx);

        figure;
        subplot(1,2,1);
    %                 plot(10*log10(r_abs_pow(row_max,:)));
        plot(10*log10(sum(r_abs_pow_ffsar,1) / max(sum(r_abs_pow_ffsar,1))));
        hold on;

        grid on;
        subplot(1,2,2);
    %                 plot(10*log10(r_abs_pow(:,col_max)));
        plot(10*log10(sum(r_abs_pow_ffsar,2) / max(sum(r_abs_pow_ffsar,2))));
        
        title('across-track PTR')
        grid on;

        save(fullfile([getenv('HOME') '/MEGA/thesis_matlab_export/'  name_run '.mat']), 'data', 'data_pseudoDD','DOM', 'lat', 'lon', 'Lx')
    end

end


%% functions
function uw_phase_detrended = plot_phase(fp)
    figure;
    subplot(1,2,1)
    uwphase = fp.uwphase;

    [t_select,t_pulse] = fp.get_look_loc_times(1);
    
    P = polyfit(t_pulse,uwphase,2);
    uwphase_fit = P(1)*t_pulse.^2 + P(2)*t_pulse + P(3);
   
    lambda_mm = fp.DDAcf.c/fp.DDAcf.fc * 1000;
    plot(t_pulse,uwphase ./ (2*pi) * (lambda_mm/2),'.'); hold on
    plot(t_pulse,uwphase_fit ./ (2*pi) * (lambda_mm/2));
    axis([min(t_select)*1.1 max(t_select)*1.1 -0.5*pi 0.5*pi])
    title(fp.mission);
    xlabel('time [s]'); ylabel('unwrapped phase [mm]');
    
    uw_phase_detrended = (uwphase - (P(2)*t_pulse + P(3))) ./ (2*pi) * (lambda_mm/2);
    
    subplot(1,2,2)
    plot(t_pulse,uw_phase_detrended,'.'); hold on
    title('detrended');
    xlabel('time [s]'); ylabel('unwrapped phase [mm]');
    axis([min(t_select)*1.1 max(t_select)*1.1 -0.5*pi 0.5*pi])
end