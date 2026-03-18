classdef TESTER_FFSAR_Processor_TransponderTest < matlab.unittest.TestCase   
    properties
        data_dir = fullfile([getenv('HOME') '/TUDTUM/ffsar_data/']);
        CONST, TR, FFSAR_default_processing_settings
    end
    
%     properties (MethodSetupParameter)
%         proc_sets = FFSAR_processing_settings;
%     end
    
    properties (TestParameter)
        % Crete transponder (CDN1): (lat,lon)=(35.3379, 23.7795)
        % Svalbard transponder: (lat,lon)=(78.23052306, 15.39376997)
        tr_file_site_dom = {
            {'cs/l1a/CS_LTA__SIR1SAR_FR_20120529T025446_20120529T025502_C001.DBL' 'svalbard' [78.0 78.8; 15.0 15.4]}
            {'s3a/l1a/crete/S3A_SR_1_SRA_A__20190325T192508_20190325T201537_20191215T170531_3029_043_013______MR1_R_NT_004.SEN3/measurement_l1a.nc' 'crete' [35.33 35.345; 23 24]}
            {'s3b/l1a/crete/S3B_SR_1_SRA_A__20190620T083403_20190620T092432_20191107T073220_3029_026_335______MR1_R_NT_004.SEN3/measurement_l1a.nc' 'crete' [35.33 35.345; 23 24]}
            {'s6a/S6A_GPP__P4__FR__1A_20210116T201350_20210116T201413_0001.NC' 'crete' [35.33 35.345; 23 24]}
            {'s6a/l1a-l1b-l2/crete/S6A_P4_1A_HR______20210901T212456_20210901T222027_20210902T134010_3331_030_018_009_EUM__OPE_ST_F03.nc' 'crete' [35.33 35.345; 23 24]}
            {'s6a/l1a-l1b-l2/crete/S6A_P4_1A_HR______20211120T051222_20211120T060838_20211121T100952_3376_038_018_009_EUM__OPE_ST_F04.nc' 'crete' [35.33 35.345; 23 24]}
            };
    end
    
    methods
        function obj = TESTER_FFSAR_Processor_TransponderTest()
            [obj.CONST,obj.TR,obj.FFSAR_default_processing_settings] = FFSAR_LoadCommonSettings;
        end
    end
    
    methods (Test)
        function testFlatPhase(testCase,tr_file_site_dom)
            FName = [testCase.data_dir tr_file_site_dom{1}];
            tr_latlon = testCase.TR.(tr_file_site_dom{2});
            DOM = tr_file_site_dom{3};
            
            is_s6 = contains(FName, 'S6A');
            
            % process only look location that is closest to the transponder
            proc_sets = testCase.FFSAR_default_processing_settings;
            proc_sets.transponder_test_mode = true;
            proc_sets.transponder_latlonalt = [tr_latlon.lat, tr_latlon.lon, tr_latlon.h];
            proc_sets.split_aperture = true;
            
            if is_s6
%                 proc_sets.integration_time = 0.5;
                proc_sets.integration_time = 2;
            end
            
            % init and set FFSAR_Processor object
            ffproc = FFSAR_Processor(FName,DOM,proc_sets);
            
            % run processing
            ffproc.setup_proc()
            ffproc.proc()
           
            if ~proc_sets.transponder_test_mode
                CS1b = ffproc.CS1b;
                r_abs_pow = abs(CS1b.SAR.data);
                r_abs_pow_opt = abs(CS1b.SAR.data_opt);
                
                figure;
                ax1 = subplot(1,2,1);
                imagesc(log10(r_abs_pow));
                title('original');
                ax2 = subplot(1,2,2);
                imagesc(log10(r_abs_pow_opt));
                title('cleansed');
                linkaxes([ax1 ax2], 'xy');
                
                [~,max_idx] = max(abs(r_abs_pow(:)));
                [row_max,col_max]=ind2sub(size(abs(r_abs_pow)),max_idx);
                
                figure;
                subplot(1,2,1);
%                 plot(10*log10(r_abs_pow(row_max,:)));
                plot(10*log10(sum(r_abs_pow,1) / max(sum(r_abs_pow,1))));
                hold on;
%                 plot(10*log10(r_abs_pow_opt(row_max,:)));
                plot(10*log10(sum(r_abs_pow_opt,1) / max(sum(r_abs_pow_opt,1))));
                title('along-track PTR')
                grid on;
                subplot(1,2,2);
%                 plot(10*log10(r_abs_pow(:,col_max)));
                plot(10*log10(sum(r_abs_pow,2) / max(sum(r_abs_pow,2))));
                hold on;
%                 plot(10*log10(r_abs_pow_opt(:,col_max)));
                plot(10*log10(sum(r_abs_pow_opt,2) / max(sum(r_abs_pow,2))));
                title('across-track PTR')
                grid on;
                
                save('/home/schlembach/TUDTUM/tmp/simulation_coastalffsar_transponder.mat', 'CS1b', 'DOM', 'proc_sets')
            else
                uw_phase_detrended = plot_phase(ffproc);
                assert(std(uw_phase_detrended)< 2.2);
            end
        end      
        
    end
end


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