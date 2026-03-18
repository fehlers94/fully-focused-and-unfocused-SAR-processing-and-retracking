function [dummy] = FFSAR_bulk_L1a_to_L1b(year,month,DOM)
    
    %YYYY = 2018:2020;
    %sprintf('\n qsub-eval -l CS_%i_%02i.log "mlab -eval \\"FFSAR_bulk_L1a_to_L1b(''%i'',''%02i'',[51.9 53; -30 30])\\""\n',repmat([reshape(repmat(YYYY,12,1),numel(YYYY)*12,1),reshape(repmat(1:12,numel(YYYY),1)',numel(YYYY)*12,1)],1,2)')

    %but for memory troubles use rather 2 processors:
    %YYYY = 2018:2020;
    %sprintf('\n qsub-eval -nn 1 -np 2 -l S3_%i_%02i.log "mlab -eval2 \\"FFSAR_bulk_L1a_to_L1b(''%i'',''%02i'',[51.9 53; -30 30])\\""\n',repmat([reshape(repmat(YYYY,12,1),numel(YYYY)*12,1),reshape(repmat(1:12,numel(YYYY),1)',numel(YYYY)*12,1)],1,2)')


    %DOM = [51.9 53; -30 30]
    %year = '2017'
    %month = '02'
    %num_workers = 2
    %##########################################################################
    LoadCommonSettings

    % read the L1a data directory
    files   = dir([PathDATA,'/ROFI/L1a/*_SR_1_SRA_A*']);
    L1A_dataset  = struct([]);

    % get only files with matching year and date
    j=int16(1);
    for i=1:length(files)
        if strcmp(files(i).name(17:20),year) & strcmp(files(i).name(21:22),month)
            L1A_dataset(j).name = files(i).name;
            L1A_dataset(j).mission = files(i).name(1:3);
            j = j+1;
        end
    end

    L1A_dataset.name

    % start L1b processing of all files
    for i=1:length(L1A_dataset)
        L1b_filepath = fullfile(PathDATA,'ROFI','L1b',[L1A_dataset(i).name '.mat']);
        if ~isfile(L1b_filepath)    
            % Load parameters and settings
            DDAcf   = DDA_ConfigFile(L1A_dataset(i).mission,'SAR');

            % determine FFSAR processing settings
            [CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
            FFSAR_processing_settings.combine_n_look_locations = 50;
            FFSAR_processing_settings.output_pseudo_delay_doppler_processing = true;
            FFSAR_processing_settings.num_coherent_bursts = 1;
            FFSAR_processing_settings.integration_time = DDAcf.BRI*180.01; % corresponding to 180 bursts as in L1b products
            %FFSAR_processing_settings.integration_time = 2;
            %FFSAR_processing_settings.along_track_oversampling = 1/(0.25*585); % corresponding to ~80 Hz
            FFSAR_processing_settings.along_track_oversampling = 1; % corresponding to FFSAR IRF speckle noise resolution
            FFSAR_processing_settings.simulate_range_walk = false;
            FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant';
            %FFSAR_processing_settings.n_cores_parfor = num_workers;

            L1a_filepath = fullfile(PathDATA,'ROFI','L1a',L1A_dataset(i).name,'measurement_l1a.nc');
            ffproc = FFSAR_Processor(L1a_filepath, DOM, FFSAR_processing_settings);
            %ffproc = FFSAR_Processor(L1a_file, DOM, FFSAR_processing_settings);
            ffproc.setup_proc();
            tic
            ffproc.proc();
            toc
            CS1b = ffproc.CS1b;

            save(L1b_filepath,'CS1b','-v7.3' )
        else
            warning(['L1b-file already existing: ', L1A_dataset(i).name, '\n'])
        end
    end
    
    dummy = 0;
end