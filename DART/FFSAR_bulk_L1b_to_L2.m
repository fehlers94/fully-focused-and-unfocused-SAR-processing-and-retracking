function [dummy] = FFSAR_bulk_L1b_to_L2(year,month,DOM)

    %YYYY = 2018:2020;
    %sprintf('\n qsub-eval -l CS_%i_%02i.log "mlab -eval \\"FFSAR_bulk_L1b_to_L2(''%i'',''%02i'',[51.9 53; -30 30])\\""\n',repmat([reshape(repmat(YYYY,12,1),numel(YYYY)*12,1),reshape(repmat(1:12,numel(YYYY),1)',numel(YYYY)*12,1)],1,2)')

    %YYYY = 2018;
    %sprintf('\n qsub-eval -nn 1 -np 2 -l L2proc_%i_%02i.log "mlab -eval2 \\"FFSAR_bulk_L1b_to_L2(''%i'',''%02i'',[51.9 53; -30 30])\\""\n',repmat([reshape(repmat(YYYY,12,1),numel(YYYY)*12,1),reshape(repmat(1:12,numel(YYYY),1)',numel(YYYY)*12,1)],1,2)')

    
    %DOM = [51.9 53; -30 30]
    %year = '2018'
    %month = '11'
    %##########################################################################
    LoadCommonSettings

    % read the L1b data directory to find files in the right month of year

    files   = dir(fullfile(PathDATA,'ROFI','L1b','*_SR_1_SRA_A*'));
    L1B_dataset  = struct([]);

    % get only files with matching year and date
    j=int16(1);
    for i=1:length(files)
        if strcmp(files(i).name(17:20),year) & strcmp(files(i).name(21:22),month)
            L1B_dataset(j).name = files(i).name;
            L1B_dataset(j).mission = files(i).name(1:3);
            j = j+1;
        end
    end

    L1B_dataset.name

    % find the fitting L2 files
    if ~isempty(L1B_dataset)
        for i = 1:length(L1B_dataset)

            % example of file names, unfortunately sensing time is not always matching!
            %'S3A_SR_2_WAT____20160513T213002_20160513T221508_20191210T103845_2706_004_114______MR1_R_NT_004.SEN3'
            %'S3A_SR_1_SRA_A__20160513T212439_20160513T221508_20190828T233304_3029_004_114______MR1_R_NT_004.SEN3'
            %'XXXOOOOOOOOOOOOOXXXXXXXXOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOXXXXXXXOOOOOOXXXXXXXXXXXXXXXXX'

            %mission = L1A_files(i).name(1:3);
            date = L1B_dataset(i).name(17:25);
            cycle_orbit = L1B_dataset(i).name(69:77);

            L2_file = dir([PathDATA,'/ROFI/L2/','*','_SR_2_WAT','*',date,'*',cycle_orbit,'*.SEN3']);

            if ~isempty(L2_file) & (length(L2_file)==1)
                fprintf(['Matching L2 file available for L1B file: ', L1B_dataset(i).name, '\n'])
                L1B_dataset(i).L2_filepath = fullfile(L2_file.folder, L2_file.name, 'enhanced_measurement.nc');
            elseif length(L2_file)>1
                warning(['Multiple Matching L2 files available for given L1B file: ', L1B_dataset(i).name, '\n'])
            else
                warning(['No Matching L2 file available for given L1B file: ', L1B_dataset(i).name, '\n'])
            end
        end
    else
        warning(['No L1b data files available in the specified time interval', '\n'])
    end


    for i= 1:length(L1B_dataset)
        L2_filepath = fullfile(PathRSLT,'ROFI_L2',[L1B_dataset(i).name]);
        if isempty(L1B_dataset(i).L2_filepath)
            warning(['L2 processing lacks corrections and cannot be performed for file: ', L1B_dataset(i).name, '\n'])
        elseif ~isfile(L2_filepath)
            fprintf(['Performing L2 processing for L1B file: ', L1B_dataset(i).name, '\n'])

            % initialise L2 results struct, which at the same time determines the sampling frequencies
            % and FF-SAR multilooking
            CS2 = [];
            CS2.('f_020Hz')=[];
            CS2.('f_040Hz')=[];
            CS2.('f_080Hz')=[];
            CS2.('f_120Hz')=[];
            CS2.('f_240Hz')=[];
            CS2.('f_020Hz_matched')=[];

            % load the L1b file
            load(fullfile(PathDATA,'ROFI','L1b',L1B_dataset(i).name))
            % afterwards, 'CS1b' has been loaded into workspace

            %% read L2 data with corrections for specific L1b file and make interpolator object for the corrections based on latitude coordinate
            [EUML2,COR] = S3_read_L2(L1B_dataset(i).L2_filepath,DOM); % EUML2 contains L2 data, whereas COR contains interpolators of the corrections dependent on latitude coordinate

            % get the first latitude value of FF-SAR and get the start time this corresponds to
            % T0 = interp1(EUML2.lat_20_ku,EUML2.time_20_ku,CS1b.GEO.LAT(1),'linear') % time of the first FF-SAR waveform (ElapsedTime = 0)

            % var = 'cog_cor_01';
            % plot(EUML2.lat_01,EUML2.(var),'ro');hold on
            % plot(CS1b.GEO.LAT,COR.(var)(CS1b.GEO.LAT),'b.');hold on

            %%
            % define the multilooking frequency
            % resolutions of 300, 150, 75, 50, 25 m
            % require about 20 Hz, 40 Hz, 80 Hz, 120 Hz, 240 Hz
            % update about:  354.7024  177.3512   88.6756   59.1171   29.5585
            % do I make this with respect to time or with respect to latitude (length)?
            % best to plot both, ground along track distance and 

            % check arclen on ref ellipsoid
            % arclen = distance(0,51,0,52,referenceEllipsoid('WGS 84'),'deg') % seems to give result in meters?
            % arclen = distance(CS1b.GEO.LON(1),CS1b.GEO.LAT(1),CS1b.GEO.LON(2),CS1b.GEO.LAT(2),referenceEllipsoid('WGS 84'),'deg') % seems to give result in meters?
            % get along track distances of neighbouring points in dataset
            % arclen = distance(CS1b_m.GEO.LON(1:end-1),CS1b_m.GEO.LAT(1:end-1),CS1b_m.GEO.LON(2:end),CS1b_m.GEO.LAT(2:end),referenceEllipsoid('WGS 84'),'deg');
            % plot(arclen)

            %% make for loop over different multilooking scenarios

            scenarios = fields(CS2);

            for j = 1:numel(scenarios)

                %%% get the time stamps/indices for multilooking of the L1b data
                freq = str2num(scenarios{j}(3:5))
                match_EUML2 = contains(scenarios{j},'matched')

                n_multilook = round(numel(CS1b.GEO.Elapsed_Time)./(CS1b.GEO.Elapsed_Time(end)-CS1b.GEO.Elapsed_Time(1))./freq)
                t_multilook = 1/(2*freq):1/freq:CS1b.GEO.Elapsed_Time(end); % this is either given by the frequency or by an array from the L2 file

                %[minValue,closestIndex] = min(abs(t_multilook' - CS1b.GEO.Elapsed_Time),[],2); % not memory efficient

                if ~match_EUML2
                    closestIndex = interp1(CS1b.GEO.Elapsed_Time, 1:numel(CS1b.GEO.Elapsed_Time), t_multilook, 'nearest', 'extrap'); % memory saving and fast
                else
                    closestIndex = interp1(CS1b.GEO.LAT, 1:numel(CS1b.GEO.LAT), EUML2.lat_20_ku, 'nearest', 'extrap'); % memory saving and fast
                    % cut off the parts that exceed the first and last index
                    closestIndex(closestIndex<=1) = [];
                    closestIndex(closestIndex>=numel(CS1b.GEO.LAT)) = [];

                    eumlat = EUML2.lat_20_ku;
                    eumind = 1:numel(EUML2.lat_20_ku);
                    idx_eumlat = (DOM(1)-0.2<eumlat)&(eumlat < DOM(3)+0.2);

                    matchedIndexEUM = interp1(eumlat(idx_eumlat), eumind(idx_eumlat), CS1b.GEO.LAT(closestIndex), 'nearest', 'extrap'); % 20 Hz index in L2 file
                end
                %%% multilook and retrack the data for FF and UF SAR
                CS1b_m = FF_SAR_Average_Waveforms(CS1b,n_multilook,L1B_dataset(i).mission,closestIndex); % multilooked data
                %imagesc(log10(CS1b_m.SAR.data)); colorbar()
                [CS2FF,~] = SAR_L1b_to_L2(L1B_dataset(i).mission,CS1b_m,DOM+[-0.01 0.01; 0 0],'SAMOSA2',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false';'do_subwaveform_fit','true'});

                CS1b_m.SAR.data = CS1b_m.SAR.data_pseudoDD;
                [CS2UF,~] = SAR_L1b_to_L2(L1B_dataset(i).mission,CS1b_m,DOM+[-0.01 0.01; 0 0],'SAMOSA2',{'SAMOSAimpl', "'DPM_v2_5_2'";'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','false';'do_subwaveform_fit','true'});
                %plot(CS2.HEI);hold on;plot(CS2.ALT-CS2.range,'ro') % compares well, HEI is uncorrected height resulting from satellite range minus raw retracked range

                %%% correct the range and SSHi data with the COR interpolator of the L2 file

                x = CS1b_m.GEO.LAT';
                CS2FF.instr_corrected_range = CS2FF.range + COR.cog_cor_01(x) + COR.mod_instr_cor_range_01_ku(x);
                CS2UF.instr_corrected_range = CS2UF.range + COR.cog_cor_01(x) + COR.mod_instr_cor_range_01_ku(x);

                CS2FF.SSHi = CS2FF.ALT' - (CS2FF.instr_corrected_range + ...
                                          COR.mod_dry_tropo_cor_meas_altitude_01(x) + ...
                                          COR.mod_wet_tropo_cor_meas_altitude_01(x) + ...
                                          COR.load_tide_sol2_01(x) + ...
                                          COR.iono_cor_gim_01_ku(x) + ...
                                          COR.solid_earth_tide_01(x) + ...
                                          0.468*COR.pole_tide_01(x)   );

                CS2UF.SSHi = CS2UF.ALT' - (CS2UF.instr_corrected_range + ...
                                          COR.mod_dry_tropo_cor_meas_altitude_01(x) + ...
                                          COR.mod_wet_tropo_cor_meas_altitude_01(x) + ...
                                          COR.load_tide_sol2_01(x) + ...
                                          COR.iono_cor_gim_01_ku(x) + ...
                                          COR.solid_earth_tide_01(x) + ...
                                          0.468*COR.pole_tide_01(x)   );

                if match_EUML2
                    % write relevant EUMETSAT L2 range data in into CS2.(scenarios{j}).EUM the Altitude/Latitude/range and SSHi
                    CS2EUM = [];
                    vars = ["time_20_ku","lat_20_ku","alt_20_ku","range_ocean_20_ku","UTC_sec_20_ku"];
                    for k =1:numel(vars)
                        CS2EUM.(vars(k)) = EUML2.(vars(k))(matchedIndexEUM);
                    end
                    CS2EUM.SSHi = CS2EUM.alt_20_ku - (CS2EUM.range_ocean_20_ku + ...
                                                      COR.mod_dry_tropo_cor_meas_altitude_01(x) + ...
                                                      COR.mod_wet_tropo_cor_meas_altitude_01(x) + ...
                                                      COR.load_tide_sol2_01(x) + ...
                                                      COR.iono_cor_gim_01_ku(x) + ...
                                                      COR.solid_earth_tide_01(x) + ...
                                                      0.468*COR.pole_tide_01(x)   );
                   CS2.(scenarios{j}).EUM = CS2EUM;
                end
                CS2FF.n_multilook = n_multilook;
                CS2.(scenarios{j}).FF = CS2FF;
                CS2.(scenarios{j}).UF = CS2UF;
            end
            %%% also save the interpolators for the corrections for later access
            CS2.COR = COR;

            %start time (human readable)
            CS2.StartTime = datestr(CS2.f_020Hz.FF.TIME(1));
            CS2.L2_filepath = L1B_dataset(1).L2_filepath;

            % save mat-file in L2 folder
            save(L2_filepath,'CS2','-v7.3' );
        else
            warning(['Self-processed L2-file already existing: ', L1B_dataset(i).name, '\n'])
        end
    end
end