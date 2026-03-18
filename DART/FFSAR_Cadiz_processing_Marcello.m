clear all
close all
clc

%% define parameters
mission = 'S3B';
orbit = '114';

%% define L1A and L2 directories

L1A_dir = ['/home/fehelers/PhD Delft/Projects/TUMTUD_CADIZ/Data/',mission,'_',orbit];
L2_dir = '/home/fehelers/PhD Delft/Projects/TUMTUD_CADIZ/Data/L2';
FF_dir = '/home/fehelers/PhD Delft/Projects/TUMTUD_CADIZ/Data/L2_FFSAR/';

log_dir = '/home/fehelers/PhD Delft/Projects/TUMTUD_CADIZ/Data/logs/';
diary_fname = [log_dir,'proc_log_',mission,'_',orbit,'.txt'];
diary(diary_fname)
fprintf(['################ PROCESSING RUN: ',char(datetime),' ################# \n'])

%% check if there is an L2 file for all L1A files in the folder

L1A_files   = dir([L1A_dir,'/**/',mission,'*_SR_1_SRA_A*']);
dataset  = struct([]);

for i = 1:length(L1A_files)

    % example of file names, unfortunately sensing time is not always matching!
    %'S3A_SR_2_WAT____20160513T213002_20160513T221508_20191210T103845_2706_004_114______MR1_R_NT_004.SEN3'
    %'S3A_SR_1_SRA_A__20160513T212439_20160513T221508_20190828T233304_3029_004_114______MR1_R_NT_004.SEN3'
    %'XXXOOOOOOOOOOOOOXXXXXXXXOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOXXXXXXXOOOOOOXXXXXXXXXXXXXXXXX'

    %mission = L1A_files(i).name(1:3);
    date = L1A_files(i).name(17:25);
    cycle_orbit = L1A_files(i).name(69:77);

    L2_file = dir([L2_dir,'/**/*',mission,'*','_SR_2_WAT','*',date,'*',cycle_orbit,'*']);

    if ~isempty(L2_file) & (length(L2_file)==1)
        fprintf('Matching L2 file available for given L1A file \n')
        j = length(dataset)+1;
        dataset(j).L1A   = [L1A_files(i).folder, '/', L1A_files(i).name, '/measurement_l1a.nc'];
        dataset(j).L2    = [L2_file.folder, '/', L2_file.name, '/enhanced_measurement.nc'];
    elseif length(L2_file)>1
        warning(['Multiple Matching L2 files available for given L1A file: ', L1A_files(i).name])
    else
        warning(['No Matching L2 file available for given L1A file: ', L1A_files(i).name])
    end
end

%% ############## PROCESSING FFSAR #################
addpath('/home/fehelers/PhD_Delft/DART/DART')

LoadCommonSettings
DDAcf   = DDA_ConfigFile(mission,'SAR');
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;

DOM = [36.7 37.3; -30 30] % actual DOM
%DOM = [36.7 36.72; -30 30] % actual DOM

%%
for i = 5:numel(dataset)
    % read the L2 file data
    LAT                    = ncread(dataset(i).L2,'lat_20_ku');
    LON                    = ncread(dataset(i).L2,'lon_20_ku');
    IDX_mask = ingeoquad(LAT,LON,DOM(1,:),DOM(2,:));
    IDX = find(IDX_mask);
    startLoc=min(IDX);
    count=max(IDX)-startLoc+1;
    
    if ~isempty(startLoc)

        container.DATA_EUML2.LAT                    = ncread(dataset(i).L2,'lat_20_ku',startLoc,count);
        container.DATA_EUML2.LON                    = ncread(dataset(i).L2,'lon_20_ku',startLoc,count);
        %container.DATA_EUML2.SWH                    = ncread(dataset(i).L2,'swh_ocean_20_ku',startLoc,count);
        %container.DATA_EUML2.SWHcor                 = ncread(dataset(i).L2,'net_instr_cor_swh_20_ku',startLoc,count);
        container.DATA_EUML2.HEI                    = ncread(dataset(i).L2,'alt_20_ku',startLoc,count) - (ncread(dataset(i).L2,'range_ocean_20_ku',startLoc,count) - 0.55590); %reverse, that cog_cor has been added to range already to reproduce approximately 'raw' retracking result
        %%

        tic
        fprintf(['Processing file ',num2str(i),' of ',num2str(numel(dataset)),'\n']);
        % L1b processing
        ffproc = FFSAR_Processor(dataset(i).L1A,DOM,FFSAR_processing_settings);
        ffproc.setup_proc();
        ffproc.proc();

        % %% save
        % save([log_dir 'ffproc.mat'], 'ffproc');
        % 
        % %% load
        % ffproc = load([log_dir 'ffproc.mat']).ffproc;
        CS1b = ffproc.CS1b;

        %% try to get the right frequency (multilooking) to match the L2 latitudes
        l2lat = container.DATA_EUML2.LAT;
        fflat = CS1b.GEO.LAT;

        [minValue,closestIndex] = min(abs(l2lat - fflat),[],2);
        % get the corrected number of multilooks necessary to agree approximately with the L2 20 Hz spacing
        N_avg = round(mean(diff(closestIndex)));

        oversample = 4; % from 20Hz to 80Hz means 4 times 'oversampling'
        closestIndex = closestIndex + round(N_avg/oversample*(0:(oversample-1)));
        closestIndex = closestIndex(:);
        closestIndex(closestIndex+ceil(N_avg) > numel(CS1b.GEO.LAT)) = [];

        % figure;
        % plot(container.DATA_EUML2.LON,container.DATA_EUML2.LAT,'b.');hold on
        % plot(wrapTo360(CS1b.GEO.LON(closestIndex)),...
        %      CS1b.GEO.LAT(closestIndex),'ro');hold on

        %% retracking of FFSAR (testwise)
        CS1ba = FF_SAR_Average_Waveforms(CS1b,round(N_avg/4),mission,closestIndex);
        CS1ba.SAR.BeamAngle = zeros(size(CS1ba.GEO.LAT));
        CS1ba.MEA.ref_range = CS1ba.MEA.win_delay*CONST.c/2;
        CS1ba.SAR.scale_factor_ku = ones(size(CS1ba.GEO.LAT));

        %% plotting averaged waveforms
        figure;
        imagesc(1:256,CS1ba.GEO.LAT,log10(CS1ba.SAR.data(:,1:end)'))
        ylabel('LAT')
        xlabel('range bin')
        title('averaged FF-SAR waveforms')

        %% retracking of FF-SAR
        [DATA_FF,CS_FF]                           = SAR_L1b_to_L2(mission,CS1ba,DOM,'SAMOSA2',{'ApplyWfClass','false'; 'UseFLAGS','false'; 'UseCOR','false'; 'UseSurfTypeSeaOnly','true'; 'UseEmptyBeamAngle','true'});

        %% plotting retracked results vs. L2 file
        figure;
        plot(container.DATA_EUML2.LAT,container.DATA_EUML2.HEI,'bo');hold on
        plot(DATA_FF.LAT,DATA_FF.HEI,'ro');hold on

        [lat,idx] = sort(DATA_FF.LAT);
        hei_sorted = movmean(DATA_FF.HEI(idx),4);
        std(detrend(hei_sorted));

        legend(['EUMETSAT L2 20Hz, detrended std=',num2str(std(detrend(container.DATA_EUML2.HEI))),' m'],...
               ['FF-SAR 20Hz average, detrended std=',num2str(std(detrend(hei_sorted))),' m'])
        grid on
        xlabel('LAT')
        ylabel('HEI')

        ffname = split(dataset(i).L1A,'/');
        ffname = ffname{end-1};
        saveas(gcf,[FF_dir,'plots/',ffname,'.png'])
        %%
    %     figure;
    %     %plot(detrend(container.DATA_EUML2.HEI),'r.'); hold on
    %     [lat,idx] = sort(DATA_FF.LAT)
    %     hei_sorted = movmean(DATA_FF.HEI(idx),4)
    %     plot(detrend(hei_sorted),'b.')
    %     std(detrend(hei_sorted))

        %% save variables to right member of dataset
        container.CS1b = CS1ba;
        container.DATA_FFL2 = DATA_FF;
        container.L1A_file = ffname;

        save([FF_dir, ffname ,'.mat'], 'container');
        close all
        clear container CS1b DATA_FF ffproc
    end
    toc
end
