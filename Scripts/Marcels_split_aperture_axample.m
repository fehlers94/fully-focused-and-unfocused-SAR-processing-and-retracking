% #########################################################################
% simple script to make L1a to L1b FFSAR processing for S6A with a split aperture
% #########################################################################

mission = 'S6A'
L1a_filepath = [getenv('HOME') '/TUDTUM/ffsar_data/' 's6a/l1a-l1b-l2/crete/S6A_P4_1A_HR______20210901T212456_20210901T222027_20210902T134010_3331_030_018_009_EUM__OPE_ST_F03.nc'];

% loading some satellite parameters
DDAcf = DDA_ConfigFile(mission,'SAR');

% loading default FF-SAR processing settings
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;

DOM = [TR.crete.lat-0.005 TR.crete.lat+0.005; -180 180]; % choose DOM symmetrically around the transponder here, otherwise the latlonbox of your choice

% now, here the along-track oversampling
FFSAR_processing_settings.integration_time = 2; % this determines the length of the complete aperture in s
FFSAR_processing_settings.along_track_oversampling = 1; % this determines the oversampling w.r.t. the resolution obtained from the complete integration time (whole aperture)
FFSAR_processing_settings.combine_n_look_locations = 50; % determines the number of consecutive waveforms that are obtained using the phasor (a constant RCMC)
FFSAR_processing_settings.split_aperture = true; % this determines that you want to split the aperture (integration times) into 4 parts

FFSAR_processing_settings.output_pseudo_delay_doppler_processing = false; % not needed for you. But possible...
FFSAR_processing_settings.simulate_range_walk = false; % also not needed...
FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant'; % crucial!

% #########
% You know, that the resolution of the radargrams from the split aperture is worse by a
% factor of 4. So setting 'along_track_oversampling' to 0.25 in order to save diskspace is probably the preferred option for you?
% However, then you should also set combine_n_look_locations to 12 or
% so, because otherwise the along-track distance over which you apply a
% constant RCMC increases from ~50*1 m to ~50*4 = 200 m, which is probably
% too coarse!
% #########

%% making the actual focussing
ffproc = FFSAR_Processor(L1a_filepath, DOM, FFSAR_processing_settings);
ffproc.setup_proc();
ffproc.proc();

CS1b = ffproc.CS1b;
%%
save(['somepath' '.mat'],'CS1b','-v7.3' )

%% plot the four different radargrams after power detection

field = {'dataIQ_050T_025T', 'dataIQ_025T_000T',  'dataIQ_000T_025T', 'dataIQ_025T_050T'};

figure;
for i = 1:4
    ax(i) = subplot(2,2,i)
    wfdata = abs(CS1b.SAR.(field{i}));
    wfdata = wfdata./max(wfdata(:));
    imagesc(log10(wfdata));
    title(field{i})
    %colormap('pink')
    colorbar()
    ylabel('range gate number')
    xlabel('waveform number')
end

linkaxes(ax,'xy')