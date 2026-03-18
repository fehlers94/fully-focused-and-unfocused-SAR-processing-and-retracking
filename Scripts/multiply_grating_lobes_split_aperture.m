range = [-3,0]
%%
field = {'dataIQ_050T_025T', 'dataIQ_025T_000T',  'dataIQ_000T_025T', 'dataIQ_025T_050T'};

% normalize all radargrams
for i = 1:4
    CS1b.SAR.(field{i}) = abs(CS1b.SAR.(field{i})).^2;
    
    norm = sum(CS1b.SAR.(field{i}));
    norm =  sum(norm(:));
    CS1b.SAR.(field{i}) = CS1b.SAR.(field{i})./norm;
end

%% add offset

for i = 1:4
    CS1b.SAR.(field{i}) = CS1b.SAR.(field{i});
end


%%

figure;
for i = 1:4
    ax(i) = subplot(6,1,i+1)
    imagesc(log10(CS1b.SAR.(field{i})))
    colormap('pink')
    caxis(range)
end

C = cat(3,CS1b.SAR.(field{1}),CS1b.SAR.(field{2}),CS1b.SAR.(field{3}),CS1b.SAR.(field{4}));
C_sorted = sort(C,3);
CS1b.SAR.filtered = sum(C_sorted(:,:,1:2),3);
CS1b.SAR.filtered = CS1b.SAR.filtered./sum(CS1b.SAR.filtered(:));

%CS1b.SAR.filtered = ( CS1b.SAR.(field{1}).*CS1b.SAR.(field{2}).*CS1b.SAR.(field{3}).*CS1b.SAR.(field{4}) ).^(1/4);

ax(5) = subplot(6,1,6)
imagesc(log10((CS1b.SAR.filtered)))
colormap('pink')
caxis(range)


ax_sum = subplot(6,1,1)
plot(sum(CS1b.SAR.filtered,1),'r-'); hold on
for i=1:4
    plot(sum(CS1b.SAR.(field{i}),1),'k-'); hold on
end

linkaxes([ax ax_sum],'x')
linkaxes(ax,'y')

%% test approach with real data:

path = '/home/fehelers/PhD Delft/Projects/FFSAR4ROFI/Data/ROFI/L1a/S3A_SR_1_SRA_A__20181129T101709_20181129T110739_20181225T013537_3029_038_279______MAR_O_NT_003.SEN3/measurement_l1a.nc'

DDAcf = DDA_ConfigFile('S3A','SAR');
[CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
FFSAR_processing_settings.combine_n_look_locations = 50;
FFSAR_processing_settings.output_pseudo_delay_doppler_processing = true;
FFSAR_processing_settings.num_coherent_bursts = 1;
FFSAR_processing_settings.simulate_range_walk = false;
FFSAR_processing_settings.focal_point_height_interpolation = 'piecewise_constant';
FFSAR_processing_settings.along_track_oversampling = 1;
FFSAR_processing_settings.split_aperture = true;
FFSAR_processing_settings.integration_time = 2;

% choose DOM symmetrically around the transponder
DOM = [52.6 52.8; -180 180];

ffproc = FFSAR_Processor_sa(path, DOM, FFSAR_processing_settings);
ffproc.setup_proc();
%%
ffproc.proc();

CS1b = ffproc.CS1b;

%% normalize the different radargrams:

field = {'dataIQ_050T_025T', 'dataIQ_025T_000T',  'dataIQ_000T_025T', 'dataIQ_025T_050T'};

% normalize all radargrams
for i = 1:4
    CS1b.SAR.(field{i}) = abs(CS1b.SAR.(field{i})).^2;
    
    norm = max(CS1b.SAR.(field{i}));
    norm =  max(norm(:));
    CS1b.SAR.(field{i}) = CS1b.SAR.(field{i})./norm;
end

%% multilook the radargram prior to filtering
CS1b_avg = FF_SAR_Average_Waveforms(CS1b,20,'S3A');


%% filter the radargrams according to the defined rule
C = cat(3,CS1b_avg.SAR.(field{1}),CS1b_avg.SAR.(field{2}),CS1b_avg.SAR.(field{3}),CS1b_avg.SAR.(field{4}));
C_sorted = sort(C,3);
CS1b_avg.SAR.filtered = mean(C_sorted(:,:,1:3),3);
%CS1b_avg.SAR.filtered = CS1b_avg.SAR.filtered./max(CS1b_avg.SAR.filtered(:));
CS1b_avg.SAR.data = CS1b_avg.SAR.data./max(CS1b_avg.SAR.data(:))


%%
figure;
for i = 1:4
    ax(i) = subplot(1,5,i)
    imagesc(log10(CS1b_avg.SAR.(field{i})))
    colormap('pink')
    caxis([-5,1])
end

ax5 = subplot(1,5,5)
imagesc(log10(CS1b_avg.SAR.filtered))
colormap('pink')
caxis([-5,1])

linkaxes([ax ax5],'xy')


%%
ax1 = subplot(1,3,1)
imagesc(log10(CS1b_avg.SAR.data))
%caxis([-5,0])

ax2 = subplot(1,3,2)
imagesc(log10(CS1b_avg.SAR.filtered))
%caxis([-5,0])
linkaxes([ax1 ax2],'xy')