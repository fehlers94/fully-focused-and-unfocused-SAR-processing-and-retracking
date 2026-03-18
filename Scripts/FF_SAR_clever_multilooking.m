%% ##########################################
%               coastal case
% ##########################################
clear all;

filename = 'S3A_SR_1_SRA_A__20211210T101721_20211210T110750_20220105T012443_3029_079_279______MAR_O_NT_004.SEN3.mat';

load(fullfile('/home/fehelers/PhD Delft/Documents/Courses/Teaching Assistant/Summer school SAR altimetry/coastal',filename));

%%
m_per_look = 1e3*deg2km(distance(CS1b.GEO.LAT(1),CS1b.GEO.LON(1),CS1b.GEO.LAT(end),CS1b.GEO.LON(end)))/numel(CS1b.GEO.LON);




%% some visualisation of the waveforms
n=5
figure;
ax1 = subplot(1,2,1)
imagesc(log10(movmean(CS1b.SAR.data,n,2)))
colormap('pink')
caxis([8,13])

ax2 = subplot(1,2,2)
imagesc(log10(CS1b.SAR.data_pseudoDD))
colormap('pink')
caxis([8,13])

linkaxes([ax1,ax2],'xy')
example = [];

%% use Kolmogorov Smirnoff test in order to judge whether single look powers are drawn from an identical distribution (maybe start with 50 samples?)
n_s = 20;
[n_r,n_w] = size(CS1b.SAR.data);

% do the test along the 128th range gate
data = CS1b.SAR.data(:,1:240000);
n_w = 240000;
data = reshape(data,[n_r n_s n_w/n_s]); %reshaping the data


%% doing 2 sample test with consecutive along-track neighbours

test = zeros(n_r,n_w/n_s);
p_test = zeros(n_r,n_w/n_s);

comp_dists = round([45 135 225] / (m_per_look*n_s) ) % these are the neighbouring indices to compare the statistics to. They have been chosen maximally distant to the replica at 90, 180, 270 m
comp_dists = [-comp_dists comp_dists]

%%
count = 0
for r = 1:256
    for s = -min(comp_dists)+1:(n_w/n_s)-max(comp_dists)-1
        %cdfplot(data(r,:,s));hold on
        %[h,p] = kstest2(data(r,:,s),data(r,:,s-1),'Alpha',0.05); % test neighbouring columns
        %[h,p] = kstest2(data(r,:,s),[data(r,:,s),data(r,:,4800)],'Alpha',0.05); % test against a manually chosen "good" column
        comp_data = [data(r,:,s-comp_dists)];
        comp_data = comp_data(:);
        [h,p] = kstest2(data(r,:,s),comp_data,'Alpha',0.05); % test against a multitude of directly neighbouring columns not affected by the same replica
        test(r,s) = h;
        p_test(r,s) = p;
    end
    count=count+1
end

% brainstorm:
% it needs to be a "local" comparison because sigma0 and waveheight and ssh
% are changing. Hence, we need a local criterion/local distribution to
% compare to. Neighbours is unsatisfacory, since it only highlights
% differences but canot handle extended contaminations. So some predefined exponential
% distribution with given mean would be good? Not sure whether it is
% accurate to do this, try out. Otherwise, include more data in reference
% comparison? Thing is, when e.g. a moving median is applied and there is
% as much contamination as open ocean, then the procedure will break down,
% failing to identify the spurious regions. So essentially, some sort of
% a-priori waveform shape input is needed?

%% plot average waveforms and test outcome
ax1 = subplot(1,3,1)
imagesc(log10(squeeze(sum(data,2))));
ax2 = subplot(1,3,2)
imagesc(test);
ax3 = subplot(1,3,3)
imagesc(p_test<0.05);

linkaxes([ax1 ax2 ax3],'xy')



%%    

%[h,p] = kstest2(x1,x2,'Alpha',0.01)




%%

% % example.d1.FF = CS1b.SAR.data(45:66,24000:25500);
% % example.d1.UF = CS1b.SAR.data_pseudoDD(45:66,24000:25500);
% 
% % example.d1.FF = CS1b.SAR.data(40:60,28000:30000);
% % example.d1.UF = CS1b.SAR.data_pseudoDD(40:60,28000:30000);
% 
% % example.d1.FF = CS1b.SAR.data(5:45,74200:75100);
% % example.d1.UF = CS1b.SAR.data_pseudoDD(5:45,74200:75100);

% best so far!!!
% example.d1.FF = CS1b.SAR.data(5:20,73200:76100);
% example.d1.UF = CS1b.SAR.data_pseudoDD(5:20,73200:76100);

% example.d1.FF = CS1b.SAR.data(5:40,73200:76100);
% example.d1.UF = CS1b.SAR.data_pseudoDD(5:40,73200:76100);

% example.d1.FF = CS1b.SAR.data(35:60,136000:138500);
% example.d1.UF = CS1b.SAR.data_pseudoDD(35:60,136000:138500);

% best number 2
% example.d1.FF = CS1b.SAR.data(15:45,117200:118400);
% example.d1.UF = CS1b.SAR.data_pseudoDD(15:45,117200:118400);

example.d1.FF = CS1b.SAR.data(25:65,112400:114000);
example.d1.UF = CS1b.SAR.data_pseudoDD(25:65,112400:114000);



x_dist = m_per_look*(1:size(example.d1.UF,2))

example.d1.FF = example.d1.FF./sum(example.d1.FF(:));
example.d1.UF = example.d1.UF./sum(example.d1.UF(:));



figure;
ax3 = subplot(3,1,1)
FF=sum(example.d1.FF)
plot(x_dist,(FF./max(FF)-0.12)./0.78);hold on
UF = sum(example.d1.UF)
plot(x_dist,(UF-min(UF))./max((UF-min(UF))))
colorbar()
ylabel('normalized power')
legend('FF-SAR','UF-SAR')

ax1 = subplot(3,1,2)
imagesc(x_dist,1:size(example.d1.FF,1),example.d1.FF)
title('FF-SAR')

colormap('pink')
%caxis([8,13])
colorbar()

ax2 = subplot(3,1,3)
imagesc(x_dist,1:size(example.d1.FF,1),example.d1.UF)
title('UF-SAR')
colormap('pink')
xlabel('along-track distance (m)')

%caxis([2e-5,5e-5])
linkaxes([ax1,ax2,ax3],'x')
colorbar()

%%#####################################################
%% start prototyping some code for removing outliers:
%%#####################################################

% illustrate UF-SAR data with obvious contamination all along the waveform:
imagesc(log10(CS1b.SAR.data_pseudoDD(:,81000:85000)));
colormap('pink')

%% define an n waveforms wide multilooking
n = 10; % multilooking prior to filtering (for robustness)
M = 1000; % along-track window for median
rb = 11; % range window for median

data = (CS1b.SAR.data(:,81000:85000));
data = movmean(data,n,2);
data = data./max(data(:));

% take the logarithm (reasonable for histogram applications, since monotonic function)
data = log10(data);
figure;
imagesc(data)


%% compute movmedian and movmad as a reference over M waveforms
data_m = medfilt2(data,[rb M]);
%data_mad = movmad(data,M,2);
figure
imagesc(data_m)

%%
mask = (data>(data_m+0.2));
data_new = data;
data_new(mask) = NaN;
figure;
imagesc(data_new.*(~mask))
colorbar()




%%
%data1 = log10(data);
%data1 = movmean(data1,10,2);
n = 300
data_movmean = log10(movmean(data,n,2));
data_movstd = movstd(data_movmean,n,1,2);

figure;
imagesc(data_movmean)
colormap('pink')

figure;
imagesc(data_movstd)
colormap('pink')

