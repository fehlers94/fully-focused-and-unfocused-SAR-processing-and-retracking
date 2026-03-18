% this script is intended to see how the SAR altimeters ground resolution
% will be altered over one orbit due to orbit parameter changes

clear all
%%
%l2name = 

% % pacific:
% l1name = '/home/fehelers/Downloads/S3A_SR_1_SRA____20191231T075312_20191231T084340_20200209T102117_3028_053_163______MR1_R_NT_004.SEN3/S3A_SR_1_SRA____20191231T075312_20191231T084340_20200209T102117_3028_053_163______MR1_R_NT_004.SEN3/measurement.nc'

% atlantic: 
%l1name = '/home/fehelers/Downloads/S3A_SR_1_SRA____20191231T230205_20191231T235234_20200209T101741_3029_053_172______MR1_R_NT_004.SEN3/S3A_SR_1_SRA____20191231T230205_20191231T235234_20200209T101741_3029_053_172______MR1_R_NT_004.SEN3/measurement.nc'

% pacific2:
l1name = '/home/fehelers/Downloads/S3B_SR_1_SRA____20191231T195106_20191231T204134_20200210T143207_3028_034_028______MR1_R_NT_004.SEN3/S3B_SR_1_SRA____20191231T195106_20191231T204134_20200210T143207_3028_034_028______MR1_R_NT_004.SEN3/measurement.nc'

[S3,CS1b] = S3_L1b_read(l1name)


CS1b.GEO.V_abs = sqrt(CS1b.GEO.V.Vx.^2+CS1b.GEO.V.Vy.^2+CS1b.GEO.V.Vz.^2)

time_diffs = diff(CS1b.GEO.TAI.secs);
time_diffs(end+1) = NaN;
CS1b.GEO.time_diffs = time_diffs;

%% calculate theoretical resolution from altitude and v_abs etc.
DDAcf = DDA_ConfigFile('S3A','SAR')

Lx = CS1b.MEA.ref_range*DDAcf.c./(2*DDAcf.fc*64*DDAcf.PRI*CS1b.GEO.V_abs);
%Lx = CS1b.GEO.H*DDAcf.c./(2*DDAcf.fc*64*DDAcf.PRI*CS1b.GEO.V_abs);

%%
subplot(1,2,1)
plot(Lx);hold on
subplot(1,2,2)
plot(time_diffs)

