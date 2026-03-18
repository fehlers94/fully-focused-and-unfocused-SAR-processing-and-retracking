addpath(genpath('/data/Projects/CryoSat2_Greenland/Software/CryoSat_Matlab_Reader_Package_v1_9/'))

defval('DOM',[52.3 53.3; 4.5 7.5])

%% L1b & L2 processing own routines - stable version
[CS1b3,run_check3]                   = FF_SAR_L1a_to_L1b('2010/08/CS_LTA__SIR1SAR_FR_20100803T233904_20100803T234014_C001.DBL',DOM);
CS1b3a                               = FF_SAR_Average_Waveforms(CS1b3,600);
CS1b3a.SAR.BeamAngle                 = cell(numel(CS1b3a.GEO.LAT),1);
CS1b3a.COR.ocean_equilibrium_tide(:) = 0;
CS1b3a.COR.ocean_longperiod_tide(:)  = 0;
CS1b3a.COR.dac(:)                    = 0;
[DATA3a,~]                           = SAR_L1b_to_L2('CS',CS1b3a,DOM,'SAMOSA2',{'ApplyWfClass','false'});

%% Compare to tide gauge water levels
PathDATA = '/data/Projects/Versatile_Hydrodynamics/Data';

DATAin = DATA3a;

%Interpolate NLGEO2018 to observation locations
NLGEO2018            = load(fullfile(PathDATA,'Gravity/Geoid_Models/NLGEO2018/NLGEO2018_Netherlands.mat'));
QG                   = griddedInterpolant(NLGEO2018.SYN.LONgrd',NLGEO2018.SYN.LATgrd',NLGEO2018.SYN.NLGEO2018','cubic','none');
DATAin.('NLGEO2018') = QG(DATAin.LON,DATAin.LAT);
clear('QG')

%Load approximate polygon of Lake IJssel
poly   = dlmread(fullfile(PathDATA,'Tide_Gauge','GriddedData','Lake_IJssel','IJsselmeer.txt'));
    
%Identify altimeter points within Lake IJssel & acquired before Jan 2019
INdata = true(size(DATAin.LAT));

%Create mesh for gridded tide gauge measurements
[lonTG,latTG] = meshgrid(min(poly(:,1)):0.01:max(poly(:,1)),min(poly(:,2)):0.02:max(poly(:,2)));

%Identify grid points within Lake IJssel
INgrd = inpolygon(lonTG,latTG,poly(:,1),poly(:,2));

%Bin observations in time domain. Note that there does not exist an
%ascii grid for all TSteps
TSteps   = floor(min(DATAin.TIME)):10/1440:ceil(max(DATAin.TIME));
FileExst = false(numel(TSteps),1);
for i=2010:2018
    FNames                    = dir(fullfile(PathDATA,'Tide_Gauge','GriddedData','Lake_IJssel',num2str(i),'*.ascii'));
    FNames                    = datenum(regexprep({FNames.name},'.ascii',''),'yyyymmddHHMMSS');
    [~,locb]                  = ismember(round(FNames,5),round(TSteps,5));
    FileExst(locb(locb ~= 0)) = true;
end
TSteps        = TSteps(FileExst);
[~,~,BIN_TSt] = histcounts(DATAin.TIME,[TSteps(1)-5/1440,TSteps(1:end-1)+.5*diff(TSteps),TSteps(end)+5/1440]);

%Compute DT corrections
CorrDTdom = nan(nnz(INdata),3);
for i=unique(BIN_TSt)'
    %Read gridded water levels at time step i
    MNWL2 = dlmread(fullfile(PathDATA,'Tide_Gauge','GriddedData','Lake_IJssel',datestr(TSteps(i),'yyyy'),fullfile(sprintf('%s.ascii',datestr(TSteps(i),'yyyymmddHHMMSS')))));
    
    %Read gridded water levels at time step i-1
    MNWL1 = dlmread(fullfile(PathDATA,'Tide_Gauge','GriddedData','Lake_IJssel',datestr(TSteps(i-1),'yyyy'),fullfile(sprintf('%s.ascii',datestr(TSteps(i-1),'yyyymmddHHMMSS')))));
    
    %Read gridded water levels at time step i+1
    MNWL3 = dlmread(fullfile(PathDATA,'Tide_Gauge','GriddedData','Lake_IJssel',datestr(TSteps(i+1),'yyyy'),fullfile(sprintf('%s.ascii',datestr(TSteps(i+1),'yyyymmddHHMMSS')))));
    
    %Identify RA/GNSS data points acquired close to time step i
    IDX     = find(BIN_TSt == i);
    
    %Interpolate gridded water levels to measurement epoch of ALL observations
    WLint    = interp1(TSteps(i-1:i+1),[MNWL1,MNWL2,MNWL3]',DATAin.TIME(IDX))';
    
    %Interpolate water levels at grid points to obs. locations
    CorrDTdom(IDX,1:2) = [DATAin.LAT(IDX)',DATAin.LON(IDX)'];
    for j=1:numel(IDX)
        F                   = scatteredInterpolant(lonTG(INgrd),latTG(INgrd),WLint(:,j),'linear','none');
        CorrDTdom(IDX(j),3) = F(DATAin.LON(IDX(j)),DATAin.LAT(IDX(j)));
        if isnan(CorrDTdom(IDX(j),3))
            F                   = scatteredInterpolant(lonTG(INgrd),latTG(INgrd),WLint(:,j),'nearest','none');
            CorrDTdom(IDX(j),3) = F(DATAin.LON(IDX(j)),DATAin.LAT(IDX(j)));
        end
    end
end

%Convert cm --> m
CorrDTdom(:,3) = CorrDTdom(:,3)./100;

IDX = DATA3a.WFc' == 0 & DATA3a.MF <= 4 & abs(DATA3a.HEI-CorrDTdom(:,3)-DATAin.NLGEO2018') <= .1;
Histogram(double(DATA3a.HEI(IDX)-CorrDTdom(IDX,3)-DATAin.NLGEO2018(IDX)'),1,20,'Differences [meters]')
figure,plot(DATA3a.LAT(IDX),DATA3a.HEI(IDX)-CorrDTdom(IDX,3)-DATAin.NLGEO2018(IDX)','.')
