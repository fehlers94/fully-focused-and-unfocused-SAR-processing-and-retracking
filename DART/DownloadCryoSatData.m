function DownloadCryoSatData(MODE,LEV,YY,MM,DOM)

%DownloadCryoSatData downloads the L1b/L2 data from the cryosat data server

% YYYY = 2010:2019;
% sprintf('qsub-eval -l DownloadL1_%i.log "mlab -eval \\"DownloadCryoSatData([],''L1'',%i)\\""\n',[YYYY;YYYY])
% sprintf('qsub-eval -l DownloadL2_%i.log "mlab -eval \\"DownloadCryoSatData([],''L2'',%i)\\""\n',[YYYY;YYYY])
% sprintf('qsub-eval -l DownloadL1_%i.log "mlab -eval \\"DownloadCryoSatData(''SAR'',''L1'',%i,[],[51 55;2 10])\\""\n',[YYYY;YYYY])

%% Settings
%General settings
LoadCommonSettings

%Computation settings
defval('MODE','SIN')                   %Identifies mode of data to be downloaded "LRM"/"SAR"/"SIN"
defval('LEV','L1')                     %Identifies level of data to be downloaded "L1"/"L2"/"L2I"
defval('YY',[])                        %Identify year for which data should be downloaded [integer]
defval('MM',[])                        %Identify month for which data should be downloaded [integer]
defval('DOM',[59 84;-74 -10])          %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]

%Set paths
PathDataMAIN = fullfile(PathDATA,'RadAlt','CryoSat');                %MAIN data folder
PathLEV      = fullfile(PathDataMAIN,sprintf('SIR_%s_%s',MODE,LEV)); %Path to L1b/2 data
PathHDR      = fullfile(PathDataMAIN,'science-pds.cryosat.esa.int',sprintf('SIR_%s_%s',MODE,LEV));

%% Download all HDR files
cd(PathDataMAIN);
if isempty(YY)
    system(sprintf('wget --reject=DBL -m --user="cryosat246" --password="tOLPiTWn" ftp://science-pds.cryosat.esa.int/%s/',sprintf('SIR_%s_%s',MODE,LEV)));
else
    if isempty(MM)
        system(sprintf('wget --reject=DBL -m --user="cryosat246" --password="tOLPiTWn" ftp://science-pds.cryosat.esa.int/%s/%i/',sprintf('SIR_%s_%s',MODE,LEV),YY));
    else
        system(sprintf('wget --reject=DBL -m --user="cryosat246" --password="tOLPiTWn" ftp://science-pds.cryosat.esa.int/%s/%i/%02i/',sprintf('SIR_%s_%s',MODE,LEV),YY,MM));
    end
end
cd(Pth_m);

%% Download all DBL files that comprise data within DOM
if isempty(YY), YY = dir(fullfile(PathHDR,'2*')); else YY = dir(fullfile(PathHDR,sprintf('%i*',YY))); end
for i=1:numel(YY)
    MM = dir(fullfile(PathHDR,YY(i).name));
    MM = MM(3:end);
    MM = MM([MM.isdir]);
    for j=1:numel(MM)
        HDRfnames = dir(fullfile(PathHDR,YY(i).name,MM(j).name,'*.HDR'));
        for k=1:numel(HDRfnames)
            %Read HDR file
            [~,Start_Lat,Start_Lon,Stop_Lat,Stop_Lon] = Cryo_HDR_read(fullfile(PathHDR,YY(i).name,MM(j).name,HDRfnames(k).name));
            [latout,lonout]                           = interpm([Start_Lat;Stop_Lat],[Start_Lon;Stop_Lon],1);
            
            %Download DBL file if any point in latout/lonout is inside DOM
            if any(ingeoquad(latout,lonout,DOM(1,:),DOM(2,:)))
                PthFTP = sprintf('ftp://science-pds.cryosat.esa.int/SIR_%s_%s/%s/%s/%s',MODE,LEV,YY(i).name,MM(j).name,regexprep(HDRfnames(k).name,'.HDR','.DBL'));
                system(sprintf('wget -P %s --user="cryosat246" --password="tOLPiTWn" %s %s',fullfile(PathLEV,YY(i).name,MM(j).name),PthFTP));
                movefile(fullfile(PathHDR,YY(i).name,MM(j).name,HDRfnames(k).name),fullfile(PathLEV,YY(i).name,MM(j).name));
            end
        end
    end
end

end