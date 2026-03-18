function [HDR,Start_Lat,Start_Lon,Stop_Lat,Stop_Lon] = Cryo_HDR_read(full_filename)

%Read HDR file
HDR = cell(1);
fid = fopen(full_filename);
IDX = 0;
while ~feof(fid)
    IDX        = IDX + 1;
    HDR{IDX,1} = strtrim(fgetl(fid));
end
fclose(fid);

%Extract Start_Lat, Start_Long, Stop_Lat, and Stop_Long
IDX       = strncmpi(HDR,'<Start_Lat',10);
Start_Lat = str2double(regexprep(HDR{IDX},{'<Start_Lat  unit="10-6 deg">','</Start_Lat>'},''))/1E6;
IDX       = strncmpi(HDR,'<Start_Long',11);
Start_Lon = str2double(regexprep(HDR{IDX},{'<Start_Long  unit="10-6 deg">','</Start_Long>'},''))/1E6;
IDX       = strncmpi(HDR,'<Stop_Lat',9);
Stop_Lat  = str2double(regexprep(HDR{IDX},{'<Stop_Lat  unit="10-6 deg">','</Stop_Lat>'},''))/1E6;
IDX       = strncmpi(HDR,'<Stop_Long',10);
Stop_Lon  = str2double(regexprep(HDR{IDX},{'<Stop_Long  unit="10-6 deg">','</Stop_Long>'},''))/1E6;

end
