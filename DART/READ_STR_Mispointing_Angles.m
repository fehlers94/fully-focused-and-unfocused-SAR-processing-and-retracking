function mpSTR = READ_STR_Mispointing_Angles(tSTR)

%Reads star tracker (STR) mispointing angles dataset available at
%https://earth.esa.int/web/guest/missions/esa-eo-missions/cryosat/str-attref
%and stores time series as mat-file.

%% Settings
LoadCommonSettings

%Path to STR data
MisPointingPath = fullfile(fullfile(PathDATA,'RadAlt','CryoSat','STR_MispointingAngles'));

%% Read/store data 
for YEAR = unique(tSTR.Year)

    %Get file names
    FNames = cell(numel(tSTR),1);
    for i = 1:numel(tSTR)
        try
            Year      = num2str(tSTR.Year(i));
            TMP       = dir(fullfile(MisPointingPath,Year,sprintf('CS_OFFL_STR_ATTREF_%4i%02i%02i*.EEF',tSTR.Year(i),tSTR.Month(i),tSTR.Day(i))));
            if isempty(TMP)
                TMP       = dir(fullfile(MisPointingPath,Year,sprintf('CS_TEST_STR_ATTREF_%4i%02i%02i*.EEF',tSTR.Year(i),tSTR.Month(i),tSTR.Day(i))));
            end
            FNames{i} = [Year,filesep,TMP.name];
        catch
            error('%s not available!',fullfile(MisPointingPath,Year,sprintf('CS_*_STR_ATTREF_%4i%02i%02i*.EEF',tSTR.Year(i),tSTR.Month(i),tSTR.Day(i))));
        end
    end

    %Read/store time series
    YMDHMS                    = deal(nan(95279,6));
    [TIME,PITCH,ROLL,YAW]     = deal(nan(95279,numel(FNames)));
    [fTIME,fPITCH,fROLL,fYAW] = deal(nan(95279,1));
    for k = 1:numel(FNames)
        [FNpth,FName,~] = fileparts(FNames{k});
        try
            load(fullfile(MisPointingPath,FNpth,[FName,'.mat']));
        catch
            %Open file, skip header, check nr of entries
            fid   = fopen(fullfile(MisPointingPath,FNames{k}),'r','b');
            count = cell2mat(textscan(fid,'<List_of_Attitude_Angles count="%f">','Headerlines',118));
            if count ~= 95279, error('Check %s',FNames{k}), end
            
            %Read blocks
            % <Attitude_Angles source="STR_2">
            % <Time ref="TAI">TAI=2011-01-01T21:56:00.000000</Time>
            % <Pitch unit="deg">0.074884</Pitch>
            % <Roll unit="deg">-0.184723</Roll>
            % <Yaw unit="deg">-0.070827</Yaw>
            % </Attitude_Angles>
            for i=1:count
                tline       = fgetl(fid); % <Attitude_Angles source="STR_2">
                YMDHMS(i,:) = fscanf(fid,'					<Time ref="TAI">TAI=%f-%f-%fT%f:%f:%f</Time>',[1 6]);
                if contains(tline,'NOMINAL')
                    fgetl(fid);fgetl(fid);fgetl(fid);
                else
                    fPITCH(i)   = fscanf(fid,'					<Pitch unit="deg">%f</Pitch>',1);
                    fROLL(i)    = fscanf(fid,'					<Roll unit="deg">%f</Roll>',1);
                    fYAW(i)     = fscanf(fid,'					<Yaw unit="deg">%f</Yaw>',1);
                end
                fgetl(fid);fgetl(fid);  % </Attitude_Angles>
            end
            fTIME = datenum(YMDHMS(:,1),YMDHMS(:,2),YMDHMS(:,3),YMDHMS(:,4),YMDHMS(:,5),YMDHMS(:,6));
            
            %Save as mat-file
            save(fullfile(MisPointingPath,FNpth,[FName,'.mat']),'fTIME','fPITCH','fROLL','fYAW');
        end
        
        %Copy vectors to arrays
        TIME(:,k) = fTIME; PITCH(:,k) = fPITCH; ROLL(:,k) = fROLL; YAW(:,k) = fYAW;
    end

    %Copy unique values to array
    [~,IDX]         = unique(TIME(:));
    mpSTR.('Time')  = TIME(IDX);
    mpSTR.('Pitch') = PITCH(IDX);
    mpSTR.('Roll')  = ROLL(IDX);
    mpSTR.('Yaw')   = YAW(IDX);
    
end

end
