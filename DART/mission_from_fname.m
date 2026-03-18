function mission = mission_from_fname(s)
%         [filepath,name,~] = fileparts(fname);
%         if contains(name, 'measurement')
%             [~,name,~] = fileparts(filepath);
%             mission=split(name, '_');
%             mission=mission{1};
%         else
%             mission=split(name, '_');
%             mission=mission{1};
%         end
    if contains(s, 'S3A')
        mission = 'S3A';
    elseif contains(s, 'S6A')
        mission = 'S6A';
    elseif contains(s, 'S3B')
        mission = 'S3B';
    elseif contains(s, 'CS')
        mission = 'CS';
    end
end
