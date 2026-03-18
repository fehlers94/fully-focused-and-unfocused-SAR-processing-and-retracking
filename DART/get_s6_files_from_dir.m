function [ncs, cycle_passes] = get_s6_files_from_dir(basedir, target_cycle, target_pass)
    filelist = dir(fullfile(basedir, '**/*.nc'));
    filelist = filelist(~[filelist.isdir]);  %remove folders from list

    ncs = {};
    cycle_passes = {};
    j = 1;
    for i=1:numel(filelist)
      fpath = fullfile(filelist(i).folder,filelist(i).name);
      [path, fname, fext] = fileparts(fpath);
      
      if contains(fname, 'measurement')
          [~, base_id] = fileparts(path);
      else
          base_id = fname;
      end
      
      parts = split(base_id, '_');

      if contains(base_id, '_1A_') || contains(base_id, '_1B_') 
          cycle = str2num(parts{14});
          pass = str2num(parts{15});
      elseif contains(base_id, '_2_') 
          cycle = str2num(parts{9});
          pass = str2num(parts{10});
      end

      if (~target_pass || (target_pass && pass == target_pass)) && (~target_cycle || (target_cycle && cycle == target_cycle)) && ~contains(fpath, 'RED')
        ncs{j} = fpath;
        cycle_passes{j} = [cycle pass];
        j = j + 1;
      end

    end
end