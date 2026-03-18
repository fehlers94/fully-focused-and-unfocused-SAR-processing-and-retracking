classdef FFSAR_Processor_sa < handle
    properties
        FName, DOM, proc_sets
        mission, is_cs, is_s3, is_s6
        nc_grp, finfo
        DDAcf, IDX_mask_cs
        
        CONST
        
        Nl, r_res_zp, W, t, tp_vec, T, Ts, BRI %preliminary data
        tm_num_burst_select, ind_first_rel_dom_burst, ind_last_rel_dom_burst
        t_burst, t_burst_corr
%       t_pulse, t_select, t_select_corr,
        LAI, VX_b, VY_b, VZ_b, Rtrk_b, lat_b, lon_b, alt_b,
        t_l, lon_l, lat_l, alt_l, Rtrk_l, h_l, h_l_interpolator, pri_l, range_ref_l, inds_nearest_l_to_b %interpolations
        R0_i
        ALLl
        qq, uwphase
        
        CS1a, CS1b
        wav, wav_select
    end
    methods      
        function obj = FFSAR_Processor_sa(FName,DOM,FFSAR_processing_settings,focus_height_interpolator)
            LoadCommonSettings
            [obj.CONST,~,FFSAR_default_processing_settings] = FFSAR_LoadCommonSettings;
            
            defval('proc_sets', FFSAR_default_processing_settings);
            obj.proc_sets = FFSAR_processing_settings;
            
            obj.FName = FName;
            defval('DOM',[])
            obj.DOM = DOM;
            
            defval('focus_height_interpolator', []);
            obj.h_l_interpolator = focus_height_interpolator;

            obj.mission = mission_from_fname(FName);
            obj.is_cs = strcmp(obj.mission,'CS');
            obj.is_s3 = contains(obj.mission,'S3');
            obj.is_s6 = contains(obj.mission,'S6');

            % TODO: introduce dynamic T? obj.T=0.886*obj.DDAcf.c*obj.DDAcf.alt_nom/(2*obj.DDAcf.fc*obj.DDAcf.v_nom*obj.DDAcf.res_a);
            if obj.proc_sets.integration_time
                obj.T = obj.proc_sets.integration_time;
            else
                if obj.is_cs
                    obj.T=1.9;
                elseif obj.is_s3
                    obj.T=2.1;
                elseif obj.is_s6
                    obj.T=3.413;
                end               
            end

            obj.nc_grp='/';
            if obj.is_s6
                obj.nc_grp='data_140/ku/';
            end
            if ~obj.is_cs
                obj.finfo = ncinfo(FName,obj.nc_grp);
            end
        end

        function read_crop_data(obj)
            buffer_time = 1 * obj.T;
            total_buffer_time = (obj.T/2 + buffer_time);
            
            %Get indices to read cropped data (for all missions except CS)
            if ~obj.is_cs
                for i=1:numel(obj.finfo.Variables)
                    if contains(obj.finfo.Variables(i).Name,'lat')
                        latvar = obj.finfo.Variables(i).Name;
                    elseif contains(obj.finfo.Variables(i).Name,'lon')
                        lonvar = obj.finfo.Variables(i).Name;
                    elseif startsWith(obj.finfo.Variables(i).Name,'time') & ~contains(obj.finfo.Variables(i).Name,'plrm')
                        timevar = obj.finfo.Variables(i).Name;
                    end
                end

                lat_nc = ncread(obj.FName,[obj.nc_grp latvar]);
                lon_nc = ncread(obj.FName,[obj.nc_grp lonvar]);
                time_nc = ncread(obj.FName,[obj.nc_grp timevar]);
                
                n_buffer_20 = ceil(1/20 / mean(diff(time_nc)));
                IDX_mask = ingeoquad(lat_nc(:),lon_nc(:),obj.DOM(1,:),obj.DOM(2,:));

                IDX = find(IDX_mask);                
                IDX = [(min(IDX) - n_buffer_20:min(IDX)-1)';IDX;(max(IDX)+1:max(IDX) + n_buffer_20)']; % pre- and append buffer to get full integer 20-hz records after averaging
                ind_first_dom_burst = min(IDX);
                count_dom_bursts = max(IDX)- ind_first_dom_burst + 1;
                
                % pre-long and extend ind_first_burst/count_all_bursts indices to accomodate T/2 at the edges + some buffer_time
                [~, ind_start_w_buffer] = min(abs(time_nc - (time_nc(ind_first_dom_burst) - total_buffer_time)));
                [~, ind_end_w_buffer] = min(abs(time_nc - (time_nc(ind_first_dom_burst + count_dom_bursts - 1) + total_buffer_time)));
                count_all_bursts = ind_end_w_buffer - ind_start_w_buffer;
                
                obj.ind_first_rel_dom_burst = ind_first_dom_burst - ind_start_w_buffer + 1;
                obj.ind_last_rel_dom_burst = obj.ind_first_rel_dom_burst + count_dom_bursts; % normally -1, but add some more look locs here
                
                clear lat_nc lon_nc time_nc
            end

            %Read data
            switch obj.mission
                case 'CS'
                    [~,obj.CS1a] = Cryo_L1b_read(obj.FName);
                case {'S3A', 'S3B', 'S6A'}
                    obj.CS1a = S3_S6_L1a_read(obj.FName,ind_start_w_buffer,count_all_bursts);
                    obj.CS1a.GEO.Start_Time = datetime(2000,1,1,0,0,obj.CS1a.TIME(obj.ind_first_rel_dom_burst), 'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
                otherwise
                    error('mission: %s not recognized',obj.mission)
            end

            %Load Parameter Settings (Reference ellipsoid, SIRAL Characteristics, etc.)
            if obj.is_cs
                obj.DDAcf = DDA_ConfigFile(obj.mission,'SAR',obj.CS1a.GEO.Baseline_ID);
            else
                obj.DDAcf = DDA_ConfigFile(obj.mission,'SAR');
            end

            %If segment crosses -180 meridian, wrap angle in degrees to [0 360]
            if obj.is_cs
                if max(abs(diff(obj.CS1a.GEO.LON(obj.CS1a.FBR.N_pulses(:) > 0)))) > 350
                    obj.CS1a.GEO.LON = wrapTo360(obj.CS1a.GEO.LON);
                end
            else
                if max(abs(diff(obj.CS1a.GEO.LON))) > 350
                    obj.CS1a.GEO.LON = wrapTo180(obj.CS1a.GEO.LON);
                end
            end

            %Correct window delay for DORIS USO drift (pp 20 CryoSat Product Handbook)
            if obj.is_cs
                if isequal(obj.CS1a.GEO.Baseline_ID,'C')
                    obj.CS1a.MEA.win_delay = obj.CS1a.MEA.win_delay.*(obj.CS1a.GEO.USO+1);
                end
            end

            %Correct window delay for internal path delay
            if obj.is_cs
                obj.CS1a.MEA.win_delay = obj.CS1a.MEA.win_delay + (2*obj.CS1a.MEA.ins_range_corr_rx_tx/obj.DDAcf.c);
            end

            %Transform acquisition time to datenum format
            if obj.is_cs
                obj.CS1a.('TIME') = datenum('2000','yyyy') + obj.CS1a.GEO.TAI.days(:) + obj.CS1a.GEO.TAI.secs(:)./86400 + obj.CS1a.GEO.TAI.microsecs(:)./1e6./86400;
            % else %do nothing here, because we have the time vector explicitely given
            % in the netCDF file
            end

            %Crop data (CS)
            if obj.is_cs
                if isempty(obj.DOM)
                    obj.IDX_mask_cs = true(numel(obj.CS1a.GEO.LAT),1);
                else
                    if ischar(obj.DOM)
                        switch obj.DOM
                            case {'JakobsHavn','jakobshavn'}
                                %Load polygon that outlines the JakobsHavn glacier
                                load(fullfile(PathDATA,'Geography','jakobshavn_basin.mat'));
                                obj.IDX_mask_cs = inpolygon(obj.CS1a.GEO.LON(:),obj.CS1a.GEO.LAT(:),polyLon,polyLat);
                            otherwise
                                error('DOMain of interest not reckognized')
                        end
                    else
                        obj.IDX_mask_cs = ingeoquad(obj.CS1a.GEO.LAT(:),obj.CS1a.GEO.LON(:),obj.DOM(1,:),obj.DOM(2,:));
                    end
                    
                    IDX = find(obj.IDX_mask_cs);
                    ind_first_dom_burst = min(IDX);
                    count_dom_bursts = max(IDX)- ind_first_dom_burst + 1;
                    
                    obj.ind_first_rel_dom_burst = ind_first_dom_burst;
                    obj.ind_last_rel_dom_burst = ind_first_dom_burst + count_dom_bursts; % normally -1, but add some more look locs here
                    
                    % pre-long and extend startLoc/count indices to accomodate T/2 at the edges + some buffer_time
                    time_nc = obj.CS1a.TIME;
                    [~, ind_start_w_buffer] = min(abs(time_nc - (time_nc(ind_first_dom_burst) - total_buffer_time)));
                    [~, ind_end_w_buffer] = min(abs(time_nc - (time_nc(ind_first_dom_burst + count_dom_bursts - 1) + total_buffer_time)));
                    
                    obj.IDX_mask_cs(ind_start_w_buffer:ind_end_w_buffer) = true(ind_end_w_buffer+ind_start_w_buffer-1, 1);
                end
            end
        end
        
        function mask_data(obj)
            if obj.is_cs
                IDXfd = (obj.CS1a.GEO.MCD_FLAG.Block_Degraded | obj.CS1a.GEO.MCD_FLAG.Blank_Block | ...
                    obj.CS1a.GEO.MCD_FLAG.Datation_Degraded     | obj.CS1a.GEO.MCD_FLAG.Orbit_Propag_Err | ...
                    obj.CS1a.GEO.MCD_FLAG.Orbit_File_Change     | obj.CS1a.GEO.MCD_FLAG.Orbit_Discontinuity | ...
                    obj.CS1a.GEO.MCD_FLAG.Echo_Saturation       | obj.CS1a.GEO.MCD_FLAG.Other_Echo_Err | ...
                    obj.CS1a.GEO.MCD_FLAG.Rx1_Err_SARin         | obj.CS1a.GEO.MCD_FLAG.Rx2_Err_SARin | ...
                    obj.CS1a.GEO.MCD_FLAG.Wind_Delay_Incon      | obj.CS1a.GEO.MCD_FLAG.AGC_Incon | ...
                    obj.CS1a.GEO.MCD_FLAG.CAL1_Corr_Miss        | obj.CS1a.GEO.MCD_FLAG.CAL1_Corr_IPF | ...
                    obj.CS1a.GEO.MCD_FLAG.DORIS_USO_Corr        | obj.CS1a.GEO.MCD_FLAG.Complex_CAL1_Corr_IPF | ...
                    obj.CS1a.GEO.MCD_FLAG.TRK_ECHO_Err          | obj.CS1a.GEO.MCD_FLAG.RX1_ECHO_Err | ...
                    obj.CS1a.GEO.MCD_FLAG.RX2_ECHO_Err          | obj.CS1a.GEO.MCD_FLAG.NPM_Incon | ...
                    obj.CS1a.GEO.MCD_FLAG.Attitude_Corr_Missing | obj.CS1a.GEO.MCD_FLAG.CAL1_Corr_Type);

                obj.CS1a.GEO.IDXfd = ~obj.IDX_mask_cs | IDXfd(:);
            else
                obj.CS1a.GEO.IDXfd = false(numel(obj.CS1a.GEO.LAT), 1);
                
                if obj.is_s6
                    obj.CS1a.GEO.IDXfd = obj.CS1a.GEO.IDXfd | logical(obj.CS1a.GEO.MCD_FLAG); %consider any set flag (value~=0) as badF
                end
            end
        end
        
        function calibrate_fbr_data(obj)
            if obj.is_cs
                obj.CS1a = Calibrate_CryoSat_FBRdata(obj.CS1a,obj.DDAcf);
            elseif obj.is_s3
                obj.CS1a = Calibrate_S3_L1Adata(obj.CS1a,obj.DDAcf);
            elseif obj.is_s6
                obj.CS1a = Calibrate_S6_L1Adata(obj.CS1a,obj.DDAcf);
            end
        end
        
        function setup_preliminaries(obj)
            % Some settings for looks
            IDX1=true(length(obj.CS1a.GEO.LAT),1);
            if obj.is_cs
                IDX1=obj.CS1a.FBR.N_pulses ~= 0;
            end
            lat=obj.CS1a.GEO.LAT(IDX1);
            lon=obj.CS1a.GEO.LON(IDX1);
            
            
            L=obj.T*median(obj.CS1a.GEO.V.V(IDX1))/2;
            H=median(obj.CS1a.GEO.H(IDX1));
            res_a=(((L^2+H^2)^0.5+obj.DDAcf.lambda0/4)^2-H^2)^0.5-L;
            Re_local=rsphere('euler',lat(1),lon(1),lat(end),lon(end),[obj.DDAcf.RefEll.SemimajorAxis obj.DDAcf.RefEll.Eccentricity]);
            obj.Nl=median(obj.CS1a.GEO.V.V(:))/(res_a*((Re_local+median(obj.CS1a.GEO.H(:)))/Re_local)); % number of looks per seconds [1/s]

            t_res=1/obj.DDAcf.B; % time resolution of the waveform [s]
            r_res=t_res*obj.DDAcf.c/2; % range resolution of the waveform [m]
            obj.r_res_zp=r_res/obj.DDAcf.os_ZP;

            if obj.DDAcf.ApplyWindow
                if strcmpi(obj.DDAcf.Window,'hamming')
                    obj.W = hamming(obj.DDAcf.Nb,'periodic')';
                elseif strcmpi(obj.DDAcf.Window,'hanning')
                    obj.W = hanning(obj.DDAcf.Nb)';
                else
                    error('Window not recognized!');
                end
            else
                obj.W = ones(1,obj.DDAcf.Nb);
            end

            % define fast-time vector
            if ~obj.is_s6
                obj.Ts=obj.DDAcf.tau_u/obj.DDAcf.Np;
            else
                obj.Ts=obj.DDAcf.B/obj.DDAcf.Bt*obj.DDAcf.tau_u/obj.DDAcf.Np;
            end
            obj.t=0:obj.Ts:(obj.DDAcf.Np-1)*obj.Ts; % fast time samples [s]
            obj.t=obj.t-max(obj.t)*(0.5);
            if obj.is_s6
                pri = obj.CS1a.MEA.PRI';
            else
                pri = obj.DDAcf.PRI;
            end
            obj.tp_vec=(-(obj.DDAcf.Nb-1)/2:1:(obj.DDAcf.Nb-1)/2)'*pri; % vector of pulse timing [s]
        end
        
        function do_interpolations(obj)
            %% Interpolate fields to ground locations
            if obj.is_cs
                obj.t_burst=obj.CS1a.GEO.Elapsed_Time(~obj.CS1a.GEO.IDXfd); % original time vector
            else
                obj.t_burst=obj.CS1a.TIME(~obj.CS1a.GEO.IDXfd);
            end
            obj.t_burst=obj.t_burst-obj.t_burst(1);

            if obj.is_s6
                % do outlier detection here, before using mean statistic
                burst_time_diffs = diff(obj.t_burst);
                BRI=median(burst_time_diffs);
                
                % look for jumps
                jump_idx = find((0.95 > burst_time_diffs/BRI)|(burst_time_diffs/BRI > 1.05)); % the median estimate is a bit coarse, but consider anything that is more than 5% of an approximated to be a jump BRI
                
                burst_time_diffs(jump_idx) = [];
                
                % now overwrite BRI using the mean for higher accuracy
                obj.BRI=mean(burst_time_diffs);
            else
                obj.BRI=obj.DDAcf.BRI;
            end
            
            if ~obj.is_s6
                obj.t_burst_corr=(0:length(obj.t_burst)-1)'*obj.BRI; % corrected time vector, but not considering any time jumps due to satellite's 'coffee breaks'

                % look for jumps in obj.t_burst
                burst_time_diffs = diff(obj.t_burst);
                jump_idx = find((0.999 > burst_time_diffs/obj.BRI)|(burst_time_diffs/obj.BRI > 1.001));

                t_burst_jump_cor = zeros(size(obj.t_burst));
                t_burst_jump_cor(jump_idx+1) = burst_time_diffs(jump_idx) - obj.BRI; % idx+1 because of forward diff operation and - BRI because we replace a BRI with an arbitrary jump
                t_burst_jump_cor = cumsum(t_burst_jump_cor);
                obj.t_burst_corr = obj.t_burst_corr + t_burst_jump_cor; % corrected time vector
            else
                obj.t_burst_corr = obj.t_burst;
            end
            
            %Construct array of valid entries of fields on a burst basis
            ma = ~obj.CS1a.GEO.IDXfd;
            ALLb = [
                obj.CS1a.GEO.LAT(ma),...
                obj.CS1a.GEO.LON(ma),...
                obj.CS1a.GEO.H(ma),...
                obj.CS1a.GEO.H_rate(ma),...
                obj.CS1a.GEO.V.Vx(ma),...
                obj.CS1a.GEO.V.Vy(ma),...
                obj.CS1a.GEO.V.Vz(ma),...
                obj.CS1a.GEO.V.V(ma),...
                obj.CS1a.GEO.TAI.days(ma),...
                obj.CS1a.GEO.TAI.secs(ma),...
                obj.CS1a.GEO.TAI.secs(ma).*1E6,...
                obj.CS1a.MEA.win_delay(ma),...
                obj.CS1a.GEO.BaseLine.X(ma),...
                obj.CS1a.GEO.BaseLine.Y(ma),...
                obj.CS1a.GEO.BaseLine.Z(ma),...
                obj.CS1a.GEO.Beam.X(ma),...
                obj.CS1a.GEO.Beam.Y(ma),...
                obj.CS1a.GEO.Beam.Z(ma),...
%                 obj.CS1a.GEO.Antenna_Bench_Pitch(ma),...
%                 obj.CS1a.GEO.Antenna_Bench_Roll(ma),...
                ];

            if ~obj.is_cs
                ALLb = [ALLb,...
                obj.CS1a.GEO.Antenna_Bench_Pitch(ma),...
                obj.CS1a.GEO.Antenna_Bench_Roll(ma),...
                ];
 
            end
            
            %Compute the piecewise polynomial forms of the cubic spline interpolants to the data values ALLb
            dt=1/(obj.Nl*obj.proc_sets.along_track_oversampling); % time between looks
            
            % time vector for look location interpolation
            obj.t_l = obj.t_burst(obj.ind_first_rel_dom_burst):dt:obj.t_burst(obj.ind_last_rel_dom_burst);
            
            % cut t_l to integers of n_combe
            if obj.proc_sets.combine_n_look_locations > 1
                obj.t_l = obj.t_l(1:floor(length(obj.t_l) / obj.proc_sets.combine_n_look_locations) * obj.proc_sets.combine_n_look_locations);
            end
            
            % get the nearest neighbors
            obj.inds_nearest_l_to_b = nearestpoint(obj.t_l, obj.t_burst, 'nearest');
            
            pp = spline(obj.t_burst,ALLb');
            obj.ALLl = ppval(pp,obj.t_l);

            % read relevant data for FF-SAR
            obj.lon_b=wrapTo180(obj.CS1a.GEO.LON(~obj.CS1a.GEO.IDXfd)); 
            obj.lat_b=obj.CS1a.GEO.LAT(~obj.CS1a.GEO.IDXfd); 
            obj.alt_b=obj.CS1a.GEO.H(~obj.CS1a.GEO.IDXfd); 
            obj.VX_b=obj.CS1a.GEO.V.Vx(~obj.CS1a.GEO.IDXfd); 
            obj.VY_b=obj.CS1a.GEO.V.Vy(~obj.CS1a.GEO.IDXfd); 
            obj.VZ_b=obj.CS1a.GEO.V.Vz(~obj.CS1a.GEO.IDXfd); 
            % roll_b=obj.CS1a.GEO.BaseLine.X(~obj.CS1a.GEO.IDXfd);
            % pitch_b=-obj.CS1a.GEO.Beam.Y(~obj.CS1a.GEO.IDXfd);
            obj.LAI=obj.CS1a.MEA.LAI(~obj.CS1a.GEO.IDXfd); 
            
            % tracking window
            if obj.is_cs
                wd_b=obj.CS1a.MEA.win_delay(~obj.CS1a.GEO.IDXfd);
                obj.Rtrk_b=wd_b*obj.DDAcf.c/2; % one-way window delay in meters
            else
                obj.Rtrk_b=obj.CS1a.MEA.ref_range;
            end
            
            % grab L1a waveforms
            obj.wav=double(obj.CS1a.FBR.data(:,:,~obj.CS1a.GEO.IDXfd));

            % interpolate burst data onto new look locations
            obj.lon_l  = polyval(polyfit(obj.t_burst,obj.lon_b,4),obj.t_l)';
            obj.lat_l  = polyval(polyfit(obj.t_burst,obj.lat_b,4),obj.t_l)';
            obj.alt_l  = polyval(polyfit(obj.t_burst,obj.alt_b,4),obj.t_l)';
            obj.Rtrk_l = polyval(polyfit(obj.t_burst,obj.Rtrk_b,4),obj.t_l)';
            
            % switch on focal point height estimation, either 'pchip' or
            % 'piecewise constant', the former might lead to slight to extreme waveheight
            % biases, depending on the mission tracker behaviour
            if strcmp(obj.proc_sets.focal_point_height_interpolation,'piecewise_constant')
                h_trk = obj.alt_b - obj.Rtrk_b; % height over ref. ell. to which the bursts are references
                h_b_const = zeros(size(h_trk));
                
                radius = 2.3; % in [m]; tolerated difference between consequitive jumps of h_trk before a new constant level is approached. Naturally, about 2 m seems to be the on board tracker threshold for S3.

                idx_start = 1;
                while idx_start < numel(h_trk)
                    piece_max = cummax(h_trk(idx_start:end));
                    piece_min = cummin(h_trk(idx_start:end));

                    idx_rel_end = find(abs(piece_max - piece_min) < radius, 1, 'last');
                    h_b_const(idx_start : idx_start+idx_rel_end-1) = 1/2*( piece_max(idx_rel_end) + piece_min(idx_rel_end) );

                    idx_start = idx_start + idx_rel_end;
                end
                % now interpolate h_l onto pulse timings and get Rtrk_l and so on
                obj.h_l = interp1(obj.t_burst,h_b_const,obj.t_l,'nearest');
                obj.range_ref_l = obj.alt_l - obj.h_l';
                obj.Rtrk_l = obj.range_ref_l;
                
            elseif strcmp(obj.proc_sets.focal_point_height_interpolation,'pchip_smooth')
                % some alternative pchip interpolation for Rtrk_l (through the center times
                % of the tracker range settings to mitigate the problems with relatively bumpy surfaces / steep slopes)
                trk_change_flag = (diff(obj.Rtrk_b) ~= 0);
                if sum(trk_change_flag) > 4
                    IDXtrk = find(trk_change_flag);
                    t_b_trk = obj.t_burst(IDXtrk);

                    % get center of tracking setting
                    % BRI/2 needs to be added for centering because of the 'righthanded' diff operation
                    if ~obj.is_s6
                        t_b_center = 0.5 * (t_b_trk(1:end-1) + t_b_trk(2:end) + obj.DDAcf.BRI);
                    else
                        t_b_center = 0.5 * (t_b_trk(1:end-1) + t_b_trk(2:end) + obj.BRI);
                    end

                    Rtrk_center = obj.Rtrk_b(IDXtrk(2:end));
                    trk_interpolator = pchip(t_b_center,Rtrk_center);
                    obj.Rtrk_l = ppval(trk_interpolator,obj.t_l)';
                end
            end
                
            % nearest neighbor interpolation to get interpolated pri values
            if obj.is_s6
                obj.pri_l = obj.CS1a.MEA.PRI(obj.inds_nearest_l_to_b);
            end
            
            % define range reference, to which we refer our look locations
            % this assumes that our tracker range Rtrk_l smoothly  
            obj.range_ref_l=obj.Rtrk_l;
            
            % estimates of elevation of look locations
            obj.h_l=obj.alt_l-obj.range_ref_l;
            
            % if an interpolator is given for h_l, then overwrite all
            % interpolations from before:
            if ~isempty(obj.h_l_interpolator)
                obj.h_l = obj.h_l_interpolator(obj.lat_l);
                obj.range_ref_l = obj.alt_l - obj.h_l;
                obj.Rtrk_l = obj.range_ref_l;
            end
           
            % satellite location at the zenith of the look locations
            % [x_l,y_l,z_l]=geodetic2enu(obj.lat_l,obj.lon_l,obj.alt_l,obj.lat_l,obj.lon_l,obj.h_l,wgs84Ellipsoid,'degrees');
            
            % compute shortest distance between track and look location
            % optimization: take directly Rtrk_l as it is the zenith distance between the satellite and the focus point
            % obj.R0_i=sqrt(x_l.^2+y_l.^2+z_l.^2);
            obj.R0_i=obj.range_ref_l;
            
            % figure; plot(obj.lat_l, obj.h_l); hold on; plot(obj.lat_b, obj.alt_b - obj.Rtrk_b); legend(['pulse-basis (pchip-interpolated), std=' num2str(std(obj.h_l))],'burst-basis'); xlabel('lat [deg]'); ylabel('h_l [m]')
        end
        
        function init_L1b_struct(obj)
            [~,obj.CS1b] = Cryo_L1b_struct('SIR_SAR_L1.DBL',numel(obj.h_l),obj.DDAcf,1);
            obj.CS1b.SAR.data = zeros(obj.DDAcf.Np*obj.DDAcf.os_ZP, numel(obj.h_l)); % required because S6 uses 256 bins per waveform, which is not considered in Cryo_L1b_struct func
            obj.CS1b.SAR.dataIQ = zeros(obj.DDAcf.Np*obj.DDAcf.os_ZP, numel(obj.h_l));
            if obj.proc_sets.split_aperture
                obj.CS1b.SAR.dataIQ_050T_025T = zeros(obj.DDAcf.Np*obj.DDAcf.os_ZP, numel(obj.h_l));
                obj.CS1b.SAR.dataIQ_025T_000T = zeros(obj.DDAcf.Np*obj.DDAcf.os_ZP, numel(obj.h_l));
                obj.CS1b.SAR.dataIQ_000T_025T = zeros(obj.DDAcf.Np*obj.DDAcf.os_ZP, numel(obj.h_l));
                obj.CS1b.SAR.dataIQ_025T_050T = zeros(obj.DDAcf.Np*obj.DDAcf.os_ZP, numel(obj.h_l));
            end
            if obj.proc_sets.output_pseudo_delay_doppler_processing
                obj.CS1b.SAR.data_pseudoDD = zeros(obj.DDAcf.Np*obj.DDAcf.os_ZP, numel(obj.h_l));
                obj.CS1b.SAR.data_pseudoDD_zeroDop = zeros(obj.DDAcf.Np*obj.DDAcf.os_ZP, numel(obj.h_l));
            end
            obj.CS1b.SAR.num_bursts = zeros(1, numel(obj.h_l));
            obj.CS1b.SAR.beam_param = squeeze(obj.CS1b.SAR.beam_param);
            obj.CS1b.SAR.BeamAngle = NaN(490, numel(obj.h_l)); % Sentinel-3 has maximum 256 doppler beams on one ground point, Sen-6 maximum 490
            obj.CS1b.SAR.stack_mask_start_stop = zeros(490, numel(obj.h_l)); % same
            
            obj.CS1b.SAR.scale_factor_ku = ones(size(obj.ALLl(1,:)));
            
            obj.CS1b.COR.Doppler_correction = zeros(size(obj.ALLl(1,:)));
            
            obj.CS1b.GEO.LAT       = obj.ALLl(1,:);
            obj.CS1b.GEO.LON       = obj.ALLl(2,:);
            obj.CS1b.GEO.H         = obj.ALLl(3,:);
            obj.CS1b.MEA.win_delay = obj.Rtrk_l'/(obj.DDAcf.c/2);
            obj.CS1b.GEO.H_rate    = obj.ALLl(4,:);
            obj.CS1b.GEO.V.Vx      = obj.ALLl(5,:);
            obj.CS1b.GEO.V.Vy      = obj.ALLl(6,:);
            obj.CS1b.GEO.V.Vz      = obj.ALLl(7,:);
            obj.CS1b.GEO.V.V       = obj.ALLl(8,:);
            obj.CS1b.GEO.TAI.days  = obj.ALLl(9,:);
            obj.CS1b.GEO.TAI.secs  = obj.ALLl(10,:);
            obj.CS1b.GEO.TAI.microsecs = obj.ALLl(11,:);
            
            % TODO: check this, writing to obj.CS1a was changed to obj.CS1b
            obj.CS1b.GEO.BaseLine.X= obj.ALLl(13,:);
            obj.CS1b.GEO.BaseLine.Y= obj.ALLl(14,:);
            obj.CS1b.GEO.BaseLine.Z= obj.ALLl(15,:);
            obj.CS1b.GEO.Beam.X    = obj.ALLl(16,:);
            obj.CS1b.GEO.Beam.Y    = obj.ALLl(17,:);
            obj.CS1b.GEO.Beam.Z    = obj.ALLl(18,:);
            
            if obj.is_cs
                %Compute roll, pitch, and yaw angles of the antenna bench (p. 16 of the CryoSat product handbook)
                obj.CS1b.GEO.Antenna_Bench_Roll  = rad2deg(obj.CS1b.GEO.BaseLine.X);
                obj.CS1b.GEO.Antenna_Bench_Yaw   = rad2deg(-obj.CS1b.GEO.BaseLine.Y); 
                obj.CS1b.GEO.Antenna_Bench_Pitch = rad2deg(-obj.CS1b.GEO.Beam.Y);
            else
                obj.CS1b.GEO.Antenna_Bench_Pitch    = obj.ALLl(19,:);
                obj.CS1b.GEO.Antenna_Bench_Roll    = obj.ALLl(20,:);
            end
            
            obj.CS1b.GEO.Start_Time = obj.CS1a.GEO.Start_Time;
            obj.CS1b.GEO.Elapsed_Time = obj.t_l - obj.t_l(1);
            obj.CS1b.GEO.Baseline_ID = 'OWN';
            obj.CS1b.GEO.BRI       = obj.BRI;
            
            if obj.is_s6
                obj.CS1b.MEA.PRI = obj.pri_l';
                obj.CS1b.MEA.p4_mode_flag = obj.CS1a.MEA.p4_mode_flag(obj.inds_nearest_l_to_b)';
            end
            
            %Set flag obj.CS1b.GEO.MCD_FLAG.Block_Degraded;
            obj.CS1b.GEO.MCD_FLAG.Block_Degraded = double(isnan(obj.CS1b.GEO.LAT));
            
            % export some other settings
            obj.CS1b.GEO.os_zp_factor = obj.DDAcf.os_ZP;
            obj.CS1b.GEO.integration_time = obj.T;
        end
        
        function setup_proc(obj)
            obj.read_crop_data()
            obj.mask_data()
            obj.calibrate_fbr_data()
            obj.setup_preliminaries()
            obj.do_interpolations()
            obj.init_L1b_struct()
        end
        
        function inds = look_loc_inds(obj)
            if isempty(obj.h_l)
                error('please run setup_proc() first')
            end
            
            if obj.proc_sets.transponder_test_mode
                tr_latlonalt = obj.proc_sets.transponder_latlonalt;
                d=distance(tr_latlonalt(1),tr_latlonalt(2),obj.lat_l,obj.lon_l,[obj.CONST.a obj.CONST.e]);
                [~,inds]=min(d);
                if inds < obj.proc_sets.combine_n_look_locations/2
                    inds = ceil(obj.proc_sets.combine_n_look_locations/2);
                end
            else
                combine_n_locs = obj.proc_sets.combine_n_look_locations;
                if mod(combine_n_locs,2) == 0
                    inds = (combine_n_locs/2:combine_n_locs:length(obj.h_l));
                else
                    inds = ((combine_n_locs+1)/2:combine_n_locs:length(obj.h_l));
                end
            end
        end
        
        function [t_select, t_pulse] = get_look_loc_times(obj, ii)
            look_loc_inds = obj.look_loc_inds();
            
            i = look_loc_inds(ii);
            q = obj.qq(:,ii);
            
            t_select = obj.t_burst_corr(q) - obj.t_l(i);
            
            if ~isvector(obj.tp_vec)
                tp_vec_select = obj.tp_vec(:,q);
            else
                tp_vec_select = obj.tp_vec;
            end
            
            t_pulse=t_select'+tp_vec_select;
            t_pulse=t_pulse(:);
        
        end
        
        function proc(obj)
            look_loc_inds = obj.look_loc_inds();

            fprintf('Started processing of %d look locations (up to index %d)...\n',numel(look_loc_inds), look_loc_inds(end));
            
            obj.qq = logical(abs(obj.t_burst-obj.t_l(look_loc_inds))<obj.T/2);
            
            % process multiple consecutive waveforms (processing-speed optimization)
            n_looks = obj.proc_sets.combine_n_look_locations;
            f_b=1/obj.T;
            if n_looks > 1
                j = (-(floor(n_looks/2)-1):floor(n_looks/2))';
            else
                j = 0;
            end
            
            is_s3 = obj.is_s3;
            is_cs = obj.is_cs;
            is_s6 = obj.is_s6;
            DDAcf = obj.DDAcf;
            W = obj.W;
            t = obj.t;
            Ts = obj.Ts;
            T = obj.T;
            r_res_zp = obj.r_res_zp;
            FName = obj.FName;
            
            % burst-based data, access via q
            t_burst = obj.t_burst;
            t_burst_corr = obj.t_burst_corr;
            VX_b = obj.VX_b;
            VY_b = obj.VY_b;
            VZ_b = obj.VZ_b;
            Rtrk_b = obj.Rtrk_b;
            wav = obj.wav;
            proc_sets = obj.proc_sets;

            tm_burst_num = [];
            if is_s6
                tm_burst_num = obj.CS1a.MEA.tm_burst_num;
            end
            
            tp_vec = obj.tp_vec;
            lat_b = obj.lat_b;
            lon_b = obj.lon_b;
            alt_b = obj.alt_b;
            
            % look-loc-based data, access via i
            lat_l = obj.lat_l;
            lon_l = obj.lon_l;
            h_l = obj.h_l;
            t_l = obj.t_l;
            pri_l = obj.pri_l;
            Rtrk_l = obj.Rtrk_l;
            R0_i = obj.R0_i;
            %% memory optimisation
%             % temp output vars
%             dataIQ = zeros([numel(look_loc_inds), obj.DDAcf.Np*obj.DDAcf.os_ZP, n_looks]);
%             if proc_sets.split_aperture
%                 dataIQ_0_T = zeros([numel(look_loc_inds), obj.DDAcf.Np*obj.DDAcf.os_ZP, n_looks]);
%                 dataIQ_T_0 = zeros([numel(look_loc_inds), obj.DDAcf.Np*obj.DDAcf.os_ZP, n_looks]);
%             end
% %             dataIQ = zeros(obj.DDAcf.Np*obj.DDAcf.os_ZP, numel(obj.h_l));
% %             data = zeros(obj.DDAcf.Np*obj.DDAcf.os_ZP, numel(obj.h_l));
%             data_pseudoDD = zeros([numel(look_loc_inds), obj.DDAcf.Np*obj.DDAcf.os_ZP, n_looks]);
%             data_pseudoDD_zeroDop = zeros([numel(look_loc_inds), obj.DDAcf.Np*obj.DDAcf.os_ZP, n_looks]);
%             win_delay = zeros(numel(look_loc_inds), n_looks);
%             Doppler_correction = zeros(numel(look_loc_inds), n_looks);
%             n_max_beamangles = 490;
%             BeamAngle = NaN(numel(look_loc_inds), n_max_beamangles, n_looks);
% 
%             num_bursts = zeros(numel(look_loc_inds), n_looks);
%             stack_mask_start_stop = zeros(numel(look_loc_inds), 490, n_looks);
            %%
            if obj.proc_sets.transponder_test_mode && numel(look_loc_inds) == 1
                uwphase = zeros(numel(look_loc_inds), sum(obj.qq)*obj.DDAcf.Nb);
            end
                        
            for ii=1:numel(look_loc_inds)
%             parpool(proc_sets.n_cores_parfor)
%             parfor (ii=1:numel(look_loc_inds))%proc_sets.n_cores_parfor)
                i = look_loc_inds(ii);
                q = obj.qq(:,ii);
                fprintf('processing nadir look location i=%d ...\n',i);

                % selection of bursts (radius of obj.T second around zenith)
                t_select = t_burst(q);
                t_select_corr = t_burst_corr(q);
                
                if is_s6
                    tm_num_burst_select = tm_burst_num(q);
                end
                
                VX_select = VX_b(q);
                VY_select = VY_b(q);
                VZ_select = VZ_b(q);
                Rtrk_select = Rtrk_b(q);
                wav_select = wav(:,:,q);
                
                % now if s6 is used, the pulse repetition can vary and so tp_vec becomes a matrix which needs to be accounted for
                if ~isvector(tp_vec)
                    tp_vec_select = tp_vec(:,q);
                else
                    tp_vec_select = tp_vec;
                end

                % satellite locations at burst in local reference frame of location 'l(i)'
                [x_select,y_select,z_select]=geodetic2enu(lat_b(q), lon_b(q), alt_b(q), lat_l(i), lon_l(i), h_l(i),wgs84Ellipsoid,'degrees');
                % since we will need the normalized local nadir vector from each of these
                % bursts in the focal reference system later on, calculate an x,y,z vector for alt_b - 1 m
                [x_sel_low,y_sel_low,z_sel_low]=geodetic2enu(lat_b(q), lon_b(q), alt_b(q) - 1, lat_l(i), lon_l(i), h_l(i),wgs84Ellipsoid,'degrees');
                x_local_nadir = x_sel_low - x_select;
                y_local_nadir = y_sel_low - y_select;
                z_local_nadir = z_sel_low - z_select; % adds up to unit length
                % convert velocity vector into local reference frame
                [vx_b,vy_b,vz_b] = ecef2enuv(VX_select,VY_select,VZ_select, lat_l(i), lon_l(i),'degrees');

                % reference the time to the look loc time
                t_select = t_select - t_l(i);
                t_select_corr = t_select_corr - t_l(i);

                % (corrected) pulse timing
                t_pulse=t_select_corr'+tp_vec_select;
                t_pulse=t_pulse(:);

                % rearrange waveforms (for all pulses/bursts consecutively)
                wav_select=transpose(reshape(wav_select, [DDAcf.Np, DDAcf.Nb*length(t_select)]));

                % Apply window (not range, but in azimuth over single bursts, corresponding to L1b processing)!
                if DDAcf.ApplyWindow
                    wav_select=wav_select.*repmat(W',size(wav_select,1)/DDAcf.Nb,1);
                end

                % compute ranges between burst locations and focal point
                R0=sqrt(x_select.^2+y_select.^2+z_select.^2);
                
                % use definition of scalar product to get looking angle
                % from each burst to the focus point
                beam_angles = real(acos( (-x_select.*x_local_nadir - y_select.*y_local_nadir - z_select.*z_local_nadir)./R0)); % angle between nadir (0,0,1) and (x_select,y_select,z_select). Always positive, now define them (like in S3) so that negative is in direction of satellite velocity and positive is looking 'back'
                % use therefore x (or y), which have zero crossing at focal point to define the sign of angle
                if x_select(1) < 0
                    beam_angles = beam_angles .* sign(x_select);
                else
                    beam_angles = beam_angles .* sign(-x_select);
                end
                
                % radial velocity variation in slow time at every burst (from target positive)
                vr_0=(vx_b.*x_select+vy_b.*y_select+vz_b.*z_select)./R0;

                % interpolate radii and radial velocity to pulse locations (this also overcomes the datation issues)
                R0_p  = polyval(polyfit(t_select,R0,4),t_pulse);
                vr0_p  = polyval(polyfit(t_select,vr_0,4),t_pulse);
                if proc_sets.simulate_range_walk % if range walk shall be simulated, then assume only a single distance to focus point for all pulses
                    R0_p_rw = repmat(R0',64,1);
                    R0_p_rw = R0_p_rw(:);
                else
                    R0_p_rw = [];
                end
                % the window delay (tracker range) is kept constant over 64 pulses
                Rtrk = repmat(Rtrk_select', DDAcf.Nb,1);
                Rtrk=Rtrk(:);

                %% Fully Focused SAR algorithm   
                % RCMC, takes care of relative motion between target and satellite
                % (Doppler and Range) (see Egido IIIA)
                f_D = 2*DDAcf.fc*vr0_p / DDAcf.c;
                if ~proc_sets.simulate_range_walk
                    RCMC = exp(-2*pi*1i*(2*DDAcf.s*(R0_p)/DDAcf.c-f_D) * (t+Ts*DDAcf.nTs));
                else
                    RCMC = exp(-2*pi*1i*(2*DDAcf.s*(R0_p_rw)/DDAcf.c-f_D) * (t+Ts*DDAcf.nTs));
                end
                % for Sentinel-3 L1b waveforms the reference bin is 44
                % (1-based), but it seems that for the l1a files it is
                % actually 45 (1-based). This needs to be corrected for
                % with a tracker shift of one range resolution, namely
                % obj.DDAcf.os_ZP*obj.r_res_zp. This is consistent with the
                % cls implementation
                % https://github.com/cls-obsnadir-dev/SMAP-FFSAR which also
                % use refbin 44 (0-based), though 43 is stated in comments
                if ~is_s3
                    CTRP = exp(2*pi*1i*2*DDAcf.s*(Rtrk)/DDAcf.c*(t+Ts*DDAcf.nTs));
                else
                    CTRP = exp(2*pi*1i*2*DDAcf.s*(Rtrk - DDAcf.os_ZP*r_res_zp)/DDAcf.c*(t+Ts*DDAcf.nTs));
                end
                %% calculate effective radargram shift on burst basis for masking
                % This is to have only full bursts contribute to the pseudo Delay
                % Doppler implementation (mitigate resolutions of e.g. 600 m when half a burst would be flagged)
                f_D_b = 2*DDAcf.fc*vr_0/DDAcf.c;
                range_shift_b = R0 - Rtrk_select - f_D_b*DDAcf.c/(2*DDAcf.s);
                bin_shift = round(range_shift_b/(r_res_zp));
                bin_mask = DDAcf.Np*DDAcf.os_ZP - bin_shift;
                
                %% compute low quality interpolations of R0 and f_D t_select_corr for masking of aliased ambiguities similar to the L1b product
                if is_s6
                    PRF = 1/pri_l(i);
                    pR0 = polyfit(t_select_corr,R0 - f_D_b*DDAcf.c/(2*DDAcf.s),2);  % parabola
                    pfD = polyfit(t_select_corr,f_D_b,1); % linear

                    % find times, around which we alias:
                    folding_times = [-PRF/(2*pfD(1)) - pfD(2)/pfD(1); PRF/(2*pfD(1)) - pfD(2)/pfD(1)];
                    t_min = min(folding_times);
                    t_max = max(folding_times);
                    t_0 = - pfD(2)/pfD(1);
                    dt_dop = t_max - t_0;
                    
                    % compute folded range history
                    R0_alias1 = polyval(pR0,-(t_select_corr - t_0 + 2*dt_dop));
                    R0_alias2 = polyval(pR0,-(t_select_corr - t_0 - 2*dt_dop));
                    R0_alias = min(R0_alias1,R0_alias2);
                    
                    % now compute difference to rcmc shift around zero doppler bin
                    range_shift_b_alias = R0_alias - Rtrk_select - f_D_b*DDAcf.c/(2*DDAcf.s);
                    range_shift_b_sym = polyval(pR0,-(t_select_corr - t_0)) - Rtrk_select - f_D_b*DDAcf.c/(2*DDAcf.s);

                    % get burst closest to nadir
                    [~,nadir_ind] = min(abs(beam_angles));
                    
                    nadir_pos = DDAcf.os_ZP*46 +1; % empirically by comparison with L1b product, this means the leading edge normally starts after range bin 47, we do not know exactly which value ESA uses, but this is not too critical for range bin masking
                    
                    nadir_tracker_offset = Rtrk_select(nadir_ind) - Rtrk_l(i); % ESA likely chooses the reference range of the burst closest to nadir as reference for this leading edge bin threshold. Hence, we account for the offset introduced by our slighlty different tracker range
                    
                    %alias_mask = round((range_shift_b_alias - range_shift_b_sym)/(obj.r_res_zp)) + nadir_pos;
                    alias_mask = round(  (range_shift_b_alias - range_shift_b_sym + nadir_tracker_offset)/(r_res_zp) + nadir_pos);
                    alias_mask = min(alias_mask,DDAcf.Np*DDAcf.os_ZP);
                    alias_mask(alias_mask<0) = DDAcf.Np*DDAcf.os_ZP;
                    
                    bin_mask = min(bin_mask,alias_mask);
                    bin_mask = round(bin_mask);
                else
                    alias_mask = [];
                end
                
                wav_select = wav_select.*RCMC.*CTRP;

                % zero-padding (set zp to 2 if you want to zero-pad)
                wav_select=padarray(wav_select,[0 (DDAcf.os_ZP-1)*DDAcf.Np/2],0,'both'); 

                % Fast-time fft, range compression (see Egido IIIB)
                wav_select=fft(wav_select,[],2);
                if ~is_s6
                    wav_select=fftshift(wav_select,2);
                end
                
                %% mask range bins that periodically entered the radargram on the wrong side
                if is_s6
                    for n = 1:length(R0)
                        if bin_shift(n)>0
                            if bin_shift(n)<DDAcf.Np*DDAcf.os_ZP
                                wav_select(1+(n-1)*DDAcf.Nb:n*DDAcf.Nb, DDAcf.Np*DDAcf.os_ZP + 1 - bin_shift(n):end) = 0;
                            else
                                wav_select(1+(n-1)*DDAcf.Nb:n*DDAcf.Nb, 1:end) = 0;
                            end
                        elseif bin_shift(n)<0
                            if -bin_shift(n)<DDAcf.Np*DDAcf.os_ZP
                                wav_select(1+(n-1)*DDAcf.Nb:n*DDAcf.Nb, 1:-bin_shift(n)) = 0;
                            else
                                wav_select(1+(n-1)*DDAcf.Nb:n*DDAcf.Nb, 1:end) = 0;
                            end
                        end
                        
                        %mask range bins that may contain aliases
                        wav_select(1+(n-1)*DDAcf.Nb:n*DDAcf.Nb, alias_mask(n):end) = 0;
                    end
                else
                    for n = 1:length(R0)
                        if bin_shift(n)>0
                            if bin_shift(n)<DDAcf.Np*DDAcf.os_ZP
                                wav_select(1+(n-1)*DDAcf.Nb:n*DDAcf.Nb, DDAcf.Np*DDAcf.os_ZP + 1 - bin_shift(n):end) = 0;
                            else
                                wav_select(1+(n-1)*DDAcf.Nb:n*DDAcf.Nb, 1:end) = 0;
                            end
                        elseif bin_shift(n)<0
                            if -bin_shift(n)<DDAcf.Np*DDAcf.os_ZP
                                wav_select(1+(n-1)*DDAcf.Nb:n*DDAcf.Nb, 1:-bin_shift(n)) = 0;
                            else
                                wav_select(1+(n-1)*DDAcf.Nb:n*DDAcf.Nb, 1:end) = 0;
                            end
                        end
                    end
                end
                %%

                % Remove 'known' phase jumps after every burst and additionally after every LAI change 
                if ~is_s6
                    diff_LAI=[0; Rtrk_select(2:end)-Rtrk_select(1:end-1)]; % absolute removed w.r.t. earlier versions
                    phase_LAI=round(cumsum(diff_LAI)'/(1/80E6*DDAcf.c/2))*1.18*pi; %(1/80E6*c/2) is basically the tracker range jumps than can be made by the instrument

                    intra_burst_corr=1;
                    if is_cs
                        intra_burst_corr=0.5;           
                    end

                    phase_jumps=intra_burst_corr*pi*(0:length(t_select)-1)+phase_LAI-mean(phase_LAI);

                    phase_jumps=repmat(phase_jumps,DDAcf.Nb,1);
                    phasor_jumps=exp(-1i*phase_jumps(:));
                    wav_select=wav_select.*phasor_jumps;
                end
                %% debug phase plot
                if is_s6 && ~contains(FName, 'GPP')
                    %phase_slope = -567.8585; % randomly fitted value over the transponder, will depend on USO frequency and typical tracker jump magnitude
                    phase_slope = -567.8585*1.002;
                    phase_LAI = Rtrk*phase_slope;
                    phase_jumps=(phase_LAI-mean(phase_LAI));

                    phasor_jumps=exp(-1i*phase_jumps(:));
                    wav_select=wav_select.*phasor_jumps;
                end
                
                % compute distances 
                n_pulses = length(t_pulse);
                Ns_zp=DDAcf.os_ZP*DDAcf.Np;
                R_i=ones(n_pulses,1)*((1:Ns_zp)-DDAcf.RefBin*DDAcf.os_ZP)*r_res_zp+R0_i(i)*ones(n_pulses,Ns_zp);
                R_p=((R0_p*ones(1,Ns_zp)).^2+(R_i.^2-R0_i(i)^2)).^0.5;

                tau_i=2*(R_p-Rtrk*ones(1,Ns_zp))/DDAcf.c;

                % RVP, takes care of a higher order effects related to the chirp (see Egido IIIC)
                if ~is_s6
                    RVP=exp(2*pi*1i*DDAcf.s/2*(tau_i).^2);
                    wav_select=wav_select.*RVP;
                end

                % along-track summing and RRP (see Egido IIID)
                AF=exp(2*pi*1i*DDAcf.fc*(tau_i));
                wav_select=wav_select.*AF;

                f_f=j*f_b./proc_sets.along_track_oversampling;
                phasor=exp(-2*pi*1i*f_f*(t_pulse(:)'));
                SW=phasor*wav_select(:,:);
                
                % splitting the aperture to form 4 distinct FF-SAR waveforms with
                % only fourth the resolution:
                if proc_sets.split_aperture
                    split_inds = [0 round(size(t_pulse,1)*0.25) round(size(t_pulse,1)*0.5) round(size(t_pulse,1)*0.75) size(t_pulse,1)];
                    
                    SW_050T_025T = phasor(:,split_inds(1)+1:split_inds(2))*wav_select(split_inds(1)+1:split_inds(2),:);
                    SW_025T_000T = phasor(:,split_inds(2)+1:split_inds(3))*wav_select(split_inds(2)+1:split_inds(3),:);
                    SW_000T_025T = phasor(:,split_inds(3)+1:split_inds(4))*wav_select(split_inds(3)+1:split_inds(4),:);
                    SW_025T_050T = phasor(:,split_inds(4)+1:split_inds(5))*wav_select(split_inds(4)+1:split_inds(5),:);
                    
                     
%                     SW_T_0 = phasor(:,1:split_ind)*wav_select(1:split_ind,:);
%                     SW_0_T = phasor(:,split_ind+1:end)*wav_select(split_ind+1:end,:);
                    %dataIQ_T_0(ii,:,:) = transpose(SW_T_0);
                    %dataIQ_0_T(ii,:,:) = transpose(SW_0_T);
                    
                    obj.CS1b.SAR.dataIQ_050T_025T(:,i+j) = transpose(SW_050T_025T);
                    obj.CS1b.SAR.dataIQ_025T_000T(:,i+j) = transpose(SW_025T_000T);
                    obj.CS1b.SAR.dataIQ_000T_025T(:,i+j) = transpose(SW_000T_025T);
                    obj.CS1b.SAR.dataIQ_025T_050T(:,i+j) = transpose(SW_025T_050T);
                    
                    %obj.CS1b.SAR.dataIQ_T_0(:,i+j) = transpose(SW_T_0);
                    %obj.CS1b.SAR.dataIQ_0_T(:,i+j) = transpose(SW_0_T);
                end
                
                % write out unwrapped phase in transpondertest case
                if obj.proc_sets.transponder_test_mode && numel(look_loc_inds) == 1
                    [~,I_wav] = max(abs(wav_select(floor(length(wav_select(:,1))/2),:)));
                    uwphase(ii,:) = transpose(unwrap(angle(wav_select(:, I_wav))-mean(unwrap(angle(wav_select(:, I_wav))))));
                end
                
                %dataIQ(ii,:,:) = transpose(SW);
                %win_delay(ii,:) = Rtrk_l(i)/(DDAcf.c/2);
                obj.CS1b.SAR.dataIQ(:,i+j) = transpose(SW);
                obj.CS1b.MEA.win_delay(:,i+j) = Rtrk_l(i)/(DDAcf.c/2);
                
                [~,nadir_idx] = min(t_pulse.^2);
                %Doppler_correction(ii,:) = DDAcf.c*f_D(nadir_idx)./(2*DDAcf.s);
                %BeamAngle(ii,:,:) = vertcat(repmat(beam_angles,1,n_looks), NaN(n_max_beamangles-numel(beam_angles),n_looks));
                %stack_mask_start_stop(ii,:,:) = vertcat(repmat(bin_mask,1,n_looks), NaN(n_max_beamangles-numel(beam_angles),n_looks));
                obj.CS1b.COR.Doppler_correction(i+j) = DDAcf.c*f_D(nadir_idx)./(2*DDAcf.s);
                n_max_beamangles = 490;
                obj.CS1b.SAR.BeamAngle(:,i+j) = vertcat(repmat(beam_angles,1,n_looks), NaN(n_max_beamangles-numel(beam_angles),n_looks));
                obj.CS1b.SAR.stack_mask_start_stop(:,i+j) = vertcat(repmat(bin_mask,1,n_looks), NaN(n_max_beamangles-numel(beam_angles),n_looks));
                
                %% pseudo DD SAR waveform implementation:
                if proc_sets.output_pseudo_delay_doppler_processing
                    num_bursts_dd = size(wav_select,1)/DDAcf.Nb;
                    num_coherent_bursts = proc_sets.num_coherent_bursts; % determines how many bursts are coherently added (1: Delay Doppler processing, all: FF-SAR)
                    num_rbins = size(wav_select,2);
                    
                    if mod(num_bursts_dd,num_coherent_bursts) == 0
                        DDphasor = reshape(phasor,[n_looks DDAcf.Nb*num_coherent_bursts num_bursts_dd/num_coherent_bursts]);
                        %size(DDphasor)
                        %test
                        %DDphasor(:,:,1) == phasor(:,1:64)

                        DDwav = reshape(wav_select', [num_rbins DDAcf.Nb*num_coherent_bursts num_bursts_dd/num_coherent_bursts]);
                        DDwav = permute(DDwav,[2 1 3]);
                        %size(DDwav)
                    else
                        %choose only an integer amount of bursts that is
                        %dividable by num_coherent_bursts
                        num_bursts_dd = floor(num_bursts_dd/num_coherent_bursts)*num_coherent_bursts;
                        
                        DDphasor = reshape(phasor(:,1:num_bursts_dd*DDAcf.Nb),[n_looks DDAcf.Nb*num_coherent_bursts num_bursts_dd/num_coherent_bursts]);
                        %size(DDphasor)
                        %test
                        %DDphasor(:,:,1) == phasor(:,1:64)

                        DDwav = reshape(wav_select(1:num_bursts_dd*DDAcf.Nb,:)', [num_rbins DDAcf.Nb*num_coherent_bursts num_bursts_dd/num_coherent_bursts]);
                        DDwav = permute(DDwav,[2 1 3]);
                        %size(DDwav)
                    end
                    
                    % make efficient use of the multiprod function not to use
                    % for loop of matrix multiplications
                    %DDdata = multiprod(DDphasor,DDwav,[1 2],[1 2]); %
                    %overflows memory for zero-padding without any obvious
                    %reason! (array size should doubles but memory usage
                    %should not diverge), do for loop instead...
                    DDdata = zeros(size(DDphasor,1),size(DDwav,2),size(DDwav,3));
                    for k = 1:size(DDwav,3)
                        DDdata(:,:,k) = DDphasor(:,:,k)*DDwav(:,:,k);
                    end
                    %size(DDdata)

                    % take absolute mean square of zero doppler beam,
                    % determined by minimal time distance
                    [~,zD_idx] = min(t_select.^2);
                    start_idx_zD = max(zD_idx-2,1);
                    end_idx_zD = min(zD_idx+2,size(DDdata,3));
                    
                    obj.CS1b.SAR.data_pseudoDD_zeroDop(:,i+j) = flip((sum(abs(DDdata(:,:,start_idx_zD:end_idx_zD)).^2,3))',2);
                    
                    % take absolute square sum of doppler beam stack (all bursts)
                    DDdata = sum(abs(DDdata).^2, 3);
                    %size(DDdata)

                    %data_pseudoDD(ii,:,:) = flip(DDdata',2);
                    %num_bursts(ii,:,:) = num_bursts_dd;
                    obj.CS1b.SAR.data_pseudoDD(:,i+j) = flip(DDdata',2);
                    obj.CS1b.SAR.num_bursts(i+j) = num_bursts_dd;
                end
                
%                 %% saving pseudo DD SAR waveform (doppler beam stack) for all ground positions (works but very slow)
%                 num_bursts = size(wav_select,1)/DDAcf.Nb;
%                 num_rbins = size(wav_select,2);
%                 DDSAR = zeros(num_rbins, n_looks);
%                 % for each ground look
%                 for n = 1:n_looks
%                     % multiply the right phasor for the ground location
%                     DD_bs = wav_select(:,:) .* phasor(n,:)';
%                     % and accumulate the pulses over bursts
%                     DD_bs = movsum(DD_bs, DDAcf.Nb);
%                     % for Nb=64, the window for the first 1:64 indices is centered around 64/2+1 = 33
%                     % so take only each 64th sample
%                     DD_bs = abs(DD_bs(DDAcf.Nb/2+1:DDAcf.Nb:end,:)).^2;
%                     % and save the summed absolute square DD waveform
%                     DDSAR(:,n) = sum(DD_bs)';
%                 end
%                 
%                 obj.CS1b.SAR.data_pseudoDD(:,i+j) = flip(DDSAR,2);
                
%                 %% saving DDSAR waveform doppler beam stack for all
%                 ground positions (does not work for memory reasons)
%                 phasor3d(1,:,:) = phasor';
%                 % Now wav_select' has dimensions (range x pulses).
%                 % Phasor has dimensions (1 x pulses x ground points).
%                 wav3d = wav_select'.*phasor3d;
%                 wav3d = permute(wav3d, [2,1,3]);
%                 % wav3d has dimensions (pulses x range x ground points).
%                 % We want to bring this to (bursts x range x ground points)
%                 % by adding over all consequtive 64 pulses in dimension 2.
%                 % This can be realised either by matrix multiplication with something 
%                 % with dimensions (bursts x pulses) using multiproc, or
%                 % just by movsum:
%                 wav3d = movsum(wav3d, 64);
%                 % the window is centered around 64/2+1 = 33
%                 wav3d = wav3d(33:64:end,:,:);
%                 
%                 % now take the absolute
%                 wav3d = abs(wav3d).^2;
%                 
%                 % and add doppler beams together and squeeze the dimension
%                 % out
%                 wav3d = squeeze(sum(wav3d,1));
%                 
%                 % exceeds memory limits easily...
            end
%             delete(gcp); % comment in when calling parpool with parfor
            
            %% memory optimisation
%             obj.CS1b.SAR.dataIQ = reshape(permute(dataIQ,[2,3,1]),size(dataIQ,2),[]);
%             if proc_sets.split_aperture
%                 obj.CS1b.SAR.dataIQ_T_0 = reshape(permute(dataIQ_T_0,[2,3,1]),size(dataIQ_T_0,2),[]);
%                 obj.CS1b.SAR.dataIQ_0_T = reshape(permute(dataIQ_0_T,[2,3,1]),size(dataIQ_0_T,2),[]);
%             end 
%             if obj.proc_sets.output_pseudo_delay_doppler_processing
%                 obj.CS1b.SAR.data_pseudoDD = reshape(permute(data_pseudoDD,[2,3,1]),size(data_pseudoDD,2),[]);
%                 obj.CS1b.SAR.data_pseudoDD_zeroDop = reshape(permute(data_pseudoDD_zeroDop,[2,3,1]),size(data_pseudoDD_zeroDop,2),[]);
%             end
%             obj.CS1b.SAR.num_bursts = reshape(num_bursts, 1, []);
%             obj.CS1b.SAR.BeamAngle = reshape(permute(BeamAngle,[2,3,1]),size(BeamAngle,2),[]);
%             obj.CS1b.SAR.stack_mask_start_stop = reshape(permute(stack_mask_start_stop,[2,3,1]),size(stack_mask_start_stop,2),[]);
% 
%             obj.CS1b.MEA.win_delay = reshape(win_delay', 1, []);
%             obj.CS1b.COR.Doppler_correction = reshape(Doppler_correction', 1, []);
            %%
            obj.CS1b.SAR.data = abs(obj.CS1b.SAR.dataIQ).^2;
            obj.CS1b.MEA.tracker_range = obj.CS1b.MEA.win_delay*(obj.DDAcf.c/2);
            obj.CS1b.MEA.ref_range = obj.CS1b.MEA.tracker_range;  %duplicate of tracker_range, TODO: remove either of both
            
            
            
            % write out unwrapped phase in transpondertest case
            if obj.proc_sets.transponder_test_mode && numel(look_loc_inds) == 1
                obj.uwphase = transpose(uwphase);
            end
            
            if obj.is_cs
                obj.CS1b = CopyAncAuxDataFields(obj.CS1a,obj.CS1b);
            end

            % Wrap angle in degrees to [-180 180]
            obj.CS1a.GEO.LON = wrapTo180(obj.CS1a.GEO.LON);
            obj.CS1b.GEO.LON = wrapTo180(obj.CS1b.GEO.LON);

            % Normalize to range [0,65534]
            if obj.is_cs
                FAC                                 = repmat(range(obj.CS1b.SAR.data),size(obj.CS1b.SAR.data,1),1,1)./(65534*(1 - (repmat(min(obj.CS1b.SAR.data),size(obj.CS1b.SAR.data,1),1,1)./obj.CS1b.SAR.data)));
                [DUM,obj.CS1b.SAR.echo_scale_power] = log2(squeeze(min(FAC)));
                obj.CS1b.SAR.echo_scaling           = DUM*1E9;
                obj.CS1b.SAR.data                   = obj.CS1b.SAR.data.*(1./min(FAC));
            end           
        end
   end
end