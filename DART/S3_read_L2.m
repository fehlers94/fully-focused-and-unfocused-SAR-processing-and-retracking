 function [EUML2,COR_Interpolator] = S3_read_L2(FName,DOM)
    EUML2 = [];
    
%     % some code for searching netcdf variables for keywords
%     finfo = ncinfo(FName)
%     %% discovering variables
%     for i = 1:length(finfo.Variables)
%         if contains(finfo.Variables(i).Name,'ssh')
%             finfo.Variables(i).Name
%         end
%     end
    
%     %in case of reading only parts (not recommended, influence of 1 Hz variable may cause severe data loss at the borders of up to ~7km)      
%     LAT20                    = ncread(FName,'lat_20_ku');
%     LON20                    = ncread(FName,'lon_20_ku');
%     IDX_mask20 = ingeoquad(LAT,LON,DOM(1,:),DOM(2,:));
%     IDX20 = find(IDX_mask20);
%     startLoc20=min(IDX20);
%     count20=max(IDX20)-startLoc20+1;

    vars = ["time_20_ku" , "UTC_sec_20_ku","time_01",...
            "lat_20_ku","lat_01",...
            "lon_20_ku","lon_01",...
            "alt_20_ku","alt_01",...
            "iono_cor_alt_20_ku",... % based on C-band pulses - ocean only
            "iono_cor_gim_01_ku",... % model-based, preferred
            "mod_dry_tropo_cor_zero_altitude_01",... 
            "mod_dry_tropo_cor_meas_altitude_01",... % preferred, 10 m altitude refer to 1 mm error
            "dop_cor_20_ku",...
            "dop_cor_l1b_20_ku",...
            "mod_wet_tropo_cor_zero_altitude_01",...
            "mod_wet_tropo_cor_meas_altitude_01",... % preferred
            "comp_wet_tropo_cor_01_ku",... % composite of model and radiometer data, details?
            "rad_wet_tropo_cor_01_ku",... % from radiometer, measured, may be bad at coasts
            "solid_earth_tide_01",...
            "ocean_tide_sol1_01",... %GOT4.10c
            "ocean_tide_sol2_01",... %FES2014b (default?)
            "load_tide_sol1_01",... %GOT4.10c
            "load_tide_sol2_01",... %FES2014b (default?)
            "sea_state_bias_01_ku",... % bad for several reasons in coastal zone, own correction required?
            "pole_tide_01",...
            "inv_bar_cor_01",...
            "hf_fluct_cor_01",...
            "cog_cor_01",...
            "mod_instr_cor_range_01_ku",...
            "range_ocean_20_ku"
            ];
    
    for i = 1:numel(vars)
        EUML2.(vars(i)) = ncread(FName,vars(i));
    end
    
    % now make a struct of gridded interpolators from the data, depending
    % on latitude
    
    vars_interp = ["iono_cor_alt_20_ku",... % based on C-band pulses - ocean only
            "iono_cor_gim_01_ku",... % model-based, preferred
            "mod_dry_tropo_cor_zero_altitude_01",... 
            "mod_dry_tropo_cor_meas_altitude_01",... % preferred, 10 m altitude refer to 1 mm error
            "dop_cor_20_ku",...
            "dop_cor_l1b_20_ku",...
            "mod_wet_tropo_cor_zero_altitude_01",...
            "mod_wet_tropo_cor_meas_altitude_01",... % preferred
            "comp_wet_tropo_cor_01_ku",... % composite of model and radiometer data, details?
            "rad_wet_tropo_cor_01_ku",... % from radiometer, measured, may be bad at coasts
            "solid_earth_tide_01",...
            "ocean_tide_sol1_01",... %GOT4.10c
            "ocean_tide_sol2_01",... %FES2014b (default?)
            "load_tide_sol1_01",... %GOT4.10c
            "load_tide_sol2_01",... %FES2014b (default?)
            "sea_state_bias_01_ku",... % bad for several reasons in coastal zone, own correction required?
            "pole_tide_01",...
            "inv_bar_cor_01",...
            "hf_fluct_cor_01",...
            "cog_cor_01",...
            "mod_instr_cor_range_01_ku"];
    
    for i = 1:numel(vars_interp)
        var = vars_interp(i);
        
        [~,sort_idx20] = sort(EUML2.lat_20_ku); % necessary because griddedInterpolant only accepts ascending coordinate vectors
        [~,sort_idx01] = sort(EUML2.lat_01);
        if contains(var,'01')
            lat = EUML2.lat_01(sort_idx01);
            val = EUML2.(var)(sort_idx01);
            idx_lat = (DOM(1)-0.2<lat)&(lat < DOM(3)+0.2);
            
            COR_Interpolator.(var)  = griddedInterpolant(lat(idx_lat),val(idx_lat),'linear');
        elseif contains(var,'20_ku')
            lat = EUML2.lat_20_ku(sort_idx20);
            val = EUML2.(var)(sort_idx20);
            idx_lat = (DOM(1)-0.2<lat)&(lat < DOM(3)+0.2); %plus minus 30 km in order to circumvent that data from 01 Hz product becomes NaN 
            
            COR_Interpolator.(var)  = griddedInterpolant(lat(idx_lat),val(idx_lat),'linear');
        else
            warning('Variable does not contain 01 or 20_ku: interpolation not possible')
        end
    end
 end
 
% IMPORTANT NOTES ON THE NECESSARY CORRECTIONS
%
% %% variables to get from L2 file
% 'iono_cor_alt_20_ku' % based on C-band pulses - ocean only
% 'iono_cor_gim_01_ku' % model-based, preferred
% 
% 'mod_dry_tropo_cor_zero_altitude_01'
% 'mod_dry_tropo_cor_meas_altitude_01' % preferred, 10 m altitude refer to 1 mm error
% 
% 'dop_cor_20_ku'
% 'dop_cor_l1b_20_ku'
% 
% 'mod_wet_tropo_cor_zero_altitude_01'
% 'mod_wet_tropo_cor_meas_altitude_01' % preferred
% 'comp_wet_tropo_cor_01_ku' % composite of model and radiometer data, details?
% 'rad_wet_tropo_cor_01_ku' % from radiometer, measured, may be bad at coasts
% 
% 'solid_earth_tide_01'
% 
% 'ocean_tide_sol1_01' %GOT4.10c
% 'ocean_tide_sol2_01' %FES2014b (default?)
% 
% 'load_tide_sol1_01' %GOT4.10c
% 'load_tide_sol2_01' %FES2014b (default?)
% 
% 'sea_state_bias_01_ku' % bad for several reasons in coastal zone, own correction required?
% 
% 'pole_tide_01'
% 
% 'inv_bar_cor_01'
% 
% 'hf_fluct_cor_01'
% 
% 'cog_cor_01'
% 
% 'mod_instr_cor_range_01_ku'
% 
% 
% 'ssha_20_ku'
%     %is calculated via:
%     altitude of satellite (alt_20_ku) 
%     - Ku band corrected ocean altimeter range (range_ocean_20_ku)
%         is hopefully calculated via:
%         tracker range 
%         [+ (USO drift correction (uso_cor_20_ku))
%         + (internal path correction (int_path_cor_20_ku))] % contained in L1a tracker range
%         + epoch 
%         + distance antenna-COG (cog_cor_01), 
%         + Doppler correction (dop_cor_20_ku), % implicitly included in L1b processing?
%         + modeled instrumental errors correction (mod_instr_cor_range_01_ku)
%         (+ system bias)
%     - altimeter ionospheric correction on Ku band (iono_cor_alt_20_ku) 
%     - model dry tropospheric correction (mod_dry_tropo_cor_zero_altitude_01) 
%     - radiometer wet tropospheric correction (rad_wet_tropo_cor_01_ku) 
%     - sea state bias correction in Ku band (sea_state_bias_01_ku) 
%     - solid earth tide height (solid_earth_tide_01) 
%     - geocentric ocean tide height solution 2 = FES (ocean_tide_sol2_01) 
%     - geocentric pole tide height (pole_tide_01) 
%     - inverted barometer height correction (inv_bar_cor_01) 
%     - high frequency fluctuations of the sea surface topography (hf_fluct_cor_01 for NTC/STC off line products only) 
%     - mean sea surface (mean_sea_surf_sol2_20_ku)
%     
%     'range_ocean_20_ku'
%         %contains: 
%         USO drift correction (uso_cor_20_ku), % included in L1a already
%         internal path correction (int_path_cor_20_ku), % included in L1a already
%         distance antenna-COG (cog_cor_01), 
%         Doppler correction (dop_cor_20_ku), % implicitly included in L1b processing?
%         modeled instrumental errors correction (mod_instr_cor_range_01_ku)
%         and system bias
% 
%     Level-1 tracker range is Distance between the altimeter reference point and the surface height associated to a range gate used as reference inside the tracking window (reference tracking point), corrected for USO frequency drift and internal path correction
% 
% This makes that the range should be given by:
% 
% retracked range R 
% + distance antenna-COG (cog_cor_01), 
% %+ Doppler correction (dop_cor_20_ku), %likely implicitly included in L1b processing
% + modeled instrumental errors correction (mod_instr_cor_range_01_ku)
% %+-? system bias
% 
% and then with h - R - all other corrections as in Dinardo 2018 it should be fully consistent, also with what Cornelis did.
% 
% 



% %%
% 
% %%
% 
% CS.COR.TIME1Hz                   = datenum('2000','yyyy') + double(ncread(fullfile(FNameL2.folder,FNameL2.name),'UTC_day_01')) + ncread(fullfile(FNameL2.folder,FNameL2.name),'UTC_sec_01')/86400;
% CS.COR.TIME20Hz                  = datenum('2000','yyyy') + double(ncread(fullfile(FNameL2.folder,FNameL2.name),'UTC_day_20_ku')) + ncread(fullfile(FNameL2.folder,FNameL2.name),'UTC_sec_20_ku')/86400;
% CS.COR.LAT1Hz                    = ncread(fullfile(FNameL2.folder,FNameL2.name),'lat_01');
% CS.COR.LON1Hz                    = ncread(fullfile(FNameL2.folder,FNameL2.name),'lon_01');
% CS.COR.LAT20Hz                   = ncread(fullfile(FNameL2.folder,FNameL2.name),'lat_20_ku');
% CS.COR.LON20Hz                   = ncread(fullfile(FNameL2.folder,FNameL2.name),'lon_20_ku');
% CS.COR.iono_cor_alt_20_ku        = ncread(fullfile(FNameL2.folder,FNameL2.name),'iono_cor_alt_20_ku');
% CS.COR.iono_cor_gim_01_ku        = ncread(fullfile(FNameL2.folder,FNameL2.name),'iono_cor_gim_01_ku');
% CS.COR.dry_trop                  = ncread(fullfile(FNameL2.folder,FNameL2.name),'mod_dry_tropo_cor_zero_altitude_01');
% CS.COR.wet_trop                  = ncread(fullfile(FNameL2.folder,FNameL2.name),'rad_wet_tropo_cor_01_ku');
% CS.COR.SSB                       = ncread(fullfile(FNameL2.folder,FNameL2.name),'sea_state_bias_01_ku');
% CS.COR.solidearth_tide           = ncread(fullfile(FNameL2.folder,FNameL2.name),'solid_earth_tide_01');
% CS.COR.ocean_equilibrium_tide    = ncread(fullfile(FNameL2.folder,FNameL2.name),'ocean_tide_sol1_01');
% CS.COR.ocean_longperiod_tide     = zeros(size(CS.COR.ocean_equilibrium_tide)); %This value has already been added to the two geocentric ocean tide height values recorded in the product (ocean_tide_sol1_01 and ocean_tide_sol2_01).
% CS.COR.ocean_loading_tide        = zeros(size(CS.COR.ocean_equilibrium_tide)); %This value has already been added to the corresponding ocean tide height value recorded in the product (ocean_tide_sol2_01)."
% CS.COR.geocentric_polar_tide     = ncread(fullfile(FNameL2.folder,FNameL2.name),'pole_tide_01');
% CS.COR.inv_bar                   = ncread(fullfile(FNameL2.folder,FNameL2.name),'inv_bar_cor_01');
% CS.COR.hf_fluct_cor              = ncread(fullfile(FNameL2.folder,FNameL2.name),'hf_fluct_cor_01');
% CS.COR.cog_cor_01                = ncread(fullfile(FNameL2.folder,FNameL2.name),'cog_cor_01');
% CS.COR.mod_instr_cor_range_01_ku = ncread(fullfile(FNameL2.folder,FNameL2.name),'mod_instr_cor_range_01_ku');
% 
% % additions
% CS.COR.wet_trop                  = ncread(fullfile(FNameL2.folder,FNameL2.name),'comp_wet_tropo_cor_01_ku');
% 
% 
% % CorrST01           = 
% % CS.COR.dry_trop
% % CS.COR.wet_trop
% % %CS.COR.ocean_equilibrium_tide
% % %CS.COR.ocean_longperiod_tide
% % CS.COR.ocean_loading_tide
% % CS.COR.solidearth_tide
% % CS.COR.geocentric_polar_tide
% % %CS.COR.inv_bar
% % %CS.COR.hf_fluct_cor
% % CS.COR.iono_cor_gim_01_ku
% % CS.COR.cog_cor_01
% % CS.COR.mod_instr_cor_range_01_ku;
% % 
% % CorrST01i          = 
% % CS.COR.dry_trop
% % CS.COR.wet_trop
% % CS.COR.ocean_loading_tide                   % is in fact read as zeroes?
% % CS.COR.iono_cor_gim_01_ku
% % CS.COR.solidearth_tide
% % 0.468*CS.COR.geocentric_polar_tide
% % CS.COR.cog_cor_01
% % CS.COR.mod_instr_cor_range_01_ku;
% % 
% % %'how to translate them:'
% % CS.COR.iono_cor_alt_20_ku        = ncread(fullfile(FNameL2.folder,FNameL2.name),'iono_cor_alt_20_ku');
% % CS.COR.iono_cor_gim_01_ku        = ncread(fullfile(FNameL2.folder,FNameL2.name),'iono_cor_gim_01_ku');
% % CS.COR.dry_trop                  = ncread(fullfile(FNameL2.folder,FNameL2.name),'mod_dry_tropo_cor_zero_altitude_01');
% % CS.COR.wet_trop                  = ncread(fullfile(FNameL2.folder,FNameL2.name),'rad_wet_tropo_cor_01_ku');
% % CS.COR.SSB                       = ncread(fullfile(FNameL2.folder,FNameL2.name),'sea_state_bias_01_ku');
% % CS.COR.solidearth_tide           = ncread(fullfile(FNameL2.folder,FNameL2.name),'solid_earth_tide_01');
% % CS.COR.ocean_equilibrium_tide    = ncread(fullfile(FNameL2.folder,FNameL2.name),'ocean_tide_sol1_01');
% % CS.COR.ocean_longperiod_tide     = zeros(size(CS.COR.ocean_equilibrium_tide)); %This value has already been added to the two geocentric ocean tide height values recorded in the product (ocean_tide_sol1_01 and ocean_tide_sol2_01).
% % CS.COR.ocean_loading_tide        = zeros(size(CS.COR.ocean_equilibrium_tide)); %This value has already been added to the corresponding ocean tide height value recorded in the product (ocean_tide_sol2_01)."
% % CS.COR.geocentric_polar_tide     = ncread(fullfile(FNameL2.folder,FNameL2.name),'pole_tide_01');
% % CS.COR.inv_bar                   = ncread(fullfile(FNameL2.folder,FNameL2.name),'inv_bar_cor_01');
% % CS.COR.hf_fluct_cor              = ncread(fullfile(FNameL2.folder,FNameL2.name),'hf_fluct_cor_01');
% % CS.COR.cog_cor_01                = ncread(fullfile(FNameL2.folder,FNameL2.name),'cog_cor_01');
% % CS.COR.mod_instr_cor_range_01_ku = ncread(fullfile(FNameL2.folder,FNameL2.name),'mod_instr_cor_range_01_ku');
