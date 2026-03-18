function [EUML2,COR_Interpolator] = S6_read_L2(FName,DOM)
%     FName = '/home/fehelers/PhD Delft/Projects/hydrology river/Data/L2/S6A_P4_2__HR______20210108T224246_20210108T233859_20220508T091227_3373_006_070_035_EUM__REP_NT_F06.SEN6/S6A_P4_2__HR______20210108T224246_20210108T233859_20220508T091227_3373_006_070_035_EUM__REP_NT_F06.SEN6/S6A_P4_2__HR_STD__NT_006_070_20210108T224246_20210108T233859_F06.nc'
%     DOM = [44.25 44.66; -30 30]
    % ncinfo(FName,'data_01').Variables.Name
    % dummy = ncread(FName,"data_01/ku/iono_cor_gim")

    vars_01 =  ["data_01/ku/iono_cor_gim",...
                "data_01/ku/model_instr_cor_range_ocean",...
                "data_01/latitude",...
                "data_01/longitude",...
                "data_01/altitude",...
                "data_01/pole_tide",...
                "data_01/solid_earth_tide",...
                "data_01/time"
                ];

    vars_20 =  ["data_20/ku/load_tide_sol1",...
                "data_20/ku/load_tide_sol2",...
                "data_20/ku/latitude",...
                "data_20/ku/longitude",...
                "data_20/ku/altitude",...
                "data_20/ku/model_dry_tropo_cor_measurement_altitude",...
                "data_20/ku/model_wet_tropo_cor_measurement_altitude",...
                "data_20/ku/time"
                ];

    for i = 1:numel(vars_01)
        [~,v] = fileparts(vars_01(i));
        EUML2.('data_01').(v) = ncread(FName,vars_01(i));
    end

    for i = 1:numel(vars_20)
        [~,v] = fileparts(vars_20(i));
        EUML2.('data_20').(v) = ncread(FName,vars_20(i));
    end

    vars_interp_01 = ["iono_cor_gim",...
                "model_instr_cor_range_ocean",...
                "pole_tide",...
                "solid_earth_tide"
                ];

    vars_interp_20 = ["load_tide_sol1",...
                "load_tide_sol2",...
                "model_dry_tropo_cor_measurement_altitude",...
                "model_wet_tropo_cor_measurement_altitude"
                ];


    for i = 1:numel(vars_interp_01)
        var = vars_interp_01(i);
        eum = EUML2.data_01;

        % necessary because griddedInterpolant only accepts ascending coordinate vectors
        [~,sort_idx01] = sort(eum.latitude);
        lat = eum.latitude(sort_idx01);
        val = eum.(var)(sort_idx01);
        idx_lat = (DOM(1)-0.2<lat)&(lat < DOM(3)+0.2);

        COR_Interpolator.(var)  = griddedInterpolant(lat(idx_lat),val(idx_lat),'linear');
    end

    for i = 1:numel(vars_interp_20)
        var = vars_interp_20(i);
        eum = EUML2.data_20;

        % necessary because griddedInterpolant only accepts ascending coordinate vectors
        [~,sort_idx20] = sort(eum.latitude);
        lat = eum.latitude(sort_idx20);
        val = eum.(var)(sort_idx20);
        idx_lat = (DOM(1)-0.2<lat)&(lat < DOM(3)+0.2);

        COR_Interpolator.(var)  = griddedInterpolant(lat(idx_lat),val(idx_lat),'linear');
    end

    % %% double check the corrections
    % figure;
    % 
    % plot(EUML2.data_01.latitude,EUML2.data_01.iono_cor_gim);hold on
    % plot(EUML2.data_20.latitude,COR_Interpolator.iono_cor_gim(EUML2.data_20.latitude));hold on
    % 
    % figure; 
    % 
    % plot(EUML2.data_20.latitude,EUML2.data_20.model_dry_tropo_cor_measurement_altitude);hold on
    % plot(EUML2.data_20.latitude,COR_Interpolator.model_dry_tropo_cor_measurement_altitude(EUML2.data_20.latitude));hold on

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
