function CS1b = CopyAncAuxDataFields(CS1a,CS1b)

%COPYANCAUXDATAFIELDS copies all ancillary and auxiliary data fields from
%CS1A to CS1b

%Input:
%CS1a: Latitudes, corrections, and surface type
%CS1b: Latitudes

%Output:
%CS1a.COR: geometry mask

%% Copy correction status and error flags
%Check whether any correction status flag has been set
CorrStatusFlag = [CS1a.COR.corr_status.dry_trop;CS1a.COR.corr_status.wet_trop;...
    CS1a.COR.corr_status.inv_bar;CS1a.COR.corr_status.dac;CS1a.COR.corr_status.gim_iono;...
    CS1a.COR.corr_status.model_iono;CS1a.COR.corr_status.ocean_equilibrium_tide;...
    CS1a.COR.corr_status.ocean_longperiod_tide;CS1a.COR.corr_status.ocean_loading_tide;...
    CS1a.COR.corr_status.solidearth_tide;CS1a.COR.corr_status.geocentric_polar_tide;...
    CS1a.COR.corr_status.surface_type];
if any(CorrStatusFlag(:)) ~= 1
    TMP = interp1(CS1a.GEO.LAT(1,:)',CorrStatusFlag',CS1b.GEO.LAT(1,:)','nearest','extrap');
    CS1b.COR.corr_status.dry_trop(:)               = TMP(:,1);  CS1b.COR.corr_status.wet_trop(:)              = TMP(:,2);
    CS1b.COR.corr_status.inv_bar(:)                = TMP(:,3);  CS1b.COR.corr_status.dac(:)                   = TMP(:,4);
    CS1b.COR.corr_status.gim_iono(:)               = TMP(:,5);  CS1b.COR.corr_status.model_iono(:)            = TMP(:,6);
    CS1b.COR.corr_status.ocean_equilibrium_tide(:) = TMP(:,7);  CS1b.COR.corr_status.ocean_longperiod_tide(:) = TMP(:,8);
    CS1b.COR.corr_status.ocean_loading_tide(:)     = TMP(:,9);  CS1b.COR.corr_status.solidearth_tide(:)       = TMP(:,10);
    CS1b.COR.corr_status.geocentric_polar_tide(:)  = TMP(:,11); CS1b.COR.corr_status.surface_type(:)          = TMP(:,12);
else
    CS1b.COR.corr_status.dry_trop(:) = 1; CS1b.COR.corr_status.wet_trop(:) = 1;
    CS1b.COR.corr_status.inv_bar(:) = 1; CS1b.COR.corr_status.dac(:) = 1;
    CS1b.COR.corr_status.gim_iono(:) = 1; CS1b.COR.corr_status.model_iono(:) = 1;
    CS1b.COR.corr_status.ocean_equilibrium_tide(:) = 1; CS1b.COR.corr_status.ocean_longperiod_tide(:) = 1;
    CS1b.COR.corr_status.ocean_loading_tide(:) = 1; CS1b.COR.corr_status.solidearth_tide(:) = 1;
    CS1b.COR.corr_status.geocentric_polar_tide(:) = 1; CS1b.COR.corr_status.surface_type(:) = 1;
end

%Check whether any correction error flag has been set
CorrErrorFlag = [CS1a.COR.corr_error.dry_trop;CS1a.COR.corr_error.wet_trop;...
    CS1a.COR.corr_error.inv_bar;CS1a.COR.corr_error.dac;CS1a.COR.corr_error.gim_iono;...
    CS1a.COR.corr_error.model_iono;CS1a.COR.corr_error.ocean_equilibrium_tide;...
    CS1a.COR.corr_error.ocean_longperiod_tide;CS1a.COR.corr_error.ocean_loading_tide;...
    CS1a.COR.corr_error.solidearth_tide;CS1a.COR.corr_error.geocentric_polar_tide;...
    CS1a.COR.corr_error.surface_type];
if any(CorrErrorFlag(:)) ~= 0
    TMP = interp1(CS1a.GEO.LAT(1,:)',CorrErrorFlag',CS1b.GEO.LAT(1,:)','nearest','extrap');
    CS1b.COR.corr_error.dry_trop(:)               = TMP(:,1);  CS1b.COR.corr_error.wet_trop(:)              = TMP(:,2);
    CS1b.COR.corr_error.inv_bar(:)                = TMP(:,3);  CS1b.COR.corr_error.dac(:)                   = TMP(:,4);
    CS1b.COR.corr_error.gim_iono(:)               = TMP(:,5);  CS1b.COR.corr_error.model_iono(:)            = TMP(:,6);
    CS1b.COR.corr_error.ocean_equilibrium_tide(:) = TMP(:,7);  CS1b.COR.corr_error.ocean_longperiod_tide(:) = TMP(:,8);
    CS1b.COR.corr_error.ocean_loading_tide(:)     = TMP(:,9);  CS1b.COR.corr_error.solidearth_tide(:)       = TMP(:,10);
    CS1b.COR.corr_error.geocentric_polar_tide(:)  = TMP(:,11); CS1b.COR.corr_error.surface_type(:)          = TMP(:,12);
else
    CS1b.COR.corr_error.dry_trop(:) = 0; CS1b.COR.corr_error.wet_trop(:) = 0;
    CS1b.COR.corr_error.inv_bar(:) = 0; CS1b.COR.corr_error.dac(:) = 0;
    CS1b.COR.corr_error.gim_iono(:) = 0; CS1b.COR.corr_error.model_iono(:) = 0;
    CS1b.COR.corr_error.ocean_equilibrium_tide(:) = 0; CS1b.COR.corr_error.ocean_longperiod_tide(:) = 0;
    CS1b.COR.corr_error.ocean_loading_tide(:) = 0; CS1b.COR.corr_error.solidearth_tide(:) = 0;
    CS1b.COR.corr_error.geocentric_polar_tide(:) = 0; CS1b.COR.corr_error.surface_type(:) = 0;
end

%% Copy geophysical corrections
%In CS1a, the correction values are provided for the first entry of each
%record. In CS1b, the values need to be computed for the first entry of
%each record as well

%Put all geophysical correction values provided in CS1a into one array
TMP = [CS1a.COR.dry_trop;CS1a.COR.wet_trop;CS1a.COR.inv_bar;CS1a.COR.dac;...
    CS1a.COR.gim_ion;CS1a.COR.model_ion;CS1a.COR.ocean_equilibrium_tide;...
    CS1a.COR.ocean_longperiod_tide;CS1a.COR.ocean_loading_tide;...
    CS1a.COR.solidearth_tide;CS1a.COR.geocentric_polar_tide;...
    CS1a.COR.TOTAL_gim;CS1a.COR.TOTAL_model];

%Conduct interpolation
TMP = interp1(CS1a.GEO.LAT(1,:)',TMP',CS1b.GEO.LAT(1,:)','linear');

%Copy to CS1b
CS1b.COR.dry_trop(:)               = TMP(:,1);
CS1b.COR.wet_trop(:)               = TMP(:,2);
CS1b.COR.inv_bar(:)                = TMP(:,3);
CS1b.COR.dac(:)                    = TMP(:,4);
CS1b.COR.gim_ion(:)                = TMP(:,5);
CS1b.COR.model_ion(:)              = TMP(:,6);
CS1b.COR.ocean_equilibrium_tide(:) = TMP(:,7);
CS1b.COR.ocean_longperiod_tide(:)  = TMP(:,8);
CS1b.COR.ocean_loading_tide(:)     = TMP(:,9);
CS1b.COR.solidearth_tide(:)        = TMP(:,10);
CS1b.COR.geocentric_polar_tide(:)  = TMP(:,11);
CS1b.COR.TOTAL_gim(:)              = TMP(:,12);
CS1b.COR.TOTAL_model(:)            = TMP(:,13);

%% Interpolate surface type to CS1b data locations
CS1b.COR.surf_type(:) = interp1(CS1a.GEO.LAT(1,:)',CS1a.COR.surf_type',CS1b.GEO.LAT(1,:)','nearest');

%% Copy fields from measurements group
CS1b.MEA.Tx_Pwr(:)    = interp1(CS1a.GEO.LAT(CS1a.MEA.Tx_Pwr(:) ~= 0),CS1a.MEA.Tx_Pwr(CS1a.MEA.Tx_Pwr(:) ~= 0),CS1b.GEO.LAT(:),'nearest');

end