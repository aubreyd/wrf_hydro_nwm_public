md5sum: e9505c9d024b0cebaa7019b609c60351  /glade/work/jamesmcc/domains/public/croton_NY/NWM/DOMAIN_LR/Fulldom_hires.nc
ncdump -h: netcdf Fulldom_hires {
dimensions:
	y = 16 ;
	x = 15 ;
variables:
	short CHANNELGRID(y, x) ;
		CHANNELGRID:grid_mapping = "crs" ;
		CHANNELGRID:coordinates = "x y" ;
		CHANNELGRID:long_name = "CHANNELGRID" ;
		CHANNELGRID:units = "Meter" ;
		CHANNELGRID:missing_value = -32768s ;
		CHANNELGRID:esri_pe_string = "PROJCS[\"Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic_2SP\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.0],UNIT[\"Meter\",1.0]];-35691800 -29075200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision" ;
	byte FLOWDIRECTION(y, x) ;
		FLOWDIRECTION:grid_mapping = "crs" ;
		FLOWDIRECTION:coordinates = "x y" ;
		FLOWDIRECTION:long_name = "FLOWDIRECTION" ;
		FLOWDIRECTION:units = "Meter" ;
		FLOWDIRECTION:missing_value = 0b ;
		FLOWDIRECTION:esri_pe_string = "PROJCS[\"Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic_2SP\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.0],UNIT[\"Meter\",1.0]];-35691800 -29075200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision" ;
	short LAKEGRID(y, x) ;
		LAKEGRID:grid_mapping = "crs" ;
		LAKEGRID:coordinates = "x y" ;
		LAKEGRID:long_name = "LAKEGRID" ;
		LAKEGRID:units = "Meter" ;
		LAKEGRID:missing_value = -32768s ;
		LAKEGRID:esri_pe_string = "PROJCS[\"Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic_2SP\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.0],UNIT[\"Meter\",1.0]];-35691800 -29075200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision" ;
	float LATITUDE(y, x) ;
		LATITUDE:grid_mapping = "crs" ;
		LATITUDE:_CoordinateSystems = "crs" ;
		LATITUDE:long_name = "latitude coordinate" ;
		LATITUDE:units = "degrees_north" ;
		LATITUDE:_CoordinateAxisType = "Lat" ;
		LATITUDE:esri_pe_string = "PROJCS[\"Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic_2SP\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.0],UNIT[\"Meter\",1.0]];-35691800 -29075200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision" ;
	double LKSATFAC(y, x) ;
		LKSATFAC:grid_mapping = "crs" ;
		LKSATFAC:coordinates = "x y" ;
		LKSATFAC:long_name = "OVROUGHRTFAC" ;
		LKSATFAC:units = "Meter" ;
		LKSATFAC:missing_value = -3.402823e+38f ;
		LKSATFAC:esri_pe_string = "PROJCS[\"Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic_2SP\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.0],UNIT[\"Meter\",1.0]];-35691800 -29075200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision" ;
	float LONGITUDE(y, x) ;
		LONGITUDE:grid_mapping = "crs" ;
		LONGITUDE:_CoordinateSystems = "crs" ;
		LONGITUDE:long_name = "longitude coordinate" ;
		LONGITUDE:units = "degrees_east" ;
		LONGITUDE:_CoordinateAxisType = "Lon" ;
		LONGITUDE:esri_pe_string = "PROJCS[\"Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic_2SP\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.0],UNIT[\"Meter\",1.0]];-35691800 -29075200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision" ;
	float OVROUGHRTFAC(y, x) ;
		OVROUGHRTFAC:grid_mapping = "crs" ;
		OVROUGHRTFAC:coordinates = "x y" ;
		OVROUGHRTFAC:long_name = "OVROUGHRTFAC" ;
		OVROUGHRTFAC:units = "Meter" ;
		OVROUGHRTFAC:missing_value = -3.402823e+38f ;
		OVROUGHRTFAC:esri_pe_string = "PROJCS[\"Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic_2SP\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.0],UNIT[\"Meter\",1.0]];-35691800 -29075200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision" ;
	float RETDEPRTFAC(y, x) ;
		RETDEPRTFAC:grid_mapping = "crs" ;
		RETDEPRTFAC:coordinates = "x y" ;
		RETDEPRTFAC:long_name = "RETDEPRTFAC" ;
		RETDEPRTFAC:units = "Meter" ;
		RETDEPRTFAC:missing_value = -3.402823e+38f ;
		RETDEPRTFAC:esri_pe_string = "PROJCS[\"Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic_2SP\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.0],UNIT[\"Meter\",1.0]];-35691800 -29075200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision" ;
	short STREAMORDER(y, x) ;
		STREAMORDER:grid_mapping = "crs" ;
		STREAMORDER:coordinates = "x y" ;
		STREAMORDER:long_name = "STREAMORDER" ;
		STREAMORDER:units = "Meter" ;
		STREAMORDER:missing_value = -32768s ;
		STREAMORDER:esri_pe_string = "PROJCS[\"Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic_2SP\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.0],UNIT[\"Meter\",1.0]];-35691800 -29075200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision" ;
	float TOPOGRAPHY(y, x) ;
		TOPOGRAPHY:long_name = "GDAL Band Number 1" ;
		TOPOGRAPHY:grid_mapping = "crs" ;
		TOPOGRAPHY:esri_pe_string = "PROJCS[\"Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic_2SP\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.0],UNIT[\"Meter\",1.0]];-35691800 -29075200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision" ;
		TOPOGRAPHY:_FillValue = -9999.f ;
	short basn_msk(y, x) ;
		basn_msk:grid_mapping = "crs" ;
		basn_msk:coordinates = "x y" ;
		basn_msk:long_name = "basn_msk" ;
		basn_msk:units = "Meter" ;
		basn_msk:missing_value = -32768s ;
		basn_msk:esri_pe_string = "PROJCS[\"Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic_2SP\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.0],UNIT[\"Meter\",1.0]];-35691800 -29075200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision" ;
	char crs ;
		crs:transform_name = "lambert_conformal_conic" ;
		crs:grid_mapping_name = "lambert_conformal_conic" ;
		crs:esri_pe_string = "PROJCS[\"Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic_2SP\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.0],UNIT[\"Meter\",1.0]];-35691800 -29075200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision" ;
		crs:spatial_ref = "PROJCS[\"Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic_2SP\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.0],UNIT[\"Meter\",1.0]];-35691800 -29075200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision" ;
		crs:long_name = "CRS definition" ;
		crs:longitude_of_prime_meridian = 0. ;
		crs:GeoTransform = "-2303999.17655 1000.0 0 1919999.66329 0 -1000.0 " ;
		crs:_CoordinateAxes = "y x" ;
		crs:_CoordinateTransformType = "Projection" ;
		crs:standard_parallel = 30., 60. ;
		crs:longitude_of_central_meridian = -97. ;
		crs:latitude_of_projection_origin = 40. ;
		crs:false_easting = 0. ;
		crs:false_northing = 0. ;
		crs:earth_radius = 6370000. ;
		crs:semi_major_axis = 6370000. ;
		crs:inverse_flattening = 0. ;
	short frxst_pts(y, x) ;
		frxst_pts:grid_mapping = "crs" ;
		frxst_pts:coordinates = "x y" ;
		frxst_pts:long_name = "frxst_pts" ;
		frxst_pts:units = "Meter" ;
		frxst_pts:missing_value = -32768s ;
		frxst_pts:esri_pe_string = "PROJCS[\"Lambert_Conformal_Conic\",GEOGCS[\"GCS_Sphere\",DATUM[\"D_Sphere\",SPHEROID[\"Sphere\",6370000.0,0.0]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Lambert_Conformal_Conic_2SP\"],PARAMETER[\"false_easting\",0.0],PARAMETER[\"false_northing\",0.0],PARAMETER[\"central_meridian\",-97.0],PARAMETER[\"standard_parallel_1\",30.0],PARAMETER[\"standard_parallel_2\",60.0],PARAMETER[\"latitude_of_origin\",40.0],UNIT[\"Meter\",1.0]];-35691800 -29075200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision" ;
	double x(x) ;
		x:standard_name = "projection_x_coordinate" ;
		x:long_name = "x coordinate of projection" ;
		x:units = "m" ;
		x:_CoordinateAxisType = "GeoX" ;
		x:resolution = 1000. ;
	double y(y) ;
		y:standard_name = "projection_y_coordinate" ;
		y:long_name = "y coordinate of projection" ;
		y:units = "m" ;
		y:_CoordinateAxisType = "GeoY" ;
		y:resolution = 1000. ;

// global attributes:
		:Conventions = "CF-1.5" ;
		:Source_Software = "WRF-Hydro GIS Pre-processor v5 (03/2018)" ;
		:history = "Thu Sep  6 17:56:43 2018: ncks -O -d x,4146,4160 -d y,1483,1498 /glade/p_old/nwc/nwmv20_finals/CONUS/DOMAIN/Fulldom_hires_netcdf_file_1km_LongRange_NWMv2.0.nc /glade/scratch/adugger/TestCases/NY_Croton/DOMAIN_NWMv2.0_LongRange//0137462010/Fulldom_hires.nc\n",
			"Tue Apr 10 15:20:44 2018: ncap2 -O -s LKSATFAC=OVROUGHRTFAC*0.0+1000.0 Fulldom_hires_netcdf_file_1km_NWMv1.1_DEFAULT.nc Fulldom_hires_netcdf_file_1km_NWMv1.1_DEFAULT.nc\n",
			"Fri Nov 11 15:48:10 2016: ncrename -v Band1,TOPOGRAPHY Fulldom_hires_netcdf_file_1km_NWMv1.1_newdem.nc Fulldom_hires_netcdf_file_1km_NWMv1.1_newdem.nc\n",
			"Fri Nov 11 15:47:37 2016: ncks -x -v TOPOGRAPHY Fulldom_hires_netcdf_file_1km_NWMv1.1_newdem.nc Fulldom_hires_netcdf_file_1km_NWMv1.1_newdem.nc\n",
			"Fri Nov 11 15:46:53 2016: ncks -A -v Band1 DEM_working/DEM1km_rclass_Fillnodata_revy.nc Fulldom_hires_netcdf_file_1km_NWMv1.1_newdem.nc\n",
			"Fri Sep 25 11:50:01 2015: ncks -A longitude.nc Fulldom_hires_netcdf_file.nc\n",
			"Fri Sep 25 11:49:59 2015: ncks -A latitude.nc Fulldom_hires_netcdf_file.nc\n",
			"Fri Sep 25 11:49:58 2015: ncks -A frxst_pts.nc Fulldom_hires_netcdf_file.nc\n",
			"Fri Sep 25 11:49:57 2015: ncks -A LAKEGRID.nc Fulldom_hires_netcdf_file.nc\n",
			"Fri Sep 25 11:49:56 2015: ncks -A ovroughrtfac.nc Fulldom_hires_netcdf_file.nc\n",
			"Fri Sep 25 11:49:55 2015: ncks -A retdeprtfac.nc Fulldom_hires_netcdf_file.nc\n",
			"Fri Sep 25 11:49:55 2015: ncks -A gw_basns.nc Fulldom_hires_netcdf_file.nc\n",
			"Fri Sep 25 11:49:54 2015: ncks -A str_order.nc Fulldom_hires_netcdf_file.nc\n",
			"Fri Sep 25 11:49:53 2015: ncks -A CHANNELGRID.nc Fulldom_hires_netcdf_file.nc\n",
			"Fri Sep 25 11:49:53 2015: ncks -A flowdirection.nc Fulldom_hires_netcdf_file.nc\n",
			"Fri Sep 25 11:49:53 2015: ncks topography.nc Fulldom_hires_netcdf_file.nc" ;
		:history_of_appended_files = "Fri Nov 11 15:46:53 2016: Appended file DEM_working/DEM1km_rclass_Fillnodata_revy.nc had following \"history\" attribute:\n",
			"Fri Nov 11 15:44:13 2016: ncpdq -O -a -y DEM1km_rclass_Fillnodata.nc DEM1km_rclass_Fillnodata_revy.nc\n",
			"Fri Nov 11 15:33:22 2016: GDAL CreateCopy( DEM1km_rclass_Fillnodata.nc, ... )\n",
			"" ;
		:proj4 = "+proj=lcc +units=m +a=6370000.0 +b=6370000.0 +lat_1=30.0 +lat_2=60.0 +lat_0=40.0 +lon_0=-97.0 +x_0=0 +y_0=0 +k_0=1.0 +nadgrids=@null +wktext  +no_defs " ;
		:processing_notes = "Created: Tue Jul 31 15:12:16 2018" ;
		:GDAL_DataType = "Generic" ;
		:NCO = "netCDF Operators version 4.7.4 (http://nco.sf.net)" ;
}
