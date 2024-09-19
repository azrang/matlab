%Originally AFRICOM_rad_calc_launch by Nicholas Barry 7 August 2015
%Updated June 2016 for use in Cameroon
%Updated July 2021 by Lisa Shay for use in Guatemala

%Use this script to launch to solar radiation calculator for one dataset
%called datalog.txt.  This script is for troubleshooting and provides all
%the information necessary to launch rad_calc.m outside of solar_sim.m. 

data_source = 'CSV';  %set input source 'CSV' (import from 'excel' no longer supported)        
plots = 'Y';          %turn plots on ('Y') or off ('N')
file='DATA11.txt'; % File to open
[G,T,t_step,t]=Guatemala_rad_calc(data_source,plots,file);