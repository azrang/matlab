function [T,V,I,month,day,hour,minutes,seconds,t_step,t,t_min,start_time] = Guatemala_data_import(file)
%Originally AFRICOM_data_import by Nicholas Barry 8 June 2016 
%Modified for use in Guatemala in July 2021 by Lisa A. Shay
%Modified columns based on data collected by Azra Rangwala

%This script loads the data from a CSV file into MATLAB. It targets a file
%called datalog.txt and reads the data produced by the solar radiation
%meters according to the rubric below.  The datalog.txt file must be in the
%same MATLAB folder and included to function. 

%This script is designed to work on one day's data at a time.  Multiple
%days can be overlayed or sorted in excel prior to import. 

%column 7: observed temperature (in F)
%column 8: measured voltage (V)
%column 9: measured current  (mA)
%column 2: Month
%column 3: Day
%column 4: Hour
%column 5: Minute
%column 6: Seconds 
file='DATA11.txt';
raw_data=readmatrix(file);         %reads data from file 

%read data files from datafile
T=raw_data(:,7);
V=raw_data(:,8);
I=raw_data(:,9);
month=raw_data(:,2);
day=raw_data(:,3);
hour=raw_data(:,4);
minutes=raw_data(:,5);
seconds=raw_data(:,6);



%calculates time data 
t_step=(minutes(2)-minutes(1))*60;  %calculates step time 
                                       %t_step must be constant for all data
                                       %calculates t_step in seconds
if (t_step == 0)                    %measurements less than 1 minute apart
    t_step=seconds(2)-seconds(1);
end                                  

t = ones(size(seconds));             %initialize t vector (seconds)
                                     %this vector is a seconds uptime
                                     %counter nessecary for later
                                     %operations                                    
%create uptime counter for incrementation later                                      
u = 0;                                       
for u = 1:1:size(seconds)            %increments t_step size for each sample
    t(u) = t_step*u;                 %fill out t vector 
end 
clear u 

t_min=min(t)                        %calcualtes the minium time value (seconds) 
                                    %this cannot be equal to zero !                                    
                                     
%set start time (daily)                                    
start_time = hour(1);     %sets start time in hours (military time)

%Temperature conversion
T=T*(5/9)+273;                       %T from F to Kelvin                    

%end of file
