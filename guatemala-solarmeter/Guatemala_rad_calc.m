function [ G,T,t_step,t ] = Guatemala_rad_calc( data_source,plots,file)

%This file is used to calculate the solar radiation given volI_mptage, current
%and temperature provided by the solar radiation sensors
%Originally by Nicholas Barry 13 July 2015 
% Modified for use in Guatemala in July 2021 by Lisa A. Shay
% Added print statements for max solar radiation and number of hours of
% data collection.
% Corrected units for "accumulated Power" (which is energy, Wh, not W).

%INPUTS: data_source = string name of source 'CSV' ('excel' no longer
%supported)
%        plots = 'Y' for plots on, O/W plots are off 
%        file = name of CSV file from which to read data

%OUTPUTS:
%        G = solar radiation vector (w/m^2)
%        T = temperature vector  (kelvin)
%        t_step = time of sample iteration (s)
%        t = time vector (s)

 run panel_info_adafruit.m;           %load panel information 

%select correct data source 
if strcmp('excel',data_source) == 1
        run data_import_xls.m;       %load measured values from excel (no longer supported)
    elseif strcmp('CSV',data_source) == 1    
        [T,V,I,month,day,hour,minutes,seconds,t_step,t,t_min,start_time] = Guatemala_data_import(file);           %load measured values from csv file 
    else  
        error('No data source selected ! ')
end
   
hf=3600;                     %seconds to hour conversion factor
P=I.*V;                      %calculates power output of panel (watts)

%reverse saturation current of the diode (constant for panel)
Io_ref=I_scT1/(exp((q*V_ocT1)/(n*k*T1))-1); 
 
% Preallocate memory to speed execution
E_g = zeros(1,size(t,1));
Io = zeros(1,size(t,1));
Id = zeros(1,size(t,1));
Iph = zeros(1,size(t,1));
G = zeros(1,size(t,1));

for h=1:1:size(t,1)         

    E_g(h)=1.16-.000702*(T(h).^2)./(T(h)-1108);          %calculates the bandgap energy 
    I_phT1=I_scT1*(G_ref/1000);
    a=(I_scT2- I_scT1)/I_scT1*1/(T2-T1);
    Io(h)=Io_ref*((T(h)/T1)^(3/n))*exp((1/T(h)-1/T1)*(-q*E_g(h))/(n*k)); %reverse breakdown current
    Vt_T1=(k*T1) / q; 
    X=Io(h)/(n*Vt_T1)*exp(V_ocT1/(n*Vt_T1)); % X is used to simplify line below 

    dvdi_Voc=-1.15/Ns /2;               %slope from I-V curve to T_ref 
    Rs=-dvdi_Voc- 1/X;                  %calculates series resistance of each cell 

    Id(h)=Io(h)*(exp((q*((V(h))/Ns+(Rs*I(h))/Ns))/(n*k*T(h)))-1);
    Iph(h)=I(h)+Id(h)+(V(h)+I(h)*Rs);    %photon current 
    G(h)=(G_ref*Iph(h))/(I_scT1*(1+a*(T(h)-T1)))*ef; %solves for insolation value
end 

time=transpose(t);  %transpose time vector for plotting 
%%
%this section uses data interpolation to creat a polynomial equation for
%the insolation verses time.
poly=6;                        %set order of polynomial approximation 

[p,mu]=polyfit(time,G,poly);   %crea    te polynomial fit to data 
O=polyval(p,time);             %evaluate polynomial over input range 

%calculate accumulated values of power and solar radiation 
accum_P=cumsum(P*t_step/3600);
accum_G=cumsum(G*t_step/3600);
total_P=accum_P(end,:);
total_G=accum_G(:,end);

fprintf('The total energy produced by the meter panel was %g Watt-hours \n',total_P)
fprintf('\n')
fprintf('The total Solar Radiation for the day was %g Wh/m^2 \n',total_G)
fprintf('\n')
fprintf('The maximum Solar Radiation for the day was %g W/m^2 \n',max(G))
fprintf('\n')
datacoll_start = hour(1)*60 + minutes(1);   % # minutes since midnight when data collection started
datacoll_end = hour(end)*60 + minutes(end); % # minutes since midnight when data collection stopped
fprintf('The total number of hours of data collection was %g hours \n', (datacoll_end - datacoll_start)/60)


    if plots == 'Y'
        Guatemala_solar_meter_plots(t,start_time,T,V,I,P,G,O,accum_P,accum_G,file);
    else 
        return
    end
end %end of file 


