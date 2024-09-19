
function Guatemala_solar_meter_plots(t,start_time,T,V,I,P,G,O,accum_P,accum_G,file)
close all 
%originally AFRICOM_solar_meter_plots.m by Nicholas Barry 8-7-2015
%Modified for use in Guatemala in July 2021 by Lisa A. Shay
%Temperature conversion corrected by Joya Debi 11 July 2021
%12 July 2021 Added filename (i.e. date) to figure titles and ensured
%Figure 4's y-axis had a minimum of 0

%this script contains the plots for the single meter solar radiation calcs.
%It is called by putting plots equal to 'Y' when using Guatemala_rad_calc.m

%Variables vs Time
figure (1)   
plot(t/3600+start_time,(T-273)/mean(T),'r')
hold on 
plot(t/3600+start_time,V,'b')         
plot(t/3600+start_time,I,'g')         
%plot(t,L,'y')
plot(t/3600+start_time,P,'c')
legend('Temperature (normalized)','Voltage (V)','Current (A)','Power(W)')
xlabel('Time (Hours)')
ylabel('Insolation Variables')
title(['Variables vs Time for ' file])
hold off 

%Temperature vs Power Output
figure (2)   
plot(((T-273.15)*9/5 + 32),P,'o')  % Convert T in Kelvin to F
xlabel('Temperature (F)')
ylabel('Power (W)')
title(['Temperature vs Power Output for ' file])

%plots solar radiation vs time
figure(3)
plot(start_time+t/3600,G,'+')
legend('Solar Radiation')
xlabel('Time (Hours)')
ylabel('Solar Radiation W/m^2')
title(['Solar Radiation for ' file])

%plot polynomial model of solar radiation
figure(4)
plot(t/3600+start_time,O,'r')
hold on 
plot(t/3600+start_time,G,'+')
legend('Model','Measurements')
xlabel('Time (hours)')
ylabel('Solar Radiation W/m^2')
title(['Solar Radiation w/ Polynomial Approximation for ' file])
y1 = ylim;
if y1(1) < 0
    ylim([0 y1(2)]);  % make sure the lower limit is 0 or greater
end
hold off 

%plots accumulation of real power 
figure (5)
plot(t/3600+start_time,accum_P)
xlabel('Time (hours)')
ylabel('Energy (Wh)')
title(['Total Energy Produced for ' file])

%plots accumulation of solar radiation
figure (6)
plot(t/3600+start_time,accum_G)
xlabel('Time (hours)')
ylabel('Solar Radiation (Wh/m^2)')
title(['Total Solar Radiation Received for ' file])

%log transformed power --- showed that it is not exponentially related
%figure(7)
%plot(((T-273.15)*9/5 + 32), log(P), 'o')
%xlabel('Temperature (F)')
%ylabel('log Power')
%title('Temperature vs log of Power Output')


