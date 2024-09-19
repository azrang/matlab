%Nicholas Barry 26 June 2015
%Panel Values: this file contains all of the values associated with the
%solar panel listed below:

%Panel Type: 2W 6V adafruit solar panel monocrystaline panel  
%weblink: http://www.adafruit.com/products/200

%Constants
k=1.38e-23;          %Boltzmann's Constant 
q=1.602e-19;         %Coulumb's Constant
n=1.3;               %Diode Ideality Factor (1<n<2)

%Panel values 
Ns=12;               %number of series cells in panel 

V_ocT1=6.15;         %open circuit voltage of panel (volts)
V_ocT1=V_ocT1/Ns;    %calculates the open circuit voltage for each cell
V_ocT2=6.0675;       %open circuit voltage of panel (volts) 
V_ocT2=V_ocT2/Ns;    %calculates the open circuit voltage for each cell

I_scT1=.375;         %short circuit current the panel amps (amps)
I_scT2=.385;         %short circuit current the panel amps (amps)

T_ref=25;            %cell refrence temperature (celcius)
T1=T_ref+273;        %converts to kelvin 
T_ref2=50;           %second temp refrence pulled from the manufactures curve 
T2=T_ref2+273;       %converts to kelvin    

G_ref=1000;          %cell light intensity refrence (w/m^2)

P_max=4.29;          %max power in watts
V_mp=6.252;          %Volatage at maxium power point 
I_mp=.349;           %current at maxium power point 

ef=.1781;            %panel efficiency 

%end of file