%This script loads the data from a CSV file into MATLAB. It targets a file
%called datalog.txt and reads the data produced by the solar radiation
%meters according to the rubric below.
%column 1: data type (1 for singular, 2 for mulitpule devices)
%column 7: observed temperature (in celcius)
%column 8: measured voltage
%column 9: measured current
%column 5: light intensity
%column 6: time stamp (need data format) (t)
%column 7: device number (1-6)
%cd '/MATLABDrive/Published' 
filename='DATA21.txt';
M= readmatrix(filename);
%data=csvread('C:\Users\Azra Rangwala\OneDrive - The Cooper Union for the Advancement of Science and Art\DATA21.TXT', 0,6); %reads data from file
%read data files from radiationdatafile.m
T=M (:,7);
V=M (:,8);
I=M (:,9);
L=M (:,13);
t=M (:,12);
%calculates time data
t_min=min(M(:,12)); %calcualtes the minium time value
%this cannot be equal to zero !
t_max=max(M(:,12)); %calcualtes the maxium time value
t_step=M(2,12)-M(1,12); %cacluatles t step time
%t step must be consistant for all data
if t_min <= 0 %checks to ensure t min is not 0
    error 't min is less than zero !'
end

%Temperature conversion
T=T+273; %T from celcuis to Kelvin
%end of file