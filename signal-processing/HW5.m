%% Designing Filters for HW 5

%1a) designing chebshev 1 bandpass filter
n = 6;
n0 = n/2;
rp = 2; 
fs = 40e3;
fhi = 12e3;
flo = 10e3;
fcrit = [flo/(fs/2), fhi/(fs/2)]; %vector 

[z,p,k] = cheby1(n0, rp, fcrit)
[b1,a1] = zp2tf(z,p,k)
figure
zplane(z,p,k);
title('Chebyshev Pole Zero Plot');

%1b) getting frequency response 
f1 = linspace(0, fs/2, 1e4);
norm_f1 = (2*pi*f1)/fs;
H = freqz(b1, a1, norm_f1);  %this is the frequency response

% graphing magnitude and phase response
mag_resp(H, f1, fs/2e3, -40, 2);
phase_resp(H, f1, fs/2e3);
Hdb = 20*log10(abs(H));

%1d) examine mag response where gain falls below -30db ouside of stopband
%freq
a=find(Hdb >= -30);
level1 = f1(a(1)-1)
level2 = f1(a(end)+1)

%4: getting frequency plots and stability values

%compute and plot frequency response 
C= 10e-9;
R= 1e3;
N= 10e4;
f2 = linspace(0, 1e6, N);
Kmax = 4;
K1 = 4- 2*sqrt(2);   %K0
K2 = 0.5*K1 + .5*Kmax;
K3 = .2*K1 + .8*Kmax;

% sys1 = tf([K1/(R*C) , 0] , [1 , (4 - K1)/(R*C) , 2 / (R^2*C^2)])
%sys2 = tf([K2/(R*C) , 0] , [1 , (4 - K2)/(R*C) , 2 / (R^2*C^2)])
%sys3 = tf([K3/(R*C) , 0] , [1 , (4 - K3)/(R*C) , 2 / (R^2*C^2)])

norm_f2 = (2*pi*f2);  %rad/sec
num_k1 = [K1/(R*C), 0];
denom_k1 = [1 , (4 - K1)/(R*C) , 2 / (R^2*C^2)];
H_k1 = freqs(num_k1, denom_k1, norm_f2);
num_k2 = [K2/(R*C) , 0];
denom_k2 = [1 , (4 - K2)/(R*C) , 2 / (R^2*C^2)];
H_k2 = freqs(num_k2, denom_k2, norm_f2);
num_k3 = [K3/(R*C) , 0];
denom_k3 = [1 , (4 - K3)/(R*C) , 2 / (R^2*C^2)];
H_k3 = freqs(num_k3, denom_k3, norm_f2);

% graphing magnitude and phase response
%xmax, ymin, ymax
mag_resp(H_k1, f2, 1e3, -35 , -5);
phase_resp(H_k1, f2, 1e3);
mag_resp(H_k2, f2, 1e3, -30 , 10);
phase_resp(H_k2, f2, 1e3);
mag_resp(H_k3, f2, 1e3, -30, 20);
phase_resp(H_k3, f2, 1e3);
    
%generate bode plots
%from 1KHz to 1MHz
ff = logspace(3,6,N);
norm_ff = (2*pi*ff);  %rad/sec
num_k1 = [K1/(R*C), 0];
denom_k1 = [1 , (4 - K1)/(R*C) , 2 / (R^2*C^2)];
H_k1_b = freqs(num_k1, denom_k1, norm_ff);
num_k2 = [K2/(R*C) , 0];
denom_k2 = [1 , (4 - K2)/(R*C) , 2 / (R^2*C^2)];
H_k2_b = freqs(num_k2, denom_k2, norm_ff);
num_k3 = [K3/(R*C) , 0];
denom_k3 = [1 , (4 - K3)/(R*C) , 2 / (R^2*C^2)];
H_k3_b = freqs(num_k3, denom_k3, norm_ff);

Hdb11 = 20*log10(abs(H_k1));
Hdb22 = 20*log10(abs(H_k2));
Hdb33 = 20*log10(abs(H_k3));

%for k1
Hdb1 = 20*log10(abs(H_k1_b));
figure;
subplot(2, 1, 1);
semilogx(ff/1e3, Hdb1);
xlabel("Frequency (kHz)");
ylabel("|H| (dB)");
title("Bode K1 Mag");
grid on;

Hph1 = rad2deg(unwrap(angle(H_k1_b)));
subplot(2, 1, 2);
semilogx(ff/1e3, Hph1);
xlabel("Frequency kHz");
ylabel("Phase (deg)");
title("Bode K1 Phase");
grid on;

%for k2
Hdb2 = 20*log10(abs(H_k2_b));
figure;
subplot(2, 1, 1);
semilogx(ff/1e3, Hdb2);
xlabel("Frequency (kHz)");
ylabel("|H| (dB)");
title("Bode K2 Mag");
grid on;

Hph2 = rad2deg(unwrap(angle(H_k2_b)));
subplot(2, 1, 2);
semilogx(ff/1e3, Hph2);
xlabel("Frequency kHz");
ylabel("Phase (deg)");
title("Bode K2 Phase");
grid on;

%for k3
Hdb3 = 20*log10(abs(H_k3_b));
figure;
subplot(2, 1, 1);
semilogx(ff/1e3, Hdb3);
xlabel("Frequency (kHz)");
ylabel("|H| (dB)");
title("Bode K3 Mag");
grid on;

Hph3 = rad2deg(unwrap(angle(H_k3_b)));
subplot(2, 1, 2);
semilogx(ff/1e3, Hph3);
xlabel("Frequency kHz");
ylabel("Phase (deg)");
title("Bode K3 Phase");
grid on;

%4c) compute wn and Q  
wn = sqrt(2)/(R*C)
Q1 = sqrt(2)/(4-K1)
Q2 = sqrt(2)/(4-K2)
Q3 = sqrt(2)/(4-K3)

%4d) frequencies at which peak occurs, and 3db below

% we need the frequency where the max gain occurs with mag responses
H1_max = max(Hdb11);
ind1 = find(Hdb11 == H1_max);
f_peak1 = f2(ind1)
H2_max = max(Hdb22);
ind2 = find(Hdb22 == H2_max);
f_peak2 = f2(ind2)
H3_max = max(Hdb33);
ind3 = find(Hdb33 == H3_max);
f_peak3 = f2(ind3)

f_lo_1 = f_lo(f2, Hdb11)
f_lo_2 = f_lo(f2, Hdb22)
f_lo_3 = f_lo(f2, Hdb33)
f_hi_1 = f_hi(f2, Hdb11)
f_hi_2 = f_hi(f2, Hdb22)
f_hi_3 = f_hi(f2, Hdb33)

%4e) checking peak freq with natural freq wn
wn1 = 2*pi*f_peak1;
e1 = isequal(wn1, wn) %these should be equal but I do not know why they are not
wn2 = 2*pi*f_peak2;
e2 = isequal(wn2, wn)
wn3 = 2*pi*f_peak3;
e3 = isequal(wn3, wn)
if ((wn1- wn)/1e3 < .05)
    if((wn2 - wn)/1e3 < .05)
        if((wn3 - wn)/1e3 < .05)
            disp("peak frequency is equal to natural frequency for all cases")
        end
    end
end


%4f) computing 3dB bandwidth
delta_f1 = f_hi_1 - f_lo_1
delta_f2 = f_hi_2 - f_lo_2
delta_f3 = f_hi_3 - f_lo_3
f_ctr1= sqrt(f_hi_1*f_lo_1)
f_ctr2= sqrt(f_hi_2*f_lo_2)
f_ctr3= sqrt(f_hi_3*f_lo_3)

%4g) 
% The peak freq should be center freq
if ((f_ctr1 - f_peak1)/1e3 < .05)
    if((f_ctr2 - f_peak2)/1e3 < .05)
        if((f_ctr3 - f_peak3)/1e3 < .05)
            disp("peak frequency is equal to center frequency for all cases")
        end
    end
end

%Q should be a measure of the spectral peak
if (Q1 - (f_ctr1/delta_f1) < .05)
    if(Q2 - (f_ctr2/delta_f2) < .05)
        if(Q3 - (f_ctr3/delta_f3) < .05)
            disp("Q is a measure of spectral peak")
        end
    end
end

function mag_resp= mag_resp(H, f, xmax, ymin, ymax)
Hdb = 20*log10(abs(H));
figure;
subplot(2, 1, 1);
mag_resp = plot(f/1e3, Hdb);
xlabel("Frequency kHz");
ylabel("|H| (dB)");
xlim([0 xmax]);
ylim([ymin ymax]);
title("Magnitude Response");
grid on;
end 

function phase_resp = phase_resp(H, f, xmax)
Hph = rad2deg(unwrap(angle(H)));
subplot(2, 1, 2);
phase_resp = plot(f/1e3, Hph);
xlabel("Frequency kHz");
ylabel("Phase (deg)");
xlim([0 xmax]);
title("Phase Response");
grid on;
end


function f_lo = f_lo(f, Hdb)
lim = max(Hdb) - 3;
%create logical array, and getting nonzero values
indd = find(Hdb < lim, 1, 'first');
vals = f(indd);
arr = Hdb < lim;
LO = find(arr ==0 , 1, 'first');
f_lo = f(LO);
end

function f_hi = f_hi(f, Hdb)
lim = max(Hdb) - 3;
%create logical array, and getting nonzero values
indd = find(Hdb < lim, 1, 'last'); %only difference is getting last indices
vals = f(indd);
arr = Hdb < lim;
HI = find(arr ==0 , 1, 'last');
f_hi = f(HI);
end 











