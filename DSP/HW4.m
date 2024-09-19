clear all;
close all;

%% HW 4 Quantization

%dig ellip filter
n_ellip = 4; %since bandstop, need halved from 8
Rp = 1.5; %passband ripple
Rs = 30; %stopband ripple
Wp_dig = [0.3 0.6]; %normalized passband freqs
fs = 40e6;
f1 = linspace(0, fs/2, 1e4);
norm_f1 = (2*pi*f1)/fs; 

[z, p, k] = ellip(n_ellip, Rp, 40, Wp_dig, 'bandpass');
[b, a] = zp2tf(z, p, k); %transfer function
h = freqz(b, a, norm_f1);

figure;
Hdb1 = 20*log10(abs(h));
plot(f1, Hdb1);
title('tf mag response');

%convert to sos form
[hsos_up, g_up] = tf2sos(b, a, 'up', inf) %up: reduces likelihood of overflow
[hsos_down, g_down] = tf2sos(b, a, 'down', inf) %down: reduces roundoff noise gain
[z_up1, p_up1, k_up1] = tf2zp(hsos_up(1, 1:3), hsos_up(1, 4:6)); % b and a in FIRST stage
[z_up2, p_up2, k_up2] = tf2zp(hsos_up(2, 1:3), hsos_up(2, 4:6)); % b and a in second stage
[z_up3, p_up3, k_up3] = tf2zp(hsos_up(3, 1:3), hsos_up(3, 4:6));
[z_up4, p_up4, k_up4] = tf2zp(hsos_up(4, 1:3), hsos_up(4, 4:6));
z_up_mags = [abs(z_up1) abs(z_up2) abs(z_up3) abs(z_up4)]
p_up_mags = [abs(p_up1) abs(p_up2) abs(p_up3) abs(p_up4)]
%yes they are are increasing
[z_down1, p_down1, k_down1] = tf2zp(hsos_down(1, 1:3), hsos_down(1, 4:6));
[z_down2, p_down2, k_down2] = tf2zp(hsos_down(2, 1:3), hsos_down(2, 4:6));
[z_down3, p_down3, k_down3] = tf2zp(hsos_down(3, 1:3), hsos_down(3, 4:6));
[z_down4, p_down4, k_down4] = tf2zp(hsos_down(4, 1:3), hsos_down(4, 4:6));
z_down_mags = [abs(z_down1) abs(z_down2) abs(z_down3) abs(z_down4)]
p_down_mags = [abs(p_down1) abs(p_down2) abs(p_down3) abs(p_down4)]
%yes they are decreasing 

H_up1 = freqz( hsos_up(1, 1:3), hsos_up(1, 4:6));
H_up2 = freqz( hsos_up(2, 1:3), hsos_up(2, 4:6));
H_up3 = freqz( hsos_up(3, 1:3), hsos_up(3, 4:6));
[H_up4, w] = freqz( hsos_up(4, 1:3), hsos_up(4, 4:6));

H_down1 = freqz( hsos_down(1, 1:3), hsos_down(1, 4:6));
H_down2 = freqz( hsos_down(2, 1:3), hsos_down(2, 4:6));
H_down3 = freqz( hsos_down(3, 1:3), hsos_down(3, 4:6));
H_down4 = freqz( hsos_down(4, 1:3), hsos_down(4, 4:6));

figure(); 
hold on;
H_up = H_up1; H_up = H_up/max(abs(H_up)); plot(w/pi, 10*log10(abs(H_up)));
H_up = H_up .* H_up2; H_up = H_up/max(abs(H_up));  plot(w/pi, 10*log10(abs(H_up)));
H_up = H_up .* H_up3; H_up = H_up/max(abs(H_up)); plot(w/pi, 10*log10(abs(H_up)));
H_up = H_up .* H_up4; H_up = H_up/max(abs(H_up)); plot(w/pi, 10*log10(abs(H_up)));
legend({'1st stage', '2nd stage', '3rd stage', '4th stage'});
xlim([0 1]);
ylim([-40 0]);
title("SOS Up Realization mag response");
xlabel("Normalized to Nyquist Frequency (rad/sec)");
ylabel("Magnitude (dB)");

figure(); 
hold on;
H_down = H_down1; H_down = H_down/max(abs(H_down)); plot(w/pi, 10*log10(abs(H_down)));
H_down = H_down .* H_down2; H_down = H_down/max(abs(H_down)); plot(w/pi, 10*log10(abs(H_down)));
H_down = H_down .* H_down3; H_down = H_down/max(abs(H_down));  plot(w/pi, 10*log10(abs(H_down)));
H_down = H_down .* H_down4; H_down = H_down/max(abs(H_down));  plot(w/pi, 10*log10(abs(H_down)));
legend({'1st stage', '2nd stage', '3rd stage', '4th stage'});
xlim([0 1]);
ylim([-40 0]);
title("SOS down Realization mag response");
xlabel("Normalized to Nyquist Frequency (rad/sec)");
ylabel("Magnitude (dB)");

%yes they are scaled 
hsos_down(2, 1:3)./hsos_up(3, 1:3);
hsos_down(4, 1:3)./hsos_up(1, 1:3); 


%allpass
y_49 = fi(-.4954, 1, 5, 3);  %get most precise representation with 5 bits
y_76 = fi(.7626, 1, 5, 3);
y_1 = fi(-1.0101, 1, 5, 3);
%original
y_13 = fi(.1336, 1, 5, 3);
y_05 = fi(.0563, 1, 5, 3);
y_15 = fi(-1.5055, 1, 5, 3);
y_12 = fi(1.2630, 1, 5, 3);
y_37 = fi(-.3778, 1, 5, 3);

og_quant_zplus1 = (y_13.data+y_05.data+y_05.data+y_13.data)/(1+y_15.data+y_12.data+y_37.data)
og_quant_zneg1 = (y_13.data-y_05.data+y_05.data-y_13.data)/(1-y_15.data+y_12.data-y_37.data)
og_ideal_zplus1 = (.1336+.0563+.0563+.1336)/(1-1.5055+1.263-.3778) 
og_ideal_zneg1 = (.1336-.0563+.0563-.1336)/(1+1.5055+1.263+.3778) 
all_quant_zplus1 = 0.5*((y_49.data+1)/(1+y_49.data)+(y_76.data+y_1.data+1)/(1+y_1.data+y_76.data))
all_quant_zneg1 = 0.5*((y_49.data-1)/(1-y_49.data)+(y_76.data-y_1.data+1)/(1-y_1.data+y_76.data))
all_ideal_zplus1 = 0.5*((-.4954+1)/(1-.4954)+(.7626+1.0101+1)/(1+.0101+.7626))
all_ideal_zneg1 = 0.5*((-.4954-1)/(1+.4954)+(.7626-1.0101+1)/(1-1.0101+.7626))

%elliptic lpf
f1 = linspace(0, pi, 1e4);
b_ellip = [.1336 .0563 .0563 .1336];
a_ellip = [1 -1.5055 1.2630 -.3778];
H_ellip = freqz(b_ellip, a_ellip, f1);  
H_quant_og = freqz([y_13.data y_05.data y_05.data y_13.data], [1 y_15.data y_12.data y_37.data], f1);
H_quant_allp = freqz(0.5*[y_49.data 1], [1 y_49.data], f1) + freqz(0.5*[y_76.data y_1.data 1], [1 y_1.data y_76.data], f1); 

max_allp = max(abs(H_ellip-H_quant_og))
max_og = max(abs(H_ellip- H_quant_allp))

Hdb1 = 10*log10(abs(H_ellip));
Hdb2 = 10*log10(abs(H_quant_og));
Hdb3 = 10*log10(abs(H_quant_allp));

figure;
plot(f1, Hdb1);
hold on;
plot(f1, Hdb2);
plot(f1, Hdb3);
xlabel("Frequency normalized");
ylabel("|H| (dB)");
title("Quantized allpass and original & elliptic mag response")
hold off;
legend('elliptic filter', 'quantized original', 'quantized allpass')
ylim([-40 5]);

%d
% The infinite precision filter overlaps completely with the quantized allpass 
% filter, while the quantized original filter is around 2 dB lower in the
% passband. 
%the local maxima in the passband for the elliptic filter and allpass
%filter occurs in the same spot. the local maxima for the quantized
%original occus at a lower dB at around the same frequency in the passband.
%the local minima in teh stopband occurs at a lower frequency than both the
%allpasss and original filter.
%the maximum gain in the stopband is highest for the quantized original
%graph, then the original graph, and then quantized allpass has the lowest
%gain.

%group delay is the derivative of the phase response 
tau1 = diff(unwrap(angle(H_ellip)));
tau2 = diff(unwrap(angle(H_quant_og))); 
tau3 = diff(unwrap(angle(H_quant_allp))); 
ang_f1 = linspace(0, pi, 1e4-1);

figure;
plot(ang_f1, tau1);
hold on;
plot(ang_f1, tau2);
plot(ang_f1, tau3);
xlabel("radians angle");
ylabel("Phase response");
title("Quantized allpass and original & elliptic phase response")
hold off;
legend('elliptic filter', 'quantized original', 'quantized allpass')
% The group delay of the quantized allpass is closer to the infinite precision filter
% so it is relatively more sensitive than the quantized original;

%making into one fraction
b1a2 = conv([y_49.data 1], [1 y_1.data y_76.data]);
b2a1 = conv([y_76.data y_1.data 1], [1 y_49.data]);
a1a2 = conv([1 y_49.data],  [1 y_1.data y_76.data]);
[z_ellip, p_ellip, k_ellip] = tf2zp(b_ellip, a_ellip)
[z_og, p_og, k_og] = tf2zp([y_13.data y_05.data y_05.data y_13.data], [1 y_15.data y_12.data y_37.data])
[z_allp, p_allp, k_allp] = tf2zp(0.5*(b1a2+b2a1), a1a2)

1- abs(p_ellip) %to check whether inside unit circle
1- abs(p_og)
1- abs(p_allp)
%yes the poles stay inside the unit circle

%the pole at z=-1 is still there for the quantized original filter, as
%well as for the allpass