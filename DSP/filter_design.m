%% DSP HW 2
close all

a1 = [6 13 6 0]; %denom coeff, positiv powers 
b1 = [-2 5 -2];   %num coeff
[z,p, k] = tf2zp(b1, a1); %does not show multiplicity so need to account for that

%% designing IIR filters for passband filters
fs = 40e6; 
Wp_dig = [9e6 12.5e6]/20e6;  %normalized with fs/2
Ws_dig = [9.5e6 12e6]/20e6;
Wp_anal = [9e6 12.5e6]*2*pi;  %multiply by 2pi to convert to rad/s
Ws_anal = [9.5e6 12e6]*2*pi;
f1 = linspace(0, fs/2, 1e4);
norm_f1 = (2*pi*f1)/fs;  %for digital
ang_f1 = 2*pi*f1; %for analog

%BUTTERWORTH
[n_butt, Wn_butt] = buttord(Wp_dig, Ws_dig, 1.5, 40); %(Rp db of passband ripple, Rs db of attenuation in stopband) 
n_butt = n_butt*2
[z, p, k] = butter(n_butt, Wn_butt, 'stop');
figure;
zplane(z, p);
title('butter digital bandstop');
[b_butt, a_butt] = butter(n_butt, Wn_butt, 'stop');
H_butt = freqz(b_butt, a_butt, norm_f1);  
mag_resp(H_butt, f1, fs/2, -50, 2);
phase_resp(H_butt, f1, fs/2);

[n_butt_anal, Wn_butt_anal] = buttord(Wp_anal, Ws_anal, 1.5, 40, 's'); %%returns Wn shoudl be angular freq
n_butt_anal = n_butt_anal *2
[z_anal, p_anal, k] = butter(n_butt_anal, Wn_butt_anal, 'stop', "s");  %Wn shoudl be angular freq
figure;
zplane(z_anal, p_anal);
title('butter analog bandstop');
[b_butt_anal, a_butt_anal] = butter(n_butt_anal, Wn_butt_anal, 'stop', "s");
H_butt_anal = freqs(b_butt_anal, a_butt_anal, ang_f1); %at angular freqs
mag_resp(H_butt_anal, f1, fs/2, -50, 2);
phase_resp(H_butt_anal, f1, fs/2);

%CHEBY1
[n_cheby1_dig, Wp_cheby1_dig] = cheb1ord(Wp_dig, Ws_dig, 1.5, 40);
n_cheby1_dig = n_cheby1_dig*2
[z, p, k] = cheby1(n_cheby1_dig, 1.5, Wp_dig, 'stop');
figure;
zplane(z, p);   %cheby1 rolls off faster than cheby2, equiripple in passband
title('cheby1 digital bandstop');
[b_cheby1, a_cheby1] = cheby1(n_cheby1_dig, 1.5, Wp_dig, 'stop');
H_cheby1 = freqz(b_cheby1, a_cheby1, norm_f1);  
mag_resp(H_cheby1, f1, fs/2, -80, 2);
phase_resp(H_cheby1, f1, fs/2);

[n_cheby1_anal, Wp_cheby1_anal] = cheb1ord(Wp_anal, Ws_anal, 1.5, 40, 's');
n_cheby1_anal = n_cheby1_anal*2
[z, p, k] = cheby1(n_cheby1_anal, 1.5, Wp_anal, 'stop', 's');
figure;
zplane(z, p);   %cheby1 rolls off faster than cheby2, equiripple in passband
title('cheby1 analog bandstop');
[b_cheby1_anal, a_cheby1_anal] = cheby1(n_cheby1_anal, 1.5, Wp_anal, 'stop', 's');
H_cheby1_anal = freqs(b_cheby1_anal, a_cheby1_anal, ang_f1);  
mag_resp(H_cheby1_anal, f1, fs/2, -80, 2);
phase_resp(H_cheby1_anal, f1, fs/2);

%CHEBY2
[n_cheby2_dig, Ws_cheby2_dig] = cheb2ord(Wp_dig, Ws_dig, 1.5, 40);
n_cheby2_dig = n_cheby2_dig*2
[z, p, k] = cheby2(n_cheby2_dig, 1.5, Ws_dig, 'stop');
figure;
zplane(z, p);   %cheby1 rolls off faster than cheby2, equiripple in passband
title('cheby2 digital bandstop');
[b_cheby2, a_cheby2] = cheby2(n_cheby2_dig, 1.5, Ws_dig, 'stop');
H_cheby2 = freqz(b_cheby2, a_cheby2, norm_f1);  
mag_resp(H_cheby2, f1, fs/2, -80, 2);
phase_resp(H_cheby2, f1, fs/2);

[n_cheby2_anal, Wp_cheby2_anal] = cheb2ord(Wp_anal, Ws_anal, 1.5, 40, 's');
n_cheby2_anal = n_cheby2_anal*2
[z, p, k] = cheby2(n_cheby2_anal, 1.5, Ws_anal, 'stop', 's');
figure;
zplane(z, p);   %cheby1 rolls off faster than cheby2, equiripple in passband
title('cheby2 analog bandstop');
[b_cheby2_anal, a_cheby2_anal] = cheby2(n_cheby2_anal, 1.5, Ws_anal, 'stop', 's');
H_cheby2_anal = freqs(b_cheby2_anal, a_cheby2_anal, ang_f1);  
mag_resp(H_cheby2_anal, f1, fs/2, -80, 2);
phase_resp(H_cheby2_anal, f1, fs/2);

%ELLIP
[n_ellip, Wn_ellip] = ellipord(Wp_dig, Ws_dig, 1.5, 40); %(Rp db of passband ripple, Rs db of attenuation in stopband) 
n_ellip = n_ellip*2
[z, p, k] = ellip(n_ellip, 1.5, 40, Wp_dig, 'stop');
figure;
zplane(z, p);
title('ellip digital bandstop');
[b_ellip, a_ellip] = ellip(n_ellip, 1.5, 40, Wp_dig, 'stop');
H_ellip = freqz(b_ellip, a_ellip, norm_f1);  
mag_resp(H_ellip, f1, fs/2, -50, 2);
phase_resp(H_ellip, f1, fs/2);

[n_ellip_anal, Wn_ellip_anal] = ellipord(Wp_anal, Ws_anal, 1.5, 40, 's');
n_ellip_anal= n_ellip_anal*2
[z, p, k] = ellip(n_ellip_anal, 1.5, 40, Wp_anal, 'stop', 's');
figure;
zplane(z, p);
title('ellip analog bandstop');
[b_ellip_anal, a_ellip_anal] = ellip(n_ellip_anal, 1.5, 40, Wp_anal, 'stop', 's');
H_ellip_anal = freqs(b_ellip_anal, a_ellip_anal, ang_f1);  
mag_resp(H_ellip_anal, f1, fs/2, -50, 2);
phase_resp(H_ellip_anal, f1, fs/2);

%for each filter, the peak passband gain looks to be at 0db for 
% each of the graphs. They approach 0db.

%at stopband edges: 9.5e6  12e6, should be less than -40db
% the gain is lower than -40db, and significantly lower for the analog filters. 
% the butterworth dig and analog filter has pretty close to stopband edge target gain.
% the cheby2 gain looks off, although it looks right on the graph *shrug*

ind_stop_1 = find(norm_f1 >= (2*pi*9.5e6)/fs, 1)
ind_stop_2 = find(norm_f1 >= (2*pi*12e6)/fs, 1)
butt_stop_gain = [20*log10(abs(H_butt(ind_stop_1))) 20*log10(abs(H_butt(ind_stop_2)))]
butt_stop_gain_anal = [20*log10(abs(H_butt_anal(ind_stop_1))) 20*log10(abs(H_butt_anal(ind_stop_2)))]
cheby1_stop_gain = [20*log10(abs(H_cheby1(ind_stop_1))) 20*log10(abs(H_cheby1(ind_stop_2)))]
cheby1_stop_gain_anal = [20*log10(abs(H_cheby1_anal(ind_stop_1))) 20*log10(abs(H_cheby1_anal(ind_stop_2)))]
cheby2_stop_gain = [20*log10(abs(H_cheby2(ind_stop_1))) 20*log10(abs(H_cheby2(ind_stop_2)))]
cheby2_stop_gain_anal = [20*log10(abs(H_cheby2_anal(ind_stop_1))) 20*log10(abs(H_cheby2_anal(ind_stop_2)))]
ellip_stop_gain = [20*log10(abs(H_ellip(ind_stop_1))) 20*log10(abs(H_ellip(ind_stop_2)))]
ellip_stop_gain_anal = [20*log10(abs(H_ellip_anal(ind_stop_1))) 20*log10(abs(H_ellip_anal(ind_stop_2)))]

%at passband edges: 9e6    12.5e6 should be -1.5
ind_pass_1 = find(norm_f1 >= (2*pi*9e6)/fs, 1);
ind_pass_2 = find(norm_f1 >= (2*pi*12.5e6)/fs, 1);
butt_pass_gain = [20*log10(abs(H_butt(ind_pass_1))) 20*log10(abs(H_butt(ind_pass_2)))];
butt_pass_gain_anal = [20*log10(abs(H_butt_anal(ind_pass_1))) 20*log10(abs(H_butt_anal(ind_pass_2)))];
cheby1_pass_gain = [20*log10(abs(H_cheby1(ind_pass_1))) 20*log10(abs(H_cheby1(ind_pass_2)))];
cheby1_pass_gain_anal = [20*log10(abs(H_cheby1_anal(ind_pass_1))) 20*log10(abs(H_cheby1_anal(ind_pass_2)))];
cheby2_pass_gain = [20*log10(abs(H_cheby2(ind_pass_1))) 20*log10(abs(H_cheby2(ind_pass_2)))];
cheby2_pass_gain_anal = [20*log10(abs(H_cheby2_anal(ind_pass_1))) 20*log10(abs(H_cheby2_anal(ind_pass_2)))];
ellip_pass_gain = [20*log10(abs(H_ellip(ind_pass_1))) 20*log10(abs(H_ellip(ind_pass_2)))];
ellip_pass_gain_anal = [20*log10(abs(H_ellip_anal(ind_pass_1))) 20*log10(abs(H_ellip_anal(ind_pass_2)))];

%% designing digital FIR filters, linear phase doubles filter order

%Equiripple
rp = 1.5;
rs= 40; 
dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20) (10^(rp/20)-1)/(10^(rp/20)+1)];   %tolerances corresponds to bands
[n, fo, ao, w] = firpmord([9e6 9.5e6 12e6 12.5e6], [1 0 1], dev, fs); %a0 denotes number of bands
b = firpm(n, fo, ao, w);
filt_length_equi = length(b);  %filterorder - 1 
figure;
stem(b);
[z, p, k] = tf2zpk(b,1);
figure;
zplane(z, p);
title('equiripple FIR digital bandstop');
H_equi = freqz(b, 1, norm_f1);
mag_resp(H_equi, f1, fs/2, -50, 2);
phase_resp(H_equi, f1, fs/2);
weights_ratio = w(1)/w(2)
tolerance_ratio = dev(2)/dev(1) %yes they are the same
max_gain_equi = max(max(20*log10(abs(H_equi(1:ind_pass_1)))), max(20*log10(abs(H_equi(ind_pass_2:length(norm_f1))))));
min_gain_equi = min(min(20*log10(abs(H_equi(1:ind_pass_1)))), min(20*log10(abs(H_equi(ind_pass_2:length(norm_f1))))));
min_max_equi = max_gain_equi-min_gain_equi
%this is about .2dB off from being 1.5dB different
max_gain_equi_stop = max(20*log10(abs(H_equi(ind_stop_1:ind_stop_2))))
%it is less than -30dB, meets specification

%kaiser window, should have greater order 
dev = [10^(rp/20) 10^(-rs/20) 10^(rp/20)];
[n, Wn, beta, ftype] = kaiserord([9e6 9.5e6 12e6 12.5e6],[1 0 1], dev, fs); %beta is sindow shape
h = fir1(n, Wn, 'stop', kaiser(n+1, beta));
filt_length_kaiser = length(h);  %filterorder - 1 
figure;
stem(h);
[z, p, k] = tf2zpk(h,1);
figure;
zplane(z, p);
title('kaiser FIR digital bandstop');
H_kaiser = freqz(h, 1, norm_f1);
mag_resp(H_kaiser, f1, fs/2, -50, 2);
phase_resp(H_kaiser, f1, fs/2);
max_gain_kaiser = max(max(20*log10(abs(H_kaiser(1:ind_pass_1)))), max(20*log10(abs(H_kaiser(ind_pass_2:length(norm_f1))))))
min_gain_kaiser = min(min(20*log10(abs(H_kaiser(1:ind_pass_1)))), min(20*log10(abs(H_kaiser(ind_pass_2:length(norm_f1))))))
min_max_kaiser = max_gain_kaiser-min_gain_kaiser
%this is not 1.5dB apart, it is a magnitude of 10 off, but I think this
%makes sense since kaiser has less ripple in the passband from the pics
max_gain_kaiser_stop = max(20*log10(abs(H_kaiser(ind_stop_1:ind_stop_2))))
%it is less than -30dB, meets specification

% e) My kaiser filter has a lower passband ripple than the specifications.
% To meet this requirement, I would increase rp to greater than 1.5dB so
% that we can acheive the 1.5dB ripple in the passband. Stopband
% attenuation is met for both filters. 


%% mag and phase responses
function mag_resp= mag_resp(H, f, xmax, ymin, ymax)
Hdb = 20*log10(abs(H));
figure;
subplot(2, 1, 1);
mag_resp = plot(f, Hdb);
xlabel("Frequency Mhz");
ylabel("|H| (dB)");
xlim([0 xmax]);
ylim([ymin ymax]);
title("Magnitude Response");
grid on;
end 

function phase_resp = phase_resp(H, f, xmax)
Hph = rad2deg(unwrap(angle(H)));
subplot(2, 1, 2);
phase_resp = plot(f, Hph);
xlabel("Frequency MHz");
ylabel("Phase (deg)");
xlim([0 xmax]);
title("Phase Response");
grid on;
end


