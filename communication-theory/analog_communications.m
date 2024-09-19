close all;
clear all;

%% FM

SNR = 10:2:30;  % 11 SNR points in dB
N =10e5;

[env_10 , ph_10] = env_phase(10, N);
[env_12 , ph_12] = env_phase(SNR(2), N);
[env_14 , ph_14] = env_phase(SNR(3), N);
[env_16 , ph_16] = env_phase(SNR(4), N);
[env_18 , ph_18] = env_phase(SNR(5), N);
[env_20 , ph_20] = env_phase(20, N);
[env_22 , ph_22] = env_phase(SNR(7), N);
[env_24 , ph_24] = env_phase(SNR(8), N);
[env_26 , ph_26] = env_phase(SNR(9), N);
[env_28 , ph_28] = env_phase(SNR(10), N);
[env_30 , ph_30] = env_phase(30, N);

ind_10 = find(abs(ph_10) > 10*pi/180);
ind_12 = find(abs(ph_12) > 10*pi/180);
ind_14 = find(abs(ph_14) > 10*pi/180);
ind_16 = find(abs(ph_16) > 10*pi/180);
ind_18 = find(abs(ph_18) > 10*pi/180);
ind_20 = find(abs(ph_20) > 10*pi/180);
ind_22 = find(abs(ph_22) > 10*pi/180);
ind_24 = find(abs(ph_24) > 10*pi/180);
ind_26 = find(abs(ph_26) > 10*pi/180);
ind_28 = find(abs(ph_28) > 10*pi/180);
ind_30 = find(abs(ph_30) > 10*pi/180);

%for 2 dB increase
[~, xph] = env_phase(SNR+2, N);
for s = 1:size(SNR+2, 2)
    n = size(find(abs(xph(:,s)) >10 *pi/180), 1);
    Pth_2(1,s) = n/N;
end


exceeds_10 =  [size(ind_10, 1) size(ind_12, 1) size(ind_14, 1) size(ind_16, 1) size(ind_18, 1) size(ind_20, 1) size(ind_22, 1) size(ind_24, 1) size(ind_26, 1) length(ind_28) length(ind_30)];
Pth_1 = exceeds_10/N;
rate = diff(log10(Pth_1))
rate_2dBplus = diff(log10(Pth_2))

figure;  %for SNR=10
h1 = histogram(env_10, 'Normalization', 'pdf');
title('SNR=10, envelope');

figure;
h2 = histogram(ph_10,'Normalization', 'pdf');
title('SNR=10, phase');
xline(10*pi/180);

figure; 
h3 = histogram(env_20,  'Normalization', 'pdf');
title('SNR=20, envelope');

figure;
h4 = histogram(ph_20, 'Normalization', 'pdf');
title('SNR=20, phase');
xline(10*pi/180);

figure;  
h5 = histogram(env_30, 'Normalization', 'pdf');
title('SNR=30, envelope');

figure;
h6 = histogram(ph_30, 'Normalization', 'pdf');
title('SNR=30, phase');
xline(10*pi/180);

figure;  %plot of SNR in dB vs log10Pth(SNR)
plot(SNR, log10(Pth_1));
title('SNR vs probability the phase>10')


function [x_env, x_phase] = env_phase(SNR, N)

var_n = 0.5*(10.^(-SNR/10));
%noise_gauss = ((sqrt(var_n).*randn(1,N)) + (sqrt(var_n)).*randn(1,N)*j).'; %to change mean, add it. and to change var, multiply it
nI = sqrt(var_n) .* randn(N,1);
nQ = sqrt(var_n) .* randn(N,1);
noise_gauss = nI + 1j.*nQ;
A=1;
signal = noise_gauss + A;
x_env = abs(signal);
x_phase = angle(signal);
end
