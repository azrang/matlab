
% HW 1 
close all

%% generating data with DFT


fs= 20e6;
f = 6e6;
N=512;
n0 = 500; %length of signal
n= 0:n0-1;  %along length of signal
k = 0:N-1;  %taking positive frequencies, freq index, peak value in vector 
w = 2*pi*f/fs;  
x = cos(w*n); 
rect_wind = (1/n0)*rectwin(n0); %scale for energy = 1    %divide by sqrt of sum squared
x_hat1 = rect_wind.' .* x;   %perform windowing
X1 = fft(x_hat1, N);  
cheby = chebwin(n0, 30).';  %window should have length of window
cheb_wind = cheby/norm(cheby);
x_hat2 = cheb_wind .* x;
X2 = fft(x_hat2, N);

Hdb1 = 20*log10(abs(X1));
%finding peak 
peak= max(Hdb1);  
k_max = find(Hdb1>= peak)-1  %indexing starts at 1

figure;
plot(k, Hdb1);
title('rect window Mag Resp');

ylim([-60 60]);  %what the y-axis shows

Hdb2 = 20*log10(abs(X2));
figure;
plot(k, Hdb2);
title('Cheby window Mag Resp');

%to graph freq response
ind = fs / N;
f_new = -fs/2+ind:ind:fs/2;   %to get length 512, from -fs/2 to fs/2

% corrupting with noise means adding
SNR = 20; 
noise = randn(size(x))*std(x)/db2mag(SNR);
r = x + noise;
x_hat1_noise = rect_wind.' .* r;   %perform windowing
X1_noise = fft(x_hat1_noise, N);  
x_hat2_noise = cheb_wind .* r;
X2_noise = fft(x_hat2_noise, N);

Hdb1_new = 20*log10(abs(X1_noise)); 
Hdb2_new = 20*log10(abs(X2_noise)); 

figure;
plot(f_new,Hdb1_new);
title('rectangular window Mag Resp with noise');
xlabel("Frequency Hz");
ylabel("|X(f)| (dB)");

figure;
plot(f_new,Hdb2_new);
title('cheby window Mag Resp with noise');
xlabel("Frequency Hz");
ylabel("|X(f)| (dB)");

%superimposed graph
k0 = k_max(1);
Hdb1_newer = 20*log10(abs(X1_noise((k0-10):(k0+10)))); 
Hdb2_newer = 20*log10(abs(X2_noise((k0-10):(k0+10)))); 
f_analog_lower = (k0-10)*fs/N;
f_analog_upper = (k0+10)*fs/N;
f_newer = f_analog_lower:ind:f_analog_upper;
figure;
plot(f_newer,Hdb1_newer);
hold on;
plot(f_newer,Hdb2_newer);
hold off;
title('Windows Mag Resp with noise k0+-10');
xlabel("Frequency Hz");
ylabel("|X(f)| (dB)");

fs= 44.1e3;
delta_f = 2; %bin spacing in Hz
M = 100; %number of blocks
N = fs/delta_f; %number of samples length
N_dft = 2^ceil(log2(N)); %dft length rounded up  

tempo = 100;
duration = 1/(tempo/60);  %to get number of seconds per beat (note)
num_samplesperbeat = duration / (1/fs);   %number of samples in a quarter note
total_samples = 0.5*N*M;  %total amount of samples if 100 blocks
%N+(M-1)L, where N is the number of bins for DFT, M is number of blocks, L=
%N/2 overlap 
num_notes =  ceil(total_samples/ num_samplesperbeat)

signal = create_signal(num_notes, N);
noverlap = N_dft/2; %to get 50% overlap
[pxx, freq] = pwelch(signal, N, noverlap, N_dft, fs);   %the second input is window 
% vector of length of the signal
Hdb_pxx = 10*log10(abs(pxx));

figure;
subplot(2,1,1);
plot(freq, Hdb_pxx);
title('periodogram hamming window');
xlabel("Frequency Hz");
ylabel("PSD dB");
subplot(2,1,2);
plot(freq, Hdb_pxx);
xlim([300 600]);   %zoomed in
title('periodogram hamming window; zoomed in');
xlabel("Frequency Hz");
ylabel("PSD dB");

figure;
spectrogram(signal, N, noverlap, N_dft, fs);

figure;
spectrogram(signal, N, noverlap, N_dft, fs);
xlim([.3 .6])

function signal = create_signal(num_notes, index)
    tones = [392 440 587.33];
    fs= 44.1e3;
    n = 0:index-1;
    signal = zeros(1, num_notes); %allocate
    for i = 1:num_notes
        rand_freq = randsample(tones, 2); %randomly choose 2 tones
        sig_1 = cos((2*pi*rand_freq(1)/fs)*n);
        sig_2 = cos((2*pi*rand_freq(2)/fs)*n);
        quarter_note = sig_1 + sig_2;
        signal = [signal ,quarter_note];
    end
end



