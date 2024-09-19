clear all;
close all;

%% HW 6 ISI

%designing the Matched Filter MF

beta = 0.3;  %rolloff
num_symb = 4;  %number of symbol
num_sps = 8; %sample per symbol period
h = rcosdesign(beta,num_symb, num_sps);
n = 0:1:length(h)-1;
g = conv(h, fliplr(h));
gprime = abs(g);
ISI = 20*log10(1/(2*sum(gprime(1:32))));
W = 325e3;   %bandwidth
f = linspace(0, W, 1e4);
G_f = freqz(g, 1, f, 2*W);   %where 2W is fs


figure;
stem(n, h ,'filled');
title('Impulse response')
xlabel('n')

figure;
stem([n 33:64], g, 'filled');
title('MF: g(n) = h[n] * h[N-n]')
xlabel('n')

figure;
plot(f, abs(G_f));
title('|G(f)|');

%% QPSKKKKK
len= 1e2; %the number of bits in data stream
%infdB_output = zeros(10, 2)
for j = 1:10
    [BER1, mean_sq_diff1] = QPSK_ber_msqdiff (inf, len,h);
    infdB_output(j,:) = [BER1 mean_sq_diff1] 
    [BER2, mean_sq_diff2] = QPSK_ber_msqdiff (-5, len,h);
    negfivedB_output(j,:) = [BER2 mean_sq_diff2]
    [BER3, mean_sq_diff3] = QPSK_ber_msqdiff (5, len,h);
    fivedB_output(j,:) = [BER3 mean_sq_diff3]
    [BER4, mean_sq_diff4] = QPSK_ber_msqdiff (10, len,h);
    tendB_output(j,:) = [BER4 mean_sq_diff4]
end

mean_infdB = mean(infdB_output, 1) 
mean_negfivedB = mean(negfivedB_output, 1) 
mean_fivedB = mean(fivedB_output, 1) 
mean_tendB = mean(tendB_output, 1) 

var_infdB = var(infdB_output, 1)
var_negfivedB = var(negfivedB_output, 1)
var_fivedB = var(fivedB_output, 1)
var_tendB = var(tendB_output, 1)

%the BER values graphed to the SNR are close to ones found on the graph
%given in the notes

function [BER, mean_sq_diff] = QPSK_ber_msqdiff (SNR, len, h)
data_in = randi([0,1], len/2, 2);  %random 1e6 length data 
Rs = 8;
n0 = 32;
%use gray coding
QPSK_symb = zeros(len/2, 1);
for k = 1:len/2
    if data_in(k, 1)==0
        if data_in(k,2)==0
            QPSK_symb(k,1)=1+1j;
        elseif data_in(k,2)==1
            QPSK_symb(k,1)=-1+1j;
        end
    elseif data_in(k,1)==1
        if data_in(k,2)==0
            QPSK_symb(k,1)=1-1j;
        elseif data_in(k,2)==1
            QPSK_symb(k,1)=-1-1j;
        end
    end
end

disp(QPSK_symb)

%add 7 zeroes in between each symbol, upsampling
upsamp_symb = reshape([QPSK_symb.';zeros(Rs-1,numel(QPSK_symb))],1,[]).';

%now need to get pulses, use conv instead of filter
p = conv(h,upsamp_symb);
hflip= fliplr(h);

%add white gaussian noise in channel
var_n = 0.5*(10.^(-SNR/10));
nI = sqrt(var_n) .* randn(length(p),1);
nQ = sqrt(var_n) .* randn(length(p),1);
noise_gauss = nI + 1j.*nQ;

received_signal = p + noise_gauss; 

%apply Matched Filter: take conjugate flip
MF_output = conv(received_signal, hflip); %why is the MATCHED FILTER these inputs

%sample starting at 32nd index and take at every 8 samples, since symbol
%rate is 8 samples/symbol
output_sampled = MF_output(n0:Rs:end-n0-1);

decoded_bits =zeros(len/2, 2); 
for k=1:length(output_sampled)
    if(real(output_sampled(k))>0)
        if(imag(output_sampled(k))>0)  %in first quadrant
            decoded_bits(k,1 ) = 0;
            decoded_bits(k,2)= 0;
        elseif(imag(output_sampled(k))<0)  %in 4th quadrant
            decoded_bits(k,1) = 1;
            decoded_bits(k,2) = 0;
        end
    elseif(real(output_sampled(k))<0)
         if(imag(output_sampled(k))>0) 
            decoded_bits(k,1 ) = 0;
            decoded_bits(k,2)= 1;
        elseif(imag(output_sampled(k))<0)  %in 4th quadrant
            decoded_bits(k,1) = 1;
            decoded_bits(k,2) = 1;
        end           
    end
end

%compute BER
error_bits = nnz(decoded_bits - data_in);  %returns number of values where nonzero
BER = error_bits / len

%compute mean square difference
mean_sq_diff = mean(abs(QPSK_symb - output_sampled).^2)

end





