clear all;
close all;

%% HW 3, MULTIRATE

%produce filter coefficients for daubechies wavelets order N
%generating db1 filters 
N1 = 1; 
wname = ['db', int2str(N1)];
[h0, h1, f0, f1] = wfilters(wname)

w = linspace(0, pi, 10e4);
h0_w = h0(1) + h0(2)*exp(-1*1i*w);
h1_w = h1(1) + h1(2)*exp(-1*1i*w);

figure;
plot(w, abs(h0_w));
hold on;
plot(w, abs(h1_w));
hold off;
legend({'H0(w)', 'H1(w)'});
title('magnitude (linear) of HAAR filters');
xlabel('freq (rad)');
ylabel('linear mag');

Pw_haar = abs(h0_w).^2 + abs(h1_w).^2;
variation = max(Pw_haar) - min(Pw_haar)   %=0, variation is very small
power = mean(Pw_haar)

%generating db5 filters
N5 = 5; 
wname = ['db', int2str(N5)];
[h0_5, h1_5, f0_5, f1_5] = wfilters(wname);
%length(dbN)= 2N

E00 = [h0_5(1) h0_5(3) h0_5(5) h0_5(7) h0_5(9)];  %evens
E01 = [h0_5(2) h0_5(4) h0_5(6) h0_5(8) h0_5(10)];  %odds
E10 = [h1_5(1) h1_5(3) h1_5(5) h1_5(7) h1_5(9)];
E11 = [h1_5(2) h1_5(4) h1_5(6) h1_5(8) h1_5(10)];

polyphase_Ez = {E00 E01; E10 E11};      %need delays of each coeff
celldisp(polyphase_Ez)

%convolution since vectors of polynomials
para_11 = conv(fliplr(E00), E00)+conv(fliplr(E10), E10);
para_12 = conv(fliplr(E00), E01)+conv(fliplr(E10), E11);
para_21 = conv(fliplr(E01), E00)+conv(fliplr(E11), E10);
para_22 = conv(fliplr(E01), E01)+conv(fliplr(E11), E11);
paraunitary_check = {para_11 para_12; para_21 para_22};
celldisp(paraunitary_check)
%there is the element 1 on the diagonal. On the off diagonal, the elements
%are almost 0 (raised to the the 10^-16) due to a small numerical error. 
%THerefore, it satisfies the condition for paraunitary

n = 10e4;
h0_5w = h0_5(1)+h0_5(2)*exp(-1*1i*w)+h0_5(3)*exp(-2*1i*w)+h0_5(4)*exp(-3*1i*w)+h0_5(5)*exp(-4*1i*w)+h0_5(6)*exp(-5*1i*w)+h0_5(7)*exp(-6*1i*w)+h0_5(8)*exp(-7*1i*w)+h0_5(9)*exp(-8*1i*w)+h0_5(10)*exp(-9*1i*w);
h1_5w = h1_5(1)+h1_5(2)*exp(-1*1i*w)+h1_5(3)*exp(-2*1i*w)+h1_5(4)*exp(-3*1i*w)+h1_5(5)*exp(-4*1i*w)+h1_5(6)*exp(-5*1i*w)+h1_5(7)*exp(-6*1i*w)+h1_5(8)*exp(-7*1i*w)+h1_5(9)*exp(-8*1i*w)+h1_5(10)*exp(-9*1i*w);

%4 stages for bank of filters
H1_z = freqz(h1_5, 1, w);
H1_z2 = freqz(h1_5, 1, 2*w);   %z^2 means 2*w
H0_z = freqz(h0_5, 1, w);
H0_z2 = freqz(h0_5, 1, 2*w);
H0_z4 = freqz(h0_5, 1, 4*w);
H1_z4 = freqz(h1_5, 1, 4*w);

figure;
plot(w, abs(h0_5w));
hold on;
plot(w, abs(h1_5w));
hold off;
legend({'H0_5(w)', 'H1_5(w)'});
title('magnitude (linear) of HAAR filters');
xlabel('freq (rad)');
ylabel('linear mag');

figure;
plot(w, abs(H1_z));
hold on;
plot(w,abs(H0_z.*H1_z2) );
plot(w, abs(H0_z.*H0_z2.*H1_z4));
plot(w, abs(H0_z.*H0_z2.*H0_z4));
legend({'stage 1', 'stage 2', 'stage 3', 'stage 4'});
title('Bank of filters');

%power complementary property
P1 = (1/2)*abs(H1_z).^2
P2 = (1/4)*abs(H0_z.*H1_z2).^2;
P3 = (1/8)*abs(H0_z.*H0_z2.*H1_z4).^2;
P4 = (1/8)*abs(H0_z.*H0_z2.*H0_z4).^2;
P_w = P1+P2+P3+P4;
P_diff = max(P_w) - min(P_w)
%this number is very small so P is constant throughout w
P_nominal = mean(P_w)
% Mk summation is 1

