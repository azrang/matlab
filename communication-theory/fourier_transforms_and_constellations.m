%% Comm Theory HW 1
close all

%%M-ary QAM


k = 2;  %4QAM
M=2^k;
mat_2 = linspace(-.5, .5, sqrt(M) )+linspace(-.5, .5, sqrt(M)).'*j;
vec_2 = reshape(mat_2, [], 1);    %each point is a different column
Eb_2 = (1/k)*(1/M)*sum(abs(vec_2).^2)
%4 points, each different columns
bits2 = k/2;
k = 4;  %16qam
M=2^k;
mat_4 = linspace(-1.5, 1.5, sqrt(M))+linspace(-1.5, 1.5, sqrt(M)).'*j;
vec_4 = reshape(mat_4, [], 1);
Eb_4 = (1/k)*(1/M)*sum(abs(vec_4).^2)
bits4 = k/2;
k=6;  %64qam
M=2^k;
mat_6 = linspace(-2.5, 2.5, sqrt(M))+linspace(-2.5, 2.5, sqrt(M)).'*j;
vec_6 = reshape(mat_6, [], 1);   
Eb_6 = (1/k)*(1/M)*sum(abs(vec_6).^2)
bits6 = k/2;
k=8;  %256QAM
M=2^k;
mat_8 = linspace(-7.5, 7.5, sqrt(M))+linspace(-7.5, 7.5, sqrt(M)).'*j;
vec_8 = reshape(mat_8, [], 1);
Eb_8 = (1/k)*(1/2^k)*sum(abs(vec_8).^2)
bits8 = k/2;

figure;
plot([2, 4, 6, 8], [Eb_2, Eb_4, Eb_6, Eb_8], 'Marker', 'o');
xlabel('k');
ylabel('Eb');
title('Eb for k values');

figure;
plot([2, 4, 6, 8], [bits2, bits4, bits6, bits8], 'Marker', 'o');
xlabel('k');
ylabel('bits per symbol per dimension');
title('bits per dimension for k values');

syms A T;
eqn = A*T*sinc(0) == 1;
A = solve(eqn, A);  % so that peak value of Y is at 0db

fs = 40e3;
f = linspace(0, fs/2, 1e8);  %large enough vector, with enough n points 

T = 0.25;
Amp = double(subs(A, T))
W = 1;

Y = Amp*exp(j*2*pi*f*.5*T)*T.*sinc(f*T).*exp((-log(2)/2)*(f/W).^2);
Hdb = 20*log10(abs(Y));
figure;
plot(f, Hdb);
xlim([0 10]); %what the x-axis shows
ylim([-200 0]);  %what the y-axis shows
xlabel("Frequency kHz");
ylabel("|Y| (dB)");
title("Output Magnitude Spectrum")

%find where mag is less than or equal to -50dB 
ind = find(Hdb<=-50, 1);
B0 = f(ind)

%use qfunc to compute and graph y(t)
stdev = sqrt(log(2))/(2*pi*W);
t = linspace(0, 2, 1e4);
y_t = Amp*(qfunc((t-T)/stdev)-qfunc((t)/stdev));
ind_t = find(abs(y_t)<= 0.1*max(abs(y_t)), 1);
T0 = t(ind_t);
B0*T0

figure;
plot(t, y_t);
xlabel("Time");
ylabel("y(t)");
title("Output Signal")






