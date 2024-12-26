%% Signals Random HW
close all

%% Heavy Tail Distribution

N1 = 1e6; %number of samples 
v=5; %number of dof
alp = 0.544;

%generating iid samples of N(0,1)
normal_dist = randn(N1, 1);
%compute amount of time abs value below 1
abs_norm = sum(-1<normal_dist & normal_dist<1);
frac_norm = abs_norm(:,1)/N1
figure; 
plot(normal_dist);
yline(1, '--')
yline(-1, '--')
title("normal distribution")
%reorganizing data and taking mean among 10 segments
seg_norm = reshape(normal_dist, [100000, 10]);
mean_norm = mean(seg_norm)

%generating Cauchy 
U = rand(N1, 1); %uniform samples
cauchy = alp*tan(pi*U);
abs_cauchy = sum(-1<cauchy & cauchy<1);
frac_cauchy = abs_cauchy(:,1)/N1
figure;
plot(cauchy);
yline(1, '--')
yline(-1, '--')
title("cauchy distribution")
seg_cauchy = reshape(cauchy, [100000, 10]);
mean_cauchy = mean(seg_cauchy)
% it is not reasonable because the averages found are all much greater than
% 0

% generating student's t dist, with variance =1
students = sqrt(3/5)*trnd(v, N1, 1);
abs_students = sum(-1<students & students<1);
frac_students = abs_students(:,1)/N1
figure;
plot(students);
yline(1, '--')
yline(-1, '--')
title("students' distribution")
seg_students = reshape(students, [100000, 10]);
mean_students = mean(seg_students)

%% ARMA and AR Models

%compute for transfer function
b1= [1 0.4 0.2];
a1= [1 -1.6 0.8];
[z, p, k]= tf2zp(b1,a1)
figure
zplane(z,p);
title('WSS Pole Zero Plot');
% a system is minimum phase if all of its poles and zeroes are inside unit
% circle. they are, and so it is minimum phase.

% generate random data and estimate the correlation 
N2 = 1e4;
var_v= 2;
v_n = sqrt(var_v)*randn(N2, 1);  %have tp multiply by standard dev, make sure to use randn to get mean 0
%apply the filter to the white noise 
y1 = filter(b1, a1, v_n); 
for m = 0:6
        y2= y1(1+m:N2);
        y3= y1(1:(N2-m));  %have to shift lag to get to fit
        r_m(m+1) = (1/N2)*dot(y2, y3)
end

figure;
stem(0:6,r_m);
hold on;
stem(-6:0, fliplr(r_m));
hold off;

%arrange correlation in toeplitx matrix
R_toeplitz1 = toeplitz(r_m);
display(R_toeplitz1)

%compute eigenvalues of R
[eigvec, eigval0] = eig(R_toeplitz1);
[eigval, idx] = sort(diag(eigval0), 'descend');
%eigenvalues in sorted order
eigvec = eigvec(:, idx);
check_psd = eigval>0 %since these values logical 1, 
%R is positive definite

%this is wrong, need it to be 7x7
%need to put in first col, first row
%y1_toeplitz= toeplitz(0:6, y1);
M= 7;   
y1_toeplitz = toeplitz(y1(M:-1:1), y1(M:N2));
R_toeplitz2= (1/(N2-M+1))*y1_toeplitz*conj(y1_toeplitz).';

R_toep_diff = norm(R_toeplitz1-R_toeplitz2) %these are not exactly the same

%% studying the PSD

[s_est,w] =pwelch(y1,hamming(512),256,512);  %s_est is psd and w is freq vector
figure;
plot(w,s_est)
%find frequency at peak
[maxvals_est,idx] = max(s_est);
max_freq = w(idx)
pole_angle1 = atan(.4/.8)
% these values are close to the max frequency, and one is the negative
% version. 
pole_angle2= atan(-.4/.8)

%% examining AR modeling

[a, varv]= aryule(y1, 4) %compute r[m], AR(4) model is found
% this is equivalent to making fourth order linear system of eqs, vector a
% has coefficient. variance is varv

display(['The variation of innovations signal is very close to the variance ' ...
    'of prediction error of AR model'])
% pass v through AR innovations filter
y4 =filter(1, a, v_n);

for m = 0:6
        y5= y4(1+m:N2);
        y6= y4(1:(N2-m));  %have to shift lag to get to fit
        r_m2(m+1) = (1/N2)*dot(y5, y6);
end

figure;
stem(0:6,r_m);
hold on;
stem(-6:0, fliplr(r_m));
hold off;
title('correlation from AR matches for abs(m)<4')

display(['these the correclation for the original signal matches with this one for ' ...
    'absolute(m)<4 as seen in the figure'])

%superimpose stem plot 
figure;
stem(y1(1:100));
hold on;
stem(y4(1:100));
hold off;
title('y1 vs y4')

display('these do not match exactly as seen in graph')



