%% HW 1: Initial Experimentation: Azra Rangwala
clc; close all;

N_iter = 1000;
M1 = 2;  %mth order linear prediction model, Mx1 weight vector
M2 = 4;
M3 = 10;
p1 = [.9 .8];
p2 = [.95 .8];
p3 = [.95 -.9];

%% Data Analysis

%% Task 1: sythesize the random signal
N_init1 = 44; %determined by max pole decay, each pole pair   
N0_1 = N_iter+M1-1;
x1 = generate_x1(p1, M1, N_iter);
x2 = generate_x1(p1, M2, N_iter);
x3 = generate_x1(p1, M3, N_iter);
x4 = generate_x1(p2, M1, N_iter);
x5 = generate_x1(p2, M2, N_iter);
x6 = generate_x1(p2, M3, N_iter);
x7 = generate_x1(p3, M1, N_iter);
x8 = generate_x1(p3, M2, N_iter);
x9 = generate_x1(p3, M3, N_iter);
X1 = generate_X1(x1, M1, N_iter);
X2 = generate_x1(x1, M2, N_iter);
X3 = generate_x1(x1, M3, N_iter);
X4 = generate_x1(x2, M1, N_iter);
X5 = generate_x1(x2, M2, N_iter);
X6 = generate_x1(x2, M3, N_iter);
X7 = generate_x1(x3, M1, N_iter);
X8 = generate_x1(x3, M2, N_iter);
X9 = generate_x1(x3, M3, N_iter);

%% Task 2:analyze the random signal

%1
r_l1 = autocorr(p1(1), p1(2));
r_l2 = autocorr(p2(1), p2(2));
r_l3 = autocorr(p3(1), p3(2));
R_l1= toeplitz(r_l1);
R_l2= toeplitz(r_l2);
R_l3= toeplitz(r_l3);

eig1 = eig(R_l1);
eig2 = eig(R_l2);
eig3 = eig(R_l3);
verify_eig_positive1 = eig1>0;
verify_eig_real1 = isreal(eig1);  %yes real and positive
verify_eig_positive2 = eig2>0;
verify_eig_real2 = isreal(eig2);  %yes real and positive
verify_eig_positive3 = eig3>0;
verify_eig_real3 = isreal(eig3);  %yes real and positive
figure;
stem(sort(eig1, 'descend'))
title('eig pole 1 plot')
figure;
stem(sort(eig2, 'descend'))
title('eig pole 2 plot')
figure;
stem(sort(eig3, 'descend'))
title('eig pole 3 plot')

%2
L=20;
eig_minmax1 = zeros(L+1, 2);
eig_minmax1(1,:) = eig(R_l1(1,1)); 
eig_minmax2 = zeros(L+1, 2);
eig_minmax2(1,:) = eig(R_l2(1,1)); 
eig_minmax3 = zeros(L+1, 2);
eig_minmax3(1,:) = eig(R_l3(1,1)); 
for i =2:L+1
    new_Rl1 = R_l1(1:i, 1:i);
    eig_sub1 = eig(new_Rl1);
    eig_minmax1(i,1)= min(eig_sub1);
    eig_minmax1(i,2)= max(eig_sub1);
    new_Rl2 = R_l2(1:i, 1:i);
    eig_sub2 = eig(new_Rl2);
    eig_minmax2(i,1)= min(eig_sub2);
    eig_minmax2(i,2)= max(eig_sub2);
    new_Rl3 = R_l3(1:i, 1:i);
    eig_sub3 = eig(new_Rl3);
    eig_minmax3(i,1)= min(eig_sub3);
    eig_minmax3(i,2)= max(eig_sub3);
end

%3 compute and plot PSD
w = linspace(-pi, pi, 1000);
S_z1 = psd_sz(p1(1), p1(2));
S_z2 = psd_sz(p2(1), p2(2));
S_z3 = psd_sz(p3(1), p3(2));
figure;
plot(w,  real(S_z1));
title('PSD of pole set 1');
xlabel('frequency')
ylabel('S(w)')
figure;
plot(w,  real(S_z2));
title('PSD of pole set 2');
xlabel('frequency')
ylabel('S(w)')
figure;
plot(w,  real(S_z3));
title('PSD of pole set 3');
xlabel('frequency')
ylabel('S(w)')
if_eigmaxlessthanpsd1 = eig_minmax1(:,2)<=max(S_z1);  %true
if_psdminlessthaneig1 = min(S_z1)<=eig_minmax1(:,1); %true
if_eigmaxlessthanpsd2 = eig_minmax2(:,2)<=max(S_z2);  %true
if_psdminlessthaneig2 = min(S_z2)<=eig_minmax2(:,1); %true
if_eigmaxlessthanpsd3 = eig_minmax3(:,2)<=max(S_z3);  %true
if_psdminlessthaneig3 = min(S_z3)<=eig_minmax3(:,1); %true
%4 compute an estimate of r from the data
x1new = generate_x1(p1, L, N_iter);
X1new = generate_X1(x1new, L, N_iter);
x2new = generate_x1(p2, L, N_iter);
X2new = generate_X1(x2new, L, N_iter);
x3new = generate_x1(p3, L, N_iter);
X3new = generate_X1(x3new, L, N_iter);
K = 1000; %makes it unbiased
Rl_1_estimate = (1/K).*(X1new*X1new.');
rl_1_estimate = Rl_1_estimate(1,:);
Rl_2_estimate = (1/K).*(X2new*X2new.');
rl_2_estimate = Rl_2_estimate(1,:);
Rl_3_estimate = (1/K).*(X3new*X3new.');
rl_3_estimate = Rl_3_estimate(1,:);

figure;
stem(r_l1)
hold on;
stem(rl_1_estimate)
title('r[l] and estimate for pole 1')
figure;
stem(r_l2)
hold on;
stem(rl_2_estimate)
title('r[l] and estimate for pole 2')
figure;
stem(r_l3)
hold on;
stem(rl_3_estimate)
title('r[l] and estimate for pole 3')

%5
%the appropriate choice for alpha= 1/sqrt(K). 

%6
delta_w = (2*pi)/(1000-1);
%they are pretty close
r_0_1_estimate = 1/(2*pi)*sum(S_z1*delta_w);
r_0_1_theoretical = rl_1_estimate(1);
r_0_2_estimate = 1/(2*pi)*sum(S_z2*delta_w);
r_0_2_theoretical = rl_2_estimate(1);
r_0_3_estimate = 1/(2*pi)*sum(S_z3*delta_w);
r_0_3_theoretical = rl_3_estimate(1);

%% LMS Algorithm
%verify stability
mew1 = .1/max(abs(S_z1));
mew2 = .1/max(abs(S_z2));
mew3 = .1/max(abs(S_z3));
stability_cond1 = 2/max(eig_minmax1(:,2)); %%yes it is stable 
stability_cond2 = 2/max(eig_minmax2(:,2)); %%yes it is stable 
stability_cond3 = 2/max(eig_minmax3(:,2)); %%yes it is stable 
K0 = 100;
iter_index = 5e3;
%M=2
[J1_p1, D1_p1] = learning_curve(mew1, iter_index, M1, p1, K0);
[J1_p2, D1_p2] = learning_curve(mew2, iter_index, M1, p2, K0);
[J1_p3, D1_p3] = learning_curve(mew3, iter_index, M1, p3, K0);
[J2_p1, D2_p1] = learning_curve(3e-5, iter_index, M1, p1, K0);
[J2_p2, D2_p2] = learning_curve(3e-5, iter_index, M1, p2, K0);
[J2_p3, D2_p3] = learning_curve(3e-5, iter_index, M1, p3, K0);
[J3_p1, D3_p1] = learning_curve(.5*stability_cond1, iter_index, M1, p1, K0);
[J3_p2, D3_p2] = learning_curve(.5*stability_cond2, iter_index, M1, p2, K0);
[J3_p3, D3_p3] = learning_curve(.5*stability_cond3, iter_index, M1, p3, K0);
[Jonerun_p1, Donerun_p1] = learning_curve(mew1, iter_index, M1, p1, 1);
[Jonerun_p2, Donerun_p2] = learning_curve(mew2, iter_index, M1, p2, 1);
[Jonerun_p3, Donerun_p3] = learning_curve(mew3, iter_index, M1, p3, 1);
%M=4
[J12_p1, D12_p1] = learning_curve(mew1, iter_index, M2, p1, K0);
[J12_p2, D12_p2] = learning_curve(mew2, iter_index, M2, p2, K0);
[J12_p3, D12_p3] = learning_curve(mew3, iter_index, M2, p3, K0);
[J22_p1, D22_p1] = learning_curve(3e-5, iter_index, M2, p1, K0);
[J22_p2, D22_p2] = learning_curve(3e-5, iter_index, M2, p2, K0);
[J22_p3, D22_p3] = learning_curve(3e-5, iter_index, M2, p3, K0);
[J32_p1, D32_p1] = learning_curve(.5*stability_cond1, iter_index, M2, p1, K0);
[J32_p2, D32_p2] = learning_curve(.5*stability_cond2, iter_index, M2, p2, K0);
[J32_p3, D32_p3] = learning_curve(.5*stability_cond3, iter_index, M2, p3, K0);
[Jonerun2_p1, Donerun2_p1] = learning_curve(mew1, iter_index, M2, p1, 1);
[Jonerun2_p2, Donerun2_p2] = learning_curve(mew2, iter_index, M2, p2, 1);
[Jonerun2_p3, Donerun2_p3] = learning_curve(mew3, iter_index, M2, p3, 1);
%M=10
[J13_p1, D13_p1] = learning_curve(mew1, iter_index, M3, p1, K0);
[J13_p2, D13_p2] = learning_curve(mew2, iter_index, M3, p2, K0);
[J13_p3, D13_p3] = learning_curve(mew3, iter_index, M3, p3, K0);
[J23_p1, D23_p1] = learning_curve(3e-5, iter_index, M3, p1, K0);
[J23_p2, D23_p2] = learning_curve(3e-5, iter_index, M3, p2, K0);
[J23_p3, D23_p3] = learning_curve(3e-5, iter_index, M3, p3, K0);
[J33_p1, D33_p1] = learning_curve(.5*stability_cond1, iter_index, M3, p1, K0);
[J33_p2, D33_p2] = learning_curve(.5*stability_cond2, iter_index, M3, p2, K0);
[J33_p3, D33_p3] = learning_curve(.5*stability_cond3, iter_index, M3, p3, K0);
[Jonerun3_p1, Donerun3_p1] = learning_curve(mew1, iter_index, M3, p1, 1);
[Jonerun3_p2, Donerun3_p2] = learning_curve(mew2, iter_index, M3, p2, 1);
[Jonerun3_p3, Donerun3_p3] = learning_curve(mew3, iter_index, M3, p3, 1);

%J infinity
J1_inf_p1 = sum(J1_p1(iter_index-49:end))/50; %for M=2
J12_inf_p1 = sum(J2_p1(iter_index-49:end))/50; %for M=2, decreases misadjustment
J1_inf_p2 = sum(J1_p2(iter_index-49:end))/50; %for M=2
J1_inf_p3 = sum(J1_p3(iter_index-49:end))/50; %for M=2
J2_inf_p1 = sum(J12_p1(iter_index-49:end))/50; %for M=4
J2_inf_p2 = sum(J12_p2(iter_index-49:end))/50; %for M=4
J2_inf_p3 = sum(J12_p3(iter_index-49:end))/50; %for M=4
J3_inf_p1 = sum(J13_p1(iter_index-49:end))/50; %for M=10
J3_inf_p2 = sum(J13_p2(iter_index-49:end))/50; %for M=10
J3_inf_p3 = sum(J13_p3(iter_index-49:end))/50; %for M=10




%just one run
%pole 1
figure;
plot(linspace(1,iter_index, iter_index), Jonerun_p1, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), Donerun_p1, 'red')
title('M=2, One RUN J and D for u = 1/max(S), pole 1');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter_index, iter_index), Jonerun2_p1, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), Donerun2_p1, 'red')
title('M=4, One RUN J and D for u = 1/max(S), pole 1');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter_index, iter_index), Jonerun3_p1, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), Donerun3_p1, 'red')
title('M=10, One RUN J and D for u = 1/max(S), pole 1');
legend('J', 'Jmin', 'D' )
%pole 2
figure;
plot(linspace(1,iter_index, iter_index), Jonerun_p2, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), Donerun_p2, 'red')
title('M=2, One RUN J and D for u = 1/max(S), pole 2');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter_index, iter_index), Jonerun2_p2, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), Donerun2_p2, 'red')
title('M=4, One RUN J and D for u = 1/max(S), pole 2');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter_index, iter_index), Jonerun3_p2, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), Donerun3_p2, 'red')
title('M=10, One RUN J and D for u = 1/max(S), pole 2');
legend('J', 'Jmin', 'D' )
%pole 3
figure;
plot(linspace(1,iter_index, iter_index), Jonerun_p3, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), Donerun_p3, 'red')
title('M=2, One RUN J and D for u = 1/max(S), pole 3');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter_index, iter_index), Jonerun2_p3, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), Donerun2_p3, 'red')
title('M=4, One RUN J and D for u = 1/max(S), pole 3');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter_index, iter_index), Jonerun3_p3, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), Donerun3_p3, 'red')
title('M=10, One RUN J and D for u = 1/max(S), pole 3');
legend('J', 'Jmin', 'D' )

%100 runs
%M=2
figure;
plot(linspace(1,iter_index, iter_index), J1_p1, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D1_p1, 'red')
title('M=2, Averaged J and D for u = 1/max(S), pole 1');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter_index, iter_index), J1_p2, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D1_p2, 'red')
title('M=2, Averaged J and D for u = 1/max(S), pole 2');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter_index, iter_index), J1_p3, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D1_p3, 'red')
title('M=2, Averaged J and D for u = 1/max(S), pole 3');
legend('J', 'Jmin', 'D' )

figure;
plot(linspace(1,iter_index, iter_index), J2_p1, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D2_p1, 'red')
title('M=2, Averaged J and D for u = 3e-5, pole 1')
legend('J', 'Jmin', 'D' )

figure;
plot(linspace(1,iter_index, iter_index), J3_p1, 'blue')
title('M=2, Averaged J')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D3_p1, 'red')
title('M=2, Averaged J and D for u = .5*stability, pole 1')
legend('J', 'Jmin', 'D' )

%M=4
figure;
plot(linspace(1,iter_index, iter_index), J12_p1, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D12_p1, 'red')
title('M=4, Averaged J and D for u = 1/max(S), pole 1');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter_index, iter_index), J12_p2, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D12_p2, 'red')
title('M=4, Averaged J and D for u = 1/max(S), pole 2');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter_index, iter_index), J12_p3, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D12_p3, 'red')
title('M=4, Averaged J and D for u = 1/max(S), pole 3');
legend('J', 'Jmin', 'D' )

figure;
plot(linspace(1,iter_index, iter_index), J22_p1, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D22_p1, 'red')
title('M=4, Averaged J and D for u = 3e-5, pole 1')
legend('J', 'Jmin', 'D' )

figure;
plot(linspace(1,iter_index, iter_index), J32_p1, 'blue')
title('M=2, Averaged J')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D32_p1, 'red')
title('M=4, Averaged J and D for u = .5*stability, pole 1')
legend('J', 'Jmin', 'D' )

%M=10
figure;
plot(linspace(1,iter_index, iter_index), J13_p1, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D13_p1, 'red')
title('M=10, Averaged J and D for u = 1/max(S), pole 1');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter_index, iter_index), J13_p2, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D13_p2, 'red')
title('M=10, Averaged J and D for u = 1/max(S), pole 2');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter_index, iter_index), J13_p3, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D13_p3, 'red')
title('M=10, Averaged J and D for u = 1/max(S), pole 3');
legend('J', 'Jmin', 'D' )

figure;
plot(linspace(1,iter_index, iter_index), J23_p1, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D23_p1, 'red')
title('M=10, Averaged J and D for u = 3e-5, pole 1')
legend('J', 'Jmin', 'D' )

figure;
plot(linspace(1,iter_index, iter_index), J33_p1, 'blue')
title('M=2, Averaged J')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter_index, iter_index), D33_p1, 'red')
title('M=10, Averaged J and D for u = .5*stability, pole 1')
legend('J', 'Jmin', 'D' )

% COMMENTS
% 1) lower u decreases rate of cvoergence and misadjustment
% 2) a higher M means that the reate of convergence decreases and
% misadjustment increases. As M increases, the evalue becomes larger, and
% so boundary condition of stability decreases
% 3) as pole gets closer to unit circle, the evalues increaes. Having u = 1/Smax makes 
% it smaller and it may reach instability

iter=50;
w_lms = lms_averaged_w(mew1, iter_index, M1, p1, K0); %it's not that close?
w_rls = rls_w_averaged(iter, M1, p1); %it's not that close?

%% RLS algorithm

iter = 50;
[J1_p1rls, D1_p1rls] = RLS(iter, M1, p1);
[J1_p2rls, D1_p2rls] = RLS(iter, M1, p2);
[J1_p3rls, D1_p3rls] = RLS(iter, M1, p3);
[J2_p1rls, D2_p1rls] = RLS(iter, M2, p1);    
[J2_p2rls, D2_p2rls] = RLS(iter, M2, p2);
[J2_p3rls, D2_p3rls] = RLS(iter, M2, p3);
[J3_p1rls, D3_p1rls] = RLS(iter, M3, p1);
[J3_p2rls, D3_p2rls] = RLS(iter, M3, p2);
[J3_p3rls, D3_p3rls] = RLS(iter, M3, p3);

%J infinity for rls
J1_inf_p1rls = sum(J1_p1rls(iter-49:end))/50; %for M=2
J1_inf_p2rls = sum(J1_p2rls(iter-49:end))/50; %for M=2
J1_inf_p3rls = sum(J1_p3rls(iter-49:end))/50; %for M=2
J2_inf_p1rls = sum(J2_p1rls(iter-49:end))/50; %for M=4
J2_inf_p2rls = sum(J2_p2rls(iter-49:end))/50; %for M=4
J2_inf_p3rls = sum(J2_p3rls(iter-49:end))/50; %for M=4
J3_inf_p1rls = sum(J3_p1rls(iter-49:end))/50; %for M=10
J3_inf_p2rls = sum(J3_p2rls(iter-49:end))/50; %for M=10
J3_inf_p3rls = sum(J3_p3rls(iter-49:end))/50; %for M=10

%the misadjustment for RLS is a lot higher than LMS


figure;
plot(linspace(1,iter, iter), J1_p1rls, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter, iter), D1_p1rls, 'red')
title('RLS M=2, J and D pole 1');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter, iter), J2_p1rls, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter, iter), D2_p1rls, 'red')
title('RLS M=4, J and D pole 1');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter, iter), J3_p1rls, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter, iter), D3_p1rls, 'red')
title('RLS M=10, J and D pole 1');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter, iter), J1_p2rls, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter, iter), D1_p2rls, 'red')
title('RLS M=2, J and D pole 2');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter, iter), J2_p2rls, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter, iter), D2_p2rls, 'red')
title('RLS M=4, J and D pole 2');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter, iter), J3_p2rls, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter, iter), D3_p2rls, 'red')
title('RLS M=10, J and D pole 2');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter, iter), J1_p3rls, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter, iter), D1_p3rls, 'red')
title('RLS M=2, J and D pole 3');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter, iter), J2_p3rls, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter, iter), D2_p3rls, 'red')
title('RLS M=4, J and D pole 3');
legend('J', 'Jmin', 'D' )
figure;
plot(linspace(1,iter, iter), J3_p3rls, 'blue')
xlabel('iterations')
hold on;
yline(1,'-','Jmin')
hold on;
plot(linspace(1,iter, iter), D3_p3rls, 'red')
title('RLS M=10, J and D pole 3');
legend('J', 'Jmin', 'D' )

function w = rls_w_averaged(iter, M1, p1)
K0 = 100;
lambda = 0.9;
J1 = zeros(iter,K0);
D1 = zeros(iter,K0);
for j = 1:K0
P1 = eye(M1);
x1_rls = generate_x1(p1, M1, iter);
X1_rls = generate_X1(x1_rls, M1, iter);
d1_rls = X1_rls(1, :);  %the desired signal is the current signal
u1_rls = X1_rls(2:M1+1, :);
k = zeros(M1,1);
s = zeros(M1,iter);
w = zeros(M1,1);
pred_error=zeros(M1,1);
R1 = (1/iter)*(u1_rls)*(u1_rls');  %autocorr
p_1 = (1/iter)*u1_rls*d1_rls';
w_opt = inv(R1)*p_1;
    for i=1:iter-1
        s(:,i) = P1(:,:)*u1_rls(:,i);
        k(:,i) = s(:,i)/(lambda+(u1_rls(:,i)')*s(:,i));
        pred_error(i) = d1_rls(i)-(w(:,:)')*u1_rls(:,i);
        w(:,:) = w(:,:)+k(:,i)*conj(pred_error(i));
        P1(:,:) = (1/lambda)*P1(:,:)-((1/lambda)*k(:,i)*(s(:,i)'));
        %check if positive definite
        check_pd = eig(P1(:,:))>0; %yes all positive
        J1(i, j) = abs(pred_error(i)).^2;
        D1(i, j) = norm(w(:,:)-w_opt).^2;
    end
w(:,:) = sum(w(1,:))/K0;
end
end

function [J_averaged, D_averaged] = RLS(iter, M1, p1)
K0 = 100;
lambda = 0.9;
J1 = zeros(iter,K0);
D1 = zeros(iter,K0);
for j = 1:K0
P1 = eye(M1);
x1_rls = generate_x1(p1, M1, iter);
X1_rls = generate_X1(x1_rls, M1, iter);
d1_rls = X1_rls(1, :);  %the desired signal is the current signal
u1_rls = X1_rls(2:M1+1, :);
k = zeros(M1,1);
s = zeros(M1,iter);
w = zeros(M1,1);
pred_error=zeros(M1,1);
R1 = (1/iter)*(u1_rls)*(u1_rls');  %autocorr
p_1 = (1/iter)*u1_rls*d1_rls';
w_opt = inv(R1)*p_1;
    for i=1:iter-1
        s(:,i) = P1(:,:)*u1_rls(:,i);
        k(:,i) = s(:,i)/(lambda+(u1_rls(:,i)')*s(:,i));
        pred_error(i) = d1_rls(i)-(w(:,:)')*u1_rls(:,i);
        w(:,:) = w(:,:)+k(:,i)*conj(pred_error(i));
        P1(:,:) = (1/lambda)*P1(:,:)-((1/lambda)*k(:,i)*(s(:,i)'));
        %check if positive definite
        check_pd = eig(P1(:,:))>0; %yes all positive
        J1(i, j) = abs(pred_error(i)).^2;
        D1(i, j) = norm(w(:,:)-w_opt).^2;
    end
    J_averaged = (1/K0)*sum(J1,2);
    D_averaged = (1/K0)*sum(D1,2);
end
end

function w = lms_averaged_w(mew, iter_index, M1, p1, K0)
error1 = zeros(iter_index,1);
w1 = zeros(M1,1);
J1 = zeros(iter_index,K0);
D1 = zeros(iter_index,K0);
for j = 1:K0
    x1 = generate_x1(p1, M1, iter_index);
    X1 = generate_X1(x1, M1, iter_index);
    d1 = X1(1, :);  %the desired signal is the current signal
    u1 = X1(2:M1+1, :); %the imput signal are the past values
    R1 = (1/iter_index)*(u1)*(u1');  %autocorr
    p_1 = (1/iter_index)*u1*d1';
    w_opt = inv(R1)*p_1;
        for i=1:iter_index-1
            error1(i,j) =  d1(i)-(w1(:,:)')*u1(:,i);
            w1(:,:) = w1(:,:)+mew*u1(:,i)*conj(error1(i,j));
            %MSE
            J1(i,j) = abs(error1(i,j)).^2;
            D1(i,j) = norm(w1(:,:)-w_opt).^2;
           
        end
        w(:,:) = sum(w1(1,:))/K0;
end
end

function [J_averaged, D_averaged] = learning_curve(mew, iter_index, M1, p1, K0)
error1 = zeros(iter_index,1);
w1 = zeros(M1,1);
J1 = zeros(iter_index,K0);
D1 = zeros(iter_index,K0);
for j = 1:K0
    x1 = generate_x1(p1, M1, iter_index);
    X1 = generate_X1(x1, M1, iter_index);
    d1 = X1(1, :);  %the desired signal is the current signal
    u1 = X1(2:M1+1, :); %the imput signal are the past values
    R1 = (1/iter_index)*(u1)*(u1');  %autocorr
    p_1 = (1/iter_index)*u1*d1';
    w_opt = inv(R1)*p_1;
        for i=1:iter_index-1
            error1(i,j) =  d1(i)-(w1(:,:)')*u1(:,i);
            w1(:,:) = w1(:,:)+mew*u1(:,i)*conj(error1(i,j));
            %MSE
            J1(i,j) = abs(error1(i,j)).^2;
            D1(i,j) = norm(w1(:,:)-w_opt).^2;
        end
%         w_averaged = sum(w1(:,:))/K0;
end 
    J_averaged = (1/K0)*sum(J1,2);
    D_averaged = (1/K0)*sum(D1,2);
%     fprintf('%d', w_averaged);
end


function x1 = generate_x1(p, M1, N_iter)
N_init1 = 44; %determined by max pole decay, each pole pair   
N0_1 = N_iter+M1-1;
v_M1_p1 = randn(1,N0_1+N_init1);  %for M=1 and pole=p1
a1 = [1 -p(1)-p(2) p(1)*p(2)];  %for pole=p1
x1_before = filter(1,a1,v_M1_p1);
x1 = x1_before(:,N_init1+1:end).';
end

function X1 = generate_X1(x1, M1, N_iter) 
N0_1 = N_iter+M1-1;
c1 = fliplr(x1(1:M1+1).');
r1 = x1(M1+1:N0_1);
X1 = toeplitz(c1,r1);
end

function r_l = autocorr(p1, p2)
L = 20;
K1 = p1/((1-p1^2)*(p1*(1+p2^2)-p2*(1+p1^2)));
K2 = p2/((1-p2^2)*(p2*(1+p2^2)-p1*(1+p2^2)));
r_l = zeros(L+1,1);
for m=1:L+1
    r_l(m) = K1*p1^abs(m-1)+K2*p2^abs(m-1);
end
end

function S_z = psd_sz(p1, p2)
    w = linspace(-pi, pi, 1000);
    z = exp(1i*w);
    S_z = 1./((1-p1.*z).*(1-p1./z).*(1-p2.*z).*(1-p2./z));
end
