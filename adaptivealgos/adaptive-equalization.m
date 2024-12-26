%% HW 3 :  Adaptive Equalization
close all;
clear;
% equalization: using the output of a decision algorithm to regenerate a
% reliable match for x[n]
%% Experiment 

N0 = 3; %channel delay between the transmitted signal and received signal. equalixer will 
%try to compensate by adjusting weights
M0 = 5;  %center tap
Mmax = 11;  %the number of taps the adaptive equalizer will have 
Niter = 20;
Niter2 = 50;
K = 10e2;
alpha = [.1 .2 .3];
noisepw = [ -30 -10];
lambda = .9;
delta = .01;
Ntrain = Mmax-M0+Niter-1; %length of training sequence
Ntrain2 = Mmax-M0+Niter2-1;
str = {'-1', '1'};
idx = randi(2,[Ntrain+M0+K 1]);
idx2 =  randi(2,[Ntrain2+M0+K 1]);
x1 = str2double(str(idx)).';  %to get iid +-1
x2 = str2double(str(idx2)).';  %to get iid +-1
w = zeros(Mmax, 1); %initial weight
w(M0) = 1;
[y11, d_alp1_Pdb1, u_alp1_Pdb1] = generate_du(noisepw(1), alpha(1), x1, Niter);
[y12, d_alp2_Pdb1, u_alp2_Pdb1] = generate_du(noisepw(1), alpha(2), x1, Niter);
[y13, d_alp3_Pdb1, u_alp3_Pdb1] = generate_du(noisepw(1), alpha(3), x1, Niter);
[y21, d_alp1_Pdb2, u_alp1_Pdb2] = generate_du(noisepw(2), alpha(1), x1, Niter);
[y22, d_alp2_Pdb2, u_alp2_Pdb2] = generate_du(noisepw(2), alpha(2), x1, Niter);
[y23, d_alp3_Pdb2, u_alp3_Pdb2] = generate_du(noisepw(2), alpha(3), x1, Niter);
[wf11, pred_error11, k11, P11] = rls(d_alp1_Pdb1, u_alp1_Pdb1, lambda, delta, w, Niter);
[wf11_qrd,pred_error11_qrd,k11_qrd,gamma11,Pch11] = invqrdrls(d_alp1_Pdb1,u_alp1_Pdb1, lambda,delta,w,Niter);
[wf12, pred_error12, k12, P12] = rls(d_alp2_Pdb1, u_alp2_Pdb1, lambda, delta, w, Niter);
[wf12_qrd,pred_error12_qrd,k12_qrd,gamma12,Pch12] = invqrdrls(d_alp2_Pdb1,u_alp2_Pdb1, lambda,delta,w,Niter);
[wf13, pred_error13, k13, P13] = rls(d_alp3_Pdb1, u_alp3_Pdb1, lambda, delta, w, Niter);
[wf13_qrd,pred_error13_qrd,k13_qrd,gamma13,Pch13] = invqrdrls(d_alp3_Pdb1,u_alp3_Pdb1, lambda,delta,w,Niter);
[wf21, pred_error21, k21, P21] = rls(d_alp1_Pdb2, u_alp1_Pdb2, lambda, delta, w, Niter);
[wf21_qrd,pred_error21_qrd,k21_qrd,gamma21,Pch21] = invqrdrls(d_alp1_Pdb2,u_alp1_Pdb2, lambda,delta,w,Niter);
[wf22, pred_error22, k22, P22] = rls(d_alp2_Pdb2, u_alp2_Pdb2, lambda, delta, w, Niter);
[wf22_qrd,pred_error22_qrd,k22_qrd,gamma22,Pch22] = invqrdrls(d_alp2_Pdb2,u_alp2_Pdb2, lambda,delta,w,Niter);
[wf23, pred_error23, k23, P23] = rls(d_alp3_Pdb2, u_alp3_Pdb2, lambda, delta, w, Niter);
[wf23_qrd,pred_error23_qrd,k23_qrd,gamma23,Pch23] = invqrdrls(d_alp3_Pdb2,u_alp3_Pdb2, lambda,delta,w,Niter);
% in all of these weight vectors, the 3rd value is greater and dominates
% all the other values. The lowest value I get in the final weight vector
% is for wf23_qrd, which is when the noise power = -10 and alpha = .10, and
% this value is .8522

%get arrays for the for loop
y_1 = [y11 y12 y13]; %for all alphas, pdb = -30
y_2 = [y21 y22 y23];
wf_1_rls = [wf11 wf12 wf13];
wf_2_rls = [wf21 wf22 wf23];
wf_2_qrd = [wf21_qrd wf22_qrd wf23_qrd];
wf_1_qrd = [wf11_qrd wf12_qrd wf13_qrd];

%SNIR THEORY
var11 = 10^(noisepw(1)/10); %linear power
var12 = 10^(noisepw(2)/10); %linear power
snirtheory_raw11 = -10*log10(4*abs(alpha(1))^2+var11);
snirtheory_raw12 = -10*log10(4*abs(alpha(2))^2+var11);
snirtheory_raw13 = -10*log10(4*abs(alpha(3))^2+var11);
snirtheory_raw21 = -10*log10(4*abs(alpha(1))^2+var12);
snirtheory_raw22 = -10*log10(4*abs(alpha(2))^2+var12);
snirtheory_raw23 = -10*log10(4*abs(alpha(3))^2+var12);
sniroptimal11 = -10*log10(var11);
sniroptimal12 = -10*log10(var12);
var11 = 10^(noisepw(1)/10); %linear power
sniroptimal11 = -10*log10(var11);
snirtheory_raw11 = -10*log10(4*abs(alpha(1))^2+var11);


% for pdb = -30
for i= 1:3
% Compute x_est
xest11 = zeros(1, K-Mmax+1);
xest11_qrd = zeros(1, K-Mmax+1);
for n = 1:K-Mmax+1
    u_n = y_1(n + Ntrain + M0:n + Ntrain + M0 - 1 + Mmax, i);
    xest11(n) = wf_1_rls(:,i)' * u_n;
    xest11_qrd(n)= wf_1_qrd(:,i)'*u_n;
end
x12 = x1(Ntrain+M0:Ntrain+M0+size(xest11, 2)-1);

% Raw SNIR
snir_raw11 = -10*log10(mean(abs(y_1(4:end,i)-x1(1:end-3)).^2));   
%to line up
% Equalized SNIR
equalized_error_rls11 = xest11 - x12.';
snir_eq11_rls = -10*log10(mean(abs(equalized_error_rls11).^2));
equalized_error_qrdrls11 = xest11_qrd - x12.';
snir_eq11_qrdrls = -10*log10(mean(abs(equalized_error_qrdrls11).^2));

% Display results
disp(['for Noise power = -30 and alpha = ', num2str(alpha(i))])
disp(['SNIR theoretical RAW: ', num2str(snirtheory_raw11)])
disp(['SNIR optimal: ',num2str(sniroptimal11)])
disp(['SNIR Raw: ', num2str(snir_raw11)]);
disp(['SNIR rls Equalized: ', num2str(snir_eq11_rls)]);
disp(['SNIR qrd rls Equalized: ', num2str(snir_eq11_qrdrls)]);
end

%for pdb = -10
for i= 1:3
% Compute x_est
xest11 = zeros(1, K-Mmax+1);
xest11_qrd = zeros(1, K-Mmax+1);
for n = 1:K-Mmax+1
    u_n = y_2(n + Ntrain + M0:n + Ntrain + M0 - 1 + Mmax, i);
    xest11(n) = wf_2_rls(:,i)' * u_n;
    xest11_qrd(n)= wf_2_qrd(:,i)'*u_n;
end
x12 = x1(Ntrain+M0:Ntrain+M0+size(xest11, 2)-1);

% Raw SNIR
snir_raw11 = -10*log10(mean(abs(y_2(4:end,i)-x1(1:end-3)).^2));   
%to line up
% Equalized SNIR
equalized_error_rls11 = xest11 - x12.';
snir_eq11_rls = -10*log10(mean(abs(equalized_error_rls11).^2));
equalized_error_qrdrls11 = xest11_qrd - x12.';
snir_eq11_qrdrls = -10*log10(mean(abs(equalized_error_qrdrls11).^2));

% Display results
disp(['for Noise power = -10 and alpha = ', num2str(alpha(i))])
disp(['SNIR theoretical RAW: ', num2str(snirtheory_raw11)])
disp(['SNIR optimal: ',num2str(sniroptimal11)])
disp(['SNIR Raw: ', num2str(snir_raw11)]);
disp(['SNIR rls Equalized: ', num2str(snir_eq11_rls)]);
disp(['SNIR qrd rls Equalized: ', num2str(snir_eq11_qrdrls)]);
end

%in general, the results are consistent in that the SNIR optimal > 
%SNIR theoretical raw > SNIR raw > SNIR equalized

%for Niter = 50. Yes in general, I do get better results. 
[y11, d_alp1_Pdb1, u_alp1_Pdb1] = generate_du(noisepw(1), alpha(1), x2, Niter2);
[y12, d_alp2_Pdb1, u_alp2_Pdb1] = generate_du(noisepw(1), alpha(2), x2, Niter2);
[y13, d_alp3_Pdb1, u_alp3_Pdb1] = generate_du(noisepw(1), alpha(3), x2, Niter2);
[y21, d_alp1_Pdb2, u_alp1_Pdb2] = generate_du(noisepw(2), alpha(1), x2, Niter2);
[y22, d_alp2_Pdb2, u_alp2_Pdb2] = generate_du(noisepw(2), alpha(2), x2, Niter2);
[y23, d_alp3_Pdb2, u_alp3_Pdb2] = generate_du(noisepw(2), alpha(3), x2, Niter2);
[wf112, pred_error112, k112, P112] = rls(d_alp1_Pdb1, u_alp1_Pdb1, lambda, delta, w, Niter2);
[wf11_qrd2,pred_error11_qrd2,k11_qrd2,gamma112,Pch112] = invqrdrls(d_alp1_Pdb1,u_alp1_Pdb1, lambda,delta,w,Niter2);
[wf122, pred_error122, k122, P122] = rls(d_alp2_Pdb1, u_alp2_Pdb1, lambda, delta, w, Niter2);
[wf12_qrd2,pred_error12_qrd2,k12_qrd2,gamma122,Pch122] = invqrdrls(d_alp2_Pdb1,u_alp2_Pdb1, lambda,delta,w,Niter2);
[wf132, pred_error132, k132, P132] = rls(d_alp3_Pdb1, u_alp3_Pdb1, lambda, delta, w, Niter2);
[wf13_qrd2,pred_error13_qrd2,k13_qrd2,gamma132,Pch132] = invqrdrls(d_alp3_Pdb1,u_alp3_Pdb1, lambda,delta,w,Niter2);
[wf212, pred_error212, k212, P212] = rls(d_alp1_Pdb2, u_alp1_Pdb2, lambda, delta, w, Niter2);
[wf21_qrd2,pred_error21_qrd2,k21_qrd2,gamma212,Pch212] = invqrdrls(d_alp1_Pdb2,u_alp1_Pdb2, lambda,delta,w,Niter2);
[wf222, pred_error222, k222, P222] = rls(d_alp2_Pdb2, u_alp2_Pdb2, lambda, delta, w, Niter2);
[wf22_qrd2,pred_error22_qrd2,k22_qrd2,gamma222,Pch222] = invqrdrls(d_alp2_Pdb2,u_alp2_Pdb2, lambda,delta,w,Niter2);
[wf232, pred_error232, k232, P232] = rls(d_alp3_Pdb2, u_alp3_Pdb2, lambda, delta, w, Niter2);
[wf23_qrd2,pred_error23_qrd2,k23_qrd2,gamma232,Pch232] = invqrdrls(d_alp3_Pdb2,u_alp3_Pdb2, lambda,delta,w,Niter2);

var21 = 10^(noisepw(1)/10); %linear power
var22 = 10^(noisepw(2)/10); %linear power
snirtheory_raw112 = -10*log10(4*abs(alpha(1))^2+var21);
snirtheory_raw122 = -10*log10(4*abs(alpha(2))^2+var21);
snirtheory_raw132 = -10*log10(4*abs(alpha(3))^2+var21);
snirtheory_raw212 = -10*log10(4*abs(alpha(1))^2+var22);
snirtheory_raw222 = -10*log10(4*abs(alpha(2))^2+var22);
snirtheory_raw232 = -10*log10(4*abs(alpha(3))^2+var22);
sniroptimal112 = -10*log10(var21);
sniroptimal122 = -10*log10(var22);

function [y1, d, u] = generate_du(noisepw, alpha, x1, Niter)
N0 = 3;
K = 10e2;
M0 = 5;  
Mmax = 11;  
Ntrain = Mmax-M0+Niter-1;
var1 = 10^(noisepw/10); %linear power
% var1 = 10^(noisepw/10); %linear power
v1 = sqrt(var1)/sqrt(2)*(randn(Ntrain+M0+K,1));
h1 = [zeros(1, N0-1), alpha, 1, -alpha];
y1 = filter(h1, 1, x1) + v1;
u = zeros(Mmax, Niter);
for k = 1:Niter
u(:,k) = flip(y1(k:Mmax+k-1)); %columns contain the input vectors
d(k) = x1(Mmax-M0+k-1);  %target signal values for training
end
end

function [w,pred_error,k,gamma,P]= invqrdrls(d, U,lambda,delta,w,Niter)
% gamma: sequence of conversion factors, Pch: final value of the cholesky
% factor of the P matrix
P = eye(size(U, 1))/delta;
k = zeros(size(U, 1),Niter);
pred_error=zeros(length(d),1);
gamma=zeros(length(d),1);
P = chol(P, 'lower');
% prearray= zeros(size(U, 1)+1, size(U, 1)+1);
for i=1:Niter-1
    prearray = [1 U(:,i)'*P/sqrt(lambda) ; zeros(size(U, 1), 1) P/sqrt(lambda)];
    postarray = (qr(prearray'))';
    a = diag(postarray); a(a<0) = -1; a(a>0) = 1; D = diag(a);
    postarray = postarray*D; 
    P = postarray(2:end, 2:end);
    k(:,i) = (postarray(1,1)^-1)*postarray(2:end,1);
    pred_error(i) = d(i)-(w(:,:)')*U(:,i);
    w(:,:) = w(:,:)+k(:,i)*conj(pred_error(i));
    gamma(i) = (postarray(1,1))^(-2);
end
end

function [w,pred_error,k,P1]= rls(d, U, lambda,delta,w,Niter)
%w is the MxNiter matrix where each column are the w tap-weight vectors for Niter
%xi is the error, k is the kalman gain vector
P1 = eye(size(U, 1))/delta;
s = zeros(size(U, 1),Niter);
k = zeros(size(U, 1),Niter);
pred_error=zeros(length(d),1);
for i=1:Niter-1
    s(:,i) = P1(:,:)*U(:,i);
    k(:,i) = s(:,i)/(lambda+(U(:,i)')*s(:,i));
    pred_error(i) = d(i)-(w(:,:)')*U(:,i);
    w(:,:) = w(:,:)+k(:,i)*conj(pred_error(i));
    P1(:,:) = (1/lambda)*P1(:,:)-((1/lambda)*k(:,i)*(s(:,i)'));
end
end


