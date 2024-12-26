%% HW 2 : Beamforming Problems
close all;

d_lambda = 0.5;
N = 100;  %number of snapshots
L = 3; %number of sources 
M = 41; % number of sensors (20 on xaxis, 20 on yaxis, 1 origin)

%initializations
% M 3-D vectors for sensor location
sensor_locations = [0 0 0];
m = [-10:-1 1:10];
for i = 1:length(m)
    sensor_locations(i+1,:) = [m(i) 0 0];
    sensor_locations(length(m)+i+1,:) = [0 m(i) 0];
end

AOA = [15 20 30; 30 40 -40];
rel_power = [0 -4 -8];
noisepow = -12;
%% 1a: generate S A V X for sensor array signal model
a_hat1 = a_hat(AOA(1, 1),AOA(2, 1));
a_hat2 = a_hat(AOA(1, 2),AOA(2, 2));
a_hat3 = a_hat(AOA(1, 3),AOA(2, 3));
s_vec1 = steering_vector(d_lambda, a_hat1, sensor_locations); %now do for all a_hat and generate S
s_vec2 = steering_vector(d_lambda, a_hat2, sensor_locations); %now do for all a_hat and generate S
s_vec3 = steering_vector(d_lambda, a_hat3, sensor_locations); %now do for all a_hat and generate S
S1 = generate_S(s_vec1, s_vec2, s_vec3);
V1 = generate_V(M,N,noisepow);
A1 = generate_A(N, rel_power);
X1 = S1*A1 + V1;
%% 1b: compute theoretical correlation matrix of each snapshot u[n]
R_U = thoer_corr_Ru(S1, rel_power, noisepow, M);
estimate_R_U = (1/N)*X1*X1';

%% 2: SVD and MUSIC / MVDR Spectra
%a 
sing_values = svd(X1);
figure;
stem(sing_values); %yes there are 3 dominant singular values 
title('singular values of X')
ratio_sing = sing_values(4)/sing_values(3);
%b
eig_R_U = eig(R_U);
figure;
stem(sort(eig_R_U,'descend')); %yes tehre are 3 dominant values
title('evals of R_U')
ratio_eigRU = eig_R_U(4)/eig_R_U(3);
%c
[U_svd, S_svd, V_svd] = svd(X1);
%the noise subsapce are the columns of U associated with the smaller singular values
Q_svd = U_svd(:, 1:L);
%projection matrix onto U noise
P_N = eye(size(Q_svd,1))-Q_svd*Q_svd';
test_noise = norm(P_N*S1); %yes around 0
%d
% Define the theta (azimuth) and phi (elevation) grid for AOA
theta_vals = (0:5:90).';   % Coarse grid for azimuth angle (0 to 90 degrees)
phi_vals = (-180:5:180).'; % Coarse grid for elevation angle (-180 to 180 degrees)
[theta_grid, phi_grid] = meshgrid(theta_vals, phi_vals);
% Reshape the grid to create 2-column vector inputs for a_hat
theta_vec = theta_grid(:);    % Column vector of all theta values
phi_vec = phi_grid(:);        % Column vector of all phi values
angles = [theta_vec, phi_vec];  % Combine into a 2-column matrix
aoa_grid = a_hat(angles(:,1), angles(:,2));
steeringvec_grid = steering_vector(d_lambda, aoa_grid, sensor_locations);
% music SPECTRUM based on SVD of data
music_1 = 1./diag((steeringvec_grid')*P_N*steeringvec_grid);
music_1_reshape = reshape(music_1, size(theta_grid));
music_markers = [ 1./((s_vec1')*P_N*s_vec1) 1./((s_vec2')*P_N*s_vec2) 1./((s_vec3')*P_N*s_vec3) ];
figure;
mesh(theta_grid,phi_grid, abs(music_1_reshape));
% Hold the plot to add markers
title('MUSIC 1')
hold on;
plot3(AOA(1,1), AOA(2,1), abs(music_markers(1)), 'o' , AOA(1,2), AOA(2,2), abs(music_markers(2)), 'o', AOA(1,3), AOA(2,3), abs(music_markers(3)), 'o');
hold off;
%mvdr SPECTRUM based on estimate of correlation
mvdr_1 = 1./diag((steeringvec_grid')*inv(estimate_R_U)*steeringvec_grid);
mvdr_1_reshape = reshape(mvdr_1, size(theta_grid));
figure;
mesh(theta_grid,phi_grid, abs(mvdr_1_reshape));
title('MVDR 1')
% Hold the plot to add markers
hold on;
mvdr_markers = [ 1./((s_vec1')*inv(estimate_R_U)*s_vec1) 1./((s_vec2')*inv(estimate_R_U)*s_vec2) 1./((s_vec3')*inv(estimate_R_U)*s_vec3) ];
plot3(AOA(1,1), AOA(2,1), abs(mvdr_markers(1)), 'o' , AOA(1,2), AOA(2,2), abs(mvdr_markers(2)), 'o', AOA(1,3), AOA(2,3), abs(mvdr_markers(3)), 'o');
hold off;
% e
twentydeg = find(theta_vals==20);
music_20deg = music_1_reshape(:, twentydeg);
figure;
plot(phi_vals, music_20deg)
title('slice 20 music')
mvdr_20deg = mvdr_1_reshape(:, twentydeg);
figure;
plot(phi_vals, mvdr_20deg)
title('slice 20 mvdr')
% f 
thirtydeg = find(phi_vals==30);
music_30deg = music_1_reshape(thirtydeg, :);
figure;
plot(theta_vals,music_30deg)
title('slice 30 music')
mvdr_30deg = mvdr_1_reshape(thirtydeg, :);
figure;
plot(theta_vals,mvdr_30deg)
title('slice 30 mvdr')
%g
%as P_N approaches identity, the theoretical minimum for MUSIC is 1.
%According to the rayleigh quotient, the minimum is the smallest eigenvalue of C
%h
music_steer1 = music_markers(1);
music_steer2 = music_markers(2);
music_steer3 = music_markers(3);
% the theoretical min is 1 so the peaks are 1000x, 500x, 200x the minimum
mvdr_steer1 = mvdr_markers(1);
mvdr_steer2 = mvdr_markers(2);
mvdr_steer3 = mvdr_markers(3);
minimum_eig_est = min(eig(estimate_R_U));
% the theoretical min of mvdr is 2.6425e-04 so the peaks are 2500x, x, 450x the min
% i
music1_min = min(abs(music_1)); % yes very close to 1
mvdr1_min = min(abs(mvdr_1)); %yes also very close to 2.6425e-04 

%% 3: Optimal Beamforming: MVDR and GSC
%MVDR
w_mvdr_source1 = (inv(estimate_R_U)*s_vec1)/(s_vec1'*inv(estimate_R_U)*s_vec1);
w_mvdr_source2 = (inv(estimate_R_U)*s_vec2)/(s_vec2'*inv(estimate_R_U)*s_vec2);
w_mvdr_source3 = (inv(estimate_R_U)*s_vec3)/(s_vec3'*inv(estimate_R_U)*s_vec3);
mvdr_sensor1 = w_mvdr_source1'*X1;  %mvdr beamformer output for sensor 1
twoDsteeringvec_grid = steering_vector(d_lambda, aoa_grid, sensor_locations);
arr_resp_sensor1 = abs(w_mvdr_source1'*twoDsteeringvec_grid)'.^2; %array response for   
mvdr_response1 = reshape(arr_resp_sensor1, size(theta_grid));
arr_resp_sensor2 = abs(w_mvdr_source2'*twoDsteeringvec_grid).^2; %array response for   
mvdr_response2 = reshape(arr_resp_sensor2, size(theta_grid));
arr_resp_sensor3 = abs(w_mvdr_source3'*twoDsteeringvec_grid).^2; %array response for   
mvdr_response3 = reshape(arr_resp_sensor3, size(theta_grid));

%GSC
g = [1; 0 ;0;];
w_gsc_source1 = S1*inv(S1'*S1)*g;
arr_resp_gsc = abs(w_gsc_source1'*twoDsteeringvec_grid)'.^2;
gsc_response = reshape(arr_resp_gsc, size(theta_grid));
gsc_20deg = gsc_response(:, twentydeg);
gsc_30deg = gsc_response(thirtydeg, :);

figure;
plot(phi_vals, 10*log10(abs(gsc_20deg)));
title('slice 20 gsc ')
figure;
plot(theta_vals, 10*log10(abs(gsc_30deg)));
title('slice 30 gsc ')

figure;
contour(theta_grid,phi_grid, 10*log10(abs(gsc_response)));
title('GSC response')
xlabel('Azimuth \theta (degrees)');
ylabel('GSC array responses');

% 2D Plot MVDR response for Source 1
figure;
contour(theta_grid,phi_grid, 10*log10(abs(mvdr_response1)));
title('MVDR response for source 1')
xlabel('Azimuth \theta (degrees)');
ylabel('mvdr array responses');
figure;
contour(theta_grid,phi_grid, 10*log10(abs(mvdr_response2)));
title('MVDR response for source 2')
xlabel('Azimuth \theta (degrees)');
ylabel('mvdr array responses');
figure;
contour(theta_grid,phi_grid, 10*log10(abs(mvdr_response3)));
title('MVDR response for source 3')
xlabel('Azimuth \theta (degrees)');
ylabel('mvdr array responses');
% here we see in the 1D plots that we have a 0dB output for the
% corresponding osensor, and nulls at the other sensors where there is the
% 0dB response. 
mvdr_sensor1_20deg = mvdr_response1(:, twentydeg);
mvdr_sensor3_20deg = mvdr_response3(:, twentydeg);
mvdr_sensor2_20deg = mvdr_response2(:, twentydeg);
figure;
plot(phi_vals, 10*log10(abs(mvdr_sensor1_20deg)));
title('slice 20 mvdr ')
hold on;
plot(phi_vals, 10*log10(abs(mvdr_sensor2_20deg)));
plot(phi_vals, 10*log10(abs(mvdr_sensor3_20deg)));
legend('sensor 1','sensor 2', 'sensor 3')
xline(30)
hold off;
mvdr_sensor1_30deg = mvdr_response1(thirtydeg, :);
mvdr_sensor2_30deg = mvdr_response2(thirtydeg, :);
mvdr_sensor3_30deg = mvdr_response3(thirtydeg, :);
figure;
plot(theta_vals, 10*log10(abs(mvdr_sensor1_30deg)));
hold on;
plot(theta_vals, 10*log10(abs(mvdr_sensor2_30deg)));
plot(theta_vals, 10*log10(abs(mvdr_sensor3_30deg)));
xline(15)
title('slice 30 mvdr ')
legend('sensor 1','sensor 2', 'sensor 3')
hold off;
source1_theta15_phi30 = [10*log10(abs(mvdr_response1(find(phi_vals==30), find(theta_vals==15)))); 10*log10(abs(mvdr_response2(find(phi_vals==30), find(theta_vals==15)))); 10*log10(abs(mvdr_response3(find(phi_vals==30), find(theta_vals==15))))];
source2_theta20_phi40 = [10*log10(abs(mvdr_response1(find(phi_vals==40), find(theta_vals==20)))); 10*log10(abs(mvdr_response2(find(phi_vals==40), find(theta_vals==20)))); 10*log10(abs(mvdr_response3(find(phi_vals==40), find(theta_vals==20))))];
source3_theta30_phineg40 = [10*log10(abs(mvdr_response1(find(phi_vals==-40), find(theta_vals==30)))); 10*log10(abs(mvdr_response2(find(phi_vals==-40), find(theta_vals==30)))); 10*log10(abs(mvdr_response3(find(phi_vals==-40), find(theta_vals==30))))];
sources = {'@source 1'; '@source 2';'@source 3'};
Null_tables = table(sources, source1_theta15_phi30, source2_theta20_phi40, source3_theta30_phineg40)

%% 4: Adaptive Beamforming: doing adaptive MVDR and GSC using RLS and LMS
lambda = 0.9;
mew = 0.3;
[lms_gsc1, lmsgscvec1, lmsgscout1] = lms(X1, S1, w_gsc_source1, mew, N);
ave_lmsgscvec1 = mean(lmsgscvec1(:,91:100), 2);
[lms_mvdr1, lmsmvdrvec1, lmsmvdrout1] = lms(X1, S1, w_mvdr_source1, mew, N);
ave_lmsmvdrvec1 = mean(lmsmvdrvec1(:,91:100), 2);
[rls_gsc1,rlsgscvec1, rlsgscout1] = rls(X1, S1, w_gsc_source1 , lambda, N);
ave_rlsgscvec1 = mean(rlsgscvec1(:,91:100), 2);
[rls_mvdr1,rlsmvdrvec1, rlsmvdrout1] = rls(X1, S1, w_mvdr_source1 , lambda, N);
ave_rlsmvdrvec1 = mean(rlsmvdrvec1(:,91:100), 2);
arr_resp_lmsgsc1 = abs(lmsgscvec1'*twoDsteeringvec_grid).^2;
%the average vector should work better than the one at one instant, since
%it has gone through multiple trials 

figure;
plot(1:N, lms_gsc1, 'b');
title('Learning Curve - LMS for GSC Algorithm steered to source 1');
xlabel('Iteration');
ylabel('Squared Error |e(i)|^2');
figure;
plot(1:N, lms_mvdr1, 'b');
title('Learning Curve - LMS for MVDR Algorithm steered to source 1');
xlabel('Iteration');
ylabel('Squared Error |e(i)|^2');

figure;
plot(1:N, rls_gsc1, 'b');
title('Learning Curve - RLS for GSC Algorithm steered to source 1');
xlabel('Iteration');
ylabel('Squared Error |e(i)|^2');
figure;
plot(1:N, rls_mvdr1, 'b');
title('Learning Curve - RLS for MVDR Algorithm steered to source 1');
xlabel('Iteration');
ylabel('Squared Error |e(i)|^2');

%% 5: Variations
%a : Negative SNR
noisepow2 = 10;
V2 = generate_V(M,N,noisepow2);
X2 = S1*A1 + V2;
R_U_2 = thoer_corr_Ru(S1, rel_power, noisepow2, M);
estimate_R_U_2 = (1/N)*X2*X2';
sing_values_2 = svd(X2);
figure;
stem(sing_values_2); %yes there are 3 dominant singular values 
title('singular values of X with negative SNR')
ratio_sing2 = sing_values_2(4)/sing_values_2(3);
eig_R_U_2 = eig(R_U_2);
[U_svd2, S_svd2, V_svd2] = svd(X2);
%the noise subsapce are the columns of U associated with the smaller singular values
Q_svd2 = U_svd2(:, 1:L);
%projection matrix onto U noise
P_N2 = eye(size(Q_svd2,1))-Q_svd2*Q_svd2';
test_noise2 = norm(P_N2*S1);  %this is higher than the other value with positive SNR
music_2 = 1./diag((steeringvec_grid')*P_N2*steeringvec_grid);
mvdr_2 = 1./diag((steeringvec_grid')*inv(estimate_R_U_2)*steeringvec_grid);
w_mvdr_source2 = (inv(estimate_R_U_2)*s_vec1)/(s_vec1'*inv(estimate_R_U_2)*s_vec1);
mvdr_sensor2 = w_mvdr_source2'*X2;  %mvdr beamformer output for sensor 1
arr_resp_sensor2 = abs(w_mvdr_source2'*twoDsteeringvec_grid)'.^2;
[lms_gsc2, lmsgscvec2, lmsgscout2] = lms(X2, S1, w_gsc_source1, mew, N);
ave_lmsgscvec2 = mean(lmsgscvec2(:,91:100), 2);
[lms_mvdr2, lmsmvdrvec2, lmsmvdrout2] = lms(X2, S1, w_mvdr_source2, mew, N);
ave_lmsmvdrvec2 = mean(lmsmvdrvec2(:,91:100), 2);
[rls_gsc2,rlsgscvec2, rlsgscout2] = rls(X2, S1, w_gsc_source1 , lambda, N);
ave_rlsgscvec2 = mean(rlsgscvec2(:,91:100), 2);
[rls_mvdr2,rlsmvdrvec2, rlsmvdrout2] = rls(X2, S1, w_mvdr_source2 , lambda, N);
ave_rlsmvdrvec2 = mean(rlsmvdrvec2(:,91:100), 2);
arr_resp_lmsgsc2 = abs(lmsgscvec2'*twoDsteeringvec_grid).^2;
%in general, the music and mvdr spectra are less pronounced from the noise
%floor in the graphs. the eigenvalues do not dominate as much. LMS and RLS
%do not end up converging reliably

%b: fewer snapshots 
N2 = 25;
N3 = 10;
V3 = generate_V(M,N2,noisepow);
A3 = generate_A(N2, rel_power);
X3 = S1*A3 + V3;
estimate_R_U_3 = (1/N2)*X3*X3';
sing_values_3 = svd(X3);
figure;
stem(sing_values_3); %yes there are 3 dominant singular values 
title('singular values of X with N=25')
[U_svd3, S_svd3, V_svd3] = svd(X3);
%the noise subsapce are the columns of U associated with the smaller singular values
Q_svd3 = U_svd3(:, 1:L);
%projection matrix onto U noise
P_N3 = eye(size(Q_svd3,1))-Q_svd3*Q_svd3';
test_noise3 = norm(P_N3*S1);  %this is higher than the other value with positive SNR
music_3 = 1./diag((steeringvec_grid')*P_N3*steeringvec_grid);
mvdr_3 = 1./diag((steeringvec_grid')*inv(estimate_R_U_3)*steeringvec_grid);
w_mvdr_source3 = (inv(estimate_R_U_3)*s_vec1)/(s_vec1'*inv(estimate_R_U_3)*s_vec1);
mvdr_sensor3 = w_mvdr_source3'*X3;  %mvdr beamformer output for sensor 1
arr_resp_sensor3 = abs(w_mvdr_source3'*twoDsteeringvec_grid)'.^2;
[lms_gsc3, lmsgscvec3, lmsgscout3] = lms(X3, S1, w_gsc_source1, mew, N2);
ave_lmsgscvec3 = mean(lmsgscvec3(:,end-9:end), 2);
[lms_mvdr3, lmsmvdrvec3, lmsmvdrout3] = lms(X3, S1, w_mvdr_source3, mew, N2);
ave_lmsmvdrvec3 = mean(lmsmvdrvec3(:,end-9:end), 2);
[rls_gsc3,rlsgscvec3, rlsgscout3] = rls(X3, S1, w_gsc_source1 , lambda, N2);
ave_rlsgscvec3 = mean(rlsgscvec3(:,end-9:end), 2);
[rls_mvdr3,rlsmvdrvec3, rlsmvdrout3] = rls(X3, S1, w_mvdr_source3 , lambda, N2);
ave_rlsmvdrvec3 = mean(rlsmvdrvec3(:,end-9:end), 2);
arr_resp_lmsgsc3 = abs(lmsgscvec3'*twoDsteeringvec_grid).^2;
%with N=25, the results are way worse but still usable. Peaks in the spetra
%are less pronounced and performance reduces. But with N<M there are
%failures in MUSIC, MVDR, LMS and RLS. Source detection fails.

%c: 
corr = [1 0.3 0 ; 0.3 1 0; 0 0 1];
variances = [1 10^(-4/10) 10^(-8/10)];
% Compute the covariance matrix
C = diag(sqrt(variances)) * corr * diag(sqrt(variances));
% I noticed that when sources are strongly correlated, it becomes very difficult to tell them apart. 
% The MUSIC spectrum struggles to identify all the sources correctly or misses some entirely. 
% MVDR does not block interference, leading to less precise results and more errors. 
% the adaptive algorithms are unpredictable.



function [rls, beamvec, beamout] = rls(X1, S1, w_q, lambda, N)
    [~, L] = size(S1);
    w_a = zeros(L, N); % Initialize adaptive weights
    P = eye(L); % Initialize P (inverse correlation matrix)
    rls = zeros(1, N);  % Preallocate squared error storage
    beamvec = zeros(size(w_q, 1), N); % Store beamformer vectors
    beamout = zeros(1, N); % Store beamformer output
    % Iterative RLS updates
    for i = 1:N
        d(i) = w_q' * X1(:, i);
        x(:,i) = S1'*X1(:,i);  %3x1
        pii(:,i) = P * x(:,i);
        k(:,i) = pii(:,i) / (lambda + x(:,i)' * pii(:,i));  %3x1
        P = (1 / lambda) * (P - k(:,i) * x(:,i)' * P);
        e(i) = d(i) - w_a(:,i)' * x(:,i);
        w_a(:,i) = w_a(:,i) + k(:,i)*conj(e(i));
        beamvec(:,i) = w_q - S1*w_a(:,i);
        beamout(i) = beamvec(:,i).'*X1(:,i);
        rls(i) = abs(e(i))^2;
    end
end


function [lms, beamvec, beamout] = lms(X1, S1, w_q, mew, N)
    w_a = zeros(size(S1, 2), 1); % Initialize adaptive weights
    lms = zeros(1, N); 
    beamvec = zeros(size(w_q, 1), N); % Store beamformer vectors
    beamout = zeros(1, N); % Store beamformer output
for i = 1:N
    d(i) = w_q'*X1(:, i);
    x(:,i) = S1'*X1(:,i);
    e(i) = d(i) - w_a'*x(:,i);  %beamformer output
    w_a = w_a+mew*x(:,i)*conj(e(i));
   % w(i) = w_gsc_source1 - S1*w_a(i);  %beamformer vector 
    lms(i) = abs(e(i))^2; % Record squared error for learning curve
    beamvec(:,i) = w_q - S1*w_a;
    beamout(i) = e(i);
end
% Plot the learning curve
end

%theoretical correlation R of u
function R = thoer_corr_Ru(S, rel_power, noise_power, M)
noise_power_linear = 10.^(noise_power/10);
R_V = (1/M)*noise_power_linear*eye(M);
power_linear= 10.^(rel_power/10);
R_A = diag(power_linear);
R = S*R_A*S'+R_V;
end

function V = generate_V(M, N, noise_power)
var = (1/(M))*10^(noise_power/10); %linear power
V = sqrt(var)/sqrt(2)*(randn(M,N)+1j*randn(M,N)); %complex gaussian
end

function A = generate_A(N, source_power) %coefficient of the signal based on the power at the source
var= 10.^(source_power/10); 
A = zeros(length(source_power), N);
for i=1:length(var)
    A(i,:) = sqrt(var(i))/sqrt(2)*(randn(N,1)+1j*randn(N,1));
end
end

function a_hat = a_hat(theta, phi)
a = zeros(length(theta), 3);
a_hat = zeros(length(theta), 3);
for i = 1:length(theta)
a_hat(i,:) = [sin(theta(i)*pi/180)*cos(phi(i)*pi/180) sin(theta(i)*pi/180).*sin(phi(i)*pi/180) cos(theta(i)*pi/180)];
%a_hat(i,:) = a(i,:)/norm(a(i,:));
end
end

function s_vec = steering_vector(d_lambda, a_hat, sensor_location) 
M = 41;
exp_raised = zeros(M,size(a_hat, 1));
s_vec = zeros(M,size(a_hat, 1));
for i =1:M
    for k = 1:size(a_hat, 1)
    exp_raised(i, k) = -1j*(2*pi*d_lambda)*dot(a_hat(k,:).', sensor_location(i,:));
%     exp_raised(i, k) = -1j*(2*pi*d_lambda)*dot(a_hat(k,:).', sensor_location(i,:));
    s_vec(i, k) = 1/sqrt(M)*exp(exp_raised(i, k));
    end
end
end

function S = generate_S(s_vec1, s_vec2, s_vec3)
S= [s_vec1 s_vec2 s_vec3];
end
