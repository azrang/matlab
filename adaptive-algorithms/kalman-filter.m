%% HW 4 :  Kalman Filter
close all;
clear;
% creating a standard kalman filter for linear systems
%% 1 : Preliminary analysis

A1 = [ -3/2 -5/2 ; 1/2 1/2];  % the state transition matrix
A2 = [ -3 -5 ; 1 1];
C = [1 0];   %output 
Qx = [ 0.1 0 ; 0 0.1];
Qy = .05;
% N =2 state variables, p =1 output

eval_A1 = eig(A1);
eval_A2 = eig(A2);  %A1 is stable since all evals are inside the unit circle
Obs_1 = [C ; C*A1];
Obs_2 = [C ; C*A2];
rank(Obs_1); rank(Obs_2); %both ranks = N = 2, so each system is observable
% to find K the steady state prediction error cov matrix, need to solve ricatti eq
K1 = idare(A1', C', Qx, Qy, zeros(2,1), eye(2));
K2 = idare(A2', C', Qx, Qy, zeros(2,1), eye(2));

%% 2 : Kalman Filter
Niter = 100;
x0 = [1; 0];
K_n_n1 = 0.1*eye(2);  %initial value for K(n, n-1)
x_n_n1 = x0;
x = zeros(2, Niter+1);
x(:,1) = x0 ; %initialization
v_x = (randn(Niter+1, 2) * chol(Qx, 'lower')).';
v_y = (randn(Niter+1, 1) * chol(Qy, 'lower')).';
K_hat_p_less = zeros(1, Niter/2);
K_hat_e_great = zeros(1, Niter/2);
y = zeros(1, Niter+1);  % Measurements
x_n_n = zeros(2, Niter+1);  % State estimate (estimated state vector)
k_diff = zeros(1,100);

for n=1:Niter+1
    if n <= 50
        A_n = A1;
        kideal = K1;
    else
        A_n = A2;
        kideal = K2;
    end
    
    y = C*x(:,n)+v_y(:,n);
    x(:,n+1) = A_n*x(:,n)+v_x(:,n);
    R_n = inv(C * K_n_n1 * C' + Qy);
    G_n = K_n_n1 * C' * R_n;
    alpha_n = y - C * x_n_n1(:,n);  %xhat(n|n-1)
    x_n_n(:,n) = x_n_n1(:,n) + G_n * alpha_n; %xhat(n|n)
    x_n_n1(:,n+1) = A_n*x_n_n(:,n);  %because they are the same in the next iteration step
    K_n_n = (eye(2)-G_n*C)*K_n_n1*(eye(2)-G_n*C)'+G_n*Qy*G_n'; %joseph form
    K_n_n1 = A_n * K_n_n * A_n' + Qx;
    k_diff(n) = norm(K_n_n1 - kideal);
    if n<=50
        normep(:,n) = norm(x(:,n) - x_n_n1(:,n));
        normee(:,n) = norm(x(:,n) - x_n_n(:,n));
    elseif n >= 51
        normep(:,n) = norm(x(:,n) - x_n_n1(:,n));
        normee(:,n) = norm(x(:,n) - x_n_n(:,n));
    end

end

%resizingcheck
x = x(:, 1:101);
x_n_n = x_n_n(:, 1:101);
x_n_n1 = x_n_n1(:, 1:101);

idx = 50;
ep_less = x(:,1:idx)-x_n_n1(:, 1:idx);
ep_great = x(:,idx+1:end)-x_n_n1(:, idx+1:end);
Kp_great = (1/idx)*ep_great*ep_great.';
Kp_less = (1/idx)*ep_less*ep_less.';
ee_less = x(:,1:idx)-x_n_n(:, 1:idx);
ee_great = x(:,idx+1:end)-x_n_n(:, idx+1:end);
Ke_great = (1/idx)*ee_great*ee_great.';
Ke_less = (1/idx)*ee_less*ee_less.';
norm(K1-Kp_less)
norm(K2-Kp_great)
deltakless = Kp_less-Ke_less;
deltakgreat = Kp_great-Ke_great;
check1 = eig(deltakless)>0;   %these are not positive definite.
% the estimated state should be better than the predicted state.
check2 = eig(deltakgreat)>0;

figure;
plot(k_diff);
title('||K(n,n-1) - Kideal||');
%the error is very close to 0 because I change the Kideal so the results match very well. 
% The transition point occurs at n=51 and it takes two time steps to converge. 

% b
figure;
hold on;
plot(x(1, :), x(2, :), 'b-', 'DisplayName', 'True trajectory');
% Predicted trajectory
plot(x_n_n1(1, :), x_n_n1(2, :), 'r--', 'DisplayName', 'Predicted trajectory');
% Estimated trajetory
plot(x_n_n(1, :), x_n_n(2, :), 'g-.', 'DisplayName', 'Estimated trajectory');
xlabel('x(1)');
ylabel('x(2)');
title('Comparison of True, Predicted, and Estimated Trajectories');
legend show;
grid on;
hold off;

%c 

figure;
hold on;
plot(x(1, 51:end), x(2, 51:end), 'b-', 'DisplayName', 'True trajectory');
title('True Trajectory for n> 50');  %yes it does blow up

figure;  %it is blowing up here
plot(normep ,'b', 'DisplayName', 'predicted difference')
hold on; 
plot(normee,'r', 'DisplayName' , 'estimated difference')
title('predicted and estimated difference with true x');
legend show;



