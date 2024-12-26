%% HW 5 :  Unscented Kalman Filter
close all;
clear;

Niter= 500;

%initial conditions
xinit = [ 6400.4; 349.14; -1.8093; -6.7967; .6932];
Qx = diag([10^(-8) 10^(-8) 2.404*10^(-5) 2.404*10^(-5) 10^(-8)]);
Qy = diag([ 1 17*10^-3]);
x_n_n1_init = [6400; 350; -2; -7; .65];
K_n_n1 = diag([10^-4 10^-4 10^-4 10^-4 1]);
dt = 0.1;
%dynamical system output measurements

rng(0X16345678,"philox");  

%generation of sigma points 
w0 = 1/3;
Nx = 5;
sigma_n = zeros(Nx, 2*Nx+1); 
w_n = zeros(1,2*Nx+1);
w_n(1)= w0;
for i=1:Nx
% e_n = zeros(2*Nx+1, 1);
e_n = zeros(Nx, 1);
e_n(i) = 1;
sigma_n(:,i+1) = sqrt(Nx/(1-w0))*e_n;
sigma_n(:, Nx+i+1) = -sqrt(Nx/(1-w0))*e_n;
w_n(i+1) = (1-w0)/(2*Nx);
w_n(Nx+i+1) = (1-w0)/(2*Nx);
end

%% Theory and Setup


%1
emp_mean = sum(w_n.*sigma_n,2); %yes it is =0
emp_cov = zeros(Nx, Nx);
for i =1:2*Nx+1
diff = sigma_n(:,i) -emp_mean;
emp_cov = emp_cov+w_n(i)*(diff*diff.');  %yes it is =1 
end
for d = 1:Nx
    X_d = sigma_n(d, :);
    sigma_d = sqrt(emp_cov(d, d)); % variance along elements of covariance matrix
    skewness_empirical(d) = sum(w_n .* ((X_d - emp_mean(d, :)).^3)) / (sigma_d^3);  %yes it is 0 for each compoennt
end
for d = 1:Nx
    X_d = sigma_n(d, :);
    sigma_d = sqrt(emp_cov(d, d)); % variance along elements of covariance matrix
    kurtosis_empirical(d) = sum(w_n .* ((X_d - emp_mean(d,:)).^4)) / (sigma_d^4);  
    Nx/(1-w0); %yes it is equal to this value for each component
end
% the kurtosis for a gaussian is =3, and so this kurtosis is higher and
% htere is an excess kurtosis which means s the probability of a large deviation 
% from the mean is greater than for the Gaussian distribution with the same variance
w0 = -.7;
w_n = zeros(1,2*Nx+1);
w_n(1)= w0;
for i=1:Nx
e_n = zeros(Nx, 1);
e_n(i) = 1;
sigma_n(:,i+1) = sqrt(Nx/(1-w0))*e_n;
sigma_n(:, Nx+i+1) = -sqrt(Nx/(1-w0))*e_n;
w_n(i+1) = (1-w0)/(2*Nx);
w_n(Nx+i+1) = (1-w0)/(2*Nx);
end
for d = 1:Nx
    X_d = sigma_n(d, :);
    sigma_d = sqrt(emp_cov(d, d)); % variance along elements of covariance matrix
    kurtosis_empirical(d) = sum(w_n .* ((X_d - emp_mean(d,:)).^4)) / (sigma_d^4);   %close to 3, so yes it would
    %need a negative w0
end

%% The Assigment: Experiment
%1
xmeas = zeros(Nx, Niter);
ymeas = zeros(2, Niter);
beta_meas = zeros(1, Niter);
beta_est = zeros(1, Niter);
ymeas(:,1) = uhlmeas(xinit,1*dt, chol(Qy, 'lower'));
xmeas(:,1) = xinit;
[x_n_n1(:,1), K_n_n1, x_n_n(:,1), K_n_n] = ukf(x_n_n1_init, K_n_n1, Qx, Qy, ymeas(:, i+1), 1, dt);
beta0 = 0.597983;
for i=1:Niter
t = i*dt;
xmeas(:,i+1) = uhlprocsim(xmeas(:,i),t,dt,'midpoint', chol(Qx, 'lower'));   %true measurements
ymeas(:, i+1) = uhlmeas(xmeas(:,i),t, chol(Qy, 'lower'));
[x_n_n1(:,i+1), K_n_n1, x_n_n(:,i+1), K_n_n] = ukf(x_n_n1(:,i), K_n_n1, Qx, Qy, ymeas(:, i), 1, dt);  %feed in last x_n_n1 and last K_n_n1
beta_meas(i) = beta0*exp(xmeas(5,i+1));
beta_est(i) = beta0*exp(x_n_n(5,i+1));
end

figure;
plot(xmeas(1, :), xmeas(2, :), 'b-', 'DisplayName', 'True trajectory');
hold on;
plot(xmeas(1,1), xmeas(2,1), 'ro', 'DisplayName', 'actual start'); % Marker at initial point
plot(x_n_n1(1,1), x_n_n1(2,1), 'ro', 'DisplayName', 'estimated start'); % Marker at initial point
%text(xmeas(1,1), xmeas(2,1), 'Start', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
% Estimated trajetory
plot(x_n_n(1, :), x_n_n(2, :), 'g-.', 'DisplayName', 'Estimated trajectory');
xlabel('x(1)');
ylabel('x(2)');
title('Comparison of True and Estimated Trajectories');
legend show;
grid on;
hold off;
%2
figure;
plot(xmeas(3, :), xmeas(4, :), 'b-', 'DisplayName', 'True velocity');
hold on;
plot(xmeas(3,1), xmeas(4,1), 'ro', 'DisplayName', 'actual start'); % Marker at initial point
plot(x_n_n1(3,1), x_n_n1(4,1), 'ro', 'DisplayName', 'estimated start'); % Marker at initial point
%text(xmeas(1,1), xmeas(2,1), 'Start', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
% Estimated trajetory
plot(x_n_n(3, :), x_n_n(4, :), 'g-.', 'DisplayName', 'Estimated velocity');
xlabel('x(1)');
ylabel('x(2)');
title('Comparison of True and Estimated velocities');
legend show;
grid on;
hold off;
%3
figure;
i = 1:1:500;
plot(i*dt, beta_meas, 'b-', 'DisplayName', 'True beta');
hold on;
plot(i*dt, beta_est, 'g-.', 'DisplayName', 'Estimated beta');
xlabel('t');
ylabel('value');
title('Comparison of True and Estimated Betas');
legend show;
grid on;
hold off;

% I ran this 5 different times using different random numbers in the random
% number generator and each time the values looked very similar!!!!


%% FUNCTIONS
%5 UKF
function [x_n_n1, K_n_n1, x_n_n, K_n_n] = ukf(x_n_n1_init, K_n_n1, Qx, Qy, ymeas, Niter, dt)
Nx = length(x_n_n1_init);
x_n_n1 = zeros(Nx, Niter);
x_n_n1(:,1) = x_n_n1_init;
x_n_n= zeros(Nx, Niter);
%n= 1; %since function wants one iteration
t = Niter*dt;
% ymeas = uhlprocsim(xinit,t,dt,'midpoint'); %getting y output
for n= 1:Niter
mu_n1(:,n) = x_n_n1(:,n);
Cchol_n_n1 = chol(K_n_n1, 'lower');
[w_n_n1, sigx_n_n1] = sigmapoints(mu_n1(:,n), Cchol_n_n1);
%apply unscented transform
t = n*dt;
for i=1:size(sigx_n_n1,2)
sigy_n_n1(:,i)= uhlmeas(sigx_n_n1(:,i),t, 0);  %obtain predicted measurements; don't want noise in sigma points
end
y_n_n1 = sigmamean(w_n_n1, sigy_n_n1);
Kxy_n_n1 = sigmacov(w_n_n1,sigx_n_n1, sigy_n_n1); %cross-covariance between predicted state and predicted measurements
Kyy_n_n1 = sigmacov(w_n_n1, sigy_n_n1)+Qy; %covariance of the predicted measurements
epsilon = 1e-6;
Kyy_n_n1 = Kyy_n_n1 + epsilon * eye(size(Kyy_n_n1));
%apply the kalman correction eqs
alpha_n = ymeas - y_n_n1;
G_n = Kxy_n_n1*inv(Kyy_n_n1);
x_n_n(:,n) = x_n_n1(:,n)+G_n*alpha_n;
K_n_n = K_n_n1 - G_n*Kyy_n_n1*G_n.';
%compute the filtered sigma points
Cchol_n_n=chol(K_n_n, 'lower');        %it says my matrix is not positive definite
[w_n_n, sigx_n_n]= sigmapoints(x_n_n(:,n), Cchol_n_n);
%apply UT to obtain predicted states
sigx_n_n1 = uhlprocsim(sigx_n_n,t,dt,'midpoint');
x_n_n1(:,n) = sigmamean(sigx_n_n1,w_n_n );
K_n_n1 = sigmacov(w_n_n,sigx_n_n1)+Qx;
end
end

%generate sigma points and weights given input mean vector and cholesky facotr
function [w_n, sig] = sigmapoints(mu, Cchol) 
w0 = 1/3;
Nx = 5;
sigma_n = zeros(Nx, 2*Nx+1); 
w_n = zeros(1,2*Nx+1);
w_n(1)= w0;
for i=1:Nx
e_n = zeros(Nx, 1);
e_n(i) = 1;
sigma_n(:,i+1) = sqrt(Nx/(1-w0))*e_n;
sigma_n(:, Nx+i+1) = -sqrt(Nx/(1-w0))*e_n;
w_n(i+1) = (1-w0)/(2*Nx);
w_n(Nx+i+1) = (1-w0)/(2*Nx);
end
sig = mu + Cchol*sigma_n;
end
%3
% computes the empirical mean from given sigma points and weight
function emp_mean = sigmamean(w_n, sig)
emp_mean = sum(w_n.*sig, 2);
end
%4
% compute the empirical cross-covariance for X
function sigcov = sigmacov(w, sigx, sigy)
Nx = size(sigx,1);
emp_meanx = sigmamean(w, sigx);
if nargin < 3
    sigcov = zeros(Nx, Nx);
    emp_meany = emp_meanx;
    sigy = sigx;
else
    emp_meany = sigmamean(w, sigy);
    Ny = size(sigy,1);
    sigcov = zeros(Nx, Ny);
end
for i =1:2*Nx+1
diffx = sigx(:,i)-emp_meanx;
diffy = sigy(:,i)-emp_meany;
sigcov = sigcov+ w(i)*(diffx*diffy.');  %yes it is =1 
end
end
% Measurement Equation
% Computes measured y
function y= uhlmeas(x,t,Qmeasch)
% Fixed parameters
% Radar Center
ctr= [6375;0];
% Default noise Cholesky factors
if nargin<3
    % Default Qmeasch
        Qmeasch= diag(sqrt([1,17e-3]));
end
y= [norm(x(1:2)-ctr);atan2(x(2)-ctr(2),x(1)-ctr(1))]+Qmeasch*randn(2,1);
end
function xdot= uhlproc(x,t,Qprocch)
N= size(x,1); % # of states
UseEuler= 1;   % if 0, use midpoint
% ixed parameters
% ballistic coefficient
beta0= 0.597983;
% Other stuff
H0= 13.406;
Gm0= 3.986e5;
R0= 6374;
% Default noise Cholesky factors
if nargin<3
    % Default Qprocch
        %Qprocch= diag(sqrt([0,0,2.4064e-5,2.4064e-5,0.005]));
        Qprocch= diag(sqrt([0,0,2.4064e-5,2.4064e-5,0]));
        % No direct disturbance on position, only velocity
        % No disturbance on logbeta, so beta is in theory
        %   (unknown) constant
end
% Drag Coefficient
beta= beta0*exp(x(5));
% Computing quantities
R= norm(x(1:2)); % distance from center of earth
V= norm(x(3:4)); % speed
D= -beta*V*exp((R0-R)/H0);  % drag
G= -Gm0/R^3;  % gravity force related term
if 0
    [D,G]
end
v= Qprocch*randn(N,1);
xdot= [x(3);x(4);D*x(3)+G*x(1);D*x(4)+G*x(2);0]+v;
end

function xupdate= uhlprocsim(x,t,dt,method,Qprocch)
Qdefault= 0;
if nargin<4
    Qdefault= 1;
    UseMidpoint= 0;
    UseEuler= 1;   % default
    % Will use default Qprocch inside uhlproc function
elseif ischar(method(1))   % this is method
    if strcmp(lower(method(1)),'e') % use Euler
        UseEuler= 1;
        UseMidpoint= 0;
    elseif strcmp(lower(method(1)),'m') % use midpoint
        UseMidpoint= 1;
        UseEuler= 0;
    else
        disp('Error: specify Euler or midpoint');
        return
    end
    if nargin<5
        Qdefault= 1;
    end
elseif nargin<5   % 4th value is Q, no method specified
    Qprocch= method;
    UseEuler=1; % default
    UseMidpoint= 0;
else % 4th value is Qprocch, 5th value is method
    tmp= method;
    method= Qprocch;
    Qprocch= tmp;   % swap method and Qprocch
    if strcmp(lower(method(1)),'e') % use Euler
        UseEuler= 1;
        UseMidpoint= 0;
    elseif strcmp(lower(method(1)),'m') % use midpoint
        UseMidpoint= 1;
        UseEuler= 0;
    else
        disp('Error: specify Euler or midpoint');
        return
    end
end

if 0
   UseEuler
   UseMidpoint
   Qdefault
end

if UseEuler
    if Qdefault
        xupdate= x+dt*uhlproc(x,t);
    else
        xupdate= x+dt*uhlproc(x,t,Qprocch);
    end
elseif UseMidpoint
    if Qdefault
        xupdate= x+dt*uhlproc(x+dt/2*uhlproc(x,t,0),t+dt/2);
    else
        xupdate= x+dt*uhlproc(x+dt/2*uhlproc(x,t,0),t+dt/2,Qprocch);
    end
else
    disp('Error: only choices are Euler or midpoint');
    return
end
end
%5
%performs one iteration of UKF

