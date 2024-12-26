%% Signals Dynamical Systems HW
close all

%simulating a continuous time system under 0-input condition

E= [3 1 ; 2 1];
A_0 = [-0.3 0.4 ; -0.4 -0.3];

[eigvec, eigval0] = eig(A_0);
[eigval0] = sort(diag(eigval0), 'descend')

A= E*A_0*inv(E);
[eigvec, eigval1] = eig(A);
[eigval1] = sort(diag(eigval1), 'descend')
% yes these have the same eigenvalues

s = sym('s'); %creates symbolic var s
% one sided laplace transform of matrix
eqn = inv(s*eye(2)-A);
phi = ilaplace(eqn);
ht = matlabFunction(phi);

x_d0 = [2 1].';
%x_d0 = [2, zeros(1, 99); 1, zeros(1, 99)];
f_s= 10;
dt = 1/f_s;

%finding xd[n]
%x_d = x_d0;
x_d(:,1) = x_d0;
for n = 1:100
    %x_d(:, n)= (ht(dt))^(n*dt)*x_d0;
    %x_d(:, n) = ht(n*dt)*x_d0;
    x_d(:, n+1) = ht(dt)*x_d(:, n);
end

%finding xde[n]
%x_de = x_d0;
x_de(:,1) = x_d0;
for n1 = 1:100
    %x_de(:, n1+1)= (eye(2)+ht(dt)*dt)*x_de(:, n1)
    x_de(:, n1+1)= (dt.* A + eye(2))*x_de(:, n1);
end

%x_dm = x_d0;
x_dm(:,1) = x_d0;
for n2 = 1:100
    x_dm(:, n2+1)= (eye(2)+dt.*A+(A^2)*(dt^2)/2)*x_dm(:, n2);
end

% plotting state space trajectories
figure;
grid on;
plot(x_d(1,:), x_d(2,:), 'DisplayName', 'x_d');
hold on;
plot(x_dm(1,:), x_dm(2,:), 'DisplayName', 'x_dm');
hold on;
plot(x_de(1,:), x_de(2,:), 'DisplayName', 'x_de');
title('state space trajectory for x_dm');
xlabel('x1(t)');
ylabel('x2(t)');
legend
hold off;

diff1 = abs(x_d - x_de);
diff2 = abs(x_d - x_dm);
max(diff1, [], 'all')
max(diff2, [], 'all')

% x_dm has a lower maximum absolute error, so I would say x_dm is the
% better method, which uses midpoint

% the part of the trajectory that has the largest error is when there is a
% a sharp turn, or a sharp change between elements in each state for x_dm. 



