close all;

%% HW 3

%channel capacity in terms of pe
x = linspace(0,1,100);
y = -x+1;
figure;
plot(x, y);
title( 'channel capacity for BEC');
xlabel('pe');
ylabel('C(bits)');
%this is different from BSC, as for the BSC the capacity is equal
%to 1-phi. This in in turn creates a parabolic looking graph. I got 
% a logarithmic graph. 

%compute D and R for a Gaussian vector

var_v = [4 3 0.8 0.1];
lambda1 = .7;
lambda2 = .9;
compute_dr(var_v, .9); %R = 1.9445
compute_dr(var_v, .7); %R = 2.4034

[D, R] = compute_dr(var_v,.87)  %THIS IS R=1.9934
display('D = [.87 .87 .8 .1], R = [ 1.1005 .8929 0 0]')

Dist_x = linspace(lambda1,lambda2, 100);
R_x = zeros(1, 100);
for k =1:100
[Di(k), R_x(k)] = compute_dr(var_v, Dist_x(k));
end

figure;
plot(Di, R_x);
title('Distortion vs. Rate');

% d) increasing lambda decreases Rate and increases Distortion


[D, R] = compute_dr(var_v, k);

function [D, R] = compute_dr(var_v, level)
    for i= 1:length(var_v)
        if(level>= var_v(i))
            D(:,i) = var_v(i);
        elseif(level <= var_v(i))
            D(:,i) = level;
        end
        R(:, i) = .5*log2(var_v(i)/D(:,i));
    end
    D= sum(D);
    R = sum(R);
end
   
