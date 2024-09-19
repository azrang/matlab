
%% HW 5
clear all;
close all; 

M = 8;
alpha = (1/(3*M))*[20 12 8 8 4 4]
dl = [2*sqrt(1) sqrt(8) 4*sqrt(1) sqrt(20) 6*sqrt(1) sqrt(40) ]
beta = (dl.^2)/(4)



%keep beta = 1, which corresponds with dmin as distance as approximation
gam = 2:.001:8;
gamma = 10.^(gam/10);
approx_smallest_beta = alpha(1)*qfunc( ((dl(1)^2.*gamma)./(4)).^0.5 );

approx_union = zeros(1, length(dl));
for i= 1:length(gamma)
  %sum for each dl in parallel
  ap1(:, i) = alpha.*qfunc( (((dl).^2.*gamma(i))./(4)).^0.5 );
  approx_union(i) = sum(ap1(:,i));
end


figure;
semilogy(gam, approx_smallest_beta);  %whenever want to plot y axis in log
title('SNR vs. probability of error');
hold on;
semilogy(gam, approx_union);
hold off;
legend({'smallest beta', 'union'});
ylabel('probability of error in dB');
xlabel('SNR in dB')

%comparing P_err
%for 2 dB
P2_err_one = alpha(1)*qfunc( ((dl(1)^2.*10.^(2/10))./(4)).^0.5 );
P2_err_all = sum(alpha.*qfunc( ((dl.^2.*10.^(2/10))./(4)).^0.5 ));
ratio = (1 - P2_err_one) / P2_err_all
%for 8 dB
P8_err_one = alpha(1)*qfunc( ((dl(1)^2.*10.^(8/10))./(4)).^0.5 );
P8_err_all = sum(alpha.*qfunc( ((dl.^2.*10.^(8/10))./(4)).^0.5 ));
ratio_8 = (1 - P2_err_one) / P2_err_all

%this ratio is the same for both 2dB and 8dB. 
