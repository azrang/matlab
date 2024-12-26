clc 
%clear 
close all

%% Signals HW 3, question 3, sampling

x = 1:7;
M = 2;
down = downsamp(M,x)
up = upsamp(M, x)
h = [3,2,1,4];
g = upsamp(M, h)   %need g to have same length as x
k = downsamp(M, conv(g, x))
p = conv(h, downsamp(M, x))

%compare k and p to get maximum absolute difference 
noble1 = max(abs(k-p)) %yes these are the same, because the difference is 0

noble2 = max(abs(upsamp(M, conv(h,x)) - conv(g, upsamp(M,x)))) %yes 
% these are the same, because the difference is 0

figure
stem(0:ceil((length(x)/M-1)), down, 'filled')

figure
stem(0:M*length(x)-M, up, 'filled')  %need same length of x and y values

function down = downsamp(M, x)
down = x(1:M:length(x));
end 

function up = upsamp(M, x)
y = x(1:1/M:length(x));
for i= 2: M*length(x)-1
    if( mod( (i-1)/M , 1)~=0)  %if remainder is not 0, then 0 element
        y(i) = 0;
    end
end
up = y(find(y,1, 'first'):find(y,1,'last')); %cuts off the extra 0s at the end of upsampled
end


