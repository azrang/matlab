close all;

%% HW 5 Multidim

hlap = (1/6)*[1 4 1 ; 4 -20 4 ; 1 4 1];
hx= (1/4)*[1 0 -1 ; 2 0 -2; 1 0 -1];
hy = (1/4)*[-1 -2 -1 ; 0 0 0 ; 1 2 1];

hlapsob = conv2(hx, hx)+conv2(hy, hy);
 % all of these values are real, so they are conjugate 
% symmetric and therefore a zero phase filter. Also 5x5 FIR

%zerophase filter: z= z* , all real 

[Hlap, f1, f2 ] =freqz2(hlap);  
[Hlapsob, f11, f22 ] =freqz2(hlapsob);  
sign(Hlap);  %yes all negative
sign(Hlapsob);  %yes all negative
Hlap(imag(Hlap)>0);
Hlapsob(imag(Hlapsob)>0);

%surface plots
figure;
surf(f1, f2, Hlap)  
figure;
surf(f11, f22, Hlapsob)   
figure;
contour(f1, f2, Hlap)  %isotropic from center, bandpass
figure;
contour(f11, f22, Hlapsob)   %isotropic, lowpass

figure;
image(Lilyx)
colormap('gray')
figure;
image(Rodanx)
colormap('gray')
figure;
image(filter2(hlap, Lilyx))    %for convolution with a matrix, without freq
%response
colormap('gray')
figure;
image(filter2(hlap, Rodanx)) 
colormap('gray')
figure;
image(filter2(hlapsob, Rodanx)) 
colormap('gray')
figure;
image(filter2(hlapsob, Lilyx))  
colormap('gray')

up_Lily = upsamp_matrix(Lilyx);
up_Rodan = upsamp_matrix(Rodanx);

figure;
image(up_Lily)
colormap('gray')
figure;
image(up_Rodan)
colormap('gray')
figure;
image(up_Lily(255:800, 255:800))
colormap('gray')
figure;
image(up_Rodan(255:800, 255:800))
colormap('gray')
figure;
image(up_Lily(400:600, 400:600))
colormap('gray')

%computing magnitude 2D DFT
dft_Lily = fft2(Lilyx);
dft_up_Lily = fft2(up_Lily);
dft_Rodan = fft2(Rodanx);
dft_up_Rodan = fft2(up_Rodan);

figure;
image(abs(fftshift(dft_Lily)))
colormap('gray')
figure;
image(abs(fftshift(dft_up_Lily)))
colormap('gray')
figure;
image(abs(fftshift(dft_Rodan)))
colormap('gray')
figure;
image(abs(fftshift(dft_up_Rodan)))
colormap('gray')

%upsampling distortion 

function v = upsamp_matrix(u)
for i = 1:1:size(u, 2)
v1 = [u(:, i); zeros(size(u, 1), 1)];
v2 = reshape(v1, [size(u, 1), 2]);
v_vert(:, i+i) = v2(:, 1);
v_vert(:, i+i+1)= v2(:, 2);
end

for k = 1:1:size(v_vert, 1)
v3 = [v_vert(k, :); zeros(1, size(v_vert, 2))];
v(k+k, :) = v3(1, :);
v(k+k+1, :)= v3(2, :);
end
end


