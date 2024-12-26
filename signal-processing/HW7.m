%% Signal Generation: 
% generating L source vectors--> to create S, then B, V, then A

M=100;
K=20;
N =200;  %N=2M
L=3;
Pdb = [0 -2 -4];  %source power
Pndb = 10;  %noise power

A_data = A(M, N,K, Pdb, Pndb );
size(A_data);

%correlation matrix
R = (1/N)*A_data*A_data.';

%% Analysis 

%singular value decomposition of A
[U, sval, V] = svd(A_data);
% correlation matrix eigen decomposition
sval = diag(sval);  % converts to a vector 
[eigvec, eigval0] = eig(R);
[eigval, idx] = sort(diag(eigval0), 'descend');
%eigenvalues in sorted order
eigvec = eigvec(:, idx);

%plotting singular values of A
figure;
stem(sval);
title('SVD of A');
%plotting sorted eigenvalues of R
figure;
stem(eigval);
title('eigenvalues of R')

%taking third largest/fourth largest ratio
dominant_sval = sval(1:3)
dominant_eigen = eigval(1:3)
ratio_sval = sval(3)/sval(4)
ratio_eig = eigval(3)/eigval(4)

% compute projection matrix to noise subspace
Ul= U(:,1:3);
P_s = Ul*Ul.'; %projection onto S
P_n = eye(M) - P_s;

%compute R inverse
R_inv = inv(R);

%% check Smusic and Smvdr spectra

%plug in 3 known source vector
S_mus1= music(L,M,K, 1, P_n)
S_mus2= music(L,M,K, 2, P_n)
S_mus3= music(L,M,K, 3, P_n)
S_mvdr1= mvdr(L,M,K,1, R_inv)
S_mvdr2= mvdr(L,M,K,2, R_inv)
S_mvdr3= mvdr(L,M,K,3, R_inv)

%random source vectors 
for i = 1:20
    a = zeros(M,1);
    ind = randperm(M, K);
    a(ind) = 1/sqrt(K);
   test(:,i) = a; 
end

%finding spectra
for i=1:20
    test_musicvec = music(20,M,K, i, P_n);
    test_music(:,i) = test_musicvec;
    test_mvdrvec = mvdr(20,M,K, i, P_n);
    test_mvdr(:,i) = test_mvdrvec;
end
avg_test_mus1= mean(test_music)
avg_test_mvdr1= mean(test_mvdr)
max_test_mus1= max(test_music)
max_test_mvdr1= max(test_mvdr)
med_test_mus1=median(test_music)
med_test_mvdr1=median(test_mvdr)

%I think that MUSIC is better at identifying source vectors because 
% based on the results I got, I did get really high numbers, almost
% infinity for the source vectors. I know that I should be getting
% appreciably larger values as compared to teh test vectors, but this is
% not happening likely due to some error in my code

%% Try again
M_new= 100;
N_new=50;
Anew_data= A(M_new, N_new, K, Pdb, Pndb);
%singular value decomposition of A
[U1, sval_new, V1] = svd(A_data);
% correlation matrix eigen decomposition
sval_new = diag(sval_new);  % converts to a vector 
figure;
stem(sval_new);
ratiosval_new = sval_new(3)/sval_new(4)
% compute projection matrix to noise subspace
Ul_new= U1(:,1:3);
P_s_new = Ul_new*Ul_new.'; %projection onto S
P_n_new = eye(M_new) - P_s_new;
%plug in 3 known source vector
S_mus1n= music(L,M_new,K, 1, P_n_new)
S_mus2n= music(L,M_new,K, 2, P_n_new)
S_mus3n= music(L,M_new,K, 3, P_n_new)

%random source vectors 
for i = 1:20
    a = zeros(M_new,1);
    ind = randperm(M_new, K);
    a(ind) = 1/sqrt(K);
   testn(:,i) = a; 
end

%finding spectra
for i=1:20
    testn_musicvec = music(20,M_new,K, i, P_n_new);
    testn_music(:,i) = testn_musicvec;
end

avg_test_mus1n= mean(testn_music)
max_test_mus1n= max(testn_music)
med_test_mus1n=median(testn_music)

%the music algorithm still works, but the performance is noticeably worse

%% One Last Thing 

display(one_last(M,K, L))
% The columns of S are linearly indeependent because the determinant of 
%its gramian is non-zero. We are then able to determine that there are 3
%source vectors. If this was completely the identity matrix, then we would
%be able to say that the columns of S form an orthonormal basis of real
%values(M). Additionally, the off-diagonal entries tell you how close 
% to orthogonal your source vectors are, if you get every diagonal 
% entry to be 0, your source vectors would be orthogonal, but if not then 
% the closer to 0 you are, the closer the inner product between the 
% two vectors/ the entries are to 0. 

function one_last= one_last(M,K,L)
for i = 1:L
    a = zeros(M,1);
    ind = randperm(M, K);
    a(ind) = 1/sqrt(K);
   S(:,i) = a; 
end
one_last=S.'*S;
end


%% functions
function S_mus = music(L,M,K, vec, P_n)
for i = 1:L
    a = zeros(M,1);
    ind = randperm(M, K);
    a(ind) = 1/sqrt(K);
   S(:,i) = a; 
end
S_mus= 1/(S(:,vec).'*P_n*S(:,vec));
end

function S_mvdr = mvdr(L,M,K, vec, R_inv)
for i = 1:L
    a = zeros(M,1);
    ind = randperm(M, K);
    a(ind) = 1/sqrt(K);
   S(:,i) = a; 
end
S_mvdr= 1/(S(:,vec).'*R_inv*S(:,vec));
end


function A = A(M, N, K, Pdb, Pndb)
% creating V with noise vectors
var_V= 10.^(Pndb/10);
V = sqrt(var_V).*randn(M, N);
%creating S matrix using 3 "special" source vectors, length M 
%Pdb should have length of L=3
L=3;
for i = 1:L
    a = zeros(M,1);
    ind = randperm(M, K);
    a(ind) = 1/sqrt(K);
   S(:,i) = a; 
end
% creating B with with b coefficients
%pdb should be a diagonal matrix
var_B= 10.^(Pdb/10).';
B = (sqrt(var_B)).*randn(L, N);

A = S*B + (1/sqrt(M)).*V;
end





