%% Example for Dynamic Problem
clear all
close all 
clc
startup_genHyBRrecycle
startup_Recycle


%% Matrix for time smoothness not added
%% Setup for the Original Problem
%% 01
n = 64;
n_t = 25;
x_true = exp_dynamic_2D_FMT(n_t, n);
%x_true = test_get_x_true(n,n_t);
x_true(x_true<0)=0;
x_true = reshape(x_true,n,n,n_t);
for i = 1:min(25,n_t)
%     subplot(5,5,i), imshow(x_true(:,:,i),[])
figure(1),subplot(5,5,i), imshow(x_true(:,:,i),[])
end





%% 02

% n = 256;
% n = 64;
% nu = 0.5;
% ell = 0.25;
% n_t = 8;
% nu_t = 0.5;
% ell_t = 0.1;
% % nu = 0.5 or 1.5 or 2.5
% x_true = get_x_true(n, nu, ell, n_t, nu_t, ell_t);
% 
% for i = 1:min(25,n_t)
% %     subplot(5,5,i), imshow(x_true(:,:,i),[])
% figure(1), subplot(5,5,i), imagesc(x_true(:,:,i))
% end

%% Setup the forward problem

prob = 6;
% prob = 4; %ns and nd = 25
% 1 - Sheprical tomography from IRtool
% 2 - raytracing
% 3 - PSF
nlevel = 0.01;
A = {};
bn = {};
R = {};
for i = 1:min(25,n_t)
    [A{i},bn{i},R{i},nlevel] = Select_Forward(prob,n,x_true(:,:,i),nlevel);
end

sig = nlevel;


%% Select Q(Matern) as Kronecker product with time matrix
% xmin = [0 0 0];           %Coordinates of left corner
% xmax = [1 1 1];             %Coordinates of right corner
% nvec = [n n n_t]; 

%Additional parameters governing length scales. Can be used to model
%nonisotropic covariance kernels

% (1) Matern kernel
nu1 = 1.5; ell1 = .1;
% nu1 = .5; ell1 = .01;
% nu1 = 0.1; ell1 = 1.5;
% nu1 = 0.1; ell1 = 0.1;
% nu1 = .5; ell1 = .25;
% nu1 = .5; ell1 = .015;
k1 = @(r) matern(r,nu1,ell1);
% nu1 = inf; ell1 = .1;
 nvec = [n n n_t];
 theta = [1 1 120];
nvec_s = [n, n];
[Q1s,k1s,Q1r] = getcovfun(nvec_s,'matern', nu1, ell1);
Q_1k = getcovfun(nvec,'matern', nu1, ell1, theta);

% time matrix D
D = @(r) matern(r,nu1,ell1); 
Dt = fullmatrix(0, 1, n_t, D, 1);
Q1 = kronMat(D,Q1s);


% (2) Kronecker product
% nvec_s = [n, n];
% %nu2 = 0.5; ell2 = 0.25; 
%  nu2 =   1; ell2 =  .15;
% % %nu2 = 0.1; ell2 = 1.5;
% % %nu2 = 0.1; ell2 = 0.1; 
% [Q2s,k2s,Q2r] = getcovfun(nvec_s,'matern', nu2, ell2);
% 
% % time matrix D
% D = diag(ones(n_t+1,1),0)-diag(ones(n_t,1),1); D = D(1:end-1,:);
% Q2t = inv(D*D'+.0001*eye(n_t));
% Q_2k = kronMat(Q2t,Q2s);

%% static Matern Kernel

xmin = [0 0];         %Coordinates of left corner
xmax = [1 1];         %Coordinates of right corner
nvec = [n, n];        %Number of points in grid

nu = 1.5; ell = .001;
k = @(r) matern(r,nu,ell);
%Additional parameters governing length scales. Can be used to model
%nonisotropic covariance kernels
theta = [1.0 1.0];      %For now set them as isotropic

% Build row/column of the Toeplitz matrix
Qr = createrow(xmin,xmax,nvec,k,theta);
Qfun = @(x)toeplitzproduct(x, Qr, nvec);
Q = funMat(Qfun,Qfun, nvec.^2);








%% (2) Select Q(Rational Quadratic) as Kronecker product
% % Build row/column of the Toeplitz matrix
% % alpha2 = 100; l2 = .5;
% %alpha2 = 0.1; l2 = 0.1;
% alpha2 = 0.01; l2 = 0.01;
% % alpha2 = 100; l2 = 0.1;
% % alpha2 = 100; l2 = 0.01; % actual values
% k2 = @(r) rational_quadratic(r,alpha2,l2);
% Qr2 = createrow(xmin,xmax,nvec,k2,theta);
% Qfun2 = @(x) toeplitzproduct(x, Qr2, nvec);
% Q_2k = funMat(Qfun2,Qfun2, [nvec(1)*nvec(2)*nvec(3) nvec(1)*nvec(2)*nvec(3)]);

% nx = n; ny = n;
% xmin = [0 0];           %Coordinates of left corner
% xmax = [1 1];           %Coordinates of right corner
% nvec = [nx, ny];        %Number of points. Here this is a 256 x 256 grid
% theta = [1 1];
% 
% % alpha = 100; l = 0.01; 
% alpha = 1; l = 0.01;
% % alpha = 0.1; l = 0.05;
% k = @(r) rational_quadratic(r,alpha,l);
% % l = .25;
% % k1 = @(r) exp(-(r.^2)/(2*l^2));
% 
% 
% Qr = createrow(xmin,xmax,nvec,k,theta);
% Qfun = @(x) toeplitzproduct(x, Qr, nvec);
% Q = funMat(Qfun,Qfun, nvec.*nvec);
% figure(2),imagesc(reshape(Qr,n,n)), title(['Parameter alpha = ',num2str(alpha),' and l = ', num2str(l),' and nLevel = ', num2str(nlevel)]);















%% Solve x with various methods
maxit = 20; NS = 0; PS = [];
% Set trunction options and matrices
trunc_options.nOuter = 15; % number of outer iterations
trunc_options.nInner = 25; % maximum storage of solution vector space
trunc_options.max_mm = 20;  % maximum number of vectors to save at compression
trunc_options.compress = 'SVD'; 
trunc_mats = [];
maxit2 = 200;

%% Optimal


for i = 1:min(25,n_t)


[PS, NS] = All_HyBRs('genHyBRrecycle1',A{i},bn{i},Q,[],R{i},'genHyBRrecycle','tikhonov','optimal',x_true(:,:,i),nlevel,maxit2,PS, NS,[],trunc_options,trunc_mats); 
 trunc.mats = PS{i,5};
    figure(20), subplot(5,5,i), imagesc(reshape(PS{i,3},n,n)),sgtitle('Optimal Reconstruction')
end

% relative error for x
Enrm_all = [];
a = reshape(x_true,n*n,n_t);
for i = 1:maxit2
   allerr = [];
   for j = 1:n_t
     allerr = [allerr;PS{j,4}.x_all{i,1}-a(:,j)];
   end
   Enrm_all(i,1) = norm(allerr)/norm(x_true(:));
end
figure(6), plot(Enrm_all,'b--','LineWidth',2)
legend('genHyBR-optimal','genHyBRrecycle-optimal')

%% UPRE


for i = 1:min(25,n_t)


[PS, NS] = All_HyBRs('genHyBRrecycle1',A{i},bn{i},Q,[],R{i},'genHyBRrecycle','tikhonov','upre',x_true(:,:,i),nlevel,maxit2,PS, NS,[],trunc_options,trunc_mats); 
% trunc.mats = PS{i,5};
    figure(21), subplot(5,5,i), imagesc(reshape(PS{i,3},n,n)),sgtitle('UPRE Reconstruction')
    
end

Enrm_all = [];
a = reshape(x_true,n*n,n_t);
for i = 1:maxit2
   allerr = [];
   for j = 1:n_t
     allerr = [allerr;PS{j,4}.x_all{i,1}-a(:,j)];
   end
   Enrm_all(i,1) = norm(allerr)/norm(x_true(:));
end
figure(7), plot(Enrm_all,'y','LineWidth',2),
legend('genHyBR-upre','genHyBRrecycle-upre')

%% WGCV


for i = 1:min(25,n_t)


[PS, NS] = All_HyBRs('genHyBRrecycle1',A{i},bn{i},Q,[],R{i},'genHyBRrecycle','tikhonov','wgcv',x_true(:,:,i),nlevel,maxit2,PS, NS,[],trunc_options,trunc_mats); 
% trunc.mats = PS{i,5};
    figure(22), subplot(5,5,i), imagesc(reshape(PS{i,3},n,n)),sgtitle('WGCV Reconstruction')
end

Enrm_all = [];
a = reshape(x_true,n*n,n_t);
for i = 1:maxit2
   allerr = [];
   for j = 1:n_t
     allerr = [allerr;PS{j,4}.x_all{i,1}-a(:,j)];
   end
   Enrm_all(i,1) = norm(allerr)/norm(x_true(:));
end
figure(8), plot(Enrm_all,'k--','LineWidth',2)
legend('genHyBR-wgcv','genHyBRrecycle-wgcv')
%% GCV


for i = 1:min(25,n_t)


[PS, NS] = All_HyBRs('genHyBRrecycle1',A{i},bn{i},Q,[],R{i},'genHyBRrecycle','tikhonov','gcv',x_true(:,:,i),nlevel,maxit2,PS, NS,[],trunc_options,trunc_mats); 
% trunc.mats = PS{i,5};
    figure(23), subplot(5,5,i), imagesc(reshape(PS{i,3},n,n)),sgtitle('GCV Reconstruction')
end

Enrm_all = [];
a = reshape(x_true,n*n,n_t);
for i = 1:maxit2
   allerr = [];
   for j = 1:n_t
     allerr = [allerr;PS{j,4}.x_all{i,1}-a(:,j)];
   end
   Enrm_all(i,1) = norm(allerr)/norm(x_true(:));
end
figure(9), plot(Enrm_all,'g--','LineWidth',2),
legend('genHyBR-gcv','genHyBRrecycle-gcv')

%% IRhybrid_flsqr Reconstruction

% for i = 1:min(25,n_t)
%     temp = x_true(:,:,i);
%     options = IRset('x_true',temp(:) , 'RegParam', 'wgcv','MaxIter', maxit2);% wgcv, gcv, dp, upre
%     if i==1
%         [X_FLSQRidp(:,i), info_FLSQRidp] = IRhybrid_flsqr(A{i},bn{i},options);
%     else 
%         [X_FLSQRidp(:,i), info_FLSQRidp] = IRhybrid_flsqr(A{i},bn{i}-bn{i-1},options);
%         X_FLSQRidp(:,i) = X_FLSQRidp(:,i) + X_FLSQRidp(:,i-1);
%     end
%         
%     
% %     semilogy(info_FLSQRidp.Enrm,'k','LineWidth',3)
%     figure(24), 
%     subplot(5,5,i), imagesc(reshape(X_FLSQRidp(:,i),n,n)),sgtitle('IRhybrid_flsqr Reconstruction')
% end

for i = 1:min(25,n_t)
    temp = x_true(:,:,i);
    options = IRset('x_true',temp(:) , 'RegParam', 'wgcv','MaxIter', maxit2);% wgcv, gcv, dp, upre

        [X_FLSQRidp(:,i), info_FLSQRidp] = IRhybrid_flsqr(A{i},bn{i},options);


        
    
%     semilogy(info_FLSQRidp.Enrm,'k','LineWidth',3)
    figure(24), 
    subplot(5,5,i), imagesc(reshape(X_FLSQRidp(:,i),n,n)),sgtitle('IRhybrid_flsqr Reconstruction')
end

%% GK hybrid Reconstruction
R = speye(size(A{1},1));
solver = 'tikhonov';
regpar = 'wgcv';
for i = 1:min(25,n_t)
    input = HyBR_lsmrset('InSolv', solver,'RegPar',regpar, 'Iter', maxit2);
    [X_l2, output_genHyBR_opt] = genHyBR(A{i},bn{i}, Q, R,input);
%     semilogy(info_FLSQRidp.Enrm,'k','LineWidth',3)
    figure(25), subplot(5,5,i), imagesc(reshape(X_l2,n,n)),sgtitle(' GK hybrid Reconstruction')
end

%% Direct reconstruction l1

for i = 1:min(25,n_t)
    lambda = 0.01;
    N = 500;
    x_l1F = l1F(A{i},bn{i},lambda,N);
%     semilogy(info_FLSQRidp.Enrm,'k','LineWidth',3)
    figure(26), subplot(5,5,i), imagesc(reshape(x_l1F,n,n)),sgtitle('l1F Direct Reconstruction')
end

%% Direct reconstruction l2

for i = 1:min(25,n_t)
    N = 100;
    lambda = 0.001;
    L = speye(size(A{1},2));
    options.lbound = 0;
    options.ubound = 1;
    x_l2CG = cglsTikConstraint(A{i},bn{i},1:N,lambda,L,options);
    %     semilogy(info_FLSQRidp.Enrm,'k','LineWidth',3)
    figure(27), subplot(5,5,i), imagesc(reshape(x_l2CG(:,end),n,n)),sgtitle('l2CG Direct Reconstruction')
end


