%% Example for Dynamic Problem
clear all
close all 
clc
startup_genHyBRrecycle
startup_Recycle
%% Setup for the Original Problem
n = 64;
n_t = 25;
x_true = dynamic_test(n_t, n);
%x_true = test_get_x_true(n,n_t);
x_true(x_true<0)=0;
x_true = reshape(x_true,n,n,n_t);
for i = 1:min(25,n_t)
%     subplot(5,5,i), imshow(x_true(:,:,i),[])
subplot(5,5,i), imagesc(x_true(:,:,i))
end

%% Setup the forward problem

% Prob - Choose forward Problem
prob = 1;
% 1 - Sheprical tomography from IRtool
% 2 - raytracing
% 3 - PSF
%nlevel = 0.01;
nlevel = 0.01;
P = [];
Aall = Select_Forward_jjh(prob,n,x_true(:,:,i),nlevel,n_t);
nrays = 90; 
for i = 1:n_t % number of sources
   Amat{i} = Aall((i-1)*nrays+1:i*nrays,:);  %((i-1)*25+1:i*25,:);
   P = blkdiag(P, Amat{i});
end
%P = kronMat(eye(n_t),Aall);
b = P*x_true(:);

% Add noise to all projection images
noise = randn(size(b(:))); 
sig = nlevel*norm(b(:)) / norm(noise(:));
bn = b(:) + sig * noise;
R = speye(length(b(:)));
%% create changing x_true
% for i = 1:n_t
%     subplot(5,5,i), imshow(x_true(:,:,i),[])
% end
% 
% figure, imagesc(reshape(bn,n/2,n/2))

%% Select Q(Matern) as Kronecker product
xmin = [0 0 0];           %Coordinates of left corner
xmax = [1 1 1];             %Coordinates of right corner
nvec = [n n n_t];
 theta = [1 1 120];
nvec_s = [n, n];
%Additional parameters governing length scales. Can be used to model
%nonisotropic covariance kernels

% (1) Matern kernel
%  nu1 = 1.5; ell1 = .1;
% %nu1 = .5; ell1 = .01;
% %nu1 = 0.1; ell1 = 1.5;
% %nu1 = 0.1; ell1 = 0.1;
% %nu1 = .5; ell1 = .25;
% %nu1 = .5; ell1 = .015;
% % k1 = @(r) matern(r,nu1,ell1);
% % nu1 = inf; ell1 = .1;
% %  nvec = [n n n_t];

% [Q1s,k1s,Q1r] = getcovfun(nvec_s,'matern', nu1, ell1);
% Q_1k = getcovfun(nvec,'matern', nu1, ell1, theta);
% 
% % time matrix D
% % D = @(r) matern(r,nu1,ell1); 
% % Dt = fullmatrix(0, 1, n_t, D, 1);
% % Q1 = kronMat(D,Q1s);
% 
% 
% %(2) Kronecker product
% nvec_s = [n, n];
% %nu2 = 0.5; ell2 = 0.25; 
%  nu2 =   1; ell2 =  .15;
% %nu2 = 0.1; ell2 = 1.5;
% %nu2 = 0.1; ell2 = 0.1; 
% [Q2s,k2s,Q2r] = getcovfun(nvec_s,'matern', nu2, ell2);
% 
% % time matrix D
% D = diag(ones(n_t+1,1),0)-diag(ones(n_t,1),1); D = D(1:end-1,:);
% Q2t = inv(D*D'+.0001*eye(n_t));
% Q_2k = kronMat(Q2t,Q2s);
% 





% (2) Select Q(Rational Quadratic) as Kronecker product
% Build row/column of the Toeplitz matrix
% alpha2 = 100; l2 = .5;
alpha2 = 0.1; l2 = 0.05;
% alpha2 = 0.01; l2 = 0.01;
% alpha2 = 100; l2 = 0.1;
% alpha2 = 100; l2 = 0.01; % actual values
k2 = @(r) rational_quadratic(r,alpha2,l2);
Qr2 = createrow(xmin,xmax,nvec,k2,theta);
Qfun2 = @(x) toeplitzproduct(x, Qr2, nvec);
Q_2k = funMat(Qfun2,Qfun2, [nvec(1)*nvec(2)*nvec(3) nvec(1)*nvec(2)*nvec(3)]);
% 


%% Solve x with various methods
maxit = 70; NS = 0; PS = [];

% LSQR
% [PS, NS] = All_HyBRs('LSQR',P,bn,Q_1k,[],R,'LSQR','tikhonov','optimal',x_true(:),nlevel,maxit,PS, NS,[],[],[]);
% [PS, NS] = All_HyBRs('LSQR',P,bn,Q_1k,[],R,'LSQR','tikhonov','gcv',x_true(:),nlevel,maxit,PS, NS,[],[],[]);
% [PS, NS] = All_HyBRs('LSQR',P,bn,Q_1k,[],R,'LSQR','tikhonov','upre',x_true(:),nlevel,maxit,PS, NS,[],[],[]);
% [PS, NS] = All_HyBRs('LSQR',P,bn,Q_1k,[],R,'LSQR','tikhonov','wgcv',x_true(:),nlevel,maxit,PS, NS,[],[],[]);


% genHyBR - Q1
[PS, NS] = All_HyBRs('genHyBR',P,bn,Q_2k,[],R,'genHyBR','tikhonov','optimal',x_true(:),nlevel,maxit,PS, NS,[],[],[]); 
[PS, NS] = All_HyBRs('genHyBR',P,bn,Q_2k,[],R,'genHyBR','tikhonov','gcv',x_true(:),nlevel,maxit,PS, NS,[],[],[]); 
[PS, NS] = All_HyBRs('genHyBR',P,bn,Q_2k,[],R,'genHyBR','tikhonov','upre',x_true(:),nlevel,maxit,PS, NS,[],[],[]);
[PS, NS] = All_HyBRs('genHyBR',P,bn,Q_2k,[],R,'genHyBR','tikhonov','wgcv',x_true(:),nlevel,maxit,PS, NS,[],[],[]); 

% % genHyBR - Q2
% [PS, NS] = All_HyBRs('genHyBR2', P,bn,Q_2k,[],R,'genHyBR','tikhonov','optimal',x_true(:),nlevel,maxit,PS, NS,[],[],[]); 
% [PS, NS] = All_HyBRs('genHyBR2', P,bn,Q_2k,[],R,'genHyBR','tikhonov','gcv',x_true(:),nlevel,maxit,PS, NS,[],[],[]); 
% [PS, NS] = All_HyBRs('genHyBR2', P,bn,Q_2k,[],R,'genHyBR','tikhonov','upre',x_true(:),nlevel,maxit,PS, NS,[],[],[]); 
% [PS, NS] = All_HyBRs('genHyBR2', P,bn,Q_2k,[],R,'genHyBR','tikhonov','wgcv',x_true(:),nlevel,maxit,PS, NS,[],[],[]); 

% % Set trunction options and matrices
% trunc_options.nOuter = 12; % number of outer iterations
% trunc_options.nInner = 20; % maximum storage of solution vector space
% trunc_options.max_mm = 15;  % maximum number of vectors to save at compression
% trunc_options.compress = 'SVD'; 
% trunc_mats = [];
maxit2 = 70;
% % genHyBR_recycle - Q1
% [PS, NS] = All_HyBRs('genHyBRrecycle1',P,bn,Q_1k,[],R,'genHyBRrecycle','tikhonov','optimal',x_true(:),nlevel,maxit2,PS, NS,[],trunc_options,trunc_mats); 
% [PS, NS] = All_HyBRs('genHyBRrecycle1',P,bn,Q_1k,[],R,'genHyBRrecycle','tikhonov','gcv',x_true(:),nlevel,maxit2,PS, NS,[],trunc_options,trunc_mats); 
% [PS, NS] = All_HyBRs('genHyBRrecycle1',P,bn,Q_1k,[],R,'genHyBRrecycle','tikhonov','upre',x_true(:),nlevel,maxit2,PS, NS,[],trunc_options,trunc_mats);
% [PS, NS] = All_HyBRs('genHyBRrecycle1',P,bn,Q_1k,[],R,'genHyBRrecycle','tikhonov','wgcv',x_true(:),nlevel,maxit2,PS, NS,[],trunc_options,trunc_mats); 
% 
% % genHyBR_recycle - Q2
% [PS, NS] = All_HyBRs('genHyBRrecycle2',P,bn,Q_2k,[],R,'genHyBRrecycle','tikhonov','optimal',x_true(:),nlevel,maxit2,PS, NS,[],trunc_options,trunc_mats); 
% [PS, NS] = All_HyBRs('genHyBRrecycle2',P,bn,Q_2k,[],R,'genHyBRrecycle','tikhonov','gcv',x_true(:),nlevel,maxit2,PS, NS,[],trunc_options,trunc_mats); 
% [PS, NS] = All_HyBRs('genHyBRrecycle2',P,bn,Q_2k,[],R,'genHyBRrecycle','tikhonov','upre',x_true(:),nlevel,maxit2,PS, NS,[],trunc_options,trunc_mats);
% [PS, NS] = All_HyBRs('genHyBRrecycle2',P,bn,Q_2k,[],R,'genHyBRrecycle','tikhonov','wgcv',x_true(:),nlevel,maxit2,PS, NS,[],trunc_options,trunc_mats); 
% 




%% Reconstruction image
for ij = 1:size(PS,1)
    PS{ij,3} = reshape(PS{ij,3},n,n,n_t);
end
for j = 1:size(PS,1)
figure, 
for i = 1:min(25,n_t)
    subplot(5,5,i), imagesc(PS{j,3}(:,:,i))
end
    sgtitle([PS{j,1} ' with ' PS{j,2}])
end
% maxit2 = PS{9,4}.iterations-1;

%% Reconstruction optimal
LABEL = {};
REG = {'optimal'};
PS_C = PS_Select(PS,LABEL,REG);

% Reconstruction
%PlotGraphs('recon',PS_C,x_true,n,maxit);
% Relative Errors
PlotGraphs('relerr',PS_C,x_true,n,maxit,maxit2);
% print(gcf,'-depsc','optimal_compare2.eps')

%% Reconstruction upre
LABEL = {};
REG = {'upre'};
PS_C = PS_Select(PS,LABEL,REG);

% Reconstruction
%PlotGraphs('recon',PS_C,x_true,n,maxit);
% Relative Errors
PlotGraphs('relerr',PS_C,x_true,n,maxit,maxit2);
% print(gcf,'-depsc','upre_compare2.eps')

%% Reconstruction wgcv
LABEL = {};
REG = {'wgcv'};
PS_C = PS_Select(PS,LABEL,REG);

% Reconstruction
%PlotGraphs('recon',PS_C,x_true,n,maxit);
% Relative Errors
PlotGraphs('relerr',PS_C,x_true,n,maxit,maxit2);
% print(gcf,'-depsc','wgcv_compare2.eps')

%% Reconstruction gcv
LABEL = {};
REG = {'gcv'};
PS_C = PS_Select(PS,LABEL,REG);

% Reconstruction
%PlotGraphs('recon',PS_C,x_true,n,maxit);
% Relative Errors
PlotGraphs('relerr',PS_C,x_true,n,maxit,maxit2);
% print(gcf,'-depsc','wgcv_compare2.eps')


























