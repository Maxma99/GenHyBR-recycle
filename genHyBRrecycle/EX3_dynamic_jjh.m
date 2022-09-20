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
nlevel = 0;
Aall = Select_Forward_jjh(prob,n,x_true(:,:,i),nlevel,n_t);
P = [];
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

% figure, imagesc(reshape(bn,n/2,n/2))


%% Select Q(Matern) as Kronecker product
xmin = [0 0 0];           %Coordinates of left corner
xmax = [1 1 1];             %Coordinates of right corner
nvec = [n n n_t];

%Additional parameters governing length scales. Can be used to model
%nonisotropic covariance kernels

% (1) Matern kernel
% nu1 = 1; ell1 = .75;
%nu1 = .5; ell1 = .01;
 nu1 = 0.1; ell1 = 1.5;
%nu1 = 0.1; ell1 = 0.1;
% nu1 = .5; ell1 = .25;
% nu1 = .5; ell1 = .015;
% k1 = @(r) matern(r,nu1,ell1);
% nu1 = inf; ell1 = .1;
%  nvec = [n n n_t];
 theta = [1 1 120];
nvec_s = [n, n];
[Q1s,k1s,Q1r] = getcovfun(nvec_s,'matern', nu1, ell1);
Q_1k = getcovfun(nvec,'matern', nu1, ell1, theta);

% time matrix D
% D = @(r) matern(r,nu1,ell1); 
% Dt = fullmatrix(0, 1, n_t, D, 1);
% Q1 = kronMat(D,Q1s);


%(2) Kronecker product
nvec_s = [n, n];
nu2 = 0.5; ell2 = 0.25; 
%nu1 = 0.1; ell1 = 1.5;
%nu2 = 0.1; ell2 = 0.1; 
[Q2s,k2s,Q2r] = getcovfun(nvec_s,'matern', nu2, ell2);

% time matrix D
D = diag(ones(n_t+1,1),0)-diag(ones(n_t,1),1); D = D(1:end-1,:);
Q2t = inv(D*D'+.0001*eye(n_t));
Q_2k = kronMat(Q2t,Q2s);






% % (2) Select Q(Rational Quadratic) as Kronecker product
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



%%
% Solve for x with variaous methods
maxit = 120; NS = 0; PS = [];

% genHyBR - Q1
 [PS, NS] = All_HyBRs('genHyBR1',P,bn,Q_1k,[],R,'genHyBR','tikhonov','optimal',x_true(:),nlevel,maxit,PS, NS,[]); 
% [PS, NS] = All_HyBRs('genHyBR1',P,bn,Q_1k,[],R,'genHyBR','tikhonov','gcv',x_true(:),nlevel,maxit,PS, NS,[]); 
% [PS, NS] = All_HyBRs('genHyBR1',P,bn,Q_1k,[],R,'genHyBR','tikhonov','upre',x_true(:),nlevel,maxit,PS, NS,[]);
% [PS, NS] = All_HyBRs('genHyBR1',P,bn,Q_1k,[],R,'genHyBR','tikhonov','wgcv',x_true(:),nlevel,maxit,PS, NS,[]); 

% % genHyBR - Q2
 [PS, NS] = All_HyBRs('genHyBR2', P,bn,Q_2k,[],R,'genHyBR','tikhonov','optimal',x_true(:),nlevel,maxit,PS, NS,[]); 
% [PS, NS] = All_HyBRs('genHyBR2', P,bn,Q_2k,[],R,'genHyBR','tikhonov','gcv',x_true(:),nlevel,maxit,PS, NS,[]); 
% [PS, NS] = All_HyBRs('genHyBR2', P,bn,Q_2k,[],R,'genHyBR','tikhonov','upre',x_true(:),nlevel,maxit,PS, NS,[]); 
%  [PS, NS] = All_HyBRs('genHyBR2', P,bn,Q_2k,[],R,'genHyBR','tikhonov','wgcv',x_true(:),nlevel,maxit,PS, NS,[]); 
% % 
% %mixHyBR - Q1,Q2
  %[PS, NS] = All_HyBRs_jjh('mixHyBR',P,bn,Q_1k,Q_2k,R,'mixHyBR',[],'optimal',x_true(:),[],maxit,PS, NS,[]); 
% [PS, NS] = All_HyBRs_jjh('mixHyBR',P,bn,Q_2k,Q_1k,R,'mixHyBR',[],'gcv',x_true(:),[],maxit,PS, NS,[]); 
% [PS, NS] = All_HyBRs_jjh('mixHyBR',P,bn,Q_2k,Q_1k,R,'mixHyBR',[],'upre',x_true(:),nlevel,maxit,PS, NS,[]);
%     [PS, NS] = All_HyBRs_jjh('mixHyBR',P,bn,Q_2k,Q_1k,R,'mixHyBR',[],'wgcv',x_true(:),nlevel,maxit,PS, NS,[]); 

% [PS, NS] = All_HyBRs_jjh('mixHyBR_j',P,bn,Q_1k,Q_2k,R,'mixHyBR',[],'optimal',x_true(:),nlevel,maxit,PS, NS,[]); 
% [PS, NS] = All_HyBRs_jjh('mixHyBR_j',P,bn,Q_1k,Q_2k,R,'mixHyBR',[],'upre',x_true(:),nlevel,maxit,PS, NS,[]);
% [PS, NS] = All_HyBRs_jjh('mixHyBR_j',P,bn,Q_1k,Q_2k,R,'mixHyBR',[],'wgcv',x_true(:),nlevel,maxit,PS, NS,[]); 

%
for ij = 1:size(PS,1)
    PS{ij,3} = reshape(PS{ij,3},n,n,n_t);
end


%
LineColor = [0, 0.4470, 0.7410;      0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 
                     0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];
                 
figure,semilogy(PS{1,4}.Enrm(2:end),'--.','color',LineColor(4,:),'linewidth',1), hold on
  semilogy(PS{2,4}.Enrm(2:end),'--+','color',LineColor(5,:),'linewidth',1), hold on
  semilogy(PS{3,4}.Enrm(2:end),'-','color',LineColor(3,:),'linewidth',1)
%    semilogy(PS{4,4}.Enrm,'-.','color',LineColor(1,:),'linewidth',1)
%  semilogy(PS{5,4}.Enrm,'-o','color',LineColor(2,:),'linewidth',1)
% %semilogy(PS{6,4}.Enrm,'*','color',LineColor(2,:),'linewidth',1)
% % 
% % %axis([1 maxit 0.45 1.1])
%  xlim([1 maxit])
% % %legend('genHyBR1-opt','genHyBR2-opt','mixHyBR-opt','mixHyBR-upre','mixHyBR-gcv','mixHyBR-wgcv')
% % legend('genHyBR1-opt','genHyBR2-opt','mixHyBR-opt')
% % ylabel('relative error','fontweight','bold')
% % xlabel('iteration','fontweight','bold')

% %% SAVE VARIABLES
% file = 'EX2_seismic_dynamic_jjh2';
% save(file);

%%

figure,
PS1 = reshape(PS{1,3},n,n,n_t);
for i = 1:min(25,n_t)
%     subplot(5,5,i), imshow(PS1(:,:,i),[])
    subplot(5,5,i), imagesc(PS1(:,:,i))
end
title('Reconstruction')

figure,
PS2 = reshape(PS{2,3},n,n,n_t);
for i = 1:min(25,n_t)
%     subplot(5,5,i), imshow(PS2(:,:,i),[])
    subplot(5,5,i), imagesc(PS2(:,:,i))
end

figure,
PS3 = reshape(PS{3,3},n,n,n_t);
for i = 1:min(25,n_t)
%     subplot(5,5,i), imshow(PS3(:,:,i),[])
    subplot(5,5,i), imagesc(PS3(:,:,i))
end

%%
idx = 10;
figure,
subplot(1,4,1), imagesc(PS1(:,:,i))
subplot(1,4,2), imagesc(PS2(:,:,i))
subplot(1,4,3), imagesc(PS3(:,:,i))
subplot(1,4,4), imagesc(x_true(:,:,i))


return
% 
%% optimal 
% mixHyBR
figure,
LABEL = {'mixHyBR'};
REG = {'wgcv'};
PS_C = PS_Select(PS,LABEL,REG);
nline = ceil(n_t/5);
for i = 1:n_t
x_recon = reshape(PS_C{1,3},n,n,n_t);
subplot(nline,5,i), imagesc(x_recon(:,:,i))
end


%% optimal 
% genHyBR1
figure,
LABEL = {'genHyBR1'};
REG = {'wgcv'};
PS_C = PS_Select(PS,LABEL,REG);
nline = ceil(n_t/5);
for i = 1:n_t
x_recon = reshape(PS_C{1,3},n,n,n_t);
subplot(nline,5,i), imagesc(x_recon(:,:,i))
end




%% Reconstruction optimal
LABEL = {};
REG = {'optimal'};
PS_C = PS_Select(PS,LABEL,REG);

% Reconstruction
%PlotGraphs('recon',PS_C,x_true,n,maxit);
% Relative Errors
PlotGraphs('relerr',PS_C,x_true,n,maxit);
% print(gcf,'-depsc','optimal_compare2.eps')

%% Reconstruction upre
LABEL = {};
REG = {'upre'};
PS_C = PS_Select(PS,LABEL,REG);

% Reconstruction
%PlotGraphs('recon',PS_C,x_true,n,maxit);
% Relative Errors
PlotGraphs('relerr',PS_C,x_true,n,maxit);
% print(gcf,'-depsc','upre_compare2.eps')

%% Reconstruction wgcv
LABEL = {};
REG = {'wgcv'};
PS_C = PS_Select(PS,LABEL,REG);

% Reconstruction
%PlotGraphs('recon',PS_C,x_true,n,maxit);
% Relative Errors
PlotGraphs('relerr',PS_C,x_true,n,maxit);
% print(gcf,'-depsc','wgcv_compare2.eps')

%% Reconstruction gcv
LABEL = {};
REG = {'gcv'};
PS_C = PS_Select(PS,LABEL,REG);

% Reconstruction
%PlotGraphs('recon',PS_C,x_true,n,maxit);
% Relative Errors
PlotGraphs('relerr',PS_C,x_true,n,maxit);
% print(gcf,'-depsc','wgcv_compare2.eps')

%% compare mix
LABEL = {'mixHyBR'};
REG = {};
PS_C = PS_Select(PS,LABEL,REG);

% Reconstruction
%PlotGraphs('recon',PS_C,x_true,n,maxit);
% Relative Errors
PlotGraphs('relerr',PS_C,x_true,n,maxit);
% print(gcf,'-depsc','wgcv_compare2.eps')