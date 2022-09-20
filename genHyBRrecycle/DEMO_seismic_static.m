%% DEMO_seismic_static.m
% larger static example
% C.Ma

%% Setup the forward problem

startup_genHyBRrecycle
clear all
close all
clc

% n = 256;
n = 256;
nu = 0.5;
ell = 0.25;
% nt = 2;
% nu_t = 2.5;
% ell_t = 0.1;
% nu = 0.5 or 1.5 or 2.5
x_true = get_x_true(n, nu, ell);

figure(1), imagesc(reshape(x_true,n,n))
title(['True Slowness Field nx 256','with nu = ',num2str(nu),' and ell = ', num2str(ell)],'fontsize',15)

% Prob - Choose forward Problem
prob = 4;
% prob = 4; %ns and nd = 25
% 1 - Sheprical tomography from IRtool
% 2 - raytracing
% 3 - PSF
nlevel = 0.01;
[A,bn,R,nlevel] = Select_Forward(prob,n,x_true,nlevel);
sigma = nlevel;

%% Setup Parameters
Regpar = 'wgcv';
solver = 'tikhonov';
maxit = 30;
maxit2 = 150;
% Set trunction options and matrices
compression = 'SVD'; % compression method
itOuter = 8; % number of outer iterations
itInner = 30; % maximum storage of solution vector space
maxvec = 20;  % maximum number of vectors to save at compression
truncmats = [];
P = [];
% Setup HyBR options for HyBRset
input = HyBRset('InSolv', solver, 'x_true', x_true(:), 'Iter', maxit2, 'RegPar', Regpar, 'nLevel', nlevel); 



%% Setup covariance matrix
% Covariance kernel of the reconstruction
nx = n; ny = n;
xmin = [0 0];           %Coordinates of left corner
xmax = [1 1];           %Coordinates of right corner
nvec = [nx, ny];        %Number of points. Here this is a 256 x 256 grid
theta = [1 1];

% alpha = 100; l = 0.01; 
% alpha = 1; l = 0.01;
alpha = 100; l = 0.3;
k = @(r) rational_quadratic(r,alpha,l);
% l = .25;
% k1 = @(r) exp(-(r.^2)/(2*l^2));

Qr = createrow(xmin,xmax,nvec,k,theta);
Qfun = @(x) toeplitzproduct(x, Qr, nvec);
Q = funMat(Qfun,Qfun, nvec.*nvec);
spy(reshape(Qr,n,n));
figure(2),imagesc(reshape(Qr,n,n)), title(['Parameter alpha = ',num2str(alpha),' and l = ', num2str(l),' and nLevel = ', num2str(nlevel)]);


%% Set and run genHyBR with Tikhonov when storage is limited 

geninput = HyBR_lsmrset('InSolv', solver,'RegPar',Regpar, 'x_true', x_true(:),'Iter', maxit, 'nLevel', nlevel);
[x_genHyBR, output_genHyBR] = genHyBR(A, bn, Q, R, geninput);

figure(4), plot(output_genHyBR.Enrm,'y--','LineWidth',2), hold on
figure(4), title(['genHyBR, HyBR recycle and genHyBR recycle with ' ,solver, ' and ', Regpar]);
figure(3), subplot(2,3,1), imagesc(reshape(x_genHyBR,size(x_true))), title('genHyBR')

%% Set and run HyBR recycle
% (1) Set trunction options and matrices
trunc_options.nOuter = itOuter; % number of outer iterations
trunc_options.nInner = itInner; % maximum storage of solution vector space
trunc_options.max_mm = maxvec;  % maximum number of vectors to save at compression
trunc_options.compress = compression; 
trunc_mats = truncmats;

removepath
startup_Recycle


% (2) Run HyBR recycle
[xhybr_recycle, toutput_recycle, ~] = HyBRrecycle(A,bn,P,input, trunc_options, trunc_mats);

figure(4), plot(toutput_recycle.Enrm,'b--','LineWidth',2), hold on
figure(3), subplot(2,3,2), imagesc(reshape(xhybr_recycle,size(x_true))), title('HyBR recycle')

removepath
startup_genHyBRrecycle

%% Set and run genHyBR recycle using TSVD for compression
% (1) Set trunction options and matrices
trunc_options.nOuter = itOuter; % number of outer iterations
trunc_options.nInner = itInner; % maximum storage of solution vector space
trunc_options.max_mm = maxvec;  % maximum number of vectors to save at compression
trunc_options.compress = compression; 
trunc_mats = truncmats;

% (2) Run genHyBR recycle

[xhybr_genrecycle, toutput_genrecycle, trunc_mats] = genHyBRrecycle(A,bn,P,Q,R,input, trunc_options, trunc_mats);

figure(4), plot(toutput_genrecycle.Enrm,'k--','LineWidth',2), hold on
figure(3), subplot(2,3,3), imagesc(reshape(xhybr_genrecycle,size(x_true))), title('genHyBR recycle')

figure(4), legend('genHyBR','HyBR-recycle','genHyBR-recycle')


%% show error image

X = [x_genHyBR(:), xhybr_recycle(:),xhybr_genrecycle(:)];
figure(3),
[E,cmin,cmax,ns,nt] = geterrorimages(X, x_true, gcf);




