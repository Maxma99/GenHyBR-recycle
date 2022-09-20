% DEMO_Deblurring.m
%
% This script sets up and computes reconstructions for a 2D deblurring 
%   problem using hybrid projection methods with recycling.
%
%      Chung, de Sturler, and Jiang. "Hybrid Projection Methods with 
%           Recycling for Inverse Problems". SISC, 2020.
%
% Note: This code requires IRTools: https://github.com/jnagy1/IRtools
%
% J. Chung, E. de Sturler, J. Jiang, C. Ma 05/2021

%% Get problem
clear all
close all

startup_genHyBRrecycle

load Grain

n = size(b,2);
gamma = [4,7];
[PSF, center] = Gauss(gamma , size(x_true,1));
PSF(center(1)+1:end,center(2)+1:end) = 0;
PSF = PSF/sum(PSF(:));
    
A = psfMatrix(PSF, 'reflexive', center);
b = A*x_true;

nLevel = 0.005; % noise level
[N, sigma] = WhiteNoise(b, nLevel, 0); % add noise 
bn = b+N;

figure(1); 
subplot(1,3,1), imshow(x_true,[]), title('True Image')
subplot(1,3,2), imshow(bn,[]), title('Noisy Blurred Image')
subplot(1,3,3), imshow(PSF,[]), title('PSF')

%% Setup Parameters
Regpar = 'optimal';
solver = 'tikhonov';
maxit = 50;
maxit2 = 350;
% Set trunction options and matrices
compression = 'SVD'; % compression method
itOuter = 16; % number of outer iterations
itInner = 50; % maximum storage of solution vector space
maxvec = 30;  % maximum number of vectors to save at compression
truncmats = [];
P = [];
% Setup HyBR options for HyBRset
input = HyBRset('InSolv', solver, 'x_true', x_true(:), 'Iter', maxit2, 'RegPar', Regpar, 'nLevel', sigma); 



%% Standard LSQR 
inputlsqr = HyBRset('InSolv', 'none', 'x_true', x_true(:), 'Iter', maxit2, 'RegPar', Regpar, 'nLevel', sigma);
[x_lsqr, output_lsqr] = HyBR(A,bn,[],inputlsqr);

figure(2), plot(output_lsqr.Enrm,'k:','LineWidth',2), hold on
figure(3), subplot(2,5,1), imshow(x_lsqr,[]), title('lsqr')

%% Standard HyBR with Tikhonov

[xhybr_tik, output_hybr_tik] = HyBR(A,bn,[],input);

figure(2), plot(output_hybr_tik.Enrm,'r--','LineWidth',2), hold on
figure(3), subplot(2,5,2), imshow(xhybr_tik,[]), title('HyBR')

%% genHyBR recycle using TSVD for compression
% (1) Set trunction options and matrices
trunc_options.nOuter = itOuter; % number of outer iterations
trunc_options.nInner = itInner; % maximum storage of solution vector space
trunc_options.max_mm = maxvec;  % maximum number of vectors to save at compression
trunc_options.compress = compression; 
trunc_mats = truncmats;

% (3) Set R and Q as covariance matrices

% (a) Select Q
%Example application of Qx in 2D
%%%%%%%%%%%%%%%%%%%%%%2d test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Consider the box (0,1)^2
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

% % Set Q as identity matrix
% Q = speye(size(A,2));

% (b) Set R as identity matrix
R = speye(size(A,1));

%% (4) Run genHyBR recycle
[xgenhybr_recycle, toutput_recycle, trunc_mats] = genHyBRrecycle(A,bn,P,Q,R,input, trunc_options, trunc_mats);
figure(2), plot(toutput_recycle.Enrm,'g','LineWidth',2), hold on
figure(3), subplot(2,5,3), imshow(reshape(xgenhybr_recycle,size(x_true)),[]), title('genHyBR-recycle')


%% genHyBR with Tikhonov when storage is limited (optimal regularization parameter)
inputgen = HyBR_lsmrset('InSolv', solver,'RegPar',Regpar, 'x_true', x_true(:),'Iter', maxit);
[x_genHyBR_opt, output_genHyBR_opt] = genHyBR(A, bn, Q, R, inputgen);
figure(2), plot(output_genHyBR_opt.Enrm,'y--','LineWidth',2), hold on
figure(3), subplot(2,5,4), imshow(reshape(xgenhybr_recycle,size(x_true)),[]), title('genHyBR')

%% (5) Run HyBR recycle

removepath
startup_Recycle
% (1) Set trunction options and matrices
trunc_options.nOuter = itOuter; % number of outer iterations
trunc_options.nInner = itInner; % maximum storage of solution vector space
trunc_options.max_mm = maxvec;  % maximum number of vectors to save at compression
trunc_options.compress = compression; 
trunc_mats = truncmats;

[xhybr_recycle, toutput_recycle, ~] = HyBRrecycle(A,bn,P,input, trunc_options, trunc_mats);

figure(2), plot(toutput_recycle.Enrm,'b--','LineWidth',2), hold on
legend('LSQR','HyBR','genHyBR-recycle','genHyBR','HyBR-recycle')
figure(3), subplot(2,5,5), imshow(reshape(xhybr_recycle,size(x_true)),[]), title('HyBR-recycle')

removepath
startup_genHyBRrecycle





%% show error image

X = [x_lsqr(:), xhybr_tik(:), xgenhybr_recycle(:), x_genHyBR_opt(:), xhybr_recycle(:)];
[E,cmin,cmax,ns,nt] = geterrorimages(X, x_true, gcf);


%% Function to construct a Gaussian blur
function [PSF, center] = Gauss( gamma , N )
%
%     [PSF, center] = Gauss( gamma , N )
%
% This function constructs the Gaussian blur PSF, and the center of the PSF
%
%   Input:
%        gamma - parameters defining the blur
%            N - size fo the PSF
%
%    Output:
%          PSF - Gaussian point spread function
%      center  - center of the PSF
%
%  J. Chung 1/24/09
%

m = N;
n = N;

%
% Set up grid points to evaluate the Gaussian function.
%
x = -fix(n/2):ceil(n/2)-1;
y = -fix(m/2):ceil(m/2)-1;
[X,Y] = meshgrid(x,y);

%
% Compute the Gaussian PSF
%
if length(gamma) == 1
  s1 = gamma(1);
  PSF = exp( -(x.^2)/(2*s1^2));
  PSFsum = sum(PSF(:));
elseif length(gamma) == 2
  s1 = gamma(1); s2 = gamma(2);
  PSF = exp( -(X.^2)/(2*s1^2) - (Y.^2)/(2*s2^2) );
  PSFsum = sum(PSF(:));
  
elseif length(gamma)==3
  s1 = gamma(1); s2 = gamma(2); s3 = gamma(3);
  num = -((X.^2)*(s1^2) + (Y.^2)*(s2^2) - 2*(s3^2)*(X.*Y));
  den = 2*(s1^2 * s2^2 - s3^4);
  PSF = exp( num / den );
  PSFsum = sum(PSF(:));
  
end
  
PSF = PSF/PSFsum;
center = [m/2+1,n/2+1];
end
