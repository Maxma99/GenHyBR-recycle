clear all
close all

graypic = double(rgb2gray(imread("VanGogh.png")))/255;
imshow(graypic);
%set blurred picture and the true picture
x_true = graypic(1:end-1,1:end-1);

gamma = [4,7];
[PSF, center] = Gauss( gamma , size(x_true,1));
PSF(center(1)+1:end,center(2)+1:end) = 0;
PSF = PSF/sum(PSF(:));
    
A = psfMatrix(PSF, 'reflexive', center);
b = A*x_true;
[m,n]=size(b);
nLevel = 0.002; % noise level
[N, sigma] = WhiteNoise(b, nLevel, 0); % add noise 
bn = b+N;

figure(1); 
subplot(1,3,1), imshow(x_true,[]), title('True Image')
subplot(1,3,2), imshow(bn,[]), title('Noisy Blurred Image')
subplot(1,3,3), imshow(PSF,[]), title('PSF')

% (1) Set trunction options and matrices
trunc_options.nOuter = 11; % number of outer iterations
trunc_options.nInner = 50; % maximum storage of solution vector space
trunc_options.max_mm = 30; % maminimum number of vectors to save at compression
trunc_options.compress = 'SVD'; 
trunc_mats = [];

% (2) Set options for hybrid method
input = HyBRset('InSolv', 'Tikhonov','Iter', 50,'RegPar', 'wgcv'); 

% (3) Set R and Q as covariance matrices

% (a) Select Q
%Example application of Qx in 2D
%%%%%%%%%%%%%%%%%%%%%%2d test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Consider the box (0,1)^2
xmin = [0 0];         %Coordinates of left corner
xmax = [1 1];         %Coordinates of right corner
nvec = [n, n];        %Number of points in grid

nu = 1.2; ell = .01;
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
P = [];
[xhybr_recycle, toutput_recycle, trunc_mats] = genHyBRrecycle(A,bn,P,Q,R,input, trunc_options, trunc_mats);
figure(2), plot(toutput_recycle.Enrm,'b','LineWidth',2), hold on
figure(3), subplot(2,4,3), imshow(reshape(xhybr_recycle,size(b)),[]), title('genHyBR-recycle-dp-tsvd')













