%this code is for test
%edit  pathdef.m
%% set startup and make rgb image to gray


clear all;
close all;
graypic = double(rgb2gray(imread("popcat1.jpg")))/255;
imshow(graypic);
%set blurred picture and the true picture
x_true = graypic;
gamma = [4,7];
[PSF, center] = Gauss( gamma , size(x_true,1));
PSF(center(1)+1:end,center(2)+1:end) = 0;
PSF = PSF/sum(PSF(:));
    
A = psfMatrix(PSF, 'reflexive', center);
b = A*x_true;

nLevel = 0.002; % noise level
[N, sigma] = WhiteNoise(b, nLevel, 0); % add noise 
bn = b+N;

figure(1); 
subplot(1,3,1), imshow(x_true,[]), title('True Image')
subplot(1,3,2), imshow(bn,[]), title('Noisy Blurred Image')
subplot(1,3,3), imshow(PSF,[]), title('PSF')

fprintf('Noise Level = %-5f\n', nLevel);
%% test data using LSQR
maxit = 250;
input = HyBRset('InSolv', 'none', 'x_true', x_true, 'Iter', maxit, 'ReOrth','on');
[x_lsqr, output_lsqr] = HyBR(A,bn,[],input);

figure(2), plot(output_lsqr.Enrm,'ko','LineWidth',2), hold on
figure(3), subplot(2,3,1), imshow(x_lsqr,[]), title('lsqr')

%% Standard HyBR with Tikhonov and dp
input = HyBRset('InSolv', 'Tikhonov', 'x_true', x_true,'Iter', 50,'RegPar','dp', 'ReOrth','on','nLevel', sigma);
[xhybr_tik, output_hybr_tik] = HyBR(A,bn,[],input);

figure(2), plot(output_hybr_tik.Enrm,'rs','LineWidth',2), hold on
figure(3), subplot(2,3,2), imshow(xhybr_tik,[]), title('HyBR')

%% HyBR recycle
% (1) Set trunction options and matrices
trunc_options.nOuter = 11; % number of outer iterations
trunc_options.nInner = 50; % maximum storage of solution vector space
trunc_options.max_mm = 30; % maminimum number of vectors to save at compression
trunc_options.compress = 'RBD'; 
trunc_mats = [];

% (2) Set options for hybrid method
input = HyBRset('InSolv', 'Tikhonov', 'x_true', x_true,'Iter', 50,'RegPar', 'dp','nLevel', sigma); 

% (3) Run HyBR recycle
[xhybr_recycle, toutput_recycle, trunc_mats] = HyBRrecycle(A,bn,[],input, trunc_options, trunc_mats);

figure(2), plot(toutput_recycle.Enrm,'b','LineWidth',2), hold on
legend('LSQR','HyBR','HyBR-recycle-dp')
figure(3), subplot(2,3,3), imshow(reshape(xhybr_recycle,size(x_true)),[]), title('HyBR-recycle-dp')

X = [x_lsqr(:), xhybr_tik(:), xhybr_recycle(:)];
[E,cmin,cmax,ns,nt] = geterrorimages(X, x_true, gcf);

fprintf('Compression Method Used: %s\n',trunc_options.compress);

