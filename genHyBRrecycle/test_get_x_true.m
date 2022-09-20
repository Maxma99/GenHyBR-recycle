function xx = test_get_x_true(n, nt)
% J. Chung 1/20/20
%

% Get a spatial image
% n = 128;
nu = .5; ell = .25;
% tic
% x_true = get_x_true(n, nu, ell);
% toc
% 
% figure, imagesc(x_true)
% title('True Slowness Field','fontsize',20)

% get a space/time image
% nt = 9; 
nu_t = .5; ell_t = .01;
tic
x_true = get_x_true(n, nu, ell, nt, nu_t, ell_t);
toc

% figure, 
% for i=1:9
%     subplot(3,3,i), imagesc(x_true(:,:,i))
% end
% title('True Slowness Field','fontsize',20)


%% Get a mixture image
x_image = dynamic(nt, 4, n); x_image = reshape(x_image, n, n, nt);
% x_image = dynamic(nt, 7, n); x_image = reshape(x_image, n, n, nt);

% figure, 
% for i=1:9
%     subplot(3,3,i), imagesc(x_image(:,:,i))
% end
% title('Image','fontsize',20)

% gamma = .5;
% gamma = 1;
gamma = 0.01;
xx = sqrt(gamma)*x_true + sqrt(1-gamma)*x_image;
% figure, 
% for i=1:9
%     subplot(3,3,i), imagesc(xx(:,:,i))
% end
% title('Image Sum','fontsize',20)