clear all
close all
grd = [64 64];

%% The flu target in grid dimensions: 3 different blobs 
fltgt_g = zeros(grd);
nblobs = 3;
fcontrast = [0.3 0.6 0.3];
blbcx = grd(1)*[0.7 0.75 0.3]; % xcentre of blobs
blbcy = grd(2)*[0.3 0.75 0.7]; % ycentre of blobs
blbrad = grd(1)*[0.1 0.1 0.1]; % radius of blobs

for i = 1:grd(1)
  for j = 1:grd(2)
      for k = 1:nblobs
        if( (i-blbcx(k))^2 + (j-blbcy(k))^2 < blbrad(k)^2)
            fltgt_g(i,j) = fltgt_g(i,j)+fcontrast(k);
        end
      end
  end
end
muatgt_g=0.05*ones(grd);
mustgt_g=1*ones(grd);
kap_homog = 1/(3*(muatgt_g(1,1)+mustgt_g(1,1))); %diffusion coefficient

figure;
subplot(2,2,1);imagesc(fltgt_g); colormap(gray);colorbar;axis square;title('target fluoresence');
subplot(2,2,2);imagesc(muatgt_g);colormap(gray);colorbar;axis square;title('Background \mu_a');
subplot(2,2,3);imagesc(mustgt_g);colormap(gray);colorbar;axis square;title('background \mu_s');

%% generate data
cons = 0.2;
u = 64; v = 64;
r_d = ones(v,2);
r_d(1:16,2) = 64;
r_d(1:16,1) = 1:4:64;
r_d(17:32,2) = 0;
r_d(17:32,1) = 1:4:64;
r_d(33:48,2) = 1:4:64;
r_d(33:48,1) = 0;
r_d(49:64,2) = 1:4:64;
r_d(49:64,1) = 64;


r_s = ones(u,2);
r_s(1:16,2) = 64;
r_s(1:16,1) = 2:4:64;
r_s(17:32,2) = 0;
r_s(17:32,1) = 2:4:64;
r_s(33:48,2) = 2:4:64;
r_s(33:48,1) = 0;
r_s(49:64,2) = 2:4:64;
r_s(49:64,1) = 64;


r_n = ones(64*64,2);    
linspec_s = {'rx','MarkerSize',8,'LineWidth',2};
linspec_d = {'b.','MarkerSize',8,'LineWidth',2};
subplot(2,2,4);plot(r_d(:,1),r_d(:,2),linspec_d{:}); 
hold on; 
plot(r_s(:,1),r_s(:,2),linspec_s{:}); 
axis([0 64 0 64]);
hold off; 

for i = 1: length(r_n(:,1))
    r_n(i,1) = mod(i-1,64)+1;
    r_n(i,2) = floor((i-1)/64)+1;
end

%% fwd

fltgt_g = fltgt_g(:);
n = grd(1)*grd(2);
k = sqrt(muatgt_g(1,1)/kap_homog);
weight = zeros(u*v,n);
projFl = zeros(u*v,1);
projExc = zeros(u*v,1);

for j = 1:u
    for l = 1:v
        cnt = (j-1)*v+l;
        projExc(cnt) = greenFunction(k,r_d(l,:),r_s(j,:));
        projFl(cnt) = 0;
        for i = 1:n
            weight(cnt,i) = greenFunction(k,r_s(j,:),r_n(i,:))*...
                greenFunction(k,r_n(i,:),r_d(l,:))/projExc(cnt);
            weight(cnt,i) = max(weight(cnt,i),0);
            projFl(cnt) = projFl(cnt)+greenFunction(k,r_s(j,:),r_n(i,:))*...
                fltgt_g(i)*greenFunction(k,r_n(i,:),r_d(l,:));
            
        end
    end
end

% add noise 
eta = 0.000001;
noise_exc = randn(size(projExc));
noise_fl = randn(size(projFl));
noise_exc = max(eta*norm(projExc)*noise_exc/norm(noise_exc),0);
noise_fl = max(eta*norm(projFl)*noise_fl/norm(noise_fl),0);

projNb = max((projFl+noise_fl)./(projExc+noise_exc),0);

%% inversion
% weight = zeros(u*v,n);
% 
% for i = 1:n
%     for j = 1:u
%         for l = 1:v
%             temp_w = greenFunction(k,r_s(j,:),r_n(i,:))*greenFunction(k,r_n(i,:),r_d(l,:))...
%                 /(greenFunction(k,r_d(l,:),r_s(j,:)));
%             weight((j-1)*v+l,i) = max(temp_w,0);
%         end
%     end
% end

lambda = 1;
L = eye(size(weight,2));
options.lbound = zeros(size(grd(1)*grd(2)));
options.ubound = ones(size(grd(1)*grd(2)));
[X] = cglsTikConstraint(weight,projNb,1:20,lambda,L,options);
%[X] = cglsAIR(weight,projFl',1:20);
figure;
for i = 1:16
    subplot(4,4,i);
    imagesc(reshape(real(X(:,i)),grd));colormap(gray),colorbar;axis square;
end

X_prime = reshape(X(:,end),grd(1),grd(2));
figure;
subplot(1,2,1);imagesc(reshape(fltgt_g,grd));colorbar;axis square;title('target fluoresence');
% subplot(2,2,2);imagesc(reshape(sim,grd));colormap(gray);colorbar;axis square;title('Image recscaling');
subplot(1,2,2);imagesc(abs(X_prime));colormap(gray);colorbar;axis square;title('reconstructed fluoresence');

%% green Function
function g = greenFunction(k,r1,r2)
    temp = norm(r1 - r2);
    g = exp(k*temp*1i)/(temp*4*pi);
end


