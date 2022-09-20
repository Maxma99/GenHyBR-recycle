bmua = imread('demo_matlab_fwd2_mua.png');
bmua = double(bmua)./255.*0.02 + 0.01;
figure;
subplot(1,2,1);
imagesc(bmua);
axis equal tight; colorbar
title('ground truth');

mua_homog = 0.2; % absorption effects
mus_homog = 1; % scattering effects
kap_homog = 1/(3*(mua_homog+mus_homog)); %diffusion coefficient
% Scale to desired parameter ranges

cons = 0.2;
u = 1; v = 16;
r_d = ones(v,2);
r_d(:,2) = 64;
r_d(:,1) = 1:4:64;

r_s = [32,-16];
r_d(end+1,:) = r_s;
r_n = ones(64*64,2);    

r = [16,16];
for i = 1: length(r_n(:,1))
    r_n(i,1) = mod(i,64);
    r_n(i,2) = floor(i-1/64)+1;
end

bmua = ones(64,64)*mua_homog;
bmus = ones(64,64)*mus_homog;
figure;
%subplot(1,2,1);
imagesc(bmua); colormap gray;
hold on;
plot(r_d(:,1), r_d(:,2),'*','linewidth',2);
axis equal tight; colorbar
title('\mu_a');
hold off;
% subplot(1,2,2);
% imagesc(bmus);
% hold on;
% plot(r_d(:,1), r_d(:,2),'*','linewidth',2);
% plot(r_s(:,1), r_s(:,2),'+','r','linewidth',2);
% axis equal tight; colorbar
% title('\mu_s');
% hold off;
n = 64*64;
k = sqrt(mua_homog/kap_homog);
%% fwd
for j = 1:u
    for l = 1:v
        p_mua = 0;
        for i = 1:n
            p_mua = p_mua+greenFunction(k,r_s,r_n(i,:)*mua_homog*greenFunction(k,r_n(i,:),r)
        end
    end
end


%% inversion
weight = zeros(u*v,n);
for i = 1:n
    for j = 1:u
        for l = 1:v
            weight(j*l,n) = a*greenFunction(k,r_s,r_n(i,:))*greenFunction(k,r_n(i,:),r)/(greenFunction(k,r,r_s));
        end
    end
end



function g = greenFunction(k,r1,r2)
    temp = norm(r1 - r2);
    g = exp(-k*temp)/(temp*4*pi);
end

