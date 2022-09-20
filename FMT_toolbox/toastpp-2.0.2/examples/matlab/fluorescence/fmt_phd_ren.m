%% MATLAB-TOAST sample script:
% Reconstruction of fluorochrome concentration from fluorescence
% measurements

clear all
close all

% A few general parameters
refind = 1.4;   % homogeneous refractive index of the medium
c0 = 0.3;       % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium

%% Create the mesh
% rad:有限元mesh半径
rad = 25;
% mkcircle(rad,6,64,4);
% mkcircle(rad,nsect,nring,nbnd)
% nsect:一个圆被分为几个section
% nring:半径被分为几个元素
% nbnd:number of high-resolution boundary rings 边缘部分 网格更密
% vtx:点坐标 3511
% idx:三角形元素顶点编号
% eltp:三角形元素背景标签，如果是同质的，即相等
[vtx,idx,eltp] = mkcircle(rad,6,32,2); %vtx :vector list  idx: Vertex of triangle element 
mesh = toastMesh(vtx,idx,eltp);
n = mesh.NodeCount;%size(vtx)

%% solution basis
grd = [64 64];
blen = prod(grd); %64*64

basis = toastBasis (mesh, grd);
nsol = basis.slen; % 基函数长度

%% Source and measurement projection vectors Source and detector locations
np = 32; %source和detector数量
for i=1:np
    qphi = (i-1)/np * 2*pi;
    qpos(i,:) = rad * [cos(qphi) sin(qphi)];
    mphi = (i-0.5)/np *2*pi;
    mpos(i,:) = rad * [cos(mphi) sin(mphi)];
end
%SetQM 方法将源位置和探测器位置附加到网格
mesh.SetQM(qpos,mpos);
% plot网格以及source和detector位置
mesh.Display
hold on
plot(qpos(:,1),qpos(:,2),'ro','MarkerFaceColor','r');
plot(mpos(:,1),mpos(:,2),'bx','MarkerFaceColor','b');

% 
dmask = mesh.DataLinkList; % 1*1024 2^np
% 从源和探测器位置，您可以使用 Qvec 和 Mvec 函数创建源和边界投影矢量
qvec = mesh.Qvec ('Neumann', 'Gaussian',0.5); %Generate a sparse matrix of source column vectors.
nQ = size(qvec,2);

mvec = mesh.Mvec ('Gaussian', 0.5, refind); %Generate a sparse matrix of measurement column vectors.
nM = size(mvec,2);

%% The flu target in grid dimensions: 3 different blobs 
fltgt_g = zeros(grd);
nblobs = 3;
fcontrast = [0.1 0.06 0.1];
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

% map to solution basis
fltgt_h   = basis.Map ('B->M', reshape(fltgt_g,[],1));
muatgt_h  = basis.Map ('B->M', reshape(muatgt_g,[],1));
mustgt_h  = basis.Map ('B->M', reshape(mustgt_g,[],1));
 
figure;
subplot(2,2,1);imagesc(fltgt_g); colormap(gray);colorbar;axis square;title('target fluoresence');
subplot(2,2,2);imagesc(muatgt_g);colormap(gray);colorbar;axis square;title('Background \mu_a');
subplot(2,2,3);imagesc(mustgt_g);colormap(gray);colorbar;axis square;title('background \mu_s');

%% modulation frequency [MHz]
freq = 0;

% refractive index
ref = ones(n,1) * refind;

%% create the simulated data

% FEM system matrix
smat = dotSysmat (mesh, muatgt_h, mustgt_h, ref, freq); %刚度矩阵

% excitation photon density
sol1=smat\qvec;
sol1=full(sol1); % mappings dont work with sparse matrices

% fluoresence photon density
sol_fl=smat\(diag(fltgt_h)*sol1);

plot_sol1eta=1;
if plot_sol1eta==1
    dr=(diag(fltgt_h)*sol1);
    figure;
    for jj=1:nQ
        subplot(6,6,jj);imagesc(reshape(basis.Map ('M->B',dr(:,jj)),grd)),colormap(gray);axis square;
    end
end

% excitation boundary data
% (size(mvec,2))*size(sol1,2) : detector个数*source个数
lgamma = reshape ((mvec.' * sol1), (size(mvec,2))*size(sol1,2), 1);
lgamma = lgamma(dmask);
xdata = [real(lgamma)];

% fluorescence boundary data
lgamma_fl = reshape ((mvec.' * sol_fl), (size(mvec,2))*size(sol1,2), 1);
lgamma_fl = lgamma_fl(dmask);
fdata = [real(lgamma_fl)];

R = sparse(diag(1./(xdata))); % scaling matrix
toastWriteVector('./data/fl_sim_mod_mat.fem',fdata)
figure;
subplot(2,2,1);imagesc(reshape(xdata,nM,nQ));colormap(gray);colorbar;axis square;title('Excitation Data RtN map');
subplot(2,2,2);imagesc(reshape(fdata,nM,nQ));colormap(gray);colorbar;axis square;title('Fluoresence Data RtN map');
subplot(2,2,3);imagesc(reshape(R*fdata,nM,nQ));colormap(gray);colorbar;axis square;title('Ratio of Fluoresence/Excitation');


%% inverse solution: reconstruct chromophore concentration

r = R*fdata;  % rescale fluorescence data by excitation. Y 通过激发重新缩放荧光数据
asol1 = smat\mvec; % adjoint 

%% Jacobian
Jr = zeros(nM*nQ,nsol); %1024*3096
for i = 1:nQ
    for j = 1:nM
        tmp_h = full(sol1(:,i).*asol1(:,j));
        tmp_s = basis.Map ('M->S',tmp_h); %B:grid M:mesh S:basis
        Jr(( i-1)*nM +j,:) = tmp_s; %每一行
    end
end
Jr = R*Jr; % rescale Jacobian 缩放雅各比矩阵

% column scaling if required S: 列归一化矩阵
c = ones(nsol,1);
for k = 1:nsol
    c(k) = norm(Jr(:,k));
end
S = sparse(diag(1./c));
S = 0.001*speye(nsol); % reset to identity (i.e. no image scaling)
Jr = Jr * S; % 进行列归一化



c = 1e5;

Algorithm=2;
lambda = 1e-4*(trace(Jr' * Jr)); % regularisation proportional to trace
if Algorithm == 1
    disp('solving using Backslash');
    tr = [r; ones(nsol,1)];
    tJ = [Jr ; lambda*eye(nsol)];
    x = tJ\tr;
else 
    disp('solving using Conjugate Gradients');
    s_tol=1e-6;
    s_iter=200;
    [x] = L1LS(Jr,r,lambda,c);
    %[x,CGflag,relres,niter] = cgs(Jr'*Jr + lambda*speye(nsol), Jr'*r,s_tol,s_iter);
   % disp(['CG iterations: ',num2str(niter)]);
end
x = S*x; %unscale
% map to solution basis
xim = basis.Map('S->B', x);
sim = basis.Map('S->B', diag(S));
 
figure;
subplot(2,2,1);imagesc(fltgt_g);colorbar;axis square;title('target fluoresence');
subplot(2,2,2);imagesc(reshape(sim,grd));colormap(gray);colorbar;axis square;title('Image recscaling');
subplot(2,2,3);imagesc(reshape(xim,grd));colormap(gray);colorbar;axis square;title('reconstructed fluoresence');
