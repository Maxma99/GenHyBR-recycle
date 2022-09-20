function [A,b] = FMT_Forward_2D(n,x_true)

%% A few general parameters
refind = 1.4;   % homogeneous refractive index of the medium
c0 = 0.3;       % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium
%% Create the mesh
% rad:有限元mesh半径
rad = 16;
% mkcircle(rad,6,64,4);
% mkcircle(rad,nsect,nring,nbnd)
% nsect:一个圆被分为几个section
% nring:半径被分为几个元素
% nbnd:number of high-resolution boundary rings 边缘部分 网格更密
% vtx:点坐标 3511
% idx:三角形元素顶点编号
% eltp:三角形元素背景标签，如果是同质的，即相等
% [vtx,idx,eltp] = mkcircle(rad,6,24,15); %vtx :vector list  idx: Vertex of triangle element 

[vtx,idx,eltp] = mkrect([[-rad -rad];[rad rad]],[n n]);



mesh = toastMesh(vtx,idx,eltp);
node = mesh.NodeCount;%size(vtx)

%% solution basis
grd = [n,n];
blen = prod(grd); %64*64

basis = toastBasis (mesh, grd);
nsol = basis.slen; % 基函数长度

%% Source and measurement projection vectors Source and detector locations
ns = 32;
nd = 320; %source和detector数量
% for i=1:ns
%     qphi = (i-1)/ns * 2*pi;
%     qpos(i,:) = rad * [cos(qphi) sin(qphi)];
% end
% for i=1:nd
%     mphi = (i-0.5)/nd *2*pi;
%     mpos(i,:) = rad * [cos(mphi) sin(mphi)];
% end
for i=1:ns
    len = 8*rad/ns;
    if i <= ns*(1/4)
        qpos(i,:) = [i*len-rad -rad];
    elseif ns*(1/4)<= i && i<= ns*(2/4)
        qpos(i,:) = [rad i*len-3*rad];
    elseif ns*(2/4)<= i && i <= ns*(3/4)
        qpos(i,:) = [5*rad-i*len rad];
    elseif i >= ns*(3/4)
        qpos(i,:) = [-rad 7*rad-i*len];
    end
end

for i=1:nd
    len = 8*rad/nd;
    if i <= nd*(1/4)
        mpos(i,:) = [i*len-rad -rad];
    elseif nd*(1/4)<= i && i<= nd*(2/4)
        mpos(i,:) = [rad i*len-3*rad];
    elseif nd*(2/4)<= i && i <= nd*(3/4)
        mpos(i,:) = [5*rad-i*len rad];
    elseif i >= nd*(3/4)
        mpos(i,:) = [-rad 7*rad-i*len];
    end
end




%SetQM 方法将源位置和探测器位置附加到网格
% mesh.SetQM(qpos,mpos);
% plot网格以及source和detector位置
figure(100),
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

muatgt_g=0.05*ones(grd);
mustgt_g=1*ones(grd);
fltgt_g = x_true;
% map to solution basis
fltgt_h   = basis.Map ('B->M', reshape(fltgt_g,[],1));
muatgt_h  = basis.Map ('B->M', reshape(muatgt_g,[],1));
mustgt_h  = basis.Map ('B->M', reshape(mustgt_g,[],1));

%% modulation frequency [MHz]
freq = 0;

% refractive index
ref = ones(node,1) * refind;

%% create the simulated data

% FEM system matrix
smat = dotSysmat (mesh, muatgt_h, mustgt_h, ref, freq); %刚度矩阵

% excitation photon density
sol1=smat\qvec;
sol1=full(sol1); % mappings dont work with sparse matrices
% fluoresence photon density
sol_fl=smat\(sol1.*(fltgt_h*ones(1,nQ)));
% plot_sol1eta=0;
% if plot_sol1eta==1
%     dr=(sol1.*(fltgt_h*ones(1,nQ)));
%     figure;
%     for jj=1:nQ
%         subplot(6,6,jj);imagesc(reshape(basis.Map ('M->B',dr(:,jj)),grd)),colormap(gray);axis square;
%     end
% end

% excitation boundary data
% (size(mvec,2))*size(sol1,2) : detector个数*source个数
lgamma = reshape ((mvec.' * sol1), (size(mvec,2))*size(sol1,2), 1);
if any(dmask)
    lgamma = lgamma(dmask);
end
xdata = [real(lgamma)];

% fluorescence boundary data
lgamma_fl = reshape ((mvec.' * sol_fl), (size(mvec,2))*size(sol1,2), 1);
if any(dmask)
    lgamma_fl = lgamma_fl(dmask);
end
fdata = [real(lgamma_fl)];

R = sparse(diag(1./(xdata))); % scaling matrix
toastWriteVector('./data/fl_sim_mod_mat.fem',fdata)
% figure;
% subplot(2,2,1);imagesc(reshape(xdata,nM,nQ));colormap(gray);colorbar;axis square;title('Excitation Data RtN map');
% subplot(2,2,2);imagesc(reshape(fdata,nM,nQ));colormap(gray);colorbar;axis square;title('Fluoresence Data RtN map');
% subplot(2,2,3);imagesc(reshape(R*fdata,nM,nQ));colormap(gray);colorbar;axis square;title('Ratio of Fluoresence/Excitation');


%% inverse solution: reconstruct chromophore concentration

r = R*fdata;  % rescale fluorescence data by excitation. Y 通过激发重新缩放荧光数据 


asol1 = smat\mvec; % adjoint 

%% Jacobian
Jr = zeros(nM*nQ,grd(1)*grd(2)); %1024*3096
for i = 1:nQ
    for j = 1:nM
        tmp_h = full(sol1(:,i).*asol1(:,j));
        tmp_s = basis.Map ('M->B',tmp_h); %B:grid M:mesh S:basis
        Jr(( i-1)*nM +j,:) = tmp_s; %每一行
    end
end
Jr = R*Jr; % rescale Jacobian 缩放雅各比矩阵

% column scaling if required S: 列归一化矩阵
% c = ones(nsol,1);
% for k = 1:nsol
%     c(k) = norm(Jr(:,k));
% end
% S = sparse(diag(1./c));
% S = 0.001*speye(nsol); % reset to identity (i.e. no image scaling)
% Jr = Jr * S; % 进行列归一化

b = r;
A = Jr;



end