% clear all
% close all
% obj = readObj('100000newpoints3D_orig.obj');
v = obj.v;
f = obj.f.v;

% img = surf2volz(v,f,0:0.1:20,0:0.1:20,0:0.1:11);
% figure;imshow3Dfull(img)

start_point = [0,0,0];
end_point = [54,54,14];
grd = end_point;
max_volume = 0.5;
nodesize = 1;

opt = 0.2;
[v1,f1] = remeshsurf(v,f,opt);

% [v1,f1] = meshcheckrepair(v,f,'deep');
% [v1,f1] = meshresample(v,f,0.2);
% figure;plotmesh(v1,f1)


[newnode,newelem] = surf2mesh(v1,f1,start_point,end_point,1,max_volume);
figure;plotmesh(newnode,newelem)
% figure;plotmesh(v1,f1)


eltp = ones(length(newelem),1)*3;

n_node = length(newnode(:,1));
n_element = length(newelem(:,1)); 

toast_mesh = toastMesh(newnode,newelem,eltp);

recon_coarseness = 1;
recon_grd = round((grd - 1) ./ recon_coarseness) + 1;
recon_basis = toastBasis(toast_mesh, recon_grd);

% figure;plotmesh(node,elem);
%% source/detector config
s_edge = 2;
source = [10, 10];
% source = length(node);
total_source = prod(source);
source_region = [min(newnode(:,1)) + s_edge , max(newnode(:,1)) - s_edge; min(newnode(:,2)) + s_edge, max(newnode(:,2)) - s_edge];
scan_x_range = [source_region(1,1) source_region(1,2)];
scan_y_range = [source_region(2,1) source_region(2,2)];
[source_x, source_y] = meshgrid(linspace(scan_x_range(1), scan_x_range(2), source(1)), linspace(scan_y_range(1), scan_y_range(2), source(2)));
s_list_xy = [reshape(source_x,[],1) reshape(source_y,[],1)];
s_list_z = min(newnode(:,3)) * ones(prod(source),1);
source_position = [s_list_xy(:,1) s_list_xy(:,2) s_list_z];
% source_position = node;

detector = [10, 10];
d_edge = 2; %mm
% d_edge = 0; %mm
detector_region = [min(newnode(:,1)) + d_edge , max(newnode(:,1)) - d_edge; min(newnode(:,2)) + d_edge, max(newnode(:,2)) - d_edge];

det_x_range = [detector_region(1,1) detector_region(1,2)];
det_y_range = [detector_region(2,1) detector_region(2,2)];
[det_x, det_y] = meshgrid(linspace(det_x_range(1), det_x_range(2), detector(1)), linspace(det_y_range(1), det_y_range(2), detector(2)));
d_list_xy = [reshape(det_x,[],1) reshape(det_y,[],1)];
% d_list_z = max(node(:,3)) * ones(prod(detector),1);
d_list_z = max(newnode(:,3)) * ones(prod(detector),1);
detector_position = [d_list_xy(:,1) d_list_xy(:,2) d_list_z];

figure;plotmesh(newnode, newelem);
hold on
plot3(source_position(:,1),source_position(:,2),source_position(:,3),'ro','MarkerFaceColor','r');
plot3(detector_position(:,1),detector_position(:,2),detector_position(:,3), 'bs','MarkerFaceColor','b');


%% opt config
ne = toast_mesh.ElementCount;
nv = toast_mesh.NodeCount;

mua = 0.001 * ones(nv,1);
mus = 1.1 * ones(nv,1);  % scattering coefficient [1/mm]
% kap = 1./(3*(mua+mus)); % diffusion coefficient [mm]
kap = 1./(3*(mus));
freq = 0;
refind = 1.4;   % refractive index
ref = refind * ones(nv,1);
c0 = 0.3;       % speed of light in vacuum [mm/ps]
c = c0/refind; % speed of light in the medium [mm/ps]

%% target config
z_node = [];
for i=1:n_node
    if newnode(i,3) <= -3
        z_node = [z_node; [newnode(i,1), newnode(i,2), newnode(i,3), i]];
    end
end

n_circle = [10, 10];
radius = 5; %mm

center = []; %origin node(s) index
circle_x_range = [min(newnode(:,1)) + radius, max(newnode(:,1)) - radius];
circle_y_range = [min(newnode(:,2)) + radius, max(newnode(:,2)) - radius];
[coordinate_x, coordinate_y] = meshgrid(linspace(circle_x_range(1), circle_x_range(2), n_circle(1)), linspace(circle_y_range(1), circle_y_range(2), n_circle(2)));
cor_xy = [reshape(coordinate_x,[],1) reshape(coordinate_y,[],1)];
cor_z = zeros(prod(n_circle),1);
circle_position = [cor_xy(:,1) cor_xy(:,2) cor_z];
cc = 45;

pattern_node = [];
for j=1:length(z_node)
    if ( (z_node(j,1) - circle_position(cc,1))^2 + (z_node(j,2) - circle_position(cc,2))^2) <= radius^2
        pattern_node = [pattern_node; z_node(j,:)];
    end
end
for k=1:length(pattern_node)
   mua(pattern_node(k,4)) = 0.1;
end

scale = 10;
show_grd = scale * (round(grd)+1) .* [1,1,1];
basis = toastBasis(toast_mesh, show_grd);
bmua = basis.Map('M->B',mua);   
bmua = reshape(bmua,show_grd);
figure;imshow3Dfull(bmua);colorbar
figure;plotmesh([newnode,mua],newelem(:,1:4))
%% source config    

toast_mesh.SetQM(source_position,detector_position);

qvec = real(toast_mesh.Qvec('Neumann','Gaussian',1));
mvec = real(toast_mesh.Mvec('Gaussian',1,refind));
noiselevel = 0;



%% generate CW data

[S,B,alpha] = dotSysmat (toast_mesh,mua,mus,ref,freq);
phi = S\qvec;
data = reshape(log(abs(mvec.'* phi)), [] , 1);
data = data + mean(data(:)).*noiselevel.*rand(size(data));

% Now build S from individual components
% S1 = toast_mesh.SysmatComponent('PFF',mua.*c,'EL');  % absorption term
S1 = toast_mesh.SysmatComponent('PFF',mua.*c);  % absorption term
S2 = toast_mesh.SysmatComponent('PDD',kap.*c);  % diffusion component
S3 = toast_mesh.SysmatComponent('BndPFF',c./(2*alpha)); % boundary term

S_FF = toast_mesh.SysmatComponent('FF');  
S1_ff = S_FF .*  mua.*c;

%% recon
recon_coarseness = 1;
recon_grd = round((grd - 1) ./ recon_coarseness) + 1;
recon_basis = toastBasis(toast_mesh, recon_grd);

mua_r = 0.001 * ones(nv,1); %initial guess
mus_r = 1.1 * ones(nv,1);

mua_r_grd_0 = recon_basis.Map('M->B',mua_r);
mus_r_0 = mus_r;

%% recon test
pj_error = [];
iteration = 30;
mua_err = [];
lambda.type = 'JtJ';
lambda.value = 10;
dmask = ones(prod(source), prod(detector));
mua_r_grd_set = {};
% J_self = GenJ_CW(toast_mesh, basis, qvec, mvec, mua_r, mus_r, ref);
% J_solus = toastJacobianCW2(toast_mesh, basis, qvec, mvec,dmask, mua_r, mus_r, ref,'none');

for it = 1:iteration
    J = GenJ_CW(toast_mesh, recon_basis, qvec, mvec, mua_r, mus_r, ref);
%     J = toastJacobianCW2(toast_mesh, basis, qvec, mvec,dmask, mua_r, mus_r, ref,'none');
%     J = toastJacobianCW(toast_mesh, basis, qvec, mvec, mua_r, mus_r, ref,'CG');

    smat = dotSysmat(toast_mesh, mua_r, mus_r, ref);
    phi = full(smat \ qvec);
    simu_mea = mvec.'* phi;
%     for n = 1:prod(source)
%         simu_mea(:,n) = simu_mea(:,n) .* cali_coe(n);
%     end
    proj = reshape(log(abs(simu_mea)),[],1);
    
    
    data_diff = (data - proj);
%     sd = ones(size(proj)) * norm(data_diff);
    pj_error = [pj_error, sum((data_diff.^2))];
    
    disp('---------------------------------');
    disp(['Iteration Number          = ' num2str(it)]);
    disp(['Projection error          = ' num2str(pj_error(end))]);
    
    if it~= 1
        p = (pj_error(end-1) - pj_error(end)) * 100 / pj_error(end-1);
        disp(['Projection error change   = ' num2str(p) '%']);
        if (~isreal(p)) || (p <= 2)
            break
        end
    end
            
    mua_r_grd = recon_basis.Map('M->B',mua_r);
%     mua_err = [mua_err, mse(bmua_r,mua_r_grd)];
    
    nn = length(mua_r_grd);
    m = length(data_diff);
%data normalization
%     for i = 1:m
%         J(i,:) = J(i,:) / sd(i);
%     end

%paramter normalisation
    for i = 1:nn
        J(:,i) = J(:,i) .* mua_r_grd(i,1);
    end
    
    if it ~= 1
        lambda.value = lambda.value ./ 10^0.25;
    end
    
    [nrow, ncol] = size(J);
   if strcmp(lambda.type, 'JtJ') 

        Hess = (J'* J);
        
        reg_mua = lambda.value * max(diag(Hess));
        %     reg_mua = lambda;
        reg = ones(ncol,1);
        reg = reg .* reg_mua;
        
        disp(['Mua Regularization        = ' num2str(reg(1,1))]);
        
        for i = 1:ncol
            Hess(i, i) = Hess(i, i) + reg(i);
        end
        
        delta_mua = mua_r_grd - mua_r_grd_0;
%         foo = Hess^(-1) * (J' * data_diff - lambda.value * delta_mua);
        foo = Hess\J'*data_diff;
   else

        Hess = (J* J');
        
        reg_mua = lambda.value * max(diag(Hess));
        reg = ones(nrow,1);
        reg = reg .* reg_mua;
        
        disp(['Mua Regularization        = ' num2str(reg(1,1))]);
        
        for i = 1:nrow
            Hess(i, i) = Hess(i, i) + reg(i);
        end
        
        foo = J' * (Hess \ data_diff);
   end


    foo = foo .* mua_r_grd;
    
    mua_r_grd = mua_r_grd + foo;
%     mua_r_grd = exp(mua_r_grd);
    mua_r_grd_set{it} = mua_r_grd;
    mua_r = recon_basis.Map('B->M',mua_r_grd);

end

r_grd = (round(grd)+1) .* [10,10,10];
r_basis = toastBasis(toast_mesh,r_grd);
mua_r_grd = r_basis.Map('M->B', mua_r);
mua_r_grd = reshape(mua_r_grd,r_grd);
norm_mua_grd = mua_r_grd ./ max(mua_r_grd(:));
figure;imshow3Dfull(norm_mua_grd);colorbar




