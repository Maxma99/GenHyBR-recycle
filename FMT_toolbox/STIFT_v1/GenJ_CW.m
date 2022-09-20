function J = GenJ_CW(mesh, basis, qvec, mvec, mua, mus, ref)
[S,~,~] = dotSysmat (mesh,mua,mus,ref,0);
dphi = S\qvec;
aphi = S\mvec;
grd = basis.Dims();
[~,nQ] = size(qvec);
[~,nM] = size(mvec);
nsol_J = prod(grd);


J = zeros(nQ*nM, nsol_J);
k = 1;
cafield = zeros(nsol_J, nM);
for i = 1:nM
    cafield(:,i) = basis.Map('M->B', aphi(:,i));
end

for i = 1:nQ
    cdfield = basis.Map('M->B', dphi(:,i));
    proj = ProjectSingle(mvec, dphi(:,i));
    for j = 1:nM
        pmdf_mua = PMDF_mua(cdfield, cafield(:,j));
        tmp_l = PMDF_log(pmdf_mua, proj(j));
        J(k,:) = tmp_l;
        k = k+1;
    end
end

sz = 1;

for i = 1:length(grd)
    sz = sz * (grd(i) / (grd(i) - 1));
end
J = J * sz;
end

function opt = PMDF_mua(dphi, aphi)
opt = -dphi .* aphi;
end

function opt = PMDF_log(pmdf, m)
opt = pmdf / m;
end

function opt = ProjectSingle(mvec, phi)
opt = zeros(size(mvec,2),1);
for i = 1: size(mvec,2)
    opt(i) = sum(phi .* mvec(:,i));
end
end